#!/usr/bin/env python3
"""Run manifest: a durable record of what a pipeline run actually did.

Before this existed there was no record of a run anywhere. A completed run could
not be reproduced - not even its own command line was recoverable - and stages
that silently produced nothing still reported success.

The manifest addresses both. It captures the resolved invocation, tool versions,
input provenance and a per-stage status, and it is rewritten after every stage so
that a killed or crashed run still leaves a usable record on disk.

Two files are written next to the run:
  run_manifest.json  - machine readable, the authoritative record
  run_manifest.txt   - the same content as a scannable summary
"""

import getpass
import hashlib
import json
import os
import platform
import re
import socket
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

SCHEMA_VERSION = 1

MANIFEST_JSON = "run_manifest.json"
MANIFEST_TXT = "run_manifest.txt"

# Hashing every input is the strongest provenance, but guard against pathological
# inputs; beyond these limits fall back to a name+size digest and say so.
_MAX_DIGEST_FILES = 5000
_MAX_DIGEST_BYTES = 4 * 1024 ** 3

_HF_SNAPSHOT_RE = re.compile(
    r"(?P<cache>.*/huggingface/hub)/datasets--(?P<org>[^/]+)--(?P<name>[^/]+)/snapshots/(?P<revision>[0-9a-f]+)"
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _package_version_and_commit(repo_dir: Optional[Path] = None) -> Dict[str, Optional[str]]:
    """Package version plus the git commit it was run from.

    Mirrors the approach already used to stamp the interactive HTML banner.
    """
    if repo_dir is None:
        repo_dir = Path(__file__).resolve().parent.parent

    version = None
    try:
        from viral_phylo import __version__ as pkg_version
        version = pkg_version
    except Exception:
        pyproject = repo_dir / "pyproject.toml"
        if pyproject.exists():
            match = re.search(r'version\s*=\s*["\']([^"\']+)["\']', pyproject.read_text(encoding="utf-8"))
            if match:
                version = match.group(1)

    commit, dirty = None, None
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            text=True, stderr=subprocess.DEVNULL, timeout=15,
        ).strip()
        status = subprocess.check_output(
            ["git", "-C", str(repo_dir), "status", "--porcelain"],
            text=True, stderr=subprocess.DEVNULL, timeout=15,
        )
        dirty = bool(status.strip())
    except Exception:
        pass

    return {"version": version, "git_commit": commit, "git_dirty": dirty}


def describe_input_directory(path: str, extensions=(".pdb", ".cif", ".mmcif")) -> Dict[str, Any]:
    """Summarise an input directory: file inventory, size, content digest, provenance.

    Detects when the inputs come from a Hugging Face snapshot and records the
    dataset and immutable revision, which pins the inputs far better than a path.
    """
    info: Dict[str, Any] = {"path": os.path.abspath(path) if path else None}
    if not path or not os.path.isdir(path):
        info["exists"] = False
        return info

    info["exists"] = True
    files: List[str] = []
    for entry in sorted(os.listdir(path)):
        if entry.lower().endswith(tuple(extensions)):
            files.append(os.path.join(path, entry))

    total = 0
    for f in files:
        try:
            total += os.path.getsize(f)
        except OSError:
            pass

    info["file_count"] = len(files)
    info["total_bytes"] = total
    info["extensions"] = sorted({os.path.splitext(f)[1].lower() for f in files})

    digest = hashlib.sha256()
    if len(files) <= _MAX_DIGEST_FILES and total <= _MAX_DIGEST_BYTES:
        for f in files:
            digest.update(os.path.basename(f).encode("utf-8"))
            try:
                with open(f, "rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(chunk)
            except OSError:
                continue
        info["digest_kind"] = "sha256-content"
    else:
        for f in files:
            try:
                digest.update(f"{os.path.basename(f)}:{os.path.getsize(f)}".encode("utf-8"))
            except OSError:
                continue
        info["digest_kind"] = "sha256-name-size"
    info["digest"] = digest.hexdigest()

    # Check the given path before resolving symlinks: a Hugging Face snapshot
    # directory is itself a tree of symlinks into a content-addressed blobs/ store,
    # so realpath() discards exactly the dataset and revision we want to record.
    abs_path = os.path.abspath(path)
    probe_paths = [abs_path]
    if files:
        probe_paths.append(os.path.abspath(files[0]))
        probe_paths.append(os.path.realpath(files[0]))
    else:
        probe_paths.append(os.path.realpath(path))

    for probe in probe_paths:
        match = _HF_SNAPSHOT_RE.match(probe)
        if match:
            tail = probe.split(match.group("revision"), 1)[-1].strip("/")
            info["huggingface"] = {
                "dataset": f"{match.group('org')}/{match.group('name')}",
                "revision": match.group("revision"),
                "subdir": (tail if probe == abs_path else tail.rsplit("/", 1)[0]) or None,
            }
            break

    if files and os.path.realpath(files[0]) != os.path.abspath(files[0]):
        info["resolved_source_dir"] = os.path.dirname(os.path.realpath(files[0]))

    return info


class RunManifest:
    """Accumulates a record of one pipeline run and writes it out incrementally."""

    def __init__(self, output_dir: str, quiet: bool = False):
        self.output_dir = output_dir
        self.quiet = quiet
        self._t0 = time.time()
        self.data: Dict[str, Any] = {
            "schema_version": SCHEMA_VERSION,
            "started_at": _utc_now(),
            "finished_at": None,
            "status": "running",
            "invocation": {},
            "environment": {},
            "versions": {},
            "input": {},
            "qc": {},
            "alignment": {},
            "clustering": {},
            "stages": [],
            "outputs": [],
            "warnings": [],
        }
        os.makedirs(output_dir, exist_ok=True)
        self.record_environment()

    # ---------------------------------------------------------------- recording

    def record_environment(self) -> None:
        env = {
            "hostname": socket.gethostname(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "python_executable": sys.executable,
            "cwd": os.getcwd(),
        }
        try:
            env["user"] = getpass.getuser()
        except Exception:
            pass
        self.data["environment"] = env
        self.data["versions"]["package"] = _package_version_and_commit()

    def record_invocation(self, argv: Optional[List[str]] = None, args: Any = None) -> None:
        """Record the command line and the fully resolved arguments.

        Resolved arguments matter as much as argv: they make every default explicit,
        so the record does not silently change meaning when a default later changes.
        """
        invocation: Dict[str, Any] = {"argv": list(argv if argv is not None else sys.argv)}
        invocation["command"] = " ".join(invocation["argv"])
        if args is not None:
            resolved = {}
            for key, value in sorted(vars(args).items()):
                if callable(value):
                    continue
                resolved[key] = value if isinstance(value, (str, int, float, bool, type(None), list)) else str(value)
            invocation["resolved_args"] = resolved
        self.data["invocation"] = invocation

    def record_tool_versions(self, overrides: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
        from viral_phylo.binaries import resolve_tool_versions

        tools = resolve_tool_versions(overrides=overrides)
        self.data["versions"]["tools"] = tools
        for name, info in tools.items():
            if not info.get("available"):
                self.add_warning(f"external tool '{name}' was not found on PATH")
        return tools

    def record_input(self, path: str, **extra: Any) -> Dict[str, Any]:
        info = describe_input_directory(path)
        info.update(extra)
        self.data["input"] = info
        return info

    def record_qc(self, excluded: List[Dict[str, Any]], retained: Optional[int] = None, **extra: Any) -> None:
        self.data["qc"] = {"excluded_count": len(excluded), "excluded": excluded, "retained": retained, **extra}

    def record_alignment(self, info: Dict[str, Any]) -> None:
        """Merge into the alignment record rather than replacing it.

        Several stages contribute here - the whole-set gate writes its verdict before
        the align stage closes and reports the aligner - so a wholesale assignment
        would silently discard whichever was written first.
        """
        self.data.setdefault("alignment", {}).update(info)

    def record_clustering(self, info: Dict[str, Any]) -> None:
        self.data["clustering"] = info

    def add_warning(self, message: str) -> None:
        if message not in self.data["warnings"]:
            self.data["warnings"].append(message)

    def record_stage(
        self,
        name: str,
        status: str,
        reason: Optional[str] = None,
        duration_s: Optional[float] = None,
        command: Optional[Any] = None,
        outputs: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        record = {
            "name": name,
            "status": status,
            "reason": reason,
            "duration_s": round(duration_s, 2) if duration_s is not None else None,
            "command": " ".join(command) if isinstance(command, (list, tuple)) else command,
            "outputs": outputs or [],
            "recorded_at": _utc_now(),
        }
        self.data["stages"].append(record)
        self.write()
        return record

    @contextmanager
    def stage(self, name: str, expect_outputs: Optional[List[str]] = None, command: Optional[Any] = None):
        """Time and record a pipeline stage.

        A stage that raises is recorded as failed and the exception propagates. A
        stage that completes but produces none of the outputs it declared is recorded
        as skipped with a reason and announced - the pipeline used to pass straight
        over exactly this case while still printing a success banner.
        """
        started = time.time()
        state: Dict[str, Any] = {"outputs": [], "reason": None, "skipped": False}
        try:
            yield state
        except Exception as exc:
            self.record_stage(
                name, "failed", reason=f"{type(exc).__name__}: {exc}",
                duration_s=time.time() - started, command=command,
            )
            self._announce(f"[!] Stage '{name}' FAILED: {type(exc).__name__}: {exc}")
            raise

        produced = [p for p in (expect_outputs or []) + list(state["outputs"]) if p and os.path.exists(p)]
        declared = (expect_outputs or []) + list(state["outputs"])

        if state["skipped"]:
            status, reason = "skipped", state["reason"] or "stage reported itself skipped"
        elif declared and not produced:
            status = "skipped"
            reason = state["reason"] or f"produced none of its {len(declared)} expected output(s)"
        else:
            status, reason = "ok", state["reason"]

        self.record_stage(
            name, status, reason=reason, duration_s=time.time() - started,
            command=command, outputs=produced,
        )
        if status == "skipped":
            self._announce(f"[!] Notice: stage '{name}' produced no output - {reason}")

    def record_outputs(self, paths: List[str]) -> None:
        inventory = []
        for path in paths:
            if path and os.path.exists(path):
                inventory.append({"path": os.path.abspath(path), "bytes": os.path.getsize(path)})
        self.data["outputs"] = inventory

    def install_exit_handlers(self) -> None:
        """Make sure the manifest never stays at "running" once the process is gone.

        An uncaught exception is recorded and the run finished as "failed"; a SIGTERM
        (``kill``) or any other exit before ``finish()`` finishes it as "interrupted".
        The killed 25 Sep follow-up run was left reading "running" indefinitely.
        """
        import atexit
        import signal
        import traceback

        previous_hook = sys.excepthook

        def _excepthook(exc_type, exc, tb):
            if self.data.get("status") == "running":
                self.data["error"] = {
                    "type": exc_type.__name__, "message": str(exc),
                    "traceback": "".join(traceback.format_exception(exc_type, exc, tb))[-4000:],
                }
                self.add_warning(f"run aborted by an uncaught {exc_type.__name__}: {exc}")
                try:
                    self.finish("failed")
                except Exception:
                    pass
            previous_hook(exc_type, exc, tb)

        def _atexit():
            if self.data.get("status") == "running":
                try:
                    self.finish("interrupted")
                except Exception:
                    pass

        def _on_sigterm(signum, frame):
            raise SystemExit(128 + signum)

        sys.excepthook = _excepthook
        atexit.register(_atexit)
        try:
            signal.signal(signal.SIGTERM, _on_sigterm)
        except (ValueError, OSError):
            pass  # not in the main thread

    def finish(self, status: str = "completed") -> str:
        self.data["status"] = status
        self.data["finished_at"] = _utc_now()
        self.data["duration_s"] = round(time.time() - self._t0, 2)
        return self.write()

    # ------------------------------------------------------------------ writing

    def _announce(self, message: str) -> None:
        if not self.quiet:
            print(f"[Manifest] {message}")

    @property
    def json_path(self) -> str:
        return os.path.join(self.output_dir, MANIFEST_JSON)

    @property
    def text_path(self) -> str:
        return os.path.join(self.output_dir, MANIFEST_TXT)

    def write(self) -> str:
        os.makedirs(self.output_dir, exist_ok=True)
        with open(self.json_path, "w", encoding="utf-8") as handle:
            json.dump(self.data, handle, indent=2, default=str)
        with open(self.text_path, "w", encoding="utf-8") as handle:
            handle.write(self.render_text())
        return self.json_path

    def render_text(self) -> str:
        d = self.data
        lines: List[str] = []
        add = lines.append

        add("=" * 78)
        add("STRUCTURAL PHYLOGENETICS - RUN MANIFEST")
        add("=" * 78)
        add(f"Status      : {d['status']}")
        add(f"Started     : {d['started_at']}")
        add(f"Finished    : {d.get('finished_at') or '(in progress)'}")
        if d.get("duration_s") is not None:
            add(f"Duration    : {d['duration_s']} s")

        pkg = d.get("versions", {}).get("package", {})
        if pkg:
            dirty = " (working tree dirty)" if pkg.get("git_dirty") else ""
            add(f"Toolkit     : v{pkg.get('version')} @ {(pkg.get('git_commit') or '?')[:12]}{dirty}")

        env = d.get("environment", {})
        if env:
            add(f"Host        : {env.get('user')}@{env.get('hostname')}  python {env.get('python')}")
            add(f"Working dir : {env.get('cwd')}")

        inv = d.get("invocation", {})
        if inv.get("command"):
            add("")
            add("-- Invocation " + "-" * 64)
            add(inv["command"])

        tools = d.get("versions", {}).get("tools", {})
        if tools:
            add("")
            add("-- External tools " + "-" * 60)
            for name, info in tools.items():
                mark = "ok " if info.get("available") else "!! "
                add(f"  {mark}{name:<12} {info.get('version') or 'NOT FOUND'}")

        src = d.get("input", {})
        if src:
            add("")
            add("-- Input " + "-" * 69)
            add(f"  path        : {src.get('path')}")
            if src.get("file_count") is not None:
                add(f"  structures  : {src.get('file_count')} ({src.get('total_bytes', 0) / 1e6:.1f} MB)")
            if src.get("digest"):
                add(f"  digest      : {src['digest'][:16]}...  [{src.get('digest_kind')}]")
            hf = src.get("huggingface")
            if hf:
                add(f"  huggingface : {hf['dataset']} @ {hf['revision'][:12]}")

        qc = d.get("qc") or {}
        if qc.get("excluded_count"):
            add("")
            add("-- QC " + "-" * 72)
            add(f"  excluded    : {qc['excluded_count']} (retained {qc.get('retained')})")
            for item in qc.get("excluded", [])[:10]:
                add(f"    - {item.get('taxon_id')}: {item.get('reason')}")
            if qc["excluded_count"] > 10:
                add(f"    ... and {qc['excluded_count'] - 10} more (see JSON)")

        stages = d.get("stages", [])
        if stages:
            add("")
            add("-- Stages " + "-" * 68)
            for stage in stages:
                mark = {"ok": "ok ", "skipped": ">> ", "failed": "!! "}.get(stage["status"], "   ")
                dur = f"{stage['duration_s']:>8.2f}s" if stage.get("duration_s") is not None else " " * 9
                add(f"  {mark}{stage['name']:<34}{dur}  {stage['status']}")
                if stage.get("reason"):
                    add(f"        reason: {stage['reason']}")

        outs = d.get("outputs", [])
        if outs:
            add("")
            add("-- Outputs " + "-" * 67)
            for item in outs:
                add(f"  {item['bytes'] / 1e6:>8.2f} MB  {item['path']}")

        warnings = d.get("warnings", [])
        if warnings:
            add("")
            add("-- Warnings " + "-" * 66)
            for warning in warnings:
                add(f"  ! {warning}")

        add("")
        add("=" * 78)
        return "\n".join(lines) + "\n"


def write_alignment_info(
    output_dir: str,
    aligner: str,
    produced: Dict[str, str],
    commands: Optional[List[Any]] = None,
    mafft_matrix: Optional[str] = None,
    tool_paths: Optional[Dict[str, str]] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> str:
    """Record which aligner really produced the files in an alignment directory.

    The MAFFT path copies its output over the ``foldmason.fasta_*.fa`` names that
    downstream code expects, so the filename alone does not say what produced it.
    This sidecar does.
    """
    from viral_phylo.binaries import tool_version

    info: Dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "written_at": _utc_now(),
        "aligner": aligner,
        "mafft_matrix": mafft_matrix,
        "produced_by": produced,
        "commands": [" ".join(c) if isinstance(c, (list, tuple)) else c for c in (commands or [])],
        "tool_versions": {},
        "note": (
            "Files named foldmason.fasta_*.fa are a legacy compatibility name used by "
            "downstream consumers; 'produced_by' records the aligner that actually wrote them."
        ),
    }
    for name in {aligner, "foldmason"}:
        try:
            info["tool_versions"][name] = tool_version(name, override=(tool_paths or {}).get(name))
        except Exception:
            info["tool_versions"][name] = None
    if extra:
        info.update(extra)

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, "alignment_info.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(info, handle, indent=2, default=str)
    return path
