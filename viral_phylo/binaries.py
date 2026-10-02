"""Binary and interpreter discovery utilities for structural phylogenetics and ML tools."""

import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional

# Tools are resolved by the same search order, so the candidate paths live as data
# rather than being copy-pasted per tool.
_EXTRA_CANDIDATES: Dict[str, List[str]] = {
    "mafft": ["/usr/local/bin/mafft"],
}

# How each tool reports its version. FoldMason, MMseqs2 and Foldseek use a bare
# `version` subcommand; IQ-TREE and MAFFT use `--version` (MAFFT writes to stderr).
# VeryFastTree is a parallel reimplementation of FastTree with the same interface;
# prefer it, then fall back to the various FastTree binary spellings.
FASTTREE_PREFERENCE = ("veryfasttree", "VeryFastTree", "FastTreeMP", "FastTree", "fasttree")

_VERSION_ARGS: Dict[str, List[str]] = {
    "foldmason": ["version"],
    "mmseqs": ["version"],
    "foldseek": ["version"],
    "iqtree": ["--version"],
    "mafft": ["--version"],
    "veryfasttree": ["-help"],
    "fasttree": ["-help"],
}

EXTERNAL_TOOLS = ("foldmason", "mafft", "iqtree", "mmseqs", "foldseek", "fasttree")


def _find_tool(name: str, override: Optional[str] = None) -> str:
    """Resolve a tool binary from an explicit override, PATH, or known conda layouts.

    Falls back to the bare name so that any failure surfaces later as a
    CalledProcessError naming the tool, rather than as a confusing path error.
    """
    if override and os.path.isfile(override) and os.access(override, os.X_OK):
        return override
    if shutil.which(name):
        return name
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    candidates = [
        os.path.join(conda_prefix, "bin", name) if conda_prefix else "",
        # The directory holding the running interpreter. When the environment's
        # python or console script is invoked by absolute path, conda is not
        # activated and CONDA_PREFIX is unset, but the tools sit right alongside it.
        os.path.join(os.path.dirname(sys.executable), name),
        os.path.expanduser(f"~/miniconda3/envs/spt/bin/{name}"),
        os.path.expanduser(f"~/miniconda3/bin/{name}"),
        f"/opt/homebrew/bin/{name}",
    ]
    candidates.extend(_EXTRA_CANDIDATES.get(name, []))
    for candidate in candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return name


def find_iqtree_bin() -> str:
    """Find iqtree binary in PATH or conda environments."""
    return _find_tool("iqtree")


def find_foldmason_bin() -> str:
    """Find foldmason binary in PATH or conda environments."""
    return _find_tool("foldmason")


def find_mafft_bin(mafft_bin: Optional[str] = None) -> str:
    """Find mafft binary in specified path, PATH, or conda environments."""
    return _find_tool("mafft", override=mafft_bin)


def find_mmseqs_bin(mmseqs_bin: Optional[str] = None) -> str:
    """Find the mmseqs binary (amino-acid sequence clustering and search)."""
    return _find_tool("mmseqs", override=mmseqs_bin)


def find_foldseek_bin(foldseek_bin: Optional[str] = None) -> str:
    """Find the foldseek binary (structural/3Di clustering and search)."""
    return _find_tool("foldseek", override=foldseek_bin)


def tool_is_available(name: str, override: Optional[str] = None) -> bool:
    """Whether a tool can actually be executed, rather than merely named."""
    resolved = _find_tool(name, override=override)
    if os.path.sep in resolved:
        return os.path.isfile(resolved) and os.access(resolved, os.X_OK)
    return shutil.which(resolved) is not None


def tool_version(name: str, override: Optional[str] = None) -> Optional[str]:
    """Return a tool's reported version string, or None if it cannot be run.

    Nothing in this package captured tool versions before, which made a completed
    run impossible to reproduce: the manifest now records these.
    """
    if name == "fasttree":
        # Report the binary that will actually run (VeryFastTree is preferred), not
        # whichever FastTree happens to be called "fasttree" on PATH.
        resolved = find_fasttree_bin(override)
        if not resolved:
            return None
        try:
            proc = subprocess.run([resolved, "-help"], capture_output=True, text=True, timeout=30, check=False)
        except (OSError, subprocess.SubprocessError):
            return None
        import re as _re
        text = (proc.stdout or "") + "\n" + (proc.stderr or "")
        match = _re.search(r"(VeryFastTree|FastTree)\s+[0-9][0-9.]*[^\n]*", text)
        return match.group(0).strip() if match else None
    if not tool_is_available(name, override=override):
        return None
    resolved = _find_tool(name, override=override)
    args = _VERSION_ARGS.get(name, ["--version"])
    try:
        proc = subprocess.run(
            [resolved, *args],
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    # MAFFT reports its version on stderr; the others use stdout.
    output = (proc.stdout or "").strip() or (proc.stderr or "").strip()
    if not output:
        return None
    return output.splitlines()[0].strip()


def resolve_tool_versions(names=None, overrides: Optional[Dict[str, str]] = None) -> Dict[str, Dict[str, object]]:
    """Resolve path, availability and version for each external tool.

    Returned shape is intended for the run manifest:
    ``{"mafft": {"path": ..., "version": ..., "available": True}, ...}``
    """
    overrides = overrides or {}
    report: Dict[str, Dict[str, object]] = {}
    for name in (names or EXTERNAL_TOOLS):
        override = overrides.get(name)
        if name == "fasttree":
            path = find_fasttree_bin(override)
            report[name] = {"path": path, "available": bool(path),
                            "version": tool_version("fasttree", override) if path else None}
            continue
        available = tool_is_available(name, override=override)
        report[name] = {
            "path": _find_tool(name, override=override),
            "available": available,
            "version": tool_version(name, override=override) if available else None,
        }
    return report


_TORCH_PROBE_CACHE: Dict[str, Dict[str, object]] = {}


def probe_torch_python(python: str, timeout: int = 180) -> Dict[str, object]:
    """Report whether an interpreter can run the embedding code, and on what device.

    Returns ``{"usable", "cuda", "torch", "transformers", "missing", "error"}``.
    Cached per interpreter, since importing torch and transformers takes seconds.
    """
    if python in _TORCH_PROBE_CACHE:
        return _TORCH_PROBE_CACHE[python]
    script = (
        "import importlib, json\n"
        "out = {'missing': []}\n"
        "for m in ('torch', 'transformers', 'scipy'):\n"
        "    try:\n"
        "        mod = importlib.import_module(m); out[m] = getattr(mod, '__version__', '?')\n"
        "    except Exception:\n"
        "        out['missing'].append(m)\n"
        "try:\n"
        "    import torch; out['cuda'] = bool(torch.cuda.is_available())\n"
        "except Exception:\n"
        "    out['cuda'] = False\n"
        "print(json.dumps(out))\n"
    )
    info: Dict[str, object] = {"usable": False, "cuda": False, "torch": None,
                               "transformers": None, "missing": [], "error": None}
    try:
        proc = subprocess.run([python, "-c", script], capture_output=True, text=True,
                              timeout=timeout, check=False)
        import json as _json
        payload = _json.loads(proc.stdout.strip().splitlines()[-1]) if proc.stdout.strip() else {}
        info.update({k: payload.get(k) for k in ("torch", "transformers")})
        info["missing"] = payload.get("missing", ["torch", "transformers", "scipy"])
        info["cuda"] = bool(payload.get("cuda"))
        info["usable"] = not info["missing"]
        if not payload:
            info["error"] = (proc.stderr or "").strip().splitlines()[-1:] or ["no output"]
    except (OSError, subprocess.SubprocessError, ValueError) as exc:
        info["error"] = f"{type(exc).__name__}: {exc}"
        info["missing"] = ["torch", "transformers", "scipy"]
    _TORCH_PROBE_CACHE[python] = info
    return info


def find_torch_python() -> str:
    """Find a Python interpreter to use for ML/embeddings.

    Candidates, in priority order:
    1. The install script environment 'spt' (e.g. envs/spt/bin/python)
    2. The currently executing Python interpreter (sys.executable)
    3. The user's active Conda environment ($CONDA_PREFIX/bin/python)
    4. Fallback candidate environments

    The chosen interpreter must actually be able to run the embedding code, so each
    candidate is probed: it has to import torch, transformers and scipy. Among those,
    one where CUDA works is preferred. Previously the active conda environment was
    returned without checking, so launching from an environment that lacked
    transformers (e.g. anaconda base) crashed the embedding stage only after every
    tree had already been built.
    """
    def is_valid_python(path: str) -> bool:
        return bool(path and os.path.isfile(path) and os.access(path, os.X_OK))

    # Candidate order is the historical priority; capability then decides.
    candidates = []
    if os.environ.get("CONDA_DEFAULT_ENV") == "spt" and os.environ.get("CONDA_PREFIX"):
        candidates.append(os.path.join(os.environ["CONDA_PREFIX"], "bin", "python"))
    conda_base_hints = [
        os.environ.get("CONDA_PREFIX", ""),
        os.path.expanduser("~/miniconda3"),
        os.path.expanduser("~/anaconda3"),
        os.path.expanduser("~/miniforge3"),
        os.path.expanduser("~/micromamba"),
        "/opt/conda",
        "/opt/homebrew/Caskroom/miniconda/base",
    ]
    for hint in conda_base_hints:
        if not hint:
            continue
        base_dir = hint.split("/envs/")[0] if "/envs/" in hint else hint
        candidates.append(os.path.join(base_dir, "envs", "spt", "bin", "python"))
    candidates.append(sys.executable)
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(os.path.join(os.environ["CONDA_PREFIX"], "bin", "python"))
    candidates += [
        os.path.expanduser("~/miniconda3/envs/nipah/bin/python"),
        os.path.expanduser("~/miniconda3/bin/python"),
        "/usr/local/bin/python3",
    ]

    seen, ordered = set(), []
    for cand in candidates:
        if is_valid_python(cand):
            real = os.path.realpath(cand)
            if real not in seen:
                seen.add(real)
                ordered.append(cand)

    probes = [(cand, probe_torch_python(cand)) for cand in ordered]
    for cand, info in probes:
        if info["usable"] and info["cuda"]:
            return cand
    for cand, info in probes:
        if info["usable"]:
            return cand
    # Nothing can run the embeddings. Return the historical choice so the caller's
    # preflight check reports exactly which packages are missing.
    return ordered[0] if ordered else sys.executable


def find_fasttree_bin(fasttree_bin: Optional[str] = None) -> Optional[str]:
    """Find a FastTree-compatible binary, preferring VeryFastTree.

    FastTree ships under several names and VeryFastTree is a drop-in parallel
    reimplementation, so try them in order of preference. Returns None when none is
    installed, rather than a bare name, so callers can report the absence properly
    instead of failing later with a confusing exec error.
    """
    if fasttree_bin and os.path.isfile(fasttree_bin) and os.access(fasttree_bin, os.X_OK):
        return fasttree_bin
    for name in FASTTREE_PREFERENCE:
        if shutil.which(name):
            return name
        candidate = os.path.join(os.path.dirname(sys.executable), name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None
