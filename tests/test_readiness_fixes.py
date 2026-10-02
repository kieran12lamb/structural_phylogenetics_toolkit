"""Regression tests for defects found in the pre-run review."""

import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path
from unittest import mock

REPO = Path(__file__).resolve().parent.parent


class TestMatrixAliasAbsolutePath(unittest.TestCase):
    """An absolute matrix path used to pass through un-aliased, so IQ-TREE -mset
    (which upper-cases its argument) would fail for every small cluster."""

    def test_absolute_path_is_aliased_relative_to_cwd(self):
        from viral_phylo.matrices import uppercase_matrix_alias
        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "store" / "matrices" / "Q.3Di.AF"
            src.parent.mkdir(parents=True)
            src.write_text("0.2\n0.1 5.4\n")
            cwd = os.getcwd()
            try:
                os.chdir(d)
                alias = uppercase_matrix_alias(str(src))
                self.assertEqual(alias, os.path.join("MATRICES", "Q.3DI.AF"))
                self.assertEqual(alias, alias.upper(), "alias must survive upper-casing unchanged")
                self.assertEqual(Path(alias).read_text(), src.read_text())
            finally:
                os.chdir(cwd)

    def test_different_matrix_with_same_name_is_not_overwritten(self):
        from viral_phylo.matrices import uppercase_matrix_alias
        with tempfile.TemporaryDirectory() as d:
            a = Path(d) / "a" / "Q.3Di.AF"; a.parent.mkdir(); a.write_text("AAA\n")
            b = Path(d) / "b" / "Q.3Di.AF"; b.parent.mkdir(); b.write_text("BBB\n")
            cwd = os.getcwd()
            try:
                os.chdir(d)
                first, second = uppercase_matrix_alias(str(a)), uppercase_matrix_alias(str(b))
                self.assertNotEqual(first, second)
                self.assertEqual(Path(first).read_text(), "AAA\n")
                self.assertEqual(Path(second).read_text(), "BBB\n")
            finally:
                os.chdir(cwd)


class TestManifestNeverLeftRunning(unittest.TestCase):
    """A killed or crashed run used to leave the manifest reading 'running'."""

    def _run(self, body):
        with tempfile.TemporaryDirectory() as d:
            script = textwrap.dedent(f"""
                import os, signal, sys
                sys.path.insert(0, {str(REPO)!r})
                from viral_phylo.manifest import RunManifest
                m = RunManifest({d!r}, quiet=True)
                m.install_exit_handlers()
            """) + textwrap.dedent(body)
            subprocess.run([sys.executable, "-c", script], capture_output=True, timeout=120)
            return json.load(open(os.path.join(d, "run_manifest.json")))

    def test_uncaught_exception_marks_failed(self):
        data = self._run("raise RuntimeError('boom')\n")
        self.assertEqual(data["status"], "failed")
        self.assertEqual(data["error"]["type"], "RuntimeError")

    def test_sigterm_marks_interrupted(self):
        data = self._run("os.kill(os.getpid(), signal.SIGTERM)\nimport time; time.sleep(5)\n")
        self.assertEqual(data["status"], "interrupted")

    def test_clean_finish_is_untouched(self):
        data = self._run("m.finish('completed')\n")
        self.assertEqual(data["status"], "completed")


class TestViewerSingleValueSupport(unittest.TestCase):
    """FastTree's single SH-like support value used to be stored and shown as UFboot."""

    def test_template_no_longer_copies_single_support_into_ufboot(self):
        js = (REPO / "viral_phylo/web/template/viewer.js").read_text()
        self.assertIn("function isSingleValueSupport(node)", js)
        self.assertIn("Single-value support", js)
        self.assertNotIn("node.ufboot = sup;", js)


class TestFastTreeProvenance(unittest.TestCase):
    def test_manifest_reports_the_fasttree_binary_that_runs(self):
        from viral_phylo.binaries import find_fasttree_bin, resolve_tool_versions
        record = resolve_tool_versions(["fasttree"])["fasttree"]
        self.assertEqual(record["path"], find_fasttree_bin())


class TestRecursionNamesAreUnique(unittest.TestCase):
    """A split parent's name used to be reused for its first child, so the recursion
    log could not be joined to the final clusters."""

    def test_split_parent_name_is_never_reused(self):
        from viral_phylo import clustering

        taxa = [f"t{i:02d}" for i in range(40)]

        def fake_extract(files, out_dir, foldmason_bin=None, threads=None):
            os.makedirs(out_dir, exist_ok=True)
            aa = os.path.join(out_dir, "unaligned_aa.fa")
            with open(aa, "w") as h:
                for t in taxa:
                    h.write(f">{t}\n{'A' * 100}\n")
            return aa, os.path.join(out_dir, "unaligned_3di.fa")

        def fake_cluster(struct_dir, work, **kwargs):
            members = sorted(os.path.splitext(n)[0] for n in os.listdir(struct_dir))
            if len(members) == 40:
                return {members[0]: members[:20], members[20]: members[20:]} if "recurse_" in work \
                    else {members[0]: members}
            return {members[0]: members}

        def fake_align(struct_dir, out_dir, **kwargs):
            n = len(os.listdir(struct_dir))
            width = 400 if n == 40 else 100          # the 40-taxon cluster is gappy
            for alph in ("3di", "aa"):
                with open(os.path.join(out_dir, f"foldmason.fasta_{alph}.fa"), "w") as h:
                    for f in sorted(os.listdir(struct_dir)):
                        h.write(f">{os.path.splitext(f)[0]}\n{'A' * 100}{'-' * (width - 100)}\n")

        with tempfile.TemporaryDirectory() as d:
            sdir = Path(d) / "structs"; sdir.mkdir()
            for t in taxa:
                (sdir / f"{t}.cif").write_text("data_x\n")
            with mock.patch("viral_phylo.alignment.extract_unaligned_sequences", fake_extract), \
                 mock.patch("viral_phylo.alignment.align_structures", fake_align), \
                 mock.patch.object(clustering, "foldseek_cluster", fake_cluster):
                summary = clustering.partition_structures_by_similarity(
                    str(sdir), str(Path(d) / "out"), mode="structural", min_cluster_size=10,
                    recurse=True, recurse_max_gap=0.55, recurse_max_expansion=2.5)

        names = [c["name"] for c in summary["clusters"]]
        split = [s for s in summary["recursion"]["splits"] if s["outcome"] == "split"]
        self.assertEqual(len(split), 1, "the gappy 40-taxon cluster should split once")
        self.assertEqual(len(names), len(set(names)), "final cluster names must be unique")
        self.assertNotIn(split[0]["cluster"], names, "a split parent's name must not be reused")
        self.assertTrue(all(c["parent"] == split[0]["cluster"] for c in summary["clusters"]),
                        "children must record their parent")


if __name__ == "__main__":
    unittest.main()


class TestThreadCap(unittest.TestCase):
    """--threads caps the whole pipeline, not just IQ-TREE."""

    def setUp(self):
        from viral_phylo.cli import THREAD_ENV_VARS
        self._saved = {v: os.environ.get(v) for v in THREAD_ENV_VARS}
        for v in THREAD_ENV_VARS:
            os.environ.pop(v, None)

    def tearDown(self):
        for v, val in self._saved.items():
            if val is None:
                os.environ.pop(v, None)
            else:
                os.environ[v] = val

    def test_cap_is_exported_to_child_thread_pools(self):
        from viral_phylo.cli import THREAD_ENV_VARS, apply_thread_cap
        self.assertEqual(apply_thread_cap("6"), 6)
        for v in THREAD_ENV_VARS:
            self.assertEqual(os.environ[v], "6")

    def test_a_lower_user_setting_is_kept(self):
        from viral_phylo.cli import apply_thread_cap
        os.environ["OMP_NUM_THREADS"] = "2"
        apply_thread_cap("8")
        self.assertEqual(os.environ["OMP_NUM_THREADS"], "2")
        self.assertEqual(os.environ["MKL_NUM_THREADS"], "8")

    def test_auto_means_no_cap(self):
        from viral_phylo.cli import apply_thread_cap
        self.assertIsNone(apply_thread_cap("AUTO"))
        self.assertNotIn("OMP_NUM_THREADS", os.environ)

    def test_aligners_are_given_the_count(self):
        from viral_phylo.alignment import _foldmason_threads, _mafft_threads
        self.assertEqual(_foldmason_threads("8"), ["--threads", "8"])
        self.assertEqual(_mafft_threads(8), ["--thread", "8"])
        self.assertEqual(_foldmason_threads("AUTO"), [])
        self.assertEqual(_mafft_threads(0), [])

    def test_every_subcommand_that_runs_tools_accepts_threads(self):
        from viral_phylo.cli import build_cli_parser
        parser = build_cli_parser()
        for argv in (["align"], ["tree", "-a", "x.fa"], ["embed", "-i", "x.fa"], ["pipeline"]):
            self.assertEqual(parser.parse_args(argv + ["--threads", "3"]).threads, "3", argv)
