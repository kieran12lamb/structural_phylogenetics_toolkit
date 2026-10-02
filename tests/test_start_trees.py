"""FastTree starting trees must never cause IQ-TREE to fail or lose a tree."""

import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from Bio import Phylo

from viral_phylo.start_trees import resolve_to_binary, validate_start_tree

REPO = Path(__file__).resolve().parent.parent


def _tree(newick):
    return Phylo.read(io.StringIO(newick), "newick")


def _binary(tree):
    root = tree.root
    return len(root.clades) == 3 and all(len(c.clades) == 2 for c in tree.find_clades() if c.clades and c is not root)


class TestRepair(unittest.TestCase):
    def test_polytomy_is_resolved_and_leaves_preserved(self):
        """FastTree joins identical sequences in a multifurcation; IQ-TREE aborts on it."""
        t = _tree("((a:0.1,b:0.0,c:0.0,d:0.0)0.9:0.2,e:0.1,(f:0.1,g:0.1)0.8:0.1);")
        n = resolve_to_binary(t, seed=1)
        self.assertEqual(n, 1)
        self.assertTrue(_binary(t))
        self.assertEqual(sorted(x.name for x in t.get_terminals()), list("abcdefg"))

    def test_resolution_is_reproducible(self):
        a, b = _tree("((a,b,c,d,e),f,g);"), _tree("((a,b,c,d,e),f,g);")
        resolve_to_binary(a, seed=7); resolve_to_binary(b, seed=7)
        out = []
        for t in (a, b):
            buf = io.StringIO(); Phylo.write(t, buf, "newick"); out.append(buf.getvalue())
        self.assertEqual(out[0], out[1])

    def test_rooted_binary_input_becomes_trifurcating(self):
        t = _tree("((a,b),(c,d));")
        resolve_to_binary(t)
        self.assertTrue(_binary(t))

    def test_labels_removed_and_lengths_positive(self):
        t = _tree("((a:0,b:-0.1)0.95:0.2,c,d);")
        resolve_to_binary(t)
        self.assertTrue(all(c.confidence is None for c in t.find_clades()))
        self.assertTrue(all(c.branch_length > 0 for c in t.find_clades() if c is not t.root))


class TestValidation(unittest.TestCase):
    def test_accepts_matching_binary_tree(self):
        t = _tree("((a,b),c,d);"); resolve_to_binary(t)
        self.assertIsNone(validate_start_tree(t, list("abcd")))

    def test_rejects_leaf_mismatch(self):
        t = _tree("((a,b),c,d);"); resolve_to_binary(t)
        self.assertIn("leaf set", validate_start_tree(t, list("abce")))

    def test_rejects_polytomy(self):
        self.assertIn("bifurcating", validate_start_tree(_tree("((a,b,c),d,e);"), list("abcde")))


# Six taxa, three of them identical - the case that made IQ-TREE reject or abort.
ALN = (">t1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKL\n>t2\nACDEFGHIKLMNPQRSTVWYACDEFGHIKL\n"
       ">t3\nACDEFGHIKLMNPQRSTVWYACDEFGHIKL\n>t4\nACDEFGHIKAMNPQRSAVWYACDEFGHIKL\n"
       ">t5\nACDEWGHIKLMNPQRSTVWYACDEFGAIKL\n>t6\nACDEFGHIKLMNPQRSTVWFACDMFGHIKL\n")


def _have(name):
    from viral_phylo.binaries import find_fasttree_bin, tool_is_available
    return find_fasttree_bin() is not None if name == "fasttree" else tool_is_available(name)


@unittest.skipUnless(_have("iqtree") and _have("fasttree"), "IQ-TREE and FastTree required")
class TestGuidedIqtree(unittest.TestCase):
    def _run(self, d):
        from viral_phylo.tree import build_tree
        aln = Path(d) / "foldmason.fasta_3di.fa"; aln.write_text(ALN)
        cwd = os.getcwd()
        try:
            os.chdir(REPO)
            build_tree(alignment_file=str(aln), method="iqtree", matrix="alphafold", rate_het="+G4",
                       bootstrap=1000, alrt=0, threads="2", output_dir=d, prefix="t", tree_type="3di",
                       guide_trees=True, guide_min_taxa=1)
        finally:
            os.chdir(cwd)
        rec = json.load(open(Path(d) / "t_alphafold.start_tree.json"))
        tree = Phylo.read(str(Path(d) / "t_alphafold.treefile"), "newick")
        return rec, sorted(x.name for x in tree.get_terminals())

    def test_guided_run_succeeds_with_identical_sequences(self):
        with tempfile.TemporaryDirectory() as d:
            rec, leaves = self._run(d)
        self.assertEqual(rec["outcome"], "guided")
        self.assertIn("-keep-ident", rec["iqtree_command"])
        self.assertEqual(leaves, [f"t{i}" for i in range(1, 7)], "no taxon may be lost")

    def test_failed_guided_run_falls_back_and_still_builds_the_tree(self):
        def broken_start_tree(alignment, out_path, alphabet, **kw):
            Path(out_path).write_text("((t1,t2),t3,not_in_alignment);\n")   # IQ-TREE will reject this
            return {"used": True, "path": out_path, "reason": None, "polytomies_resolved": 0,
                    "alphabet": alphabet, "seed": 1, "command": "stub"}
        with tempfile.TemporaryDirectory() as d, \
                mock.patch("viral_phylo.start_trees.make_start_tree", broken_start_tree):
            rec, leaves = self._run(d)
        self.assertEqual(rec["outcome"], "fell_back")
        self.assertNotIn("-t ", rec["iqtree_command"])
        self.assertEqual(leaves, [f"t{i}" for i in range(1, 7)])


if __name__ == "__main__":
    unittest.main()


class TestGuideTreesDefault(unittest.TestCase):
    def test_on_by_default_with_opt_out(self):
        from viral_phylo.cli import build_cli_parser
        parser = build_cli_parser()
        self.assertTrue(parser.parse_args(["pipeline"]).guide_trees)
        self.assertFalse(parser.parse_args(["pipeline", "--no-guide-trees"]).guide_trees)
        self.assertEqual(parser.parse_args(["pipeline"]).guide_tree_min_taxa, 80)


class TestPreflightWithoutFastTree(unittest.TestCase):
    """Starting trees are an optimisation, so a missing FastTree must not block a run
    that otherwise uses only IQ-TREE - but must still block one that needs FastTree
    for the whole-set tree."""

    def _preflight(self, master_method):
        from viral_phylo.cli import build_cli_parser, pipeline_preflight
        from viral_phylo.manifest import RunManifest
        with tempfile.TemporaryDirectory() as d:
            args = build_cli_parser().parse_args(
                ["pipeline", "--input-folder", str(REPO / "tests/fixtures/structures"),
                 "--master-method", master_method, "--method", "iqtree"])
            manifest = RunManifest(d, quiet=True)
            with mock.patch("viral_phylo.binaries.find_fasttree_bin", return_value=None):
                try:
                    pipeline_preflight(args, manifest)
                    return "ran", manifest.data["preflight"]
                except SystemExit:
                    return "blocked", manifest.data["preflight"]

    def test_iqtree_only_run_is_warned_not_blocked(self):
        outcome, pre = self._preflight("iqtree")
        self.assertEqual(outcome, "ran")
        self.assertTrue(any("starting trees are disabled" in w for w in pre["warnings"]))

    def test_fasttree_master_tree_still_blocks(self):
        outcome, pre = self._preflight("fasttree")
        self.assertEqual(outcome, "blocked")
        self.assertTrue(any("FastTree" in b for b in pre["blockers"]))
