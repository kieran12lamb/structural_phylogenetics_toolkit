"""FastTree -trans support for the 3Di substitution matrices.

The conversion conventions were established against FastTree's own input checks
and confirmed by likelihood: on a fixed 32-taxon 3Di tree, FastTree with the
converted Q.3Di.AF matched IQ-TREE -m Q.3Di.AF to within 0.003 log-lik units.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from viral_phylo.matrices import PAML_AA_ORDER, paml_to_fasttree_trans, read_paml_matrix

REPO = Path(__file__).resolve().parent.parent
MATRIX = REPO / "matrices" / "Q.3Di.AF"


def _read_trans(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        parts = line.split("\t")
        rows[parts[0]] = [float(x) for x in parts[1:]]
    return header, rows


class TestConversion(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = paml_to_fasttree_trans(str(MATRIX), os.path.join(self.tmp.name, "q.fasttree"))
        self.header, self.rows = _read_trans(self.out)
        _, self.pi = read_paml_matrix(str(MATRIX))

    def tearDown(self):
        self.tmp.cleanup()

    def test_header_has_no_leading_cell(self):
        self.assertEqual(self.header, list(PAML_AA_ORDER) + ["*"])

    def test_one_row_per_letter_in_order(self):
        self.assertEqual(list(self.rows), list(PAML_AA_ORDER))

    def test_columns_sum_to_zero(self):
        """FastTree's convention: row i, column j is the rate j -> i."""
        for j in range(20):
            self.assertAlmostEqual(sum(self.rows[a][j] for a in PAML_AA_ORDER), 0.0, places=10)

    def test_diagonal_negative_and_frequencies_carried_over(self):
        for i, a in enumerate(PAML_AA_ORDER):
            self.assertLess(self.rows[a][i], 0.0)
            self.assertAlmostEqual(self.rows[a][20], self.pi[i], places=12)

    def test_reversible_and_normalised(self):
        Q = [[self.rows[PAML_AA_ORDER[j]][i] for j in range(20)] for i in range(20)]  # un-transpose
        for i in range(20):
            for j in range(20):
                self.assertAlmostEqual(self.pi[i] * Q[i][j], self.pi[j] * Q[j][i], places=12)
        self.assertAlmostEqual(-sum(self.pi[i] * Q[i][i] for i in range(20)), 1.0, places=10)

    def test_rejects_a_truncated_matrix(self):
        bad = Path(self.tmp.name) / "bad.paml"
        bad.write_text("0.1 0.2 0.3\n")
        with self.assertRaises(ValueError):
            paml_to_fasttree_trans(str(bad), os.path.join(self.tmp.name, "x"))


FASTTREE = shutil.which("FastTree") or shutil.which("VeryFastTree")
# The environment's bin/ is not on PATH unless conda is activated; look beside Python too.
FASTTREE_BESIDE_PYTHON = os.path.join(os.path.dirname(os.path.realpath(sys.executable)), "FastTree")


@unittest.skipUnless(FASTTREE or os.path.isfile(FASTTREE_BESIDE_PYTHON), "FastTree not installed")
class TestFastTreeUsesTheMatrix(unittest.TestCase):
    """The matrix must actually change the model, not merely be accepted."""

    ALN = ">a\nACDEFGHIKLMNPQRSTVWY\n>b\nACDEFGHIKLMNPQRSTVWA\n>c\nACDEAGHIKLMNPQRSAVWY\n" \
          ">d\nACDEFGHIKLWNPQRSTVWY\n>e\nACDEFGHAKLMNPQRSTVWY\n"

    def test_build_tree_records_the_3di_matrix(self):
        from viral_phylo.tree import build_tree
        with tempfile.TemporaryDirectory() as d:
            aln = Path(d) / "foldmason.fasta_3di.fa"; aln.write_text(self.ALN)
            cwd = os.getcwd()
            try:
                os.chdir(REPO)                      # matrices/ resolves from the repo root
                build_tree(alignment_file=str(aln), method="fasttree", matrix="alphafold", rate_het="auto",
                           criterion="BIC", bootstrap=0, alrt=0, threads="2", output_dir=d, prefix="t",
                           tree_type="3di")
            finally:
                os.chdir(cwd)
            info = json.load(open(Path(d) / "t_3di_fasttree.treefile.info.json"))
            self.assertIn("-trans", info["command"])
            self.assertIn("Q.3Di.AF", info["model"])
            self.assertIsNone(info["note"])

    def test_trans_changes_the_likelihood_versus_lg(self):
        binary = FASTTREE or FASTTREE_BESIDE_PYTHON
        with tempfile.TemporaryDirectory() as d:
            aln = Path(d) / "a.fa"; aln.write_text(self.ALN)
            trans = paml_to_fasttree_trans(str(MATRIX), os.path.join(d, "q.fasttree"))

            def loglk(*model):
                err = subprocess.run([binary, *model, "-nocat", str(aln)], capture_output=True, text=True).stderr
                return [ln for ln in err.splitlines() if "LogLk" in ln][-1]
            self.assertNotEqual(loglk("-lg"), loglk("-trans", trans))


if __name__ == "__main__":
    unittest.main()
