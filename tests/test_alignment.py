import unittest
from pathlib import Path

VALID_3DI_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY-")

def strip_all_gap_columns(alignment_dict):
    """Utility function matching the browser gap stripping engine."""
    if not alignment_dict:
        return {}
    seq_len = len(next(iter(alignment_dict.values())))
    keep_indices = []
    for c in range(seq_len):
        has_residue = any(seq[c] != "-" for seq in alignment_dict.values())
        if has_residue:
            keep_indices.append(c)
    
    stripped = {}
    for tid, seq in alignment_dict.items():
        stripped[tid] = "".join(seq[i] for i in keep_indices)
    return stripped


class TestAlignment(unittest.TestCase):
    def setUp(self):
        self.sample_3di = {
            "Taxon_A": "DPCDVPDDPDPCDVPD",
            "Taxon_B": "DPCDVPDDPDP--VPD",
            "Taxon_C": "DP-DVPDDPDP--VPD",
            "Taxon_D": "DP-DVPDDPDP--VPD"
        }
        self.sample_aa = {
            "Taxon_A": "MKVLLLLALLSTLTLL",
            "Taxon_B": "MKVLLLLALLST-TLL",
            "Taxon_C": "MK-LLLLALLST-TLL",
            "Taxon_D": "MK-LLLLALLST-TLL"
        }
        self.repo_root = Path(__file__).resolve().parent.parent

    def test_3di_alphabet_validity(self):
        """Verify that 3Di records only contain valid structural states or gaps."""
        for tid, seq in self.sample_3di.items():
            chars = set(seq.upper())
            invalid = chars - VALID_3DI_ALPHABET
            self.assertEqual(len(invalid), 0, f"Found invalid 3Di characters {invalid} in {tid}")

    def test_strip_all_gap_columns(self):
        """Test that columns with all gaps are stripped properly."""
        sample_gapped = {
            "Taxon_1": "M-K--L-V",
            "Taxon_2": "A-K--I-V",
            "Taxon_3": "V-R--L-I"
        }
        stripped = strip_all_gap_columns(sample_gapped)
        self.assertEqual(len(stripped), 3)
        for tid, seq in stripped.items():
            self.assertEqual(len(seq), 4)
        self.assertEqual(stripped["Taxon_1"], "MKLV")
        self.assertEqual(stripped["Taxon_2"], "AKIV")
        self.assertEqual(stripped["Taxon_3"], "VRLI")

    def test_mat3di_matrix_exists(self):
        """Ensure mat3di.out exists in matrices/ directory and has content."""
        matrix_path = self.repo_root / "matrices/mat3di.out"
        self.assertTrue(matrix_path.exists(), "matrices/mat3di.out must exist")
        self.assertGreater(matrix_path.stat().st_size, 0, "matrices/mat3di.out must not be empty")

    def test_iqtree_3di_matrices_exist(self):
        """Ensure Q.3Di.AF and Q.3Di.LLM matrices exist in matrices/."""
        af_matrix = self.repo_root / "matrices/Q.3Di.AF"
        llm_matrix = self.repo_root / "matrices/Q.3Di.LLM"
        self.assertTrue(af_matrix.exists(), "matrices/Q.3Di.AF must exist")
        self.assertTrue(llm_matrix.exists(), "matrices/Q.3Di.LLM must exist")
        self.assertGreater(af_matrix.stat().st_size, 100)
        self.assertGreater(llm_matrix.stat().st_size, 100)


    def test_compute_alignment_coverage(self):
        """Verify alignment coverage calculation across sequences and columns."""
        from scripts.viral_phylogenetics import compute_alignment_coverage

        test_aln = {
            "Seq1": "AAAA----",  # 4 non-gaps / 8 = 0.50
            "Seq2": "AAAAAA--",  # 6 non-gaps / 8 = 0.75
            "Seq3": "AAAAAAAA"   # 8 non-gaps / 8 = 1.00
        }
        res = compute_alignment_coverage(test_aln)
        self.assertEqual(res["aln_length"], 8)
        self.assertEqual(res["taxa"]["Seq1"]["sequence_coverage"], 0.50)
        self.assertEqual(res["taxa"]["Seq2"]["sequence_coverage"], 0.75)
        self.assertEqual(res["taxa"]["Seq3"]["sequence_coverage"], 1.00)

        # Check column occupancies
        self.assertEqual(res["column_occupancy"][0], 1.0)
        self.assertEqual(res["column_occupancy"][4], round(2/3, 4))
        self.assertEqual(res["column_occupancy"][6], round(1/3, 4))

    def test_filter_alignment_by_coverage(self):
        """Verify filtering sequences below the coverage threshold."""
        from scripts.viral_phylogenetics import filter_alignment_by_coverage, write_alignment_fasta
        import tempfile

        test_aln = {
            "Core1": "AAAAAAAA",  # 1.00
            "Core2": "AAAAAA--",  # 0.75
            "Frag1": "AA------",  # 0.25
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            f_3di = Path(tmpdir) / "test_3di.fa"
            write_alignment_fasta(test_aln, str(f_3di))

            # Filter at 70% threshold
            f_res = filter_alignment_by_coverage(str(f_3di), min_coverage=0.70, output_dir=tmpdir)
            self.assertIn("Core1", f_res["passed_taxa"])
            self.assertIn("Core2", f_res["passed_taxa"])
            self.assertIn("Frag1", f_res["removed_taxa"])
            self.assertEqual(len(f_res["passed_taxa"]), 2)
            self.assertEqual(len(f_res["removed_taxa"]), 1)

            # Check that output file exists and has stripped gap columns
            self.assertTrue(Path(f_res["out_3di"]).exists())


if __name__ == "__main__":
    unittest.main()
