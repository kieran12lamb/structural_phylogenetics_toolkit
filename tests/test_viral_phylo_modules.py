"""Unit tests for the modular viral_phylo package architecture."""

import unittest
from pathlib import Path

import viral_phylo
from viral_phylo.binaries import (
    find_iqtree_bin,
    find_foldmason_bin,
    find_mafft_bin,
    find_torch_python,
)
from viral_phylo.matrices import (
    ensure_3di_matrix,
    ensure_matrix_file,
    EDMOND_MATRICES,
)
from viral_phylo.alignment import (
    parse_alignment_fasta,
    write_alignment_fasta,
    strip_all_gap_columns_dict,
    compute_alignment_coverage,
    filter_alignment_by_coverage,
)
from viral_phylo.metadata import (
    classify_columns,
    choose_default_color_column,
    EXPANDED_PALETTE,
    KNOWN_VALUE_COLORS,
)
from viral_phylo.web.builder import assemble_html
from viral_phylo.web.alignments import parse_fasta, compute_cov_dict


class TestViralPhyloPackage(unittest.TestCase):
    def setUp(self):
        self.repo_root = Path(__file__).resolve().parent.parent

    def test_package_exports_and_version(self):
        """Verify package version and top-level symbols."""
        self.assertEqual(viral_phylo.__version__, "1.0.0")
        self.assertTrue(len(viral_phylo.__all__) >= 20)
        self.assertIn("find_iqtree_bin", viral_phylo.__all__)
        self.assertIn("build_tree", viral_phylo.__all__)
        self.assertIn("compute_umap_projection", viral_phylo.__all__)

    def test_binaries_discovery(self):
        """Verify binary discovery functions return non-empty strings."""
        iq = find_iqtree_bin()
        fm = find_foldmason_bin()
        ma = find_mafft_bin()
        py = find_torch_python()

        self.assertIsInstance(iq, str)
        self.assertGreater(len(iq), 0)
        self.assertIsInstance(fm, str)
        self.assertGreater(len(fm), 0)
        self.assertIsInstance(ma, str)
        self.assertGreater(len(ma), 0)
        self.assertIsInstance(py, str)
        self.assertGreater(len(py), 0)

    def test_matrices_metadata(self):
        """Verify Edmond matrix metadata definitions."""
        self.assertIn("alphafold", EDMOND_MATRICES)
        self.assertIn("esmfold", EDMOND_MATRICES)
        self.assertEqual(EDMOND_MATRICES["alphafold"]["filename"], "Q.3Di.AF")
        self.assertEqual(EDMOND_MATRICES["esmfold"]["filename"], "Q.3Di.LLM")

    def test_alignment_core_routines(self):
        """Verify alignment parsing, gap stripping, and coverage computations."""
        mock_aln = {
            "taxa1": "ACDEF-GHIK",
            "taxa2": "A-DEF-GHIK",
            "taxa3": "---EF-GHIK",
        }
        stripped = strip_all_gap_columns_dict(mock_aln)
        self.assertEqual(len(stripped["taxa1"]), 9)

        cov = compute_alignment_coverage(mock_aln)
        self.assertIn("taxa", cov)
        self.assertEqual(cov["aln_length"], 10)
        self.assertEqual(cov["taxa"]["taxa1"]["sequence_coverage"], 0.9)
        self.assertEqual(cov["taxa"]["taxa2"]["sequence_coverage"], 0.8)
        self.assertEqual(cov["taxa"]["taxa3"]["sequence_coverage"], 0.6)

    def test_web_template_assembly(self):
        """Verify that assemble_html merges index, styles, scripts, and JSON cleanly."""
        dummy_json = '{"test": true}'
        dummy_version = "v1.0.0-test"
        dummy_options = '<option value="test">Test</option>'
        html = assemble_html(dummy_json, dummy_version, dummy_options)

        self.assertIn("<!DOCTYPE html>", html)
        self.assertIn("<style>", html)
        self.assertIn("v1.0.0-test", html)
        self.assertIn('const DATASETS = {"test": true};', html)
        self.assertIn("function renderTree()", html)
        self.assertIn("function renderUmapScatter(", html)


if __name__ == "__main__":
    unittest.main()
