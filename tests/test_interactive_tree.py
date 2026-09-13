import unittest
from pathlib import Path


class TestInteractiveTree(unittest.TestCase):
    def setUp(self):
        self.repo_root = Path(__file__).resolve().parent.parent
        self.html_path = self.repo_root / "interactive_tree.html"

    def test_interactive_tree_exists_and_valid(self):
        self.assertTrue(self.html_path.exists(), "interactive_tree.html must exist at repo root")
        content = self.html_path.read_text(encoding="utf-8")
        
        # Verify file size (standalone with embedded precomputations is > 2 MB)
        self.assertGreater(len(content), 1_000_000)

        # Verify critical script assets and fallbacks
        self.assertIn("results/ca_500_structures.js", content)
        self.assertIn("results/alignments_data.js", content)

        # Verify cohorts embedded
        self.assertIn('"1193":', content)
        self.assertIn('"500":', content)
        self.assertIn('"6":', content)

        # Verify ZIP exporter and reproducibility components
        self.assertIn("function exportSubcladePackage()", content)
        self.assertIn("function createZipArchive(files)", content)
        self.assertIn("REPRODUCE.sh", content)
        self.assertIn("REPRODUCIBILITY.md", content)

    def test_example_pdbs_exist(self):
        """Verify that the 6 benchmark example PDB structures exist in results/glycoprotein_workflow/structures."""
        structures_dir = self.repo_root / "results/glycoprotein_workflow/structures"
        self.assertTrue(structures_dir.exists(), "results/glycoprotein_workflow/structures must exist")
        pdbs = list(structures_dir.glob("*.pdb"))
        self.assertEqual(len(pdbs), 6, f"Expected exactly 6 example PDB structures, found {len(pdbs)}")
        for pdb in pdbs:
            self.assertGreater(pdb.stat().st_size, 10_000, f"{pdb.name} is too small to be a valid PDB")


if __name__ == "__main__":
    unittest.main()
