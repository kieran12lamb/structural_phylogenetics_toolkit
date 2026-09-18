import unittest
from pathlib import Path
from scripts.viral_phylogenetics import build_cli_parser, ensure_matrix_file

class TestCLI(unittest.TestCase):
    def setUp(self):
        self.parser = build_cli_parser()
        self.repo_root = Path(__file__).resolve().parent.parent

    def test_cli_parser_fetch(self):
        args = self.parser.parse_args(["fetch", "-q", "glycoprotein", "-m", "10", "-o", "structures"])
        self.assertEqual(args.subcommand, "fetch")
        self.assertEqual(args.qualifier, "glycoprotein")
        self.assertEqual(args.max_sequences, 10)
        self.assertEqual(args.output_dir, "structures")

    def test_cli_parser_align(self):
        args = self.parser.parse_args(["align", "-i", "structures", "-o", "alignments", "--aligner", "foldmason"])
        self.assertEqual(args.subcommand, "align")
        self.assertEqual(args.folder, "structures")
        self.assertEqual(args.aligner, "foldmason")

    def test_cli_parser_tree(self):
        args = self.parser.parse_args([
            "tree",
            "-a", "alignments/3di.fa",
            "--alignment-aa", "alignments/aa.fa",
            "--tree-type", "both",
            "--matrix", "alphafold",
            "--rate-heterogeneity", "auto"
        ])
        self.assertEqual(args.subcommand, "tree")
        self.assertEqual(args.tree_type, "both")
        self.assertEqual(args.matrix, "alphafold")

    def test_cli_parser_embed(self):
        args = self.parser.parse_args([
            "embed",
            "--fasta", "alignments/aa.fa",
            "--model", "esm2",
            "--metric", "cosine"
        ])
        self.assertEqual(args.subcommand, "embed")
        self.assertEqual(args.model, "esm2")
        self.assertEqual(args.metric, "cosine")

    def test_cli_parser_pipeline(self):
        args = self.parser.parse_args([
            "pipeline",
            "--input-folder", "my_structures",
            "--tree-type", "both",
            "--matrix", "both"
        ])
        self.assertEqual(args.subcommand, "pipeline")
        self.assertEqual(args.input_folder, "my_structures")
        self.assertEqual(args.matrix, "both")

    def test_cli_parser_fetch_alphafold(self):
        args = self.parser.parse_args([
            "fetch",
            "--source", "alphafold",
            "-u", "P00520,P04637",
            "--format", "pdb",
            "--download-pae",
            "-o", "afdb_dir"
        ])
        self.assertEqual(args.subcommand, "fetch")
        self.assertEqual(args.source, "alphafold")
        self.assertEqual(args.uniprot, "P00520,P04637")
        self.assertEqual(args.format, "pdb")
        self.assertTrue(args.download_pae)
        self.assertEqual(args.output_dir, "afdb_dir")

    def test_cli_parser_pipeline_alphafold(self):
        args = self.parser.parse_args([
            "pipeline",
            "--source", "alphafold",
            "-u", "P00520,P04637",
            "--tree-type", "3di",
            "--download-pae",
            "-o", "afdb_pipeline_results"
        ])
        self.assertEqual(args.subcommand, "pipeline")
        self.assertEqual(args.source, "alphafold")
        self.assertEqual(args.uniprot, "P00520,P04637")
        self.assertEqual(args.tree_type, "3di")
        self.assertTrue(args.download_pae)
        self.assertEqual(args.output_dir, "afdb_pipeline_results")


    def test_matrix_discovery(self):
        matrices_dir = str(self.repo_root / "matrices")
        af_path = ensure_matrix_file("alphafold", matrices_dir=matrices_dir)
        self.assertIsNotNone(af_path)
        self.assertIn("Q.3Di.AF", af_path)

    def test_find_torch_python_resolution(self):
        import os
        from scripts.viral_phylogenetics import find_torch_python
        from scripts.embed_and_cluster import find_torch_python as find_torch_embed
        
        py1 = find_torch_python()
        self.assertIsNotNone(py1)
        self.assertTrue(os.path.isfile(py1))
        self.assertTrue(os.access(py1, os.X_OK))

        py2 = find_torch_embed()
        self.assertIsNotNone(py2)
        self.assertTrue(os.path.isfile(py2))
        self.assertTrue(os.access(py2, os.X_OK))


    def test_cli_parser_alignment_coverage_and_multi_alignment(self):
        """Verify --min-coverage and --multi-alignment in align and pipeline subcommands."""
        # Align subcommand with coverage options
        args_align = self.parser.parse_args([
            "align", "-i", "structures", "-o", "alignments",
            "--min-coverage", "0.70", "--multi-alignment", "--filter-coverage"
        ])
        self.assertEqual(args_align.min_coverage, 0.70)
        self.assertTrue(args_align.multi_alignment)
        self.assertTrue(args_align.filter_coverage)

        # Pipeline subcommand with coverage options
        args_pipe = self.parser.parse_args([
            "pipeline", "--input-folder", "structures",
            "--min-coverage", "0.75", "--multi-alignment"
        ])
        self.assertEqual(args_pipe.min_coverage, 0.75)
        self.assertTrue(args_pipe.multi_alignment)


if __name__ == "__main__":
    unittest.main()
