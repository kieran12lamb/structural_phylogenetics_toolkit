import json
import re
import tempfile
import unittest
import zipfile
from pathlib import Path

from viral_phylo.cli import build_cli_parser
from viral_phylo.web.report import build_tree_report, newick_leaf_names, package_run_report


class TestTreeReport(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.tree = self.tmp / "example.treefile"
        self.tree.write_text("((A_1:0.1,B_2:0.2)95/100:0.05,(C_3:0.3,D_4:0.1):0.02);\n")

    def _embedded_dataset(self, html):
        match = re.search(r"const DATASETS = (.*?);\n", html)
        self.assertIsNotNone(match)
        return json.loads(match.group(1).replace("<\\/", "</"))["custom"]

    def test_leaf_names_newick_and_nexus(self):
        self.assertEqual(newick_leaf_names(self.tree.read_text()), ["A_1", "B_2", "C_3", "D_4"])
        nexus = "#NEXUS\nbegin trees;\n translate 1 'Taxon A', 2 'it''s';\n tree t = [&R] (1:0.1,2:0.2);\nend;\n"
        self.assertEqual(newick_leaf_names(nexus), ["Taxon A", "it's"])

    def test_report_without_metadata(self):
        out = build_tree_report(str(self.tree))
        self.assertEqual(out, self.tmp / "example_report.html")
        html = out.read_text()
        self.assertIn('let currentScale = "custom";', html)
        ds = self._embedded_dataset(html)
        self.assertTrue(ds["custom_tree"])
        self.assertIsNone(ds["raw_metadata"])
        self.assertIn("A_1:0.1", ds["tree_text"])

    def test_report_with_metadata(self):
        meta = self.tmp / "meta.csv"
        meta.write_text("sample,host\nA_1,Bat\nB_2,Human\nC_3,Bat\n")
        out = build_tree_report(str(self.tree), metadata_path=str(meta), output_path=str(self.tmp / "r.html"), title="My tree")
        ds = self._embedded_dataset(out.read_text())
        self.assertEqual(ds["title"], "📂 My tree")
        self.assertEqual(ds["metadata_filename"], "meta.csv")
        self.assertEqual(set(ds["raw_metadata"]["taxa"]), {"A_1", "B_2", "C_3"})
        self.assertEqual(ds["raw_metadata"]["columns"][0]["key"], "host")

    def test_script_breakout_is_escaped(self):
        self.tree.write_text("(A:1,'</script><b>':2);")
        html = build_tree_report(str(self.tree)).read_text()
        self.assertNotIn("'</script><b>'", html)

    def test_rejects_non_tree_file(self):
        bad = self.tmp / "bad.nwk"
        bad.write_text("not a tree")
        with self.assertRaises(ValueError):
            build_tree_report(str(bad))

    def test_cli_parser_report(self):
        args = build_cli_parser().parse_args(["report", "-t", "t.nwk", "-meta", "m.csv", "-o", "out.html"])
        self.assertEqual((args.subcommand, args.tree, args.metadata, args.output), ("report", "t.nwk", "m.csv", "out.html"))
        args_run = build_cli_parser().parse_args(["report", "--run-dir", "results/my_run", "-o", "my_run.zip"])
        self.assertEqual((args_run.subcommand, args_run.run_dir, args_run.output), ("report", "results/my_run", "my_run.zip"))

    def test_cli_parser_pipeline_zip_flags(self):
        p = build_cli_parser()
        a_def = p.parse_args(["pipeline"])
        self.assertTrue(a_def.zip)
        self.assertIsNone(a_def.zip_path)
        a_no = p.parse_args(["pipeline", "--no-zip"])
        self.assertFalse(a_no.zip)
        a_custom = p.parse_args(["pipeline", "--zip-path", "/tmp/bundle.zip"])
        self.assertEqual(a_custom.zip_path, "/tmp/bundle.zip")

    def test_package_run_report_creates_self_contained_zip(self):
        run_dir = self.tmp / "demo_run"
        aln_dir = run_dir / "alignment"
        phy_dir = run_dir / "phylogeny"
        c1_aln = aln_dir / "cluster_1_e1.0"
        c1_phy = phy_dir / "cluster_1_e1.0"
        scratch_dir = aln_dir / "tmpFolder"
        for d in (c1_aln, c1_phy, scratch_dir):
            d.mkdir(parents=True, exist_ok=True)

        # Write alignments, trees, metadata, manifest, ca_structures, and scratch files
        (aln_dir / "foldmason.fasta_3di.fa").write_text(">A_1\nACDEF\n>B_2\nACDEF\n")
        (aln_dir / "foldmason.fasta_aa.fa").write_text(">A_1\nMKVLV\n>B_2\nMKVLV\n")
        (c1_aln / "foldmason.fasta_3di.fa").write_text(">A_1\nACDEF\n>B_2\nACDEF\n")
        (c1_aln / "foldmason.fasta_aa.fa").write_text(">A_1\nMKVLV\n>B_2\nMKVLV\n")
        (c1_aln / "foldmason.fasta.nw").write_text("(A_1:0.1,B_2:0.2);\n")
        (scratch_dir / "temp_db.bin").write_text("scratch")
        (aln_dir / "sDB.index").write_text("scratch")

        (phy_dir / "demo_run_both.treefile").write_text("(A_1:0.1,B_2:0.2);\n")
        (phy_dir / "demo_run_aa.treefile").write_text("(A_1:0.15,B_2:0.25);\n")
        (c1_phy / "cluster_1_e1.0_both.treefile").write_text("(A_1:0.1,B_2:0.2);\n")
        (run_dir / "taxa_metadata.json").write_text(json.dumps({
            "taxa": {"A_1": {"taxon_id": "A_1", "family": "Paramyxo"}, "B_2": {"taxon_id": "B_2", "family": "Paramyxo"}},
            "columns": [{"key": "family", "label": "Family", "type": "categorical", "values": ["Paramyxo"], "colors": {"Paramyxo": "#38bdf8"}}],
            "default_color_col": "family"
        }))
        (run_dir / "run_manifest.json").write_text(json.dumps({"status": "completed", "stages": [], "warnings": []}))
        (run_dir / "run_manifest.txt").write_text("status: completed\n")
        (run_dir / "ca_structures.js").write_text("window.CA_STRUCTURES_BY_DATASET = {};\n")

        zip_out = package_run_report(run_dir)
        self.assertEqual(zip_out, run_dir / "demo_run_report.zip")
        self.assertTrue(zip_out.is_file())
        self.assertTrue((run_dir / "interactive_tree.html").is_file())
        self.assertTrue((run_dir / "alignments_data.js").is_file())

        with zipfile.ZipFile(zip_out, "r") as zf:
            names = set(zf.namelist())
            self.assertIn("demo_run_report/interactive_tree.html", names)
            self.assertIn("demo_run_report/index.html", names)
            self.assertIn("demo_run_report/ca_structures.js", names)
            self.assertIn("demo_run_report/alignments_data.js", names)
            self.assertIn("demo_run_report/taxa_metadata.json", names)
            self.assertIn("demo_run_report/run_manifest.json", names)
            self.assertIn("demo_run_report/phylogeny/demo_run_both.treefile", names)
            self.assertIn("demo_run_report/phylogeny/cluster_1_e1.0/cluster_1_e1.0_both.treefile", names)
            self.assertIn("demo_run_report/alignment/foldmason.fasta_3di.fa", names)
            self.assertIn("demo_run_report/alignment/cluster_1_e1.0/foldmason.fasta.nw", names)
            # Ensure scratch files are excluded
            self.assertFalse(any("tmpFolder" in n for n in names))
            self.assertFalse(any("sDB" in n for n in names))


if __name__ == "__main__":
    unittest.main()
