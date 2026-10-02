import json
import re
import tempfile
import unittest
from pathlib import Path

from viral_phylo.cli import build_cli_parser
from viral_phylo.web.report import build_tree_report, newick_leaf_names


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


if __name__ == "__main__":
    unittest.main()
