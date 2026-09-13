import unittest
from scripts.metadata_handler import classify_columns, choose_default_color_column


class TestMetadata(unittest.TestCase):
    def test_classify_columns(self):
        records = [
            {"id": f"tax{i}", "Family": "Rhabdoviridae" if i % 2 == 0 else "Flaviviridae", "pLDDT": f"{60.0 + i * 3.5}", "seq_len": f"{400 + i * 15}"}
            for i in range(12)
        ]
        
        col_specs = classify_columns(records, id_col="id")
        spec_map = {s["key"]: s for s in col_specs}
        
        # Family should be classified as categorical
        self.assertIn("Family", spec_map)
        self.assertEqual(spec_map["Family"]["type"], "categorical")
        self.assertTrue("colors" in spec_map["Family"] or "palette" in spec_map["Family"])
        color_dict = spec_map["Family"].get("colors", {})
        self.assertGreaterEqual(len(color_dict), 2)

        # pLDDT should be classified as continuous
        self.assertIn("pLDDT", spec_map)
        self.assertEqual(spec_map["pLDDT"]["type"], "continuous")
        self.assertLessEqual(spec_map["pLDDT"]["min"], 65.0)
        self.assertGreaterEqual(spec_map["pLDDT"]["max"], 90.0)

        # seq_len should be classified as continuous
        self.assertIn("seq_len", spec_map)
        self.assertEqual(spec_map["seq_len"]["type"], "continuous")
        self.assertLessEqual(spec_map["seq_len"]["min"], 455.0)
        self.assertGreaterEqual(spec_map["seq_len"]["max"], 505.0)

    def test_choose_default_color_column(self):
        col_specs = [
            {"key": "notes", "type": "categorical"},
            {"key": "Family", "type": "categorical"},
            {"key": "pLDDT", "type": "continuous"}
        ]
        best_col = choose_default_color_column(col_specs)
        self.assertIn(best_col, ("Family", "pLDDT"))


if __name__ == "__main__":
    unittest.main()
