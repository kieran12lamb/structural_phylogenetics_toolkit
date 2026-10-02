import tempfile
import unittest
from pathlib import Path

from viral_phylo.structures import (
    extract_ca_dict,
    parse_structure_ca,
    write_ca_bundle,
)

# mmCIF in the style FoldMason/ESMFold emit here: whitespace-delimited rows whose
# column order is declared by the _atom_site loop header. Note that the varying
# per-residue confidence (pLDDT) sits in `occupancy` and `B_iso_or_equiv` is pinned
# at 1.00 - the reverse of the AlphaFold convention.
CIF_PLDDT_IN_OCCUPANCY = """data_ESMFold
#
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.id
ATOM N N . SER A 1 1 . 1 SER A N 79.95 1.00 10.299 -6.259 -20.907 1 1
ATOM C CA . SER A 1 1 . 1 SER A CA 82.15 1.00 9.152 -5.455 -20.498 1 2
ATOM C C . SER A 1 1 . 1 SER A C 83.45 1.00 8.435 -6.082 -19.308 1 3
ATOM C CA . MET A 1 2 . 2 MET A CA 88.24 1.00 7.601 -5.912 -17.043 1 4
ATOM C CA . GLY A 1 3 . 3 GLY A CA 92.80 1.00 4.298 -7.771 -16.949 1 5
#
"""

# The conventional layout: occupancy pinned at 1.00, pLDDT in the B-factor.
CIF_PLDDT_IN_BFACTOR = CIF_PLDDT_IN_OCCUPANCY.replace(
    "A N 79.95 1.00 ", "A N 1.00 79.95 "
).replace(
    "A CA 82.15 1.00 ", "A CA 1.00 82.15 "
).replace(
    "A C 83.45 1.00 ", "A C 1.00 83.45 "
).replace(
    "A CA 88.24 1.00 ", "A CA 1.00 88.24 "
).replace(
    "A CA 92.80 1.00 ", "A CA 1.00 92.80 "
)


class TestStructureParsing(unittest.TestCase):
    def setUp(self):
        self.repo_root = Path(__file__).resolve().parent.parent
        self.fixtures = self.repo_root / "tests/fixtures/structures"
        self.tmp = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, name, text):
        path = self.tmp_path / name
        path.write_text(text, encoding="utf-8")
        return str(path)

    def test_cif_extracts_only_ca_atoms(self):
        trace = parse_structure_ca(self._write("demo.cif", CIF_PLDDT_IN_OCCUPANCY))
        # 5 ATOM rows, of which 3 are CA
        self.assertEqual(len(trace), 3)
        for residue in trace:
            self.assertEqual(len(residue), 6)

    def test_cif_maps_fields_by_header_not_fixed_columns(self):
        trace = parse_structure_ca(self._write("demo.cif", CIF_PLDDT_IN_OCCUPANCY))
        x, y, z, conf, resnum, resname = trace[0]
        self.assertEqual([x, y, z], [9.2, -5.5, -20.5])
        self.assertEqual(resnum, 1)
        self.assertEqual(resname, "SER")
        self.assertAlmostEqual(conf, 82.2, places=1)

    def test_confidence_picked_from_whichever_column_varies(self):
        """pLDDT must be found whether it sits in occupancy or in the B-factor."""
        from_occupancy = parse_structure_ca(self._write("a.cif", CIF_PLDDT_IN_OCCUPANCY))
        from_bfactor = parse_structure_ca(self._write("b.cif", CIF_PLDDT_IN_BFACTOR))
        self.assertEqual(
            [r[3] for r in from_occupancy],
            [r[3] for r in from_bfactor],
            "same pLDDT values should be recovered from either column layout",
        )
        self.assertEqual([r[3] for r in from_occupancy], [82.2, 88.2, 92.8])

    def test_pdb_fixture_parses(self):
        pdbs = sorted(self.fixtures.glob("*.pdb"))
        self.assertTrue(pdbs, "expected PDB fixtures in tests/fixtures/structures")
        trace = parse_structure_ca(str(pdbs[0]))
        self.assertGreater(len(trace), 50)
        confidences = [r[3] for r in trace]
        self.assertGreater(len(set(confidences)), 1, "per-residue confidence should vary")

    def test_malformed_and_empty_inputs_return_empty(self):
        self.assertEqual(parse_structure_ca(self._write("empty.cif", "")), [])
        self.assertEqual(parse_structure_ca(self._write("empty.pdb", "")), [])
        self.assertEqual(parse_structure_ca(self._write("junk.cif", "not a structure\n")), [])

    def test_extract_ca_dict_over_fixture_directory(self):
        ca = extract_ca_dict(str(self.fixtures), verbose=False)
        self.assertEqual(len(ca), 6, "all 6 benchmark structures should yield a backbone")
        for taxon_id, trace in ca.items():
            self.assertNotIn(".pdb", taxon_id)
            self.assertGreater(len(trace), 0)

    def test_write_ca_bundle_emits_loadable_js(self):
        ca = {"demo": [[1.0, 2.0, 3.0, 90.0, 1, "SER"]]}
        out = self.tmp_path / "ca_structures.js"
        self.assertIsNotNone(write_ca_bundle(ca, str(out)))
        text = out.read_text(encoding="utf-8")
        self.assertTrue(text.startswith("window.CA_STRUCTURES = Object.assign("))
        self.assertIn('"demo"', text)

    def test_write_ca_bundle_skips_empty(self):
        out = self.tmp_path / "none.js"
        self.assertIsNone(write_ca_bundle({}, str(out)))
        self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()


class TestAAAlignmentResolution(unittest.TestCase):
    """Regression: with --tree-type both, the FastTree path never looked for the AA
    alignment beside the 3Di one, so the AA tree was silently skipped."""

    def test_finds_paired_aa_alignment_in_same_directory(self):
        from viral_phylo.tree import resolve_aa_alignment
        with tempfile.TemporaryDirectory() as d:
            tdi = Path(d) / "foldmason.fasta_3di.fa"; tdi.write_text(">a\nAC\n")
            aa = Path(d) / "foldmason.fasta_aa.fa"; aa.write_text(">a\nAC\n")
            self.assertEqual(resolve_aa_alignment(str(tdi)), str(aa))

    def test_prefers_mafft_output_when_present(self):
        from viral_phylo.tree import resolve_aa_alignment
        with tempfile.TemporaryDirectory() as d:
            tdi = Path(d) / "foldmason.fasta_3di.fa"; tdi.write_text(">a\nAC\n")
            (Path(d) / "foldmason.fasta_aa.fa").write_text(">a\nAC\n")
            mafft = Path(d) / "mafft.fasta_aa.fa"; mafft.write_text(">a\nAC\n")
            self.assertEqual(resolve_aa_alignment(str(tdi)), str(mafft))

    def test_never_returns_the_3di_file_itself(self):
        from viral_phylo.tree import resolve_aa_alignment
        with tempfile.TemporaryDirectory() as d:
            only = Path(d) / "alignment.fa"; only.write_text(">a\nAC\n")
            self.assertIsNone(resolve_aa_alignment(str(only)))
