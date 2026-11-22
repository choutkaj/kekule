"""Offline checks of the dictionary charge table generator."""
import gzip
import tempfile
import unittest
from pathlib import Path

from charges import blocks, render, rows, tokens

DICTIONARY = """\
data_TRS
_chem_comp.id TRS
_chem_comp.name "2-AMINO-2-HYDROXYMETHYL-PROPANE-1,3-DIOL"
_chem_comp.pdbx_synonyms
;TRIS BUFFER
loop_
_chem_comp.id NOT
;
_chem_comp.pdbx_modified_date 2020-06-17
#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
TRS C  C  C 0
TRS N  N  N 1
TRS O1 O1 O 0
TRS HN1 1HN H 0
#
data_GOL
_chem_comp.id GOL
_chem_comp.pdbx_modified_date 2026-11-13
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
GOL C1 C1 C 0
GOL O1 O1 O 0
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
GOL C1 O1
data_2ES
_chem_comp.id 2ES
_chem_comp.pdbx_modified_date 2026-09-16
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
2ES AU AU AU ?
2ES "C1'" 'C1"' C 0
2ES CL CL CL -1
data_PT
_chem_comp.id PT
_chem_comp.pdbx_modified_date 2011-06-04
_chem_comp_atom.comp_id PT
_chem_comp_atom.atom_id PT
_chem_comp_atom.alt_atom_id PT
_chem_comp_atom.type_symbol PT
_chem_comp_atom.charge 2
"""


class Tokens(unittest.TestCase):
    def test_quoted_values_end_at_a_quote_before_a_space(self):
        self.assertEqual(tokens("2ES \"C1'\" 'C1\"' C 0"), ["2ES", "C1'", 'C1"', "C", "0"])
        self.assertEqual(tokens("X 'O1'1' O"), ["X", "O1'1", "O"])


class Generator(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.source = Path(directory.name) / "components.cif.gz"
        with gzip.open(self.source, "wt", encoding="utf-8") as text:
            text.write(DICTIONARY)

    def test_blocks_skip_text_fields_and_read_single_atoms(self):
        found = list(blocks(DICTIONARY.splitlines()))
        self.assertEqual([block["id"] for block in found], ["TRS", "GOL", "2ES", "PT"])
        self.assertEqual(len(found[0]["atoms"]), 4)
        self.assertEqual(found[3]["single"]["charge"], "2")

    def test_rows_list_charged_and_unknown_components_only(self):
        table, release = rows(self.source)
        self.assertEqual(release, "2026-11-13")
        self.assertEqual(
            table,
            [
                ("2ES", "AU/AU=?", "C1'", "CL/CL=-1"),
                ("PT", "PT/PT=+2"),
                ("TRS", "C", "N=+1", "O1"),
            ],
        )

    def test_render_states_the_release(self):
        text = render(*rows(self.source))
        self.assertIn("# Release: entries modified up to 2026-11-13.\n", text)
        self.assertTrue(text.endswith("TRS\tC\tN=+1\tO1\n"))


if __name__ == "__main__":
    unittest.main()
