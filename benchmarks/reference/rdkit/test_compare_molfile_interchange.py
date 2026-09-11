"""Regression checks for the full Molfile interchange comparison contract."""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from rdkit import Chem

import compare_molfile_interchange as comparison


class InterchangeContract(unittest.TestCase):
    def test_unlabelled_atom_and_bond_loss_are_differences(self):
        molecule = Chem.MolFromSmiles("CCO")
        expected = {"graph": comparison.chemical_graph(molecule), "labels": [], "groups": []}
        self.assertEqual(comparison.compare_output(expected, molecule, False)["status"], "match")
        for changed in ["CC", "C.CO", "C=CO", "CCN", "[13CH3]CO", "CC[O-]"]:
            with self.subTest(changed=changed):
                result = comparison.compare_output(expected, Chem.MolFromSmiles(changed), False)
                self.assertEqual(result["status"], "difference")
                self.assertIn("graph", result["differences"])

    def test_clearing_double_labels_preserves_axial_labels(self):
        molecule = Chem.MolFromSmiles("FC=CCl")
        with patch.object(comparison, "labels", return_value=[["bond", 0, "M"], ["bond", 1, "E"]]):
            self.assertEqual(comparison.expected_labels(molecule, True), [["bond", 0, "M"]])

    def test_existing_evidence_is_not_overwritten(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "existing.json"
            output.write_text("prior evidence", encoding="utf-8")
            with patch("sys.argv", ["check", "--probe", "unused", "--output", str(output)]), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
                comparison.main()
            self.assertEqual(raised.exception.code, 2)
            self.assertEqual(output.read_text(encoding="utf-8"), "prior evidence")

    def test_unpinned_rdkit_is_rejected_before_reading_inputs(self):
        with patch.object(comparison.rdBase, "rdkitVersion", "2026.03.5"), patch("sys.argv", ["check", "--probe", "unused", "--output", "unused.json"]), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
            comparison.main()
        self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
