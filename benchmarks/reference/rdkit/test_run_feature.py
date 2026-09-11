"""Focused regressions for the independently defined writer reference contract."""

import unittest

from rdkit import Chem

import run_feature as reference


class IsomericReferenceTests(unittest.TestCase):
    def test_required_charge_brackets_fix_total_hydrogens_without_mutating_source(self):
        # Model chemical normalization introducing charge on an inferred-H atom.
        # This is element independent: an emitted charged nitrogen must also
        # carry its hydrogen count in bracket syntax.
        source = Chem.MolFromSmiles("N")
        source.GetAtomWithIdx(0).SetFormalCharge(1)
        source.GetAtomWithIdx(0).UpdatePropertyCache(strict=True)
        self.assertEqual(source.GetAtomWithIdx(0).GetNumImplicitHs(), 4)
        projected, changes = reference.required_smiles_bracket_declarations(source)
        atom = projected.GetAtomWithIdx(0)
        self.assertTrue(atom.GetNoImplicit())
        self.assertEqual(atom.GetNumExplicitHs(), 4)
        self.assertEqual(atom.GetNumImplicitHs(), 0)
        self.assertFalse(source.GetAtomWithIdx(0).GetNoImplicit())
        self.assertEqual(source.GetAtomWithIdx(0).GetNumExplicitHs(), 0)
        self.assertEqual(len(changes), 1)
        self.assertFalse(changes[0]["source"]["no_implicit_hydrogens"])
        self.assertTrue(changes[0]["emitted"]["no_implicit_hydrogens"])
        self.assertEqual(Chem.MolToSmiles(source), Chem.MolToSmiles(projected))

    def test_uncharged_metal_neighbor_keeps_its_source_hydrogen_declaration(self):
        source = Chem.MolFromSmiles("O[Fe]=O")
        before = reference.smiles_perceived_semantic_record(source)
        projected, changes = reference.required_smiles_bracket_declarations(source)
        self.assertEqual(changes, [])
        self.assertEqual(reference.smiles_perceived_semantic_record(projected), before)
        self.assertFalse(projected.GetAtomWithIdx(0).GetNoImplicit())
        self.assertEqual(projected.GetAtomWithIdx(0).GetNumImplicitHs(), 1)

    def test_record_retains_full_source_and_reference_emission_evidence(self):
        source = Chem.MolFromSmiles("OCl(=O)(O)O")
        source_projection = reference.smiles_perceived_semantic_record(source)
        evidence = []
        record = {"record_index": 0, "status": "ok", "title": None,
                  "smiles": "OCl(=O)(O)O", "mol": source}
        expected = reference.isomeric_smiles_record(record, evidence)
        self.assertEqual(evidence[0]["decoded_source"]["normalized_perceived"], source_projection)
        changes = evidence[0]["mandatory_bracket_declarations"]
        self.assertEqual([(item["source"]["symbol"], item["source"]["formal_charge"])
                          for item in changes], [("Cl", 1), ("O", -1)])
        self.assertTrue(evidence[0]["whole_graph_check"]["equal"])
        emitted = Chem.MolFromSmiles(evidence[0]["rdkit_emitted_smiles"])
        self.assertEqual(expected["normalized_perceived"], reference.smiles_perceived_semantic_record(emitted))
        self.assertEqual(source_projection, reference.smiles_perceived_semantic_record(source))


if __name__ == "__main__":
    unittest.main()
