"""Regression checks for the optional differential runner's comparison contract."""

import json
from pathlib import Path
import tempfile
import unittest

from rdkit import Chem

import compare_cip


class ComparisonContract(unittest.TestCase):
    def test_complete_maps_include_absent_labels_and_graph_counts(self):
        expected = {"status": "ok", "atom_count": 2, "bond_count": 1,
                    "labels": {"atoms": [[0, "LowerR"]], "bonds": []}}
        self.assertEqual(compare_cip.comparison(expected, expected), "match")
        for changed in ({"atom_count": 3}, {"bond_count": 2},
                        {"labels": {"atoms": [], "bonds": []}},
                        {"labels": {"atoms": [[0, "R"]], "bonds": []}}):
            self.assertEqual(compare_cip.comparison(expected, expected | changed), "difference")

    def test_reference_failure_is_never_a_match(self):
        failed = {"status": "error", "stage": "cip"}
        self.assertEqual(compare_cip.comparison(failed, failed), "reference_error")
        self.assertEqual(compare_cip.summarize([
            {"comparison": "reference_error"}, {"comparison": "reference_only"}
        ])["reference_error"], 1)

    def test_old_json_expectations_are_recomputed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "old.json"
            path.write_text(json.dumps({"cases": [
                {"smiles": "CCO", "expected": "stale"},
                {"smiles": "CCO", "expected": "different"},
                {"smiles": "[Na+].[2H][C@](F)(Cl)Br"},
            ]}), encoding="utf-8")
            cases = compare_cip.read_cases([path], "json", False)
            self.assertEqual([case["smiles"] for case in cases],
                             ["CCO", "[Na+].[2H][C@](F)(Cl)Br"])
            self.assertTrue(all("expected" not in case for case in cases))

    def test_reference_retains_explicit_hydrogen_and_component_indices(self):
        Chem.SetUseLegacyStereoPerception(False)
        result = compare_cip.reference("[Na+].[2H][C@](F)(Cl)Br", "assertions", 1_000_000)
        self.assertEqual(result["status"], "ok")
        self.assertEqual(result["atom_count"], 6)
        self.assertEqual(result["bond_count"], 4)
        self.assertEqual(result["labels"], {"atoms": [[2, "S"]], "bonds": []})

    def test_reference_records_chirality_removed_by_sanitization(self):
        source = "CC=1C=CC=2[N@@]3CC=4C=C(C=CC4[N@](CC2C1)C3)C"
        cleaned = compare_cip.reference(source, "sanitized", 1_000_000)
        retained = compare_cip.reference(source, "assertions", 1_000_000)
        self.assertTrue(cleaned["evidence"]["parsed_tags"])
        self.assertEqual(cleaned["evidence"]["sanitized_tags"], [])
        self.assertEqual(retained["evidence"]["parsed_tags"], retained["evidence"]["sanitized_tags"])

if __name__ == "__main__":
    unittest.main()
