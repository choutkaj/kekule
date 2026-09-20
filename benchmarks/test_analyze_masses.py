import copy
import gzip
import json
import math
from pathlib import Path
import tempfile
import unittest

import analyze_masses as analysis


def record(mass=12.0, isotope=None, count=1, charge=0):
    return {"record_index": 0, "status": "ok", "title": "example",
            "formula": {"formal_charge": charge, "terms": [
                {"element": "C", "isotope": isotope, "count": count}]},
            "average_mass_da": mass, "monoisotopic_mass_da": mass}


class ReconstructionTests(unittest.TestCase):
    def test_isotope_and_ion_contributions_are_separate(self):
        electron = analysis.NATIVE_ELECTRON_MASS
        actual = record(26.0 + electron, isotope=13, count=2, charge=-1)
        expected = record(26.02, isotope=13, count=2, charge=-1)
        expected["monoisotopic_mass_da"] += 0.00055
        original = copy.deepcopy((expected, actual))

        def native(symbol, isotope):
            self.assertEqual((symbol, isotope), ("C", 13))
            return (13.0, 13.0)

        result = analysis.analyze_record(expected, actual, native, lambda *_: (13.01, 13.01), 0.00055)
        self.assertEqual(result["classification"], "mass_models_verified")
        self.assertAlmostEqual(result["deltas"]["average_mass_da"]["atomic_data_contribution_da"], -0.02)
        self.assertEqual(result["deltas"]["average_mass_da"]["charge_convention_contribution_da"], electron)
        self.assertEqual(result["deltas"]["monoisotopic_mass_da"]["charge_convention_contribution_da"], electron - 0.00055)
        self.assertEqual((expected, actual), original)

    def test_matching_formula_does_not_explain_a_corrupt_mass(self):
        expected, actual = record(), record()
        actual["average_mass_da"] += 0.01
        result = analysis.analyze_record(expected, actual, lambda *_: (12.0, 12.0), lambda *_: (12.0, 12.0), 0.00055)
        self.assertEqual(result["classification"], "unexplained_mass")
        self.assertGreater(result["native"]["average_mass_da"]["residual_da"], 0.009)

    def test_summation_roundoff_is_reported(self):
        actual = record(math.nextafter(12.0, math.inf))
        result = analysis.analyze_record(record(), actual, lambda *_: (12.0, 12.0), lambda *_: (12.0, 12.0), 0.00055)
        self.assertEqual(result["classification"], "mass_models_verified")
        self.assertGreater(result["native"]["average_mass_da"]["residual_da"], 0)

    def test_formula_and_identity_changes_remain_unexplained(self):
        expected = record()
        for actual in (record(count=2), record(isotope=13), record(charge=1)):
            result = analysis.analyze_record(expected, actual, None, None, 0.00055)
            self.assertEqual(result["classification"], "formula_difference")
        actual = record()
        actual["title"] = "other"
        self.assertEqual(analysis.analyze_record(expected, actual, None, None, 0.00055)["classification"], "record_identity_difference")
        expected["status"] = actual["status"] = "unsupported"
        self.assertEqual(analysis.analyze_record(expected, actual, None, None, 0.00055)["classification"], "invalid_record_status")

    def test_invalid_or_missing_atomic_data_is_not_verified(self):
        for constants in (lambda *_: (None, 12.0), lambda *_: (0.0, 12.0), lambda *_: (math.inf, 12.0)):
            result = analysis.analyze_record(record(), record(), constants, constants, 0.00055)
            self.assertEqual(result["classification"], "unavailable_reconstruction")
        actual = record(math.nan)
        result = analysis.analyze_record(record(), actual, lambda *_: (12.0, 12.0), lambda *_: (12.0, 12.0), 0.00055)
        self.assertEqual(result["classification"], "invalid_mass_observation")

    def test_invalid_formula_counts_and_duplicate_terms_fail(self):
        for amount in (0, -1, True, 1.5, 2**53 + 1):
            with self.assertRaises(ValueError):
                analysis.reconstruct(record(count=amount)["formula"], lambda *_: (12.0, 12.0), (0.0, 0.0))
        formula = record()["formula"]
        formula["terms"] *= 2
        with self.assertRaises(ValueError):
            analysis.reconstruct(formula, lambda *_: (12.0, 12.0), (0.0, 0.0))


class ReportTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.cases = self.root / "cases.jsonl.gz"
        case = {"dataset": "unit", "feature": "descriptor.molecular", "id": "example",
                "fixture": "focused-unit-fixture", "record_index": 0, "status": "agrees",
                "expected": {"value": {"records": [record()]}},
                "actual": {"value": {"records": [record()]}}}
        with gzip.open(self.cases, "wt", encoding="utf-8") as stream:
            for row in (case, {**case, "status": "error"}, {**case, "status": "not_applicable"}):
                stream.write(json.dumps(row) + "\n")
        self.report = {"complete": True, "error": None, "implementation": {"revision": "test"},
                       "cases": str(self.cases), "results": [{"feature": "descriptor.molecular",
                       "cases": 2, "not_applicable": 1,
                       "golden": {"reference": {"tool": "rdkit", "version": "test"}}}]}
        self.report_path = self.root / "report.json"
        self.output = self.root / "analysis.json"

    def analyze(self):
        self.report_path.write_text(json.dumps(self.report), encoding="utf-8")
        return analysis.analyze_report(self.report_path, self.output, lambda *_: (12.0, 12.0),
                                       lambda *_: (12.0, 12.0), 0.00055, {"rdkit_version": "test"}, {})

    def test_streaming_report_retains_errors_and_input_fingerprints(self):
        original_hash = analysis.sha256(self.cases)
        result = self.analyze()
        self.assertEqual(result["case_outcomes"], {"agrees": 1, "error": 1, "not_applicable": 1})
        self.assertEqual(result["component_classifications"], {"mass_models_verified": 1})
        self.assertEqual(result["cases_sha256"], original_hash)
        self.assertEqual(analysis.sha256(self.cases), original_hash)
        with gzip.open(result["details"], "rt", encoding="utf-8") as stream:
            self.assertEqual(len(list(stream)), 1)
        with self.assertRaisesRegex(ValueError, "already exists"):
            self.analyze()

    def test_incomplete_reports_and_wrong_reference_versions_are_rejected(self):
        self.report["complete"] = False
        with self.assertRaisesRegex(ValueError, "complete report"):
            self.analyze()
        self.report["complete"] = True
        self.report["results"][0]["golden"]["reference"]["version"] = "other"
        with self.assertRaisesRegex(ValueError, "installed RDKit"):
            self.analyze()
        self.assertFalse(self.output.exists())

    def test_truncated_case_files_do_not_publish_a_complete_analysis(self):
        self.report["results"][0]["cases"] += 1
        with self.assertRaisesRegex(ValueError, "row count"):
            self.analyze()
        self.assertFalse(self.output.exists())


if __name__ == "__main__":
    unittest.main()
