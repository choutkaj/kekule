"""Offline audit-contract tests; no scientific dependencies or network required."""
import json
import gzip
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import audit
import hashlib
from paths import historical, historical_bytes


class AuditTests(unittest.TestCase):
    def test_source_bytes_are_pinned(self):
        lock = audit.verify_sources()
        self.assertEqual(len(lock["sources"]), 11)
        for source in lock["sources"]:
            self.assertIn(source["revision"], source["url"])
            self.assertEqual(len(source["revision"]), 40)
        with patch.object(audit, "digest", return_value="changed"):
            with self.assertRaisesRegex(ValueError, "checksum"):
                audit.verify_sources()

    def test_inventory_keeps_every_parameter_and_section_setting(self):
        observed = audit.inventory(audit.HERE / "fixtures/rosemary.offxml")
        self.assertEqual(observed["pattern_count"], 371)
        sections = {section["name"]: section for section in observed["sections"]}
        self.assertEqual({name: len(s["parameters"]) for name, s in sections.items()},
                         dict(Constraints=1, Bonds=90, Angles=42, ProperTorsions=187,
                              ImproperTorsions=7, vdW=35, Electrostatics=0,
                              LibraryCharges=9, NAGLCharges=0))
        self.assertEqual(sections["NAGLCharges"]["attributes"]["model_file_hash"], audit.MODEL_HASH)
        self.assertEqual(sections["Bonds"]["attributes"]["fractional_bondorder_method"], "AM1-Wiberg")
        self.assertEqual(observed["bondorder_parameters"], [])
        self.assertNotIn("distance", sections["Constraints"]["parameters"][0])
        self.assertEqual(sections["vdW"]["attributes"]["scale14"], "0.5")
        self.assertEqual(sections["Electrostatics"]["attributes"]["scale14"], "0.8333333333")

    def test_selected_cases_have_external_provenance(self):
        cases = audit.cases()
        self.assertEqual(len(cases), 23)
        self.assertEqual(len({c["id"] for c in cases}), len(cases))
        self.assertTrue(all(c["source"] for c in cases))
        with self.assertRaisesRegex(ValueError, "not present"):
            audit.literal("test_lookups.py", "invented molecule")

    def test_equal_errors_are_not_agreements(self):
        error = dict(status="error", message="failure")
        self.assertEqual(audit.classify(error, error), "reference_error")
        self.assertEqual(audit.classify(dict(status="ok"), error), "implementation_error")

    def test_mapping_order_features_and_duplicates_are_asserted(self):
        reference = dict(status="ok", tagged_matches=[[1, 2], [2, 1]], rings=[True, False])
        self.assertEqual(audit.classify(reference, reference), "equal")
        for changed in [dict(reference, tagged_matches=[[1, 2]]),
                        dict(reference, tagged_matches=[[1, 2], [1, 2]]),
                        dict(reference, rings=[False, True])]:
            self.assertEqual(audit.classify(reference, changed), "mismatch")

    def test_incomplete_native_stream_is_rejected(self):
        with patch.object(audit.subprocess, "run") as run:
            run.return_value.stdout = '{}\n'
            with self.assertRaisesRegex(ValueError, "Incomplete"):
                audit.native("observer", [{}, {}])

    def test_publication_rejects_nonfinite_observations(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            audit.publish(path, dict(complete=False))
            with self.assertRaises(ValueError):
                audit.publish(path, dict(charge=float("nan")))
            self.assertEqual(json.loads(path.read_text()), dict(complete=False))

    def test_configuration_preserves_class_identity_without_generic_string_fallback(self):
        self.assertEqual(audit.describe_configuration(str), {"python_type": "builtins.str"})
        with self.assertRaises(TypeError):
            audit.describe_configuration(object())

    def test_system_errors_and_permutation_disagreements_fail_audit(self):
        record = dict(input={"id": "regression"}, status="ok", graph={"comparison": "equal"},
                      charges={"status": "ok"}, system={"status": "ok"}, domain=[True, ""])
        report = dict(stage="reference", complete=True, records=[record], smarts=[],
                      lookup_contract=[{"passed": True}])
        self.assertFalse(audit.audit_failed(report))
        record["system"] = {"status": "error"}
        self.assertTrue(audit.audit_failed(report))
        record["system"] = {"status": "ok"}
        record["charge_permutation"] = {"passed": False}
        self.assertTrue(audit.audit_failed(report))
        self.assertEqual(audit.summarize(report)["charge_permutation_failures"], ["regression"])

    def test_domain_rejection_requires_successful_system_assignment(self):
        record = dict(input={"id": "ion"}, status="ok", graph={"comparison": "equal"},
                      charges={"status": "error"}, system={"status": "ok"}, domain=[False, "domain"])
        report = dict(stage="reference", complete=True, records=[record], smarts=[],
                      lookup_contract=[{"passed": True}])
        self.assertFalse(audit.audit_failed(report))
        record["domain"] = [True, ""]
        self.assertTrue(audit.audit_failed(report))

    def test_archived_reports_are_complete_and_checksum_bound(self):
        manifest = historical("reports.lock.json")
        for entry in manifest["reports"]:
            data = historical_bytes(entry["path"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), entry["sha256"])
            report = historical(entry["path"])
            self.assertTrue(report["complete"])
            self.assertEqual(report["generation"], entry["generation"])
            self.assertEqual(len(report["records"]), 23)
            self.assertEqual(len(report["smarts"]), 8533)
            self.assertEqual(report["summary"]["graph"], {"equal": 23})
            self.assertEqual(report["summary"]["smarts"], {"equal": 8533})
            for row in report["smarts"]:
                self.assertEqual(row["expected"], row["actual"])

    def test_reference_findings_remain_visible(self):
        report = historical("reference.json.gz")
        self.assertFalse(audit.audit_failed(report))
        self.assertEqual(report["lookup_table_sizes"], {"am1bcc_charges": 13944})
        self.assertTrue(report["summary"]["lookup_contract_passed"])
        self.assertEqual(report["summary"]["charge_permutation_failures"],
                         ["lookup-entry-sulfide", "lookup-reordered-sulfide", "lookup-miss-ethane"])
        self.assertEqual(audit.summarize(report)["charge_permutation_failures"], [])
        sodium = next(r for r in report["records"] if r["input"]["id"] == "cid_923")
        self.assertEqual(sodium["charges"]["status"], "error")
        self.assertEqual(sodium["domain"][0], False)
        self.assertEqual(sodium["system"]["charges"], [{"atoms": [0], "charge": 1.0}])

    def test_charge_tolerance_boundary(self):
        self.assertFalse(audit.permutation_failed({"max_abs_error": 5e-5}))
        self.assertTrue(audit.permutation_failed({"max_abs_error": 5.001e-5}))
        self.assertTrue(audit.permutation_failed({"max_abs_error": float("nan")}))

    def test_runtime_assets_and_extended_reference_provenance(self):
        root = audit.HERE.parents[1]
        crate = root / "crates/kekule-openff"
        for name in ("rosemary.offxml", "LICENSE-forcefields", "LICENSE-nagl", "LICENSE-models"):
            self.assertEqual((crate / "data" / name).read_bytes(), (audit.HERE / "fixtures" / name).read_bytes())
        self.assertEqual((crate / "tests/fixtures/audit.json.gz").read_bytes(), historical_bytes("reference.json.gz"))
        self.assertEqual(audit.digest(crate / "data/normalizations.json"), "06dcbb613e5cfe75441c09b8b3b6fdc10f47b78430928a0ac835487d3d3be22c")
        sources = json.loads((audit.HERE / "fixtures/supplementary-sources.lock.json").read_text())["sources"]
        self.assertEqual(len(sources), 10)
        for entry in sources:
            self.assertEqual(audit.digest(audit.HERE / entry["path"]), entry["sha256"])
        manifest = historical("implementation-reports.lock.json")
        for entry in manifest["reports"]:
            self.assertEqual(hashlib.sha256(historical_bytes(entry["path"])).hexdigest(), entry["sha256"])
        extended = historical("native-reference.json.gz")
        self.assertEqual(len(extended["records"]), 66)
        self.assertTrue(all(r["passed"] for r in extended["records"]))
        lookup = historical("lookup-identity.json.gz")
        self.assertEqual(lookup["total"], 13944)
        self.assertEqual(lookup["passed"] + len(lookup["disagreements"]), lookup["total"])
        self.assertEqual(lookup["classification"], {"reference-rejected": 689, "reference-and-native-agree-stored-key-differs": 19, "native-rejected-reference-accepted": 1})


if __name__ == "__main__":
    unittest.main()
