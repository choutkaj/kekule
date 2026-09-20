import unittest
from smarts_conformance import classify, query_rows
from pathlib import Path
import hashlib
import json


class SmartsConformanceTests(unittest.TestCase):
    def test_complete_rdkit_extraction_and_provenance(self):
        root = Path(__file__).parent / "smarts-fixtures/rdkit-queries"
        lock = json.loads((root / "sources.lock.json").read_text())
        total = 0
        for source in lock["upstream"]:
            path = root/source["path"]
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), source["sha256"])
            if source["layout"] == "license":
                continue
            expected = []
            for line in path.read_text().splitlines():
                if not line.strip() or line.lstrip().startswith(("#", "//")):
                    continue
                if source["layout"] == "label-tabs-query":
                    expected.append([field.strip() for field in line.split("\t") if field.strip()][1])
                else:
                    expected.append(line.split()[0])
            actual = list(query_rows([root/"data"/(path.stem+".smarts")]))
            self.assertEqual([row["smarts"] for row in actual], expected)
            total += len(actual)
        self.assertEqual(total, 518)
        for pack in lock["packs"]:
            self.assertEqual(hashlib.sha256((root/pack["path"]).read_bytes()).hexdigest(), pack["sha256"])

    def test_comparison_retains_mapping_order_and_stereo_observations(self):
        expected = dict(status="ok", matches=[[0, 1], [1, 0]], tagged_matches=[[0, 1], [1, 0]])
        self.assertEqual(classify(expected, expected), "equal")
        changed = dict(expected, matches=[[0, 1]])
        self.assertEqual(classify(expected, changed), "mismatch")
        changed = dict(expected, tagged_matches=[[0, 1], [0, 1]])
        self.assertEqual(classify(expected, changed), "mismatch")

    def test_invalid_reference_is_not_a_dialect_exclusion(self):
        actual = dict(status="parse_error", kind="Unsupported")
        self.assertEqual(classify(dict(status="parse_error"), actual), "both_rejected")
        self.assertEqual(classify(dict(status="ok"), actual), "dialect_exclusion")
        self.assertEqual(classify(dict(status="target_error"), actual), "reference_error")

    def test_complete_openff_source_rows_and_provenance(self):
        root = Path(__file__).parent / "smarts-fixtures/openff-smarts"
        lock = json.loads((root / "sources.lock.json").read_text())
        for source in lock["sources"]:
            self.assertEqual(hashlib.sha256((root/source["path"]).read_bytes()).hexdigest(), source["sha256"])
        rows = list(query_rows(sorted(root.glob("*.offxml"))))
        self.assertEqual(len(rows), 521)
        self.assertTrue(any(row["smarts"] == "[N:1](H:2)(H:3)" for row in rows))
        self.assertTrue({"VirtualSites", "LibraryCharges", "ChargeIncrementModel", "ImproperTorsions"} <= {row["section"] for row in rows})


if __name__ == "__main__":
    unittest.main()
