"""Offline checks that the bundled component templates match their lock."""
import json
import re
import unittest

from fetch import DATA, HERE, ROOT, digest, identifiers


class BundledComponentProvenance(unittest.TestCase):
    def setUp(self):
        self.lock = json.loads((HERE / "sources.lock.json").read_text(encoding="utf-8"))

    def test_shipped_entries_match_their_lock(self):
        for entry in self.lock["components"]:
            path = ROOT / entry["path"]
            self.assertEqual(path.stat().st_size, entry["bytes"], entry["id"])
            self.assertEqual(digest(path), entry["sha256"], entry["id"])

    def test_selection_is_complete_and_exact(self):
        locked = [entry["id"] for entry in self.lock["components"]]
        self.assertEqual(locked, identifiers())
        self.assertEqual(self.lock["groups"], json.loads(json.dumps(self.lock["groups"])))
        shipped = sorted(path.stem for path in DATA.glob("*.cif"))
        self.assertEqual(shipped, sorted(locked))

    def test_rust_bundle_lists_the_locked_selection(self):
        source = (ROOT / "crates/kekule/src/io/components/bundled.rs").read_text(encoding="utf-8")
        listed = re.search(r"bundled_components!\[(.*?)\]", source, re.DOTALL)
        self.assertIsNotNone(listed)
        rust = re.findall(r'"([A-Z0-9]+)"', listed.group(1))
        self.assertEqual(rust, identifiers())


class DictionaryCodeProvenance(unittest.TestCase):
    def setUp(self):
        self.lock = json.loads((HERE / "codes.lock.json").read_text(encoding="utf-8"))

    def test_shipped_table_matches_its_lock(self):
        path = ROOT / self.lock["output"]
        self.assertEqual(digest(path), self.lock["output_sha256"])

    def test_bundled_polymer_entries_are_in_the_table(self):
        table = (ROOT / self.lock["output"]).read_text(encoding="utf-8")
        identifiers = {line.split("	")[0] for line in table.splitlines() if not line.startswith("#")}
        lock = json.loads((HERE / "sources.lock.json").read_text(encoding="utf-8"))
        for identifier in lock["groups"]["amino acids"] + lock["groups"]["nucleotides"]:
            self.assertIn(identifier, identifiers)


if __name__ == "__main__":
    unittest.main()
