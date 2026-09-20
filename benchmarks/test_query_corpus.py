"""Verify that the bundled queries preserve every row of their pinned sources."""
import hashlib
import json
from pathlib import Path
import unittest


class QueryCorpusTests(unittest.TestCase):
    def test_complete_tables_are_extracted_without_query_rewriting_or_filtering(self):
        root = Path(__file__).resolve().parent / "corpora/rdkit-queries"
        lock = json.loads((root / "sources.lock.json").read_text(encoding="utf-8"))
        entries = {entry["id"]: entry for entry in lock["entries"]}
        self.assertEqual(len(entries), len(lock["entries"]))
        self.assertEqual(lock["upstream_package"]["version"], "2026.03.3")
        self.assertEqual(lock["upstream_package"]["license"], "BSD-3-Clause")
        self.assertEqual(len(lock["upstream_package"]["sha256"]), 64)
        seen = set()
        counts = []
        for source in lock["upstream"]:
            payload = (root / source["path"]).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), source["sha256"])
            if source["layout"] == "license":
                self.assertIn(b"BSD 3-Clause License", payload)
                continue
            pack = next(pack for pack in lock["packs"] if pack["upstream_file"] == source["path"])
            expected_lines, expected_ids = [], []
            for number, raw in enumerate(payload.decode("utf-8").splitlines(), 1):
                line = raw.strip()
                if not line or line.startswith(("#", "//")):
                    continue
                if source["layout"] == "label-tabs-query":
                    fields = [part.strip() for part in line.split("\t") if part.strip()]
                    label, query = fields[:2]
                else:
                    self.assertEqual(source["layout"], "query-whitespace-label")
                    parts = line.split(None, 1)
                    query, label = parts[0], parts[1] if len(parts) == 2 else ""
                source_id = f"RDKit:{Path(source['path']).stem}:{number}"
                self.assertEqual(entries[source_id]["upstream_file"], source["path"])
                self.assertEqual(entries[source_id]["upstream_line"], number)
                expected_ids.append(source_id)
                expected_lines.append(query + ("\t" + label if label else "") + "\n")
            actual = (root / pack["path"]).read_bytes()
            self.assertEqual(actual, "".join(expected_lines).encode("utf-8"))
            self.assertEqual(hashlib.sha256(actual).hexdigest(), pack["sha256"])
            self.assertEqual(pack["members"], expected_ids)
            self.assertFalse(seen.intersection(expected_ids))
            seen.update(expected_ids)
            counts.append(len(expected_ids))
        self.assertEqual(counts, [38, 52, 428])
        self.assertEqual(seen, entries.keys())


if __name__ == "__main__":
    unittest.main()
