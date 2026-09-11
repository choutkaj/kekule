import gzip
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import migrate_v2 as migration


class MigrationTests(unittest.TestCase):
    def test_tetrahedral_permutation_preserves_configuration(self):
        element = {"type": "tetrahedral", "index": 3, "center_atom_index": 2,
                   "carriers": [{"atom_index": 4}, {"atom_index": 1},
                                {"atom_index": 5}, {"implicit_hydrogen": True}],
                   "orientation": "clockwise", "source": "smiles", "specifiedness": "specified"}
        result = migration.canonical_element(element, None)
        self.assertEqual(result["orientation"], "counter_clockwise")
        self.assertEqual(result["carriers"], [{"atom_index": 1}, {"atom_index": 4},
                                              {"atom_index": 5}, {"implicit_hydrogen": True}])
        self.assertEqual(element["orientation"], "clockwise")
        self.assertEqual(migration.canonical_element(result, None), result)
        document = migration.source_document("source", "smiles", [], [result], [element["source"]])
        self.assertEqual(document["stereo_sources"], [{"element_index": 3, "source": "smiles", "specifiedness": "specified"}])

    def test_archive_dry_run_is_read_only_and_write_keeps_original_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original = root / "original.json.gz"
            content = gzip.compress(json.dumps({"schema_version": 1}).encode(), mtime=0)
            original.write_bytes(content)
            with patch.object(migration, "ARCHIVE", root / "archive"):
                self.assertEqual(migration.archive(original, "corpus", "feature"), original)
                self.assertFalse((root / "archive").exists())
                archived = migration.archive(original, "corpus", "feature", write=True)
                self.assertEqual(archived.read_bytes(), content)
                original.write_bytes(gzip.compress(json.dumps({"schema_version": 2}).encode()))
                self.assertEqual(migration.archive(original, "corpus", "feature"), archived)
                self.assertEqual(archived.read_bytes(), content)

    def test_scoped_publication_preserves_other_corpus_provenance(self):
        previous = {"schema_version": 2, "rdkit": migration.VERSION,
                    "fixtures": [{"path": "a", "old_sha256": "unchanged"}, {"path": "b", "new_sha256": "old"}]}
        self.assertEqual(migration.merged_provenance(previous, [{"path": "b", "new_sha256": "new"}]),
                         [{"path": "a", "old_sha256": "unchanged"}, {"path": "b", "new_sha256": "new"}])
        with self.assertRaises(ValueError):
            migration.merged_provenance({**previous, "rdkit": "another version"}, [])

    def test_overlap_audit_rejects_a_dropped_reference_field(self):
        with self.assertRaises(ValueError):
            migration.changed_paths({"stereo": {"descriptor": "R"}}, {"stereo": {}})
        self.assertEqual(migration.changed_paths({"descriptor": "R"}, {"descriptor": "S"}), ["$.descriptor"])


if __name__ == "__main__":
    unittest.main()
