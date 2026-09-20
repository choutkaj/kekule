"""Check byte provenance and complete input selection for the structure corpus."""
import hashlib
import json
from pathlib import Path, PurePosixPath
import unittest


class StructureCorpusTests(unittest.TestCase):
    def test_pinned_structure_inputs_are_complete_and_unmodified(self):
        root = Path(__file__).resolve().parent / "corpora/rdkit-structures"
        lock = json.loads((root / "sources.lock.json").read_text(encoding="utf-8"))
        repository = lock["upstream_repository"]
        self.assertEqual(repository["commit"], "e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f")
        self.assertEqual(repository["tag"], "Release_2026_03_3")
        self.assertEqual(repository["license"], "BSD-3-Clause")
        for source in lock["upstream"]:
            self.assertEqual(hashlib.sha256((root / source["path"]).read_bytes()).hexdigest(),
                             source["sha256"])
        commit = json.loads((root / "upstream/commit.json").read_bytes())
        self.assertEqual(commit["sha"], repository["commit"])
        tree = json.loads((root / "upstream/tree.json").read_bytes())
        self.assertEqual(tree["sha"], repository["test_data_tree_sha1"])
        self.assertFalse(tree["truncated"])
        families = ["BMS-986142", "JDQ443", "Mrtx1719", "RP-6306", "Sotorasib", "ZM374979"]
        self.assertEqual(lock["compound_families"], families)
        selected = {}
        for item in tree["tree"]:
            path = PurePosixPath(item["path"])
            if item["type"] != "blob" or ".expected" in path.name:
                continue
            drug = path.parent.as_posix() == "atropisomers" and any(
                path.name == family + ".sdf" or path.name.startswith(family + "_")
                for family in families)
            chebi = path.parent.as_posix() == "." and path.name.startswith("chebi_")
            if (drug and path.suffix == ".sdf") or (chebi and path.suffix == ".mol"):
                selected[path.as_posix()] = item
        self.assertEqual(len(selected), 50)
        self.assertTrue(any("Bad" in path for path in selected))
        self.assertTrue(any("_3d" in path for path in selected))
        entries = lock["entries"]
        self.assertEqual(len(entries), len(selected))
        self.assertEqual({entry["id"] for entry in entries}, {"RDKit:" + p for p in selected})
        observed = set()
        for entry in entries:
            self.assertEqual(len(entry["files"]), 1)
            source = entry["files"][0]
            relative = source["upstream_path"].removeprefix(repository["path"] + "/")
            item = selected[relative]
            self.assertEqual(entry["id"], "RDKit:" + relative)
            self.assertEqual(source["upstream_path"], repository["path"] + "/" + relative)
            self.assertEqual(source["path"], "data/" + str(PurePosixPath(relative).with_suffix(".mol")))
            data = (root / source["path"]).read_bytes()
            if relative.endswith(".sdf"):
                self.assertEqual(data.splitlines()[-1].strip(), b"M  END")
                self.assertFalse(any(line.strip() == b"$$$$" or line.lstrip().startswith(b">")
                                     for line in data.splitlines()))
            self.assertEqual(len(data), item["size"])
            self.assertEqual(source["size"], item["size"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), source["sha256"])
            self.assertEqual(hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest(),
                             item["sha"])
            self.assertEqual(source["git_blob_sha1"], item["sha"])
            observed.add(source["path"])
        self.assertEqual(observed, {p.relative_to(root).as_posix()
                                    for p in (root / "data").rglob("*") if p.is_file()})


if __name__ == "__main__":
    unittest.main()
