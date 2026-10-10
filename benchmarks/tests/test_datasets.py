"""Dataset manifests, curated files and snapshot archives (stdlib only)."""
import gzip
import hashlib
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kbench import datasets  # noqa: E402
from kbench.datasets import DatasetError, Record, SourceFile, extract, pack, parse_record, query_rows  # noqa: E402


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


class QueryRows(unittest.TestCase):
    def test_label_tab_tables_skip_only_comments_and_blank_lines(self):
        text = (
            "// Label\tSMARTS\tNotes\n"
            "\n"
            "-C(=O)O\t*-C(=O)[O;D1]\tcarboxylic acids\n"
            "#N\t   *#[N;D1]\t\t\tnitriles\n"
            "  AcidChloride.Aromatic\t[$(C-!@[a])](=O)(Cl)\tAromatic\n"
        )
        self.assertEqual(
            query_rows(text, "label-tabs-query"),
            [
                (3, "*-C(=O)[O;D1]", "-C(=O)O"),
                (4, "*#[N;D1]", "#N"),
                (5, "[$(C-!@[a])](=O)(Cl)", "AcidChloride.Aromatic"),
            ],
        )

    def test_whitespace_tables_split_query_from_label(self):
        text = "# header\n#\n[Br,Cl,I][CX4;CH,CH2] Reactive_alkyl_halides\nP\n"
        self.assertEqual(
            query_rows(text, "query-whitespace-label"),
            [(3, "[Br,Cl,I][CX4;CH,CH2]", "Reactive_alkyl_halides"), (4, "P", "")],
        )


class Records(unittest.TestCase):
    def file(self, **changes):
        value = {"path": "chembl/CHEMBL25.mol", "format": "mol", "sha256": "0" * 64, "url": "https://x"}
        value.update(changes)
        return value

    def record(self, *files):
        return {"id": "chembl:CHEMBL25", "stratum": "s", "source": "chembl", "source_id": "CHEMBL25",
                "files": list(files)}

    def test_valid_record_round_trips(self):
        value = self.record(self.file(role="2d"))
        self.assertEqual(parse_record(value).to_json(), value)

    def test_rejects_malformed_files(self):
        for bad in (
            self.file(sha256="ABC"),
            self.file(format="pdb"),
            self.file(path="../escape.mol"),
            self.file(path="/abs.mol"),
            self.file(path="c:/x.mol"),
            self.file(format="smarts"),
            self.file(format="smarts", row=3, layout="unknown"),
        ):
            with self.subTest(bad=bad), self.assertRaises(DatasetError):
                parse_record(self.record(bad))


class Snapshot(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.source = self.root / "staging"
        contents = {"a/one.smi": b"CCO\n", "b/two.mol": b"\n  test\n\n  0  0\nM  END\n"}
        self.records = []
        for index, (path, data) in enumerate(sorted(contents.items())):
            (self.source / path).parent.mkdir(parents=True, exist_ok=True)
            (self.source / path).write_bytes(data)
            file = SourceFile(path, "smiles" if path.endswith(".smi") else "mol", sha(data), "https://x")
            self.records.append(Record(f"r{index}", "s", "test", str(index), (file,)))

    def tearDown(self):
        self.temporary.cleanup()

    def test_pack_is_deterministic_and_extracts_exact_files(self):
        first = pack(self.source, self.records, self.root / "one.tar.xz")
        second = pack(self.source, self.records, self.root / "two.tar.xz")
        self.assertEqual(first, second)
        self.assertEqual((self.root / "one.tar.xz").read_bytes(), (self.root / "two.tar.xz").read_bytes())
        extract(self.root / "one.tar.xz", self.root / "out", self.records)
        self.assertEqual((self.root / "out/a/one.smi").read_bytes(), b"CCO\n")
        self.assertEqual(sorted(p.name for p in (self.root / "out").rglob("*") if p.is_file()), ["one.smi", "two.mol"])

    def test_extras_travel_with_the_archive_and_are_verified(self):
        notice = b"CC0\n"
        pack(self.source, self.records, self.root / "n.tar.xz", {"NOTICES/LICENSES.md": notice})
        extract(self.root / "n.tar.xz", self.root / "out", self.records, {"NOTICES/LICENSES.md": sha(notice)})
        self.assertEqual((self.root / "out/NOTICES/LICENSES.md").read_bytes(), notice)
        with self.assertRaises(DatasetError):
            extract(self.root / "n.tar.xz", self.root / "other", self.records, {"NOTICES/LICENSES.md": sha(b"BY\n")})

    def test_pack_refuses_files_that_differ_from_the_manifest(self):
        (self.source / "a/one.smi").write_bytes(b"CCN\n")
        with self.assertRaises(DatasetError):
            pack(self.source, self.records, self.root / "bad.tar.xz")

    def archive(self, members):
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w:xz") as tar:
            for name, data in members:
                info = tarfile.TarInfo(name)
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
        path = self.root / "crafted.tar.xz"
        path.write_bytes(buffer.getvalue())
        return path

    def test_extract_rejects_tampered_extra_missing_and_unsafe_members(self):
        good = [("a/one.smi", b"CCO\n"), ("b/two.mol", b"\n  test\n\n  0  0\nM  END\n")]
        for members in (
            [("a/one.smi", b"CCN\n"), good[1]],
            good + [("c/extra.smi", b"C\n")],
            good[:1],
            good + [("../escape.smi", b"C\n")],
        ):
            with self.subTest(members=[name for name, _ in members]), self.assertRaises(DatasetError):
                extract(self.archive(members), self.root / "out", self.records)
            self.assertFalse((self.root / "out").exists())


class CommittedDatasets(unittest.TestCase):
    def test_manifests_quotas_and_curated_files_are_consistent(self):
        names = sorted(path.parent.name for path in datasets.DATASETS.glob("*/manifest.jsonl"))
        self.assertEqual(names, ["bio", "small"])
        for name in names:
            with self.subTest(dataset=name):
                dataset = datasets.Dataset(name)
                counts = {}
                for record in dataset.records.values():
                    counts[record.stratum] = counts.get(record.stratum, 0) + 1
                for stratum in dataset.selection["stratum"]:
                    if "quota" in stratum:
                        self.assertEqual(counts.get(stratum["id"]), stratum["quota"], stratum["id"])
                    if stratum.get("curated") == "all":
                        members = {r.id for r in dataset.records.values() if r.stratum == stratum["id"]}
                        self.assertTrue(members <= set(dataset.curated_ids), stratum["id"])
                dataset.files("curated")

    def test_committed_curated_references_describe_the_curated_inputs(self):
        for name in ("bio", "small"):
            with self.subTest(dataset=name):
                dataset = datasets.Dataset(name)
                data = gzip.decompress((dataset.directory / "references-curated.jsonl.gz").read_bytes())
                header = json.loads(data.split(b"\n", 1)[0])
                self.assertEqual((header["dataset"], header["version"]), (name, dataset.version))
                self.assertEqual(header["inputs"], dataset.inputs_digest("curated"))
                self.assertNotEqual(dataset.inputs_digest("curated"), dataset.inputs_digest("full"))


if __name__ == "__main__":
    unittest.main()
