"""Dataset manifests, snapshot archives and the local data cache.

A dataset is a directory under `benchmarks/datasets/` holding:

- `selection.toml`: hand-written strata, quotas, seeds and upstream queries;
- `selection.lock.json`: what selection resolved (upstream releases, pool
  sizes and hashes) and the published snapshot archive's URL and SHA-256;
- `manifest.jsonl`: one record per line, each listing its files with
  upstream URL and SHA-256;
- `curated.txt` and `curated/`: the in-git tier, byte-identical to the
  corresponding snapshot members.

The full tier is read from a hash-locked snapshot archive extracted into the
cache. Nothing here contacts an upstream database; only `select/` does.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import re
import shutil
import tarfile
import tomllib
import urllib.request
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "datasets"
MANIFEST_SCHEMA = 1
FORMATS = {"smiles", "mol", "sdf", "mmcif", "smarts"}
SHA256 = re.compile(r"[0-9a-f]{64}")
QUERY_LAYOUTS = {"label-tabs-query", "query-whitespace-label"}


def cache_root() -> Path:
    """Cache directory; `KEKULE_BENCH_CACHE` survives `cargo clean`."""
    override = os.environ.get("KEKULE_BENCH_CACHE")
    return Path(override) if override else ROOT.parent / "target" / "kekule-bench"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class SourceFile:
    """One input file of a record, as fetched from upstream."""

    path: str
    format: str
    sha256: str
    url: str
    role: str | None = None
    row: int | None = None
    layout: str | None = None

    @property
    def gzip(self) -> bool:
        return self.path.endswith(".gz")

    def to_json(self) -> dict:
        value = {"path": self.path, "format": self.format, "sha256": self.sha256, "url": self.url}
        if self.role is not None:
            value["role"] = self.role
        if self.row is not None:
            value["row"] = self.row
        if self.layout is not None:
            value["layout"] = self.layout
        return value


@dataclass(frozen=True)
class Record:
    """One molecule, structure or query, possibly in several representations."""

    id: str
    stratum: str
    source: str
    source_id: str
    files: tuple[SourceFile, ...]

    def to_json(self) -> dict:
        return {
            "id": self.id,
            "stratum": self.stratum,
            "source": self.source,
            "source_id": self.source_id,
            "files": [file.to_json() for file in self.files],
        }


class DatasetError(ValueError):
    """A manifest, selection or file that violates the dataset contract."""


def parse_record(value: dict) -> Record:
    try:
        files = tuple(
            SourceFile(
                path=file["path"],
                format=file["format"],
                sha256=file["sha256"],
                url=file["url"],
                role=file.get("role"),
                row=file.get("row"),
                layout=file.get("layout"),
            )
            for file in value["files"]
        )
        record = Record(value["id"], value["stratum"], value["source"], value["source_id"], files)
    except (KeyError, TypeError) as error:
        raise DatasetError(f"malformed record {value.get('id', '?')}: {error!r}") from None
    for file in record.files:
        if file.format not in FORMATS:
            raise DatasetError(f"{record.id}: unknown format {file.format!r}")
        if not SHA256.fullmatch(file.sha256):
            raise DatasetError(f"{record.id}: invalid sha256 for {file.path}")
        if not safe_member(file.path):
            raise DatasetError(f"{record.id}: unsafe path {file.path!r}")
        if file.format == "smarts" and (file.row is None or file.layout not in QUERY_LAYOUTS):
            raise DatasetError(f"{record.id}: SMARTS table rows need a row and a known layout")
    if not record.files:
        raise DatasetError(f"{record.id}: no files")
    return record


def safe_member(path: str) -> bool:
    parts = path.split("/")
    return bool(path) and not path.startswith("/") and ":" not in parts[0] and all(
        part not in ("", ".", "..") for part in parts
    )


class Dataset:
    """A dataset directory with its selection, manifest and curated tier."""

    def __init__(self, name: str, directory: Path | None = None):
        self.name = name
        self.directory = directory or DATASETS / name
        with open(self.directory / "selection.toml", "rb") as handle:
            self.selection = tomllib.load(handle)
        self.lock = json.loads((self.directory / "selection.lock.json").read_text(encoding="utf-8"))
        self.records = self._load_manifest()
        self.curated_ids = self._load_curated()

    @property
    def version(self) -> int:
        return int(self.selection["dataset"]["version"])

    def _load_manifest(self) -> dict[str, Record]:
        records: dict[str, Record] = {}
        paths: set[str] = set()
        with open(self.directory / "manifest.jsonl", encoding="utf-8") as handle:
            header = json.loads(handle.readline())
            if header != {"schema": MANIFEST_SCHEMA, "dataset": self.name}:
                raise DatasetError(f"unexpected manifest header {header}")
            for line in handle:
                record = parse_record(json.loads(line))
                if record.id in records:
                    raise DatasetError(f"duplicate record {record.id}")
                for file in record.files:
                    if file.path in paths and file.format != "smarts":
                        raise DatasetError(f"{record.id}: duplicate path {file.path}")
                    paths.add(file.path)
                records[record.id] = record
        strata = {stratum["id"] for stratum in self.selection["stratum"]}
        unknown = {record.stratum for record in records.values()} - strata
        if unknown:
            raise DatasetError(f"records in undeclared strata: {sorted(unknown)}")
        return records

    def _load_curated(self) -> list[str]:
        ids = [
            line.strip()
            for line in (self.directory / "curated.txt").read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        ]
        missing = [record_id for record_id in ids if record_id not in self.records]
        if missing:
            raise DatasetError(f"curated ids not in the manifest: {missing[:5]}")
        return ids

    def tier(self, tier: str) -> list[Record]:
        if tier == "curated":
            return [self.records[record_id] for record_id in self.curated_ids]
        if tier == "full":
            return list(self.records.values())
        raise DatasetError(f"unknown tier {tier!r}")

    def files(self, tier: str) -> Path:
        """Directory holding the tier's files, verified against the manifest."""
        if tier == "curated":
            directory = self.directory / "curated"
        else:
            directory = self.extracted()
        for record in self.tier(tier):
            for file in record.files:
                verify_file(directory / file.path, file.sha256)
        return directory

    def snapshot(self) -> dict:
        snapshot = self.lock.get("snapshot")
        if not snapshot:
            raise DatasetError(f"{self.name} has no published snapshot yet")
        return snapshot

    def extracted(self) -> Path:
        snapshot = self.snapshot()
        directory = cache_root() / "data" / self.name / f"v{self.version}"
        marker = directory / ".snapshot-sha256"
        if marker.exists() and marker.read_text().strip() == snapshot["sha256"]:
            return directory
        archive = self.download()
        if directory.exists():
            shutil.rmtree(directory)
        extract(archive, directory, self.records.values(), snapshot.get("extras", {}))
        marker.write_text(snapshot["sha256"] + "\n")
        return directory

    def download(self) -> Path:
        snapshot = self.snapshot()
        target = cache_root() / "archives" / snapshot["name"]
        if target.exists() and sha256_file(target) == snapshot["sha256"]:
            return target
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_suffix(target.suffix + ".part")
        request = urllib.request.Request(snapshot["url"], headers={"User-Agent": "kekule-bench"})
        with urllib.request.urlopen(request, timeout=600) as response, open(partial, "wb") as out:
            shutil.copyfileobj(response, out, 1 << 20)
        actual = sha256_file(partial)
        if actual != snapshot["sha256"]:
            partial.unlink()
            raise DatasetError(f"{snapshot['name']}: sha256 {actual} != {snapshot['sha256']}")
        partial.replace(target)
        return target


def verify_file(path: Path, expected: str) -> None:
    if not path.is_file():
        raise DatasetError(f"missing {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise DatasetError(f"{path}: sha256 {actual} != {expected}")




def member_paths(records) -> list[str]:
    return sorted({file.path for record in records for file in record.files})


def pack(source: Path, records, archive: Path, extras: dict[str, bytes] | None = None) -> str:
    """Write a deterministic `.tar.xz` of every manifest file; return its sha256.

    `extras` are further members by name, such as licence notices and the
    stored reference observations.
    """
    archive.parent.mkdir(parents=True, exist_ok=True)
    expected = {file.path: file.sha256 for record in records for file in record.files}
    members = []
    for path in member_paths(records):
        data = (source / path).read_bytes()
        if sha256_bytes(data) != expected[path]:
            raise DatasetError(f"{path}: content does not match the manifest")
        members.append((path, data))
    members += sorted((extras or {}).items())
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w:xz", preset=9, format=tarfile.PAX_FORMAT) as tar:
        for path, data in members:
            info = tarfile.TarInfo(path)
            info.size = len(data)
            info.mtime = 0
            info.mode = 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    archive.write_bytes(buffer.getvalue())
    return sha256_bytes(buffer.getvalue())


def extract(archive: Path, directory: Path, records, extras: dict[str, str] | None = None) -> None:
    """Extract exactly the manifest's files and the extras, verifying every member."""
    expected = {file.path: file.sha256 for record in records for file in record.files}
    expected.update(extras or {})
    seen = set()
    staging = directory.with_name(directory.name + ".partial")
    if staging.exists():
        shutil.rmtree(staging)
    with tarfile.open(archive, mode="r:xz") as tar:
        for member in tar:
            if not member.isfile() or not safe_member(member.name):
                raise DatasetError(f"unexpected archive member {member.name!r}")
            if member.name not in expected or member.name in seen:
                raise DatasetError(f"archive member not in the manifest: {member.name}")
            data = tar.extractfile(member).read()
            if sha256_bytes(data) != expected[member.name]:
                raise DatasetError(f"{member.name}: sha256 mismatch inside the archive")
            target = staging / member.name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
            seen.add(member.name)
    missing = set(expected) - seen
    if missing:
        shutil.rmtree(staging)
        raise DatasetError(f"archive lacks {len(missing)} manifest files, e.g. {sorted(missing)[:3]}")
    staging.replace(directory)


def query_rows(text: str, layout: str) -> list[tuple[int, str, str]]:
    """Every query row of an upstream SMARTS table as (line, SMARTS, label).

    Lines are numbered from 1. Only blank and comment lines are skipped;
    surrounding field whitespace is removed and the query is otherwise verbatim.
    """
    rows = []
    for number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if layout == "label-tabs-query":
            if not stripped or stripped.startswith("//"):
                continue
            fields = [field.strip() for field in line.split("	") if field.strip()]
            label, query = fields[0], fields[1]
        elif layout == "query-whitespace-label":
            if not stripped or stripped.startswith("#"):
                continue
            query, _, label = stripped.partition(" ")
            label = label.strip()
        else:
            raise DatasetError(f"unknown query layout {layout!r}")
        rows.append((number, query, label))
    return rows


def read_query(directory: Path, file: SourceFile) -> tuple[str, str]:
    """The (SMARTS, label) of one query record's table row."""
    text = (directory / file.path).read_text(encoding="utf-8")
    for number, query, label in query_rows(text, file.layout):
        if number == file.row:
            return query, label
    raise DatasetError(f"{file.path}: line {file.row} is not a query row")
