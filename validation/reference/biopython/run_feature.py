#!/usr/bin/env python3
"""Generate Biopython-backed golden data for macromolecular validation features."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import warnings
from pathlib import Path
from typing import Any


SUPPORTED_FEATURES = {
    "io.mmcif.parse",
    "bio.hierarchy.smcra",
}

ATOM_SITE_FIELDS = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "auth_atom_id",
    "label_alt_id",
    "label_comp_id",
    "auth_comp_id",
    "label_asym_id",
    "auth_asym_id",
    "label_seq_id",
    "auth_seq_id",
    "pdbx_PDB_ins_code",
    "occupancy",
    "B_iso_or_equiv",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "pdbx_PDB_model_num",
]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate normalized JSON golden data with Biopython."
    )
    parser.add_argument("--feature", required=True, choices=sorted(SUPPORTED_FEATURES))
    parser.add_argument("--corpus", default="tiny")
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="Repository root. Defaults to the script's containing checkout.",
    )
    parser.add_argument(
        "--fixture",
        action="append",
        help="Fixture path from the selected corpus manifest. May be repeated.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for JSON output. Defaults to validation/corpora/<corpus>/golden/<feature>.",
    )
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help="Only check that Biopython imports and print its version.",
    )
    args = parser.parse_args()

    biopython = import_biopython()
    if args.check_deps:
        print(f"Biopython {biopython['version']}")
        return 0

    repo_root = args.repo_root.resolve()
    corpus_dir = repo_root / "validation" / "corpora" / args.corpus
    manifest_path = corpus_dir / "features" / f"{args.feature}.toml"
    manifest = read_manifest(manifest_path)
    if manifest.get("corpus_id") != args.corpus:
        raise SystemExit(
            f"{manifest_path} declares corpus_id {manifest.get('corpus_id')!r}, "
            f"expected {args.corpus!r}"
        )
    fixtures = selected_fixtures(manifest, args.fixture)
    output_dir = (args.output_dir or corpus_dir / "golden" / args.feature).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for fixture in fixtures:
        fixture_path = (corpus_dir / fixture).resolve()
        if not fixture_path.exists():
            raise SystemExit(f"{manifest_path} references missing fixture: {fixture}")
        document = generate_document(
            args.feature, args.corpus, fixture, fixture_path, biopython
        )
        output_path = output_dir / f"{slugify_fixture(fixture)}.json.gz"
        write_json(output_path, document)
        print(output_path)
    return 0


def import_biopython() -> dict[str, Any]:
    try:
        import Bio
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        from Bio.PDB.MMCIFParser import MMCIFParser
    except ImportError as error:
        raise SystemExit(
            "Biopython is not importable. Create the environment from "
            "validation/reference/biopython/environment.yml before generating goldens."
        ) from error
    return {
        "version": Bio.__version__,
        "MMCIF2Dict": MMCIF2Dict,
        "MMCIFParser": MMCIFParser,
    }


def read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"missing validation manifest: {path}")
    manifest = parse_simple_manifest(path.read_text(encoding="utf-8"))
    fixtures = manifest.get("fixtures")
    if not isinstance(fixtures, list) or not all(isinstance(item, str) for item in fixtures):
        raise SystemExit(f"{path} must define fixtures as a string array")
    return manifest


def parse_simple_manifest(text: str) -> dict[str, Any]:
    manifest: dict[str, Any] = {}
    lines = iter(text.splitlines())
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = [part.strip() for part in line.split("=", 1)]
        if value == "[":
            items: list[str] = []
            for array_line in lines:
                array_line = array_line.strip()
                if array_line == "]":
                    break
                item = array_line.rstrip(",").strip()
                if item.startswith('"') and item.endswith('"'):
                    items.append(item[1:-1])
            manifest[key] = items
        elif value.startswith('"') and value.endswith('"'):
            manifest[key] = value[1:-1]
        else:
            manifest[key] = value
    return manifest


def selected_fixtures(manifest: dict[str, Any], requested: list[str] | None) -> list[str]:
    fixtures = list(manifest["fixtures"])
    if not requested:
        return fixtures
    unknown = sorted(set(requested) - set(fixtures))
    if unknown:
        raise SystemExit(f"requested fixture(s) not present in manifest: {', '.join(unknown)}")
    return [fixture for fixture in fixtures if fixture in requested]


def generate_document(
    feature_id: str,
    corpus_id: str,
    fixture: str,
    fixture_path: Path,
    biopython: dict[str, Any],
) -> dict[str, Any]:
    atom_site = atom_site_table(fixture_path, biopython["MMCIF2Dict"])
    structure = structure_summary(fixture_path, biopython["MMCIFParser"])
    return {
        "schema_version": 1,
        "feature_id": feature_id,
        "corpus_id": corpus_id,
        "fixture_id": slugify_fixture(fixture),
        "fixture_path": fixture,
        "input_sha256": sha256_file(fixture_path),
        "reference": {
            "tool": "biopython",
            "version": biopython["version"],
            "runtime_dependency": False,
        },
        "expected": {
            "atom_site_rows": atom_site,
            "structure": structure,
        },
    }


def atom_site_table(fixture_path: Path, MMCIF2Dict: Any) -> dict[str, Any]:
    try:
        raw = MMCIF2Dict(str(fixture_path))
    except Exception as error:
        return {"status": "parse_error", "error": error.__class__.__name__, "rows": []}

    values_by_field: dict[str, list[str | None]] = {}
    row_count = 0
    for field in ATOM_SITE_FIELDS:
        raw_value = raw.get(f"_atom_site.{field}")
        values = normalize_mmcif_column(raw_value)
        if values:
            row_count = max(row_count, len(values))
        values_by_field[field] = values

    rows = []
    for index in range(row_count):
        row = {
            field: normalize_missing(values[index]) if index < len(values) else None
            for field, values in values_by_field.items()
        }
        rows.append(row)

    return {
        "status": "ok",
        "row_count": row_count,
        "rows": rows,
    }


def structure_summary(fixture_path: Path, MMCIFParser: Any) -> dict[str, Any]:
    parser = MMCIFParser(QUIET=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = parser.get_structure(fixture_path.stem, str(fixture_path))
    except Exception as error:
        return {"status": "parse_error", "error": error.__class__.__name__, "models": []}

    models = []
    for model in structure:
        chains = []
        for chain in model:
            residues = []
            for residue in chain:
                residue_name = residue.get_resname().strip()
                hetflag, sequence_id, insertion_code = residue.id
                atoms = []
                for atom in residue.get_unpacked_list():
                    atoms.append(
                        {
                            "name": atom.get_name(),
                            "full_name": atom.get_fullname().strip(),
                            "altloc": normalize_missing(atom.get_altloc().strip()),
                            "element": normalize_missing(atom.element.strip()),
                            "occupancy": atom.get_occupancy(),
                            "bfactor": atom.get_bfactor(),
                            "coord": [round(float(value), 6) for value in atom.get_coord()],
                        }
                    )
                residues.append(
                    {
                        "name": residue_name,
                        "hetflag": normalize_missing(hetflag.strip()),
                        "sequence_id": sequence_id,
                        "insertion_code": normalize_missing(str(insertion_code).strip()),
                        "atoms": atoms,
                    }
                )
            chains.append({"id": chain.id, "residues": residues})
        models.append({"id": model.id, "chains": chains})

    return {"status": "ok", "models": models}


def normalize_mmcif_column(value: Any) -> list[str | None]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def normalize_missing(value: str | None) -> str | None:
    if value is None or value in {"", ".", "?"}:
        return None
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify_fixture(fixture: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", fixture).strip("._-")


def write_json(path: Path, document: dict[str, Any]) -> None:
    payload = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode("utf-8")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(payload)


if __name__ == "__main__":
    raise SystemExit(main())
