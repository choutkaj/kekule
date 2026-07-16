#!/usr/bin/env python3
"""Generate Biopython-backed golden data for macromolecular validation features."""

from __future__ import annotations

import argparse
import concurrent.futures
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import Any


SUPPORTED_FEATURES = {
    "bio.secondary-structure.dssp",
    "io.mmcif.parse",
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
    parser.add_argument("--corpus", default="pdb-100")
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
    parser.add_argument(
        "--jobs",
        type=positive_int,
        default=1,
        help="Number of independent fixture generators to run concurrently.",
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

    tasks = [
        (args.feature, args.corpus, fixture, str(corpus_dir), str(output_dir))
        for fixture in fixtures
    ]
    if args.jobs == 1:
        output_paths = [generate_fixture(task, biopython) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
            output_paths = list(executor.map(generate_fixture, tasks))
    for output_path in output_paths:
        print(output_path)
    return 0


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def generate_fixture(
    task: tuple[str, str, str, str, str],
    biopython: dict[str, Any] | None = None,
) -> Path:
    feature, corpus, fixture, corpus_dir_text, output_dir_text = task
    corpus_dir = Path(corpus_dir_text)
    output_dir = Path(output_dir_text)
    fixture_path = (corpus_dir / fixture).resolve()
    manifest_path = corpus_dir / "features" / f"{feature}.toml"
    if not fixture_path.exists():
        raise SystemExit(f"{manifest_path} references missing fixture: {fixture}")
    document = generate_document(
        feature,
        corpus,
        fixture,
        fixture_path,
        biopython or import_biopython(),
    )
    output_path = output_dir / f"{slugify_fixture(fixture)}.json.gz"
    write_json(output_path, document)
    return output_path


def import_biopython() -> dict[str, Any]:
    try:
        import Bio
        from Bio.PDB.DSSP import DSSP
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        from Bio.PDB.mmcifio import MMCIFIO
        from Bio.PDB.MMCIFParser import MMCIFParser
    except ImportError as error:
        raise SystemExit(
            "Biopython is not importable. Create the environment from "
            "validation/reference/biopython/environment.yml before generating goldens."
        ) from error
    return {
        "version": Bio.__version__,
        "DSSP": DSSP,
        "MMCIF2Dict": MMCIF2Dict,
        "MMCIFIO": MMCIFIO,
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
    if feature_id == "bio.secondary-structure.dssp":
        expected = dssp_summary(
            fixture_path,
            biopython["MMCIFParser"],
            biopython["MMCIF2Dict"],
            biopython["MMCIFIO"],
            biopython["DSSP"],
        )
        reference = dssp_reference(biopython["version"])
        return {
            "schema_version": 1,
            "feature_id": feature_id,
            "corpus_id": corpus_id,
            "fixture_id": slugify_fixture(fixture),
            "fixture_path": fixture,
            "input_sha256": sha256_file(fixture_path),
            "reference": reference,
            "expected": expected,
        }
    atom_site = atom_site_table(fixture_path, biopython["MMCIF2Dict"])
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
        "expected": {"atom_site_rows": atom_site},
    }


def dssp_reference(biopython_version: str) -> dict[str, Any]:
    executable = shutil.which("mkdssp")
    if executable is None:
        raise SystemExit(
            "mkdssp is not available. Recreate the environment from "
            "validation/reference/biopython/environment.yml."
        )
    version = subprocess.run(
        [executable, "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    return {
        "tool": "biopython",
        "version": f"Biopython {biopython_version} / {version}",
        "biopython_version": biopython_version,
        "dssp_version": version,
        "dssp_executable_sha256": sha256_file(Path(executable)),
        "command": (
            "Bio.PDB.DSSP.DSSP(model, highest_occupancy_snapshot, "
            "dssp='mkdssp', file_type='MMCIF')"
        ),
        "extended_command": (
            "mkdssp --output-format=mmcif --quiet "
            "highest_occupancy_snapshot annotated.cif"
        ),
        "runtime_dependency": False,
    }


def dssp_summary(
    fixture_path: Path,
    MMCIFParser: Any,
    MMCIF2Dict: Any,
    MMCIFIO: Any,
    DSSP: Any,
) -> dict[str, Any]:
    parser = MMCIFParser(QUIET=True)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = parser.get_structure(fixture_path.stem, str(fixture_path))[0]
            with tempfile.TemporaryDirectory(prefix="molecules-dssp-") as temp_dir:
                selected_path = Path(temp_dir) / fixture_path.name
                write_highest_occupancy_snapshot(
                    fixture_path, selected_path, MMCIF2Dict, MMCIFIO
                )
                assignments = DSSP(
                    model,
                    str(selected_path),
                    dssp="mkdssp",
                    file_type="MMCIF",
                )
                annotated_path = Path(temp_dir) / "annotated.cif"
                subprocess.run(
                    [
                        "mkdssp",
                        "--output-format=mmcif",
                        "--quiet",
                        str(selected_path),
                        str(annotated_path),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                extended_rows = dssp_extended_rows(MMCIF2Dict(str(annotated_path)))
    except Exception as error:
        return {
            "status": "reference_error",
            "error": error.__class__.__name__,
            "message": str(error),
            "residues": [],
        }
    if len(assignments) == 0:
        return {"status": "no_analyzable_residues", "residues": []}

    keys = list(assignments.keys())
    if len(keys) != len(extended_rows):
        raise RuntimeError(
            "Biopython legacy DSSP and DSSP4 mmCIF output contain different "
            f"residue counts: {len(keys)} versus {len(extended_rows)}"
        )
    keys_by_dssp_index = {int(assignments[key][0]): key for key in keys}
    residues = []
    previous_dssp_index = None
    previous_label_chain = None
    for key, extended in zip(keys, extended_rows, strict=True):
        chain_id, residue_id = key
        _, sequence_id, insertion_code = residue_id
        value = assignments[key]
        dssp_index = int(value[0])
        if previous_dssp_index is None:
            chain_break = "new_chain"
        elif dssp_index != previous_dssp_index + 1:
            chain_break = (
                "new_chain"
                if extended["label_chain_id"] != previous_label_chain
                else "gap"
            )
        else:
            chain_break = "none"
        residues.append(
            {
                "chain_id": chain_id,
                "sequence_id": sequence_id,
                "insertion_code": normalize_missing(str(insertion_code).strip()),
                "label_chain_id": extended["label_chain_id"],
                "label_sequence_id": extended["label_sequence_id"],
                "residue_name": extended["residue_name"],
                "residue_one_letter": value[1],
                "secondary_structure": " " if value[2] == "-" else value[2],
                "chain_break": chain_break,
                "phi_degrees": dssp_optional_angle(value[4]),
                "psi_degrees": dssp_optional_angle(value[5]),
                "tco": extended["tco"],
                "kappa_degrees": extended["kappa_degrees"],
                "alpha_degrees": extended["alpha_degrees"],
                "helix_positions": extended["helix_positions"],
                "sheet": extended["sheet"],
                "strand": extended["strand"],
                "ladders": extended["ladders"],
                "beta_parallel": extended["beta_parallel"],
                "acceptors": [
                    dssp_bond(value[6], value[7], dssp_index, keys_by_dssp_index),
                    dssp_bond(value[10], value[11], dssp_index, keys_by_dssp_index),
                ],
                "donors": [
                    dssp_bond(value[8], value[9], dssp_index, keys_by_dssp_index),
                    dssp_bond(value[12], value[13], dssp_index, keys_by_dssp_index),
                ],
            }
        )
        previous_dssp_index = dssp_index
        previous_label_chain = extended["label_chain_id"]
    return {"status": "ok", "residues": residues}


def dssp_extended_rows(document: dict[str, Any]) -> list[dict[str, Any]]:
    prefix = "_dssp_struct_summary."

    def column(name: str) -> list[str | None]:
        return normalize_mmcif_column(document.get(f"{prefix}{name}"))

    fields = (
        "label_asym_id",
        "label_seq_id",
        "label_comp_id",
        "helix_3_10",
        "helix_alpha",
        "helix_pi",
        "helix_pp",
        "sheet",
        "strand",
        "ladder_1",
        "ladder_2",
        "TCO",
        "kappa",
        "alpha",
    )
    columns = {field: column(field) for field in fields}
    row_count = len(columns["label_asym_id"])

    ladder_ids = normalize_mmcif_column(document.get("_dssp_struct_ladder.id"))
    ladder_types = normalize_mmcif_column(document.get("_dssp_struct_ladder.type"))
    type_by_ladder = {
        ladder_id: ladder_type
        for ladder_id, ladder_type in zip(ladder_ids, ladder_types, strict=True)
    }

    rows = []
    for row in range(row_count):
        ladder_tokens = [
            normalize_missing(columns["ladder_1"][row]),
            normalize_missing(columns["ladder_2"][row]),
        ]
        rows.append(
            {
                "label_chain_id": columns["label_asym_id"][row],
                "label_sequence_id": int(columns["label_seq_id"][row]),
                "residue_name": columns["label_comp_id"][row],
                "tco": dssp_optional_number(columns["TCO"][row]),
                "kappa_degrees": dssp_optional_number(columns["kappa"][row]),
                "alpha_degrees": dssp_optional_number(columns["alpha"][row]),
                "helix_positions": [
                    dssp_helix_position(columns["helix_3_10"][row]),
                    dssp_helix_position(columns["helix_alpha"][row]),
                    dssp_helix_position(columns["helix_pi"][row]),
                    dssp_helix_position(columns["helix_pp"][row]),
                ],
                "sheet": dssp_identifier(columns["sheet"][row], one_based=True),
                "strand": dssp_identifier(columns["strand"][row], one_based=True),
                "ladders": [
                    dssp_identifier(token, one_based=False) for token in ladder_tokens
                ],
                "beta_parallel": [
                    None
                    if token is None
                    else type_by_ladder[token].lower() == "parallel"
                    for token in ladder_tokens
                ],
            }
        )
    return rows


def dssp_optional_number(value: Any) -> float | None:
    normalized = normalize_missing(None if value is None else str(value))
    return None if normalized is None else float(normalized)


def dssp_helix_position(value: Any) -> str:
    normalized = normalize_missing(None if value is None else str(value))
    if normalized is None:
        return "none"
    return {">": "start", "<": "end", "X": "start_and_end"}.get(
        normalized, "middle"
    )


def dssp_identifier(value: Any, *, one_based: bool) -> int | None:
    normalized = normalize_missing(None if value is None else str(value))
    if normalized is None:
        return None
    if normalized.isdigit():
        identifier = int(normalized)
    else:
        identifier = 0
        for character in normalized.upper():
            if not "A" <= character <= "Z":
                raise ValueError(f"unsupported DSSP identifier {normalized!r}")
            identifier = identifier * 26 + ord(character) - ord("A") + 1
        identifier -= 1
    return identifier + 1 if one_based else identifier


def write_highest_occupancy_snapshot(
    input_path: Path, output_path: Path, MMCIF2Dict: Any, MMCIFIO: Any
) -> None:
    """Write the same explicit altloc snapshot used by mmcif::interpret.

    DSSP itself keeps the last atom-site row for duplicate backbone names. The
    feature consumes an already-selected Model instead, so reference goldens
    must run mkdssp on that same coordinate choice. Highest occupancy wins;
    ties prefer a missing altloc and then the lexicographically first label.
    """
    document = MMCIF2Dict(str(input_path))
    # DSSP validates the whole input dictionary even though this descriptive
    # archive-link category is unrelated to coordinates or polymer topology.
    # Some otherwise usable PDB entries contain duplicate keys there.
    for key in [key for key in document if key.startswith("_pdbx_database_related.")]:
        del document[key]
    atom_columns = {
        key: normalize_mmcif_column(value)
        for key, value in document.items()
        if key.startswith("_atom_site.")
    }
    row_count = len(atom_columns.get("_atom_site.id", []))
    if row_count == 0:
        raise RuntimeError(f"{input_path} has no atom_site rows")

    identity_fields = (
        "group_PDB",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "label_comp_id",
        "auth_asym_id",
        "auth_seq_id",
        "auth_comp_id",
        "label_atom_id",
        "auth_atom_id",
        "pdbx_PDB_model_num",
    )

    def field(row: int, name: str) -> str | None:
        values = atom_columns.get(f"_atom_site.{name}", [])
        return values[row] if row < len(values) else None

    def rank(row: int) -> tuple[float, int, str]:
        occupancy = field(row, "occupancy")
        try:
            occupancy_value = float(occupancy) if occupancy not in (None, ".", "?") else 0.0
        except ValueError:
            occupancy_value = 0.0
        alt_id = normalize_missing(field(row, "label_alt_id"))
        return (-occupancy_value, 0 if alt_id is None else 1, alt_id or "")

    selected_by_identity: dict[tuple[str | None, ...], int] = {}
    for row in range(row_count):
        identity = tuple(field(row, name) for name in identity_fields)
        current = selected_by_identity.get(identity)
        if current is None or rank(row) < rank(current):
            selected_by_identity[identity] = row
    selected_rows = sorted(selected_by_identity.values())

    for key, values in atom_columns.items():
        document[key] = [values[row] for row in selected_rows]
    writer = MMCIFIO()
    writer.set_dict(document)
    writer.save(str(output_path))


def dssp_optional_angle(value: Any) -> float | None:
    angle = float(value)
    return None if angle == 360.0 else angle


def dssp_bond(
    relative_index: Any,
    energy: Any,
    dssp_index: int,
    keys_by_dssp_index: dict[int, Any],
) -> dict[str, Any] | None:
    relative_index = int(relative_index)
    energy = float(energy)
    if relative_index == 0 and energy == 0.0:
        return None
    partner_key = keys_by_dssp_index.get(dssp_index + relative_index)
    if partner_key is None:
        raise RuntimeError(
            f"DSSP bond from index {dssp_index} points to missing relative index "
            f"{relative_index}"
        )
    partner_chain, partner_residue = partner_key
    _, partner_sequence, partner_insertion = partner_residue
    return {
        "partner_chain_id": partner_chain,
        "partner_sequence_id": partner_sequence,
        "partner_insertion_code": normalize_missing(str(partner_insertion).strip()),
        "energy_kcal_per_mol": energy,
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
