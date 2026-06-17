#!/usr/bin/env python3
"""Generate RDKit-backed golden data for small-molecule validation features."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path
from typing import Any


SUPPORTED_FEATURES = {
    "io.sdf.v2000.parse",
    "algo.rings.fast",
    "algo.aromaticity.rdkit-like-basic",
}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate normalized JSON golden data with RDKit."
    )
    parser.add_argument("--feature", required=True, choices=sorted(SUPPORTED_FEATURES))
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[3],
        help="Repository root. Defaults to the script's containing checkout.",
    )
    parser.add_argument(
        "--fixture",
        action="append",
        help="Fixture path from validation.toml to generate. May be repeated.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for JSON output. Defaults to validation/features/<feature>/golden.",
    )
    parser.add_argument(
        "--check-deps",
        action="store_true",
        help="Only check that RDKit imports and print its version.",
    )
    args = parser.parse_args()

    rdkit = import_rdkit()
    if args.check_deps:
        print(f"RDKit {rdkit['version']}")
        return 0

    repo_root = args.repo_root.resolve()
    feature_dir = repo_root / "validation" / "features" / args.feature
    manifest_path = feature_dir / "validation.toml"
    manifest = read_manifest(manifest_path)
    fixtures = selected_fixtures(manifest, args.fixture)
    output_dir = (args.output_dir or feature_dir / "golden").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    for fixture in fixtures:
        fixture_path = (feature_dir / fixture).resolve()
        if not fixture_path.exists():
            raise SystemExit(f"{manifest_path} references missing fixture: {fixture}")
        document = generate_document(args.feature, fixture, fixture_path, rdkit)
        output_path = output_dir / f"{slugify_fixture(fixture)}.json"
        write_json(output_path, document)
        print(output_path)
    return 0


def import_rdkit() -> dict[str, Any]:
    try:
        from rdkit import Chem, rdBase
    except ImportError as error:
        raise SystemExit(
            "RDKit is not importable. Create the environment from "
            "validation/reference/rdkit/environment.yml before generating goldens."
        ) from error
    return {"Chem": Chem, "version": rdBase.rdkitVersion}


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
    fixture: str,
    fixture_path: Path,
    rdkit: dict[str, Any],
) -> dict[str, Any]:
    records = read_sdf_records(fixture_path, rdkit["Chem"])
    if feature_id == "io.sdf.v2000.parse":
        expected = {"records": [sdf_record(record) for record in records]}
    elif feature_id == "algo.rings.fast":
        expected = {"records": [ring_record(record) for record in records]}
    elif feature_id == "algo.aromaticity.rdkit-like-basic":
        expected = {"records": [aromaticity_record(record) for record in records]}
    else:
        raise SystemExit(f"unsupported feature for RDKit generator: {feature_id}")

    return {
        "schema_version": 1,
        "feature_id": feature_id,
        "fixture_id": slugify_fixture(fixture),
        "fixture_path": fixture,
        "input_sha256": sha256_file(fixture_path),
        "reference": {
            "tool": "rdkit",
            "version": rdkit["version"],
            "runtime_dependency": False,
        },
        "expected": expected,
    }


def read_sdf_records(fixture_path: Path, Chem: Any) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    supplier = Chem.SDMolSupplier(
        str(fixture_path),
        sanitize=False,
        removeHs=False,
        strictParsing=False,
    )
    for index, mol in enumerate(supplier):
        if mol is None:
            records.append(
                {
                    "record_index": index,
                    "status": "parse_error",
                    "title": None,
                    "mol": None,
                }
            )
            continue
        records.append(
            {
                "record_index": index,
                "status": "ok",
                "title": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
                "mol": mol,
            }
        )
    return records


def sdf_record(record: dict[str, Any]) -> dict[str, Any]:
    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": [atom_json(atom) for atom in mol.GetAtoms()],
        "bonds": [bond_json(bond) for bond in mol.GetBonds()],
        "properties": molecule_properties(mol),
    }


def ring_record(record: dict[str, Any]) -> dict[str, Any]:
    from rdkit import Chem

    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    Chem.GetSymmSSSR(mol)
    rings = mol.GetRingInfo()
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_in_ring": [atom.IsInRing() for atom in mol.GetAtoms()],
        "bond_in_ring": [rings.NumBondRings(bond.GetIdx()) > 0 for bond in mol.GetBonds()],
    }


def aromaticity_record(record: dict[str, Any]) -> dict[str, Any]:
    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    sanitized = clone_and_sanitize(mol)
    if sanitized is None:
        return {
            "record_index": record["record_index"],
            "status": "sanitize_error",
            "title": record["title"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_aromatic": [atom.GetIsAromatic() for atom in sanitized.GetAtoms()],
        "bond_aromatic": [bond.GetIsAromatic() for bond in sanitized.GetBonds()],
    }


def clone_and_sanitize(mol: Any) -> Any | None:
    from rdkit import Chem

    cloned = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(cloned)
    except Exception:
        return None
    return cloned


def atom_json(atom: Any) -> dict[str, Any]:
    return {
        "index": atom.GetIdx(),
        "atomic_number": atom.GetAtomicNum(),
        "symbol": atom.GetSymbol(),
        "formal_charge": atom.GetFormalCharge(),
        "isotope": atom.GetIsotope() or None,
        "explicit_hydrogens": atom.GetNumExplicitHs(),
        "atom_map": atom.GetAtomMapNum() or None,
        "aromatic": atom.GetIsAromatic(),
    }


def bond_json(bond: Any) -> dict[str, Any]:
    return {
        "index": bond.GetIdx(),
        "begin_atom_index": bond.GetBeginAtomIdx(),
        "end_atom_index": bond.GetEndAtomIdx(),
        "bond_type": str(bond.GetBondType()),
        "is_aromatic": bond.GetIsAromatic(),
        "stereo": str(bond.GetStereo()),
    }


def molecule_properties(mol: Any) -> dict[str, str]:
    props: dict[str, str] = {}
    for name in sorted(mol.GetPropNames(includePrivate=False, includeComputed=False)):
        props[name] = mol.GetProp(name)
    return props


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify_fixture(fixture: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", fixture).strip("._-")


def write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
