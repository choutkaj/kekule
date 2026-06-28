#!/usr/bin/env python3
"""Generate RDKit-backed golden data for small-molecule validation features."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path
from typing import Any


SUPPORTED_FEATURES = {
    "algo.aromaticity.rdkit-like",
    "algo.rings.fast",
    "algo.rings.sssr",
    "algo.valence.rdkit-like",
    "chem.sanitize.rdkit-like",
    "core.conformers",
    "io.mol.v2000.parse",
    "io.mol.v2000.write",
    "io.mol.v3000.parse",
    "io.mol.v3000.write",
    "io.sdf.v2000.parse",
    "io.sdf.v2000.write",
    "io.smiles.parse",
    "io.smiles.write",
    "io.smiles.canonical",
}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate normalized JSON golden data with RDKit."
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
        help="Only check that RDKit imports and print its version.",
    )
    args = parser.parse_args()

    rdkit = import_rdkit()
    if args.check_deps:
        print(f"RDKit {rdkit['version']}")
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
        document = generate_document(args.feature, args.corpus, fixture, fixture_path, rdkit)
        output_path = output_dir / f"{slugify_fixture(fixture)}.json.gz"
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
    corpus_id: str,
    fixture: str,
    fixture_path: Path,
    rdkit: dict[str, Any],
) -> dict[str, Any]:
    if feature_id == "io.sdf.v2000.parse":
        records = read_sdf_records(fixture_path, rdkit["Chem"])
        expected = {"records": [sdf_record(record) for record in records]}
    elif feature_id == "io.sdf.v2000.write":
        records = read_records_by_suffix(fixture_path, rdkit["Chem"])
        expected = {"records": [sdf_record(record) for record in records]}
    elif feature_id == "io.mol.v2000.write":
        records = read_records_by_suffix(fixture_path, rdkit["Chem"])
        expected = {"records": [mol_record(record) for record in records]}
    elif feature_id == "io.mol.v2000.parse":
        records = read_records_by_suffix(fixture_path, rdkit["Chem"])
        expected = {"records": [mol_parse_record(record) for record in records]}
    elif feature_id == "io.mol.v3000.write":
        records = read_records_by_suffix(fixture_path, rdkit["Chem"])
        expected = {"records": [mol_record(record) for record in records]}
    elif feature_id == "io.mol.v3000.parse":
        records = read_records_by_suffix(fixture_path, rdkit["Chem"])
        expected = {"records": [mol_parse_record(record) for record in records]}
    elif feature_id == "core.conformers":
        records = read_records_by_suffix(fixture_path, rdkit["Chem"])
        expected = {"records": [conformer_record(record) for record in records]}
    elif feature_id == "io.smiles.parse":
        records = read_smiles_records(fixture_path, rdkit["Chem"], sanitize=False)
        expected = {"records": [smiles_parse_record(record) for record in records]}
    elif feature_id == "io.smiles.write":
        records = read_smiles_records(fixture_path, rdkit["Chem"], sanitize=False)
        expected = {"records": [smiles_write_record(record) for record in records]}
    elif feature_id == "io.smiles.canonical":
        records = read_canonical_smiles_records(fixture_path, rdkit["Chem"], sanitize=True)
        expected = {
            "records": [
                canonical_smiles_record(record, exact_smiles=corpus_id == "tiny")
                for record in records
            ]
        }
    elif feature_id == "algo.rings.fast":
        records = read_sdf_records(fixture_path, rdkit["Chem"])
        expected = {"records": [ring_record(record) for record in records]}
    elif feature_id == "algo.rings.sssr":
        records = read_sdf_records(fixture_path, rdkit["Chem"])
        expected = {"records": [ring_set_record(record) for record in records]}
    elif feature_id == "algo.valence.rdkit-like":
        records = read_sdf_records(fixture_path, rdkit["Chem"])
        expected = {"records": [valence_record(record) for record in records]}
    elif feature_id == "chem.sanitize.rdkit-like":
        records = read_sdf_records(fixture_path, rdkit["Chem"])
        expected = {"records": [sanitized_atom_record(record) for record in records]}
    elif feature_id == "algo.aromaticity.rdkit-like":
        records = read_sdf_records(fixture_path, rdkit["Chem"])
        expected = {"records": [aromaticity_record(record) for record in records]}
    else:
        raise SystemExit(f"unsupported feature for RDKit generator: {feature_id}")

    return {
        "schema_version": 1,
        "feature_id": feature_id,
        "corpus_id": corpus_id,
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
    raw_records = read_sdf_blocks(fixture_path)
    supplier = Chem.SDMolSupplier(
        str(fixture_path),
        sanitize=False,
        removeHs=False,
        strictParsing=False,
    )
    for index, mol in enumerate(supplier):
        raw_block = raw_records[index] if index < len(raw_records) else ""
        if mol is None:
            records.append(
                {
                    "record_index": index,
                    "status": "parse_error",
                    "title": None,
                    "mol": None,
                    "radicals": parse_mdl_radicals(raw_block),
                    "bond_stereo": parse_mdl_bond_stereo(raw_block),
                }
            )
            continue
        records.append(
            {
                "record_index": index,
                "status": "ok",
                "title": mol.GetProp("_Name") if mol.HasProp("_Name") else "",
                "mol": mol,
                "radicals": parse_mdl_radicals(raw_block),
                "bond_stereo": parse_mdl_bond_stereo(raw_block),
            }
        )
    return records


def read_records_by_suffix(fixture_path: Path, Chem: Any) -> list[dict[str, Any]]:
    if fixture_path.suffix.lower() in {".mol", ".mdl"}:
        raw_block = fixture_path.read_text(encoding="utf-8", errors="replace")
        mol = Chem.MolFromMolFile(
            str(fixture_path),
            sanitize=False,
            removeHs=False,
            strictParsing=False,
        )
        return [
            {
                "record_index": 0,
                "status": "ok" if mol is not None else "parse_error",
                "title": mol.GetProp("_Name") if mol is not None and mol.HasProp("_Name") else "",
                "mol": mol,
                "radicals": parse_mdl_radicals(raw_block),
                "bond_stereo": parse_mdl_bond_stereo(raw_block),
            }
        ]
    return read_sdf_records(fixture_path, Chem)


def read_smiles_records(fixture_path: Path, Chem: Any, sanitize: bool) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index, raw_line in enumerate(fixture_path.read_text(encoding="utf-8").splitlines()):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        smiles = parts[0]
        title = parts[1] if len(parts) > 1 else ""
        unsupported = smiles_unsupported_subset_reason(smiles)
        if unsupported is not None:
            records.append(
                {
                    "record_index": index,
                    "status": "unsupported",
                    "title": title,
                    "smiles": smiles,
                    "mol": None,
                    "radicals": {},
                    "bond_stereo": {},
                }
            )
            continue
        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
        records.append(
            {
                "record_index": index,
                "status": "ok" if mol is not None else "parse_error",
                "title": title,
                "smiles": smiles,
                "mol": mol,
                "radicals": {},
                "bond_stereo": {},
            }
        )
    return records


def read_canonical_smiles_records(
    fixture_path: Path, Chem: Any, sanitize: bool
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index, raw_line in enumerate(fixture_path.read_text(encoding="utf-8").splitlines()):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        smiles = parts[0]
        title = parts[1].strip() if len(parts) > 1 else ""
        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
        records.append(
            {
                "record_index": index,
                "status": "ok" if mol is not None else "parse_error",
                "title": title,
                "smiles": smiles,
                "mol": mol,
                "radicals": {},
                "bond_stereo": {},
            }
        )
    return records


def smiles_unsupported_subset_reason(smiles: str) -> str | None:
    if any(ch in smiles for ch in ("@", "/", "\\", "*")):
        return "unsupported"
    return None


def read_sdf_blocks(fixture_path: Path) -> list[str]:
    text = fixture_path.read_text(encoding="utf-8", errors="replace")
    blocks: list[str] = []
    current: list[str] = []
    for line in text.splitlines():
        if line == "$$$$":
            blocks.append("\n".join(current) + "\n")
            current = []
        else:
            current.append(line)
    if current:
        blocks.append("\n".join(current) + "\n")
    return blocks


def parse_mdl_radicals(block: str) -> dict[int, str]:
    radicals: dict[int, str] = {}
    code_to_radical = {1: "SINGLET", 2: "DOUBLET", 3: "TRIPLET"}
    for raw_line in block.splitlines():
        if not raw_line.startswith("M  RAD"):
            continue
        fields = raw_line.split()
        if len(fields) < 4:
            continue
        try:
            pair_count = int(fields[2])
            values = [int(field) for field in fields[3:]]
        except ValueError:
            continue
        for offset in range(0, min(len(values), pair_count * 2), 2):
            if offset + 1 >= len(values):
                break
            atom_index = values[offset] - 1
            radical = code_to_radical.get(values[offset + 1])
            if atom_index >= 0 and radical is not None:
                radicals[atom_index] = radical
    return radicals


def parse_mdl_bond_stereo(block: str) -> dict[int, dict[str, str]]:
    lines = block.splitlines()
    if len(lines) < 4 or "V2000" not in lines[3]:
        return {}
    try:
        atom_count = int(lines[3][0:3])
        bond_count = int(lines[3][3:6])
    except ValueError:
        count_fields = lines[3].split()
        if len(count_fields) < 2:
            return {}
        try:
            atom_count = int(count_fields[0])
            bond_count = int(count_fields[1])
        except ValueError:
            return {}
    bond_start = 4 + atom_count
    overrides: dict[int, dict[str, str]] = {}
    for bond_index, line in enumerate(lines[bond_start : bond_start + bond_count]):
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            order = int(fields[2])
            stereo_code = int(fields[3])
        except ValueError:
            continue
        stereo = "STEREONONE"
        direction = "NONE"
        if order == 1 and stereo_code == 1:
            direction = "BEGINWEDGE"
        elif order == 1 and stereo_code == 4:
            direction = "UNKNOWN"
        elif order == 1 and stereo_code == 6:
            direction = "BEGINDASH"
        elif order == 2 and stereo_code == 3:
            stereo = "STEREOANY"
        overrides[bond_index] = {"stereo": stereo, "bond_direction": direction}
    return overrides


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
        "atoms": [
            atom_json(atom, record["radicals"].get(atom.GetIdx()))
            for atom in mol.GetAtoms()
        ],
        "bonds": [
            bond_json(bond, record["bond_stereo"].get(bond.GetIdx()))
            for bond in mol.GetBonds()
        ],
        "properties": molecule_properties(mol),
    }


def sdf_record_basic(record: dict[str, Any]) -> dict[str, Any]:
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
        "atoms": [basic_atom_json(atom) for atom in mol.GetAtoms()],
        "bonds": [basic_bond_json(bond) for bond in mol.GetBonds()],
        "properties": molecule_properties(mol),
    }


def mol_record(record: dict[str, Any]) -> dict[str, Any]:
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
        "atoms": [
            atom_json(atom, record["radicals"].get(atom.GetIdx()))
            for atom in mol.GetAtoms()
        ],
        "bonds": [
            bond_json(bond, record["bond_stereo"].get(bond.GetIdx()))
            for bond in mol.GetBonds()
        ],
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


def ring_set_record(record: dict[str, Any]) -> dict[str, Any]:
    from rdkit import Chem

    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    rings = [list(ring) for ring in Chem.GetSymmSSSR(mol)]
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "rings": rings,
    }


def conformer_record(record: dict[str, Any]) -> dict[str, Any]:
    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    conformers = []
    for conformer in mol.GetConformers():
        conformers.append(
            [
                {
                    "atom_index": index,
                    "x": conformer.GetAtomPosition(index).x,
                    "y": conformer.GetAtomPosition(index).y,
                    "z": conformer.GetAtomPosition(index).z,
                }
                for index in range(mol.GetNumAtoms())
            ]
        )
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_count": mol.GetNumAtoms(),
        "conformers": conformers,
        "atoms": [conformer_atom_json(atom) for atom in mol.GetAtoms()],
    }


def mol_parse_record(record: dict[str, Any]) -> dict[str, Any]:
    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    conformers = conformers_json(mol)
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_count": mol.GetNumAtoms(),
        "conformers": conformers,
        "atoms": [
            atom_json(atom, record["radicals"].get(atom.GetIdx()))
            for atom in mol.GetAtoms()
        ],
    }


def conformers_json(mol: Any) -> list[list[dict[str, Any]]]:
    conformers = []
    for conformer in mol.GetConformers():
        conformers.append(
            [
                {
                    "atom_index": index,
                    "x": conformer.GetAtomPosition(index).x,
                    "y": conformer.GetAtomPosition(index).y,
                    "z": conformer.GetAtomPosition(index).z,
                }
                for index in range(mol.GetNumAtoms())
            ]
        )
    return conformers


def conformer_atom_json(atom: Any) -> dict[str, Any]:
    return basic_atom_json(atom)


def basic_atom_json(atom: Any) -> dict[str, Any]:
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


def sanitized_atom_record(record: dict[str, Any]) -> dict[str, Any]:
    sanitized = clone_and_sanitize(record["mol"]) if record["mol"] is not None else None
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
        "atoms": [basic_atom_json(atom) for atom in sanitized.GetAtoms()],
    }


def valence_record(record: dict[str, Any]) -> dict[str, Any]:
    from rdkit import Chem

    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
        }
    prepared = Chem.Mol(mol)
    try:
        prepared.UpdatePropertyCache(strict=False)
    except Exception:
        return {
            "record_index": record["record_index"],
            "status": "valence_error",
            "title": record["title"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atoms": [valence_atom_json(atom) for atom in prepared.GetAtoms()],
    }


def smiles_parse_record(record: dict[str, Any]) -> dict[str, Any]:
    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
            "title": record["title"],
            "input_smiles": record["smiles"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "input_smiles": record["smiles"],
        "raw": smiles_raw_semantic_record(mol),
        "sanitized": smiles_sanitized_semantic_record(mol),
        "write_round_trip": smiles_sanitized_semantic_record(mol),
    }


def smiles_write_record(record: dict[str, Any]) -> dict[str, Any]:
    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
            "title": record["title"],
            "input_smiles": record["smiles"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "input_smiles": record["smiles"],
        "sanitized": smiles_sanitized_semantic_record(mol),
    }


def canonical_smiles_record(record: dict[str, Any], exact_smiles: bool) -> dict[str, Any]:
    from rdkit import Chem

    mol = record["mol"]
    if mol is None:
        return {
            "record_index": record["record_index"],
            "status": record["status"],
            "title": record["title"],
            "input_smiles": record["smiles"],
        }
    canonical = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=False,
    )
    canonical_mol = Chem.MolFromSmiles(canonical, sanitize=False)
    item = {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "input_smiles": record["smiles"],
        "sanitized": smiles_sanitized_semantic_record(canonical_mol),
    }
    if exact_smiles:
        item["canonical_smiles"] = canonical
    return item


def smiles_raw_semantic_record(mol: Any) -> dict[str, Any]:
    return {
        "atom_count": mol.GetNumAtoms(),
        "bond_count": mol.GetNumBonds(),
        "atoms": [basic_atom_json(atom) for atom in mol.GetAtoms()],
        "bonds": [basic_bond_json(bond) for bond in mol.GetBonds()],
    }


def smiles_sanitized_semantic_record(mol: Any) -> dict[str, Any]:
    sanitized = clone_and_sanitize(mol)
    if sanitized is None:
        return {"status": "sanitize_error"}
    return {
        "status": "ok",
        "atom_count": sanitized.GetNumAtoms(),
        "bond_count": sanitized.GetNumBonds(),
        "atoms": smiles_sanitized_atoms_json(sanitized),
        "bonds": smiles_sanitized_bonds_json(sanitized),
    }


def smiles_sanitized_atoms_json(mol: Any) -> list[dict[str, Any]]:
    atoms = []
    for atom in mol.GetAtoms():
        item = valence_atom_json(atom)
        item.pop("index", None)
        item["isotope"] = atom.GetIsotope() or None
        item["atom_map"] = atom.GetAtomMapNum() or None
        item["aromatic"] = atom.GetIsAromatic()
        item["no_implicit_hydrogens"] = atom.GetNoImplicit()
        neighbors = []
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            neighbors.append(
                {
                    "atom": smiles_sanitized_atom_key(neighbor),
                    "bond_type": smiles_semantic_bond_type(bond),
                    "is_aromatic": bond.GetIsAromatic(),
                }
            )
        neighbors.sort(key=lambda neighbor: json.dumps(neighbor, sort_keys=True))
        item["neighbors"] = neighbors
        atoms.append((smiles_sanitized_atom_sort_key(item), item))
    atoms.sort(key=lambda item: (item[0], json.dumps(item[1], sort_keys=True)))
    return [item for _, item in atoms]


def smiles_sanitized_atom_sort_key(atom: dict[str, Any]) -> str:
    no_implicit = str(atom["no_implicit_hydrogens"]).lower()
    aromatic = str(atom["aromatic"]).lower()
    return (
        f"{atom['atomic_number']:03}|{atom['symbol']}|{atom['formal_charge']}|"
        f"{atom['isotope'] or 0}|"
        f"{atom['explicit_hydrogens']}|{atom['implicit_hydrogens']}|"
        f"{no_implicit}|{atom['explicit_valence']}|{atom['atom_map'] or 0}|{aromatic}"
    )


def smiles_sanitized_atom_key(atom: Any) -> str:
    item = valence_atom_json(atom)
    item["isotope"] = atom.GetIsotope() or None
    item["atom_map"] = atom.GetAtomMapNum() or None
    item["aromatic"] = atom.GetIsAromatic()
    item["no_implicit_hydrogens"] = atom.GetNoImplicit()
    return smiles_sanitized_atom_sort_key(item)


def smiles_sanitized_bonds_json(mol: Any) -> list[dict[str, Any]]:
    bonds = [smiles_sanitized_bond_json(bond) for bond in mol.GetBonds()]
    bonds.sort(key=lambda item: json.dumps(item, sort_keys=True))
    return bonds


def smiles_sanitized_bond_json(bond: Any) -> dict[str, Any]:
    endpoints = sorted(
        [
            smiles_sanitized_atom_key(bond.GetBeginAtom()),
            smiles_sanitized_atom_key(bond.GetEndAtom()),
        ]
    )
    return {
        "endpoint_atoms": endpoints,
        "bond_type": smiles_semantic_bond_type(bond),
        "is_aromatic": bond.GetIsAromatic(),
    }


def smiles_semantic_bond_type(bond: Any) -> str:
    if bond.GetIsAromatic():
        return "AROMATIC"
    return str(bond.GetBondType())


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


def atom_json(atom: Any, radical_override: str | None = None) -> dict[str, Any]:
    radical, unpaired_electrons = radical_json(atom, radical_override)
    return {
        "index": atom.GetIdx(),
        "atomic_number": atom.GetAtomicNum(),
        "symbol": atom.GetSymbol(),
        "formal_charge": atom.GetFormalCharge(),
        "isotope": atom.GetIsotope() or None,
        "explicit_hydrogens": atom.GetNumExplicitHs(),
        "atom_map": atom.GetAtomMapNum() or None,
        "radical": radical,
        "unpaired_electrons": unpaired_electrons,
        "aromatic": atom.GetIsAromatic(),
    }


def radical_json(atom: Any, radical_override: str | None) -> tuple[str | None, int]:
    if radical_override is not None:
        return radical_override, {"SINGLET": 0, "DOUBLET": 1, "TRIPLET": 2}[radical_override]
    unpaired_electrons = atom.GetNumRadicalElectrons()
    if unpaired_electrons == 0:
        return None, 0
    if unpaired_electrons == 1:
        return "DOUBLET", 1
    if unpaired_electrons == 2:
        return "TRIPLET", 2
    return None, unpaired_electrons


def valence_atom_json(atom: Any) -> dict[str, Any]:
    from rdkit import Chem

    return {
        "index": atom.GetIdx(),
        "atomic_number": atom.GetAtomicNum(),
        "symbol": atom.GetSymbol(),
        "formal_charge": atom.GetFormalCharge(),
        "explicit_hydrogens": atom.GetNumExplicitHs(),
        "implicit_hydrogens": atom.GetNumImplicitHs(),
        "explicit_valence": atom.GetValence(Chem.rdchem.ValenceType.EXPLICIT),
    }


def bond_json(bond: Any, stereo_override: dict[str, str] | None = None) -> dict[str, Any]:
    stereo = stereo_override["stereo"] if stereo_override else str(bond.GetStereo())
    direction = (
        stereo_override["bond_direction"] if stereo_override else str(bond.GetBondDir())
    )
    if direction == "EITHERDOUBLE":
        direction = "NONE"
    return {
        "index": bond.GetIdx(),
        "begin_atom_index": bond.GetBeginAtomIdx(),
        "end_atom_index": bond.GetEndAtomIdx(),
        "bond_type": str(bond.GetBondType()),
        "is_aromatic": bond.GetIsAromatic(),
        "stereo": stereo,
        "bond_direction": direction,
    }


def basic_bond_json(bond: Any) -> dict[str, Any]:
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
    payload = (json.dumps(document, indent=2, sort_keys=True) + "\n").encode("utf-8")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(payload)


if __name__ == "__main__":
    raise SystemExit(main())
