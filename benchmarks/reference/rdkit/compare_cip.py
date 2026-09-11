#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["rdkit==2026.3.6"]
# ///
"""Compare complete CIP descriptor maps on externally supplied SMILES inputs.

This optional scientific audit is not a runtime dependency or release gate.
Input checksums, reference version, sanitization mode, tag changes, failures,
and the probe binary checksum are retained in the output. Existing reports are
never overwritten. See docs/stereo-validation.md for the comparison contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger, rdBase


RDKIT_VERSION = "2026.03.6"
DESCRIPTORS = {
    "r": "LowerR", "s": "LowerS", "m": "LowerM", "p": "LowerP",
    "z": "SeqCis", "e": "SeqTrans", "seqCis": "SeqCis", "seqTrans": "SeqTrans",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atom_tags(molecule: Any) -> list[list[Any]]:
    return [
        [atom.GetIdx(), str(atom.GetChiralTag())]
        for atom in molecule.GetAtoms()
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    ]


def reference(smiles: str, mode: str, max_iterations: int) -> dict[str, Any]:
    stage = "parse"
    evidence: dict[str, Any] = {}
    try:
        # sanitize=False also retains explicit hydrogen vertices and atom indices.
        molecule = Chem.MolFromSmiles(smiles, sanitize=False)
        if molecule is None:
            raise ValueError("RDKit rejected SMILES syntax")
        evidence["parsed_tags"] = atom_tags(molecule)
        stage = "sanitize"
        operations = Chem.SanitizeFlags.SANITIZE_ALL
        if mode == "assertions":
            operations &= ~Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
        Chem.SanitizeMol(molecule, sanitizeOps=operations)
        evidence["sanitized_tags"] = atom_tags(molecule)
        stage = "stereo_perception"
        Chem.AssignStereochemistry(molecule, cleanIt=False, force=True)
        evidence["prepared_tags"] = atom_tags(molecule)
        for item in [*molecule.GetAtoms(), *molecule.GetBonds()]:
            if item.HasProp("_CIPCode"):
                item.ClearProp("_CIPCode")
        stage = "cip"
        Chem.AssignCIPLabels(molecule, maxRecursiveIterations=max_iterations)

        def label(item: Any) -> str:
            value = item.GetProp("_CIPCode")
            return DESCRIPTORS.get(value, value)

        atoms = sorted(
            [atom.GetIdx(), label(atom)]
            for atom in molecule.GetAtoms() if atom.HasProp("_CIPCode")
        )
        bonds = sorted(
            [*sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())), label(bond)]
            for bond in molecule.GetBonds() if bond.HasProp("_CIPCode")
        )
        return {
            "status": "ok", "atom_count": molecule.GetNumAtoms(),
            "bond_count": molecule.GetNumBonds(),
            "labels": {"atoms": atoms, "bonds": bonds}, "evidence": evidence,
        }
    except Exception as error:
        return {
            "status": "error", "stage": stage,
            "message": str(error), "evidence": evidence,
        }


def read_cases(paths: list[Path], input_format: str, stereo_only: bool) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    seen: set[str] = set()
    for path in paths:
        if input_format == "json":
            data = json.loads(path.read_text(encoding="utf-8"))
            rows = [
                {"smiles": row["smiles"], "source": row.get("source", str(path))}
                for row in data["cases"]
            ]
        else:
            rows = []
            for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.split("\t") if input_format == "suite" else line.split()
                row = {"smiles": fields[0], "source": f"{path}:{number}"}
                if input_format == "suite":
                    row["case_id"] = fields[1] if len(fields) > 1 else ""
                    row["published_labels"] = fields[2] if len(fields) > 2 else ""
                rows.append(row)
        for row in rows:
            smiles = row["smiles"]
            if smiles in seen or (stereo_only and not any(mark in smiles for mark in "@/\\")):
                continue
            seen.add(smiles)
            cases.append(row)
    return cases


def run_probe(probe: Path, cases: list[dict[str, Any]], timeout: float) -> list[dict[str, Any]]:
    try:
        process = subprocess.run(
            [str(probe)], input="\n".join(case["smiles"] for case in cases) + "\n",
            capture_output=True, text=True, encoding="utf-8", timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        if len(cases) > 1:
            return [run_probe(probe, [case], timeout)[0] for case in cases]
        return [{"status": "timeout", "seconds": timeout}]
    except OSError as error:
        return [{"status": "process_error", "message": str(error)} for _ in cases]
    lines = process.stdout.splitlines()
    if process.returncode or len(lines) != len(cases):
        if len(cases) > 1:
            return [run_probe(probe, [case], timeout)[0] for case in cases]
        return [{"status": "process_error", "exit_code": process.returncode,
                 "stdout": process.stdout, "stderr": process.stderr}]
    try:
        values = [json.loads(line) for line in lines]
        if not all(isinstance(value, dict) and "status" in value for value in values):
            raise ValueError("probe response must contain an object with status")
        return values
    except (ValueError, TypeError) as error:
        return [{"status": "protocol_error", "message": str(error)} for _ in cases]


def comparison(reference_value: dict[str, Any], actual: dict[str, Any]) -> str:
    if reference_value["status"] != "ok":
        return "reference_error"
    if actual["status"] != "ok":
        return "implementation_error"
    fields = ("atom_count", "bond_count", "labels")
    return "match" if all(reference_value[field] == actual.get(field) for field in fields) else "difference"


def summarize(cases: list[dict[str, Any]]) -> dict[str, int]:
    return {
        key: sum(case.get("comparison") == key for case in cases)
        for key in ("match", "difference", "reference_only", "reference_error", "implementation_error")
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", type=Path, default=[])
    parser.add_argument("--corpus", help="Existing benchmarks/corpora/<name>/data/packs/*.smi inputs")
    parser.add_argument("--input-format", choices=("smiles", "suite", "json"), default="smiles")
    parser.add_argument("--stereo-only", action="store_true")
    parser.add_argument("--mode", choices=("sanitized", "assertions"), required=True)
    parser.add_argument("--probe", type=Path, help="Compiled xtask cip_probe example; omit for references only")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-iterations", type=int, default=1_000_000)
    parser.add_argument("--timeout", type=float, default=30.0)
    args = parser.parse_args()
    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"expected RDKit {RDKIT_VERSION}, found {rdBase.rdkitVersion}")
    if args.max_iterations <= 0 or args.timeout <= 0:
        parser.error("iteration and time limits must be positive")
    repo = Path(__file__).resolve().parents[3]
    paths = list(args.input)
    if args.corpus:
        paths.extend(sorted((repo / "benchmarks" / "corpora" / args.corpus / "data" / "packs").glob("*.smi")))
    if not paths:
        parser.error("provide existing --input files or an installed --corpus")
    if args.output.exists():
        parser.error("output exists; choose a new filename to preserve prior evidence")
    if args.probe:
        args.probe = args.probe.resolve(strict=True)
    Chem.SetUseLegacyStereoPerception(False)
    RDLogger.DisableLog("rdApp.*")
    cases = read_cases(paths, args.input_format, args.stereo_only)
    if not cases:
        parser.error("no input records remain after filtering")
    for case in cases:
        case["reference"] = reference(case["smiles"], args.mode, args.max_iterations)
        case["comparison"] = (
            "reference_only" if case["reference"]["status"] == "ok" else "reference_error"
        )
    if args.probe:
        for index in range(0, len(cases), 32):
            batch = cases[index:index + 32]
            for case, actual in zip(batch, run_probe(args.probe, batch, args.timeout)):
                case["actual"] = actual
                case["comparison"] = comparison(case["reference"], actual)
    counts = summarize(cases)
    result = {
        "schema_version": 1, "rdkit": rdBase.rdkitVersion, "python": platform.python_version(),
        "mode": args.mode, "stereo_only": args.stereo_only, "deduplication": "identical input SMILES",
        "max_recursive_iterations": args.max_iterations, "probe_timeout_seconds": args.timeout,
        "inputs": [{"path": str(path), "sha256": sha256(path)} for path in paths],
        "probe": {"path": str(args.probe), "sha256": sha256(args.probe)} if args.probe else None,
        "case_count": len(cases), "counts": counts, "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("rdkit", "mode", "case_count", "counts")}))
    return 1 if counts["difference"] or counts["implementation_error"] or counts["reference_error"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
