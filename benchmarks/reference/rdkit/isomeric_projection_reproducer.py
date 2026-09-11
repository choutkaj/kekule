#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11"
# dependencies = ["rdkit==2026.3.6"]
# ///
"""Distinguish optional metal-neighbor and required charge brackets in SMILES.

Build the optional implementation adapter with:
    cargo build -p xtask --example smiles_write_probe --locked
Run from the repository root:
    uv run --python 3.13 benchmarks/reference/rdkit/isomeric_projection_reproducer.py \
        --probe target/debug/examples/smiles_write_probe
Append .exe on Windows. Omitting --probe checks only the external reference.
"""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

from rdkit import Chem, rdBase

from run_feature import (
    smiles_isomeric_stereo_semantic_record,
    smiles_perceived_semantic_record,
)


RDKIT_VERSION = "2026.03.6"
RDKIT_SOURCE = (
    "https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_6/"
    "Code/GraphMol/SmilesParse/SmilesWrite.cpp"
)
RDKIT_SOURCE_SHA256 = "b3ff6bff6a933b7dd1891122cbea180cfcd841b06820b79513cb934775b49efb"


def parse(smiles):
    molecule = Chem.MolFromSmiles(smiles, sanitize=False)
    if molecule is None:
        raise ValueError(f"invalid fixture or emitted SMILES: {smiles}")
    Chem.SanitizeMol(molecule)
    Chem.AssignStereochemistry(molecule, cleanIt=False, force=True)
    return molecule


def projection(molecule):
    return {
        "normalized_perceived": smiles_perceived_semantic_record(molecule),
        "stereo": smiles_isomeric_stereo_semantic_record(molecule),
    }


def graph_key(molecule):
    # This is an independent external check over the complete graph, beyond
    # the benchmark's sorted atom/neighbor projections. It does not invoke
    # Kekule's canonical writer or require identical emitted SMILES strings.
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def required_charge_projection(molecule):
    projected = Chem.Mol(molecule)
    changes = []
    for atom in projected.GetAtoms():
        if atom.GetFormalCharge() and not atom.GetNoImplicit():
            # Formal charge requires brackets in SMILES. The syntax cannot
            # retain inferred-H policy there, so encode the independently
            # perceived source total. Optional metal-neighbor brackets do not
            # enter this transformation.
            assert atom.GetSymbol() in ("Cl", "O")
            assert atom.GetTotalNumHs() == 0
            changes.append({"atom": atom.GetIdx(), "symbol": atom.GetSymbol(),
                            "formal_charge": atom.GetFormalCharge(),
                            "field": "no_implicit_hydrogens", "before": False, "after": True})
            atom.SetNoImplicit(True)
    assert {change["symbol"] for change in changes} == {"Cl", "O"}
    projected.UpdatePropertyCache(strict=True)
    assert graph_key(projected) == graph_key(molecule)
    return projected, changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"requires RDKit {RDKIT_VERSION}, found {rdBase.rdkitVersion}")
    if args.output and args.output.exists():
        parser.error("output exists; preserve it and choose a new filename")
    Chem.SetUseLegacyStereoPerception(False)
    fixture = Path(__file__).with_name("fixtures") / "metal-writer-projection.json"
    sources = json.loads(fixture.read_text(encoding="utf-8"))["cases"]
    outputs = None
    if args.probe:
        args.probe = args.probe.resolve(strict=True)
        process = subprocess.run(
            [str(args.probe)], input="\n".join(row["smiles"] for row in sources) + "\n",
            text=True, encoding="utf-8", capture_output=True, check=True, timeout=30,
        )
        outputs = [json.loads(line) for line in process.stdout.splitlines()]
        if len(outputs) != len(sources):
            raise ValueError("adapter must return exactly one record per input")
    cases = []
    for index, source in enumerate(sources):
        original = parse(source["smiles"])
        source_projection = projection(original)
        changes = []
        target = original
        if source["declaration_contract"] == "required_charge_brackets":
            target, changes = required_charge_projection(original)
        expected = projection(target)
        emitted = Chem.MolToSmiles(original, canonical=False, isomericSmiles=True)
        reparsed = parse(emitted)
        reference = projection(reparsed)
        same_reference_graph = graph_key(original) == graph_key(reparsed)
        # Retain the complete raw source and reference emission; the minimum
        # legal target is independently derived from source chemistry only.
        assert same_reference_graph, source
        assert source_projection != reference, source
        if changes:
            assert expected == reference, source
        case = {
            **source,
            "source_projection": source_projection,
            "required_declaration_changes": changes,
            "minimum_legal_source_projection": expected,
            "reference_emission": {
                "smiles": emitted, "projection": reference,
                "whole_graph_preserved": same_reference_graph,
            },
            "comparison": "reference_only",
        }
        if outputs is not None:
            output = outputs[index]
            if output.get("status") != "ok":
                raise ValueError(f"adapter failed for CID {source['cid']}: {output}")
            actual = parse(output["smiles"])
            actual_projection = projection(actual)
            same_graph = graph_key(original) == graph_key(actual)
            same_declarations = expected == actual_projection
            case["implementation_emission"] = {
                "smiles": output["smiles"], "projection": actual_projection,
                "whole_graph_preserved": same_graph,
                "minimum_legal_source_projection_preserved": same_declarations,
            }
            case["comparison"] = "match" if same_graph and same_declarations else "difference"
        cases.append(case)
    report = {
        "schema_version": 1, "rdkit": rdBase.rdkitVersion,
        "fixture_sha256": hashlib.sha256(fixture.read_bytes()).hexdigest(),
        "reference_source": {"url": RDKIT_SOURCE, "sha256": RDKIT_SOURCE_SHA256},
        "probe_sha256": hashlib.sha256(args.probe.read_bytes()).hexdigest() if args.probe else None,
        "cases": cases,
    }
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(json.dumps({"rdkit": RDKIT_VERSION, "cases": len(cases),
                          "comparisons": [case["comparison"] for case in cases]}))
    else:
        print(text, end="")
    return int(any(case["comparison"] == "difference" for case in cases))


if __name__ == "__main__":
    raise SystemExit(main())
