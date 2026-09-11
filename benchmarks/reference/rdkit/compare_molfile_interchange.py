"""Optional RDKit Molfile interchange check using externally supplied corpora.

Build the xtask molfile_interchange_probe example first. This is an independent
cross-tool check, not a golden generator or a runtime dependency.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess

from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler

from compare_cip import DESCRIPTORS, RDKIT_VERSION, sha256


def labels(molecule):
    for item in list(molecule.GetAtoms()) + list(molecule.GetBonds()):
        if item.HasProp("_CIPCode"):
            item.ClearProp("_CIPCode")
    rdCIPLabeler.AssignCIPLabels(molecule)
    return sorted(
        [kind, item.GetIdx(), DESCRIPTORS.get(item.GetProp("_CIPCode"), item.GetProp("_CIPCode"))]
        for kind, items in [("atom", molecule.GetAtoms()), ("bond", molecule.GetBonds())]
        for item in items if item.HasProp("_CIPCode")
    )


def groups(molecule):
    return sorted(
        [str(group.GetGroupType()), sorted(atom.GetIdx() for atom in group.GetAtoms()),
         sorted(bond.GetIdx() for bond in group.GetBonds())]
        for group in molecule.GetStereoGroups()
    )


def chemical_graph(molecule):
    """Complete indexed graph after the same RDKit sanitization on both sides."""
    return {
        "atoms": [
            {"index": atom.GetIdx(), "element": atom.GetAtomicNum(),
             "isotope": atom.GetIsotope(), "charge": atom.GetFormalCharge(),
             "total_hydrogens": atom.GetTotalNumHs(includeNeighbors=True),
             "radical_electrons": atom.GetNumRadicalElectrons()}
            for atom in molecule.GetAtoms()
        ],
        "bonds": [
            {"index": bond.GetIdx(),
             "endpoints": sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())),
             "order": str(bond.GetBondType())}
            for bond in molecule.GetBonds()
        ],
    }


def expected_labels(molecule, clear_double_stereo):
    result = labels(molecule)
    if clear_double_stereo:
        result = [item for item in result if not (
            item[0] == "bond" and molecule.GetBondWithIdx(item[1]).GetBondType() == Chem.BondType.DOUBLE
        )]
    return result


def compare_output(expected, molecule, clear_double_stereo):
    if molecule is None:
        return {"status": "parse_error"}
    actual = {"graph": chemical_graph(molecule), "labels": labels(molecule), "groups": groups(molecule)}
    differences = [key for key in expected if actual[key] != expected[key]]
    if clear_double_stereo and any(
        bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetStereo() in (
            Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
            Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS,
        ) for bond in molecule.GetBonds()
    ):
        differences.append("cleared_double_stereo")
    return {"status": "difference" if differences else "match", "differences": differences, "actual": actual}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"expected RDKit {RDKIT_VERSION}, found {rdBase.rdkitVersion}")
    if args.output.exists():
        parser.error("output exists; choose a new filename to preserve prior evidence")
    args.probe = args.probe.resolve(strict=True)
    inputs = []
    cases = []
    tetrahedron = None
    tetrahedron_origin = None
    for path in sorted((args.repo / "benchmarks/corpora/pubchem-1k/data/packs").glob("*.sdf")):
        inputs.append({"path": str(path.resolve()), "sha256": sha256(path)})
        for index, record in enumerate(path.read_text().split("$$$$")):
            if not record.strip():
                continue
            source = record.lstrip("\r\n").split("M  END", 1)[0] + "M  END\n"
            molecule = Chem.MolFromMolBlock(source, removeHs=False)
            if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
                continue
            expected = labels(molecule)
            origin = {"path": str(path.resolve()), "record_index": index, "title": molecule.GetProp("_Name")}
            if tetrahedron is None and len(expected) == 1 and expected[0][0] == "atom":
                tetrahedron = Chem.RemoveHs(molecule)
                tetrahedron_origin = origin
            if any(item[0] == "bond" and item[2] in ("E", "Z") for item in expected) and sum(case["kind"] == "drawn" for case in cases) < 32:
                cases.append(dict(kind="drawn", name=f"{path.name}:{index}", molecule=molecule, origin=origin))
    assert cases and tetrahedron is not None, "external corpus must supply alkene and tetrahedral records"
    source_axis = args.repo / "benchmarks/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop1.mol"
    inputs.append({"path": str(source_axis.resolve()), "sha256": sha256(source_axis)})
    axis_molecule = Chem.MolFromMolBlock(source_axis.read_text(), removeHs=False)
    axis = next(bond.GetIdx() for bond in axis_molecule.GetBonds() if bond.GetStereo() in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW))
    for kind in [Chem.StereoGroupType.STEREO_ABSOLUTE, Chem.StereoGroupType.STEREO_AND, Chem.StereoGroupType.STEREO_OR]:
        grouped = Chem.RWMol(axis_molecule)
        grouped.SetStereoGroups([Chem.CreateStereoGroup(kind, grouped, [], [axis], 1)])
        cases.append(dict(kind="axis_group", name=str(kind), molecule=grouped, origin={"path": str(source_axis.resolve())}))

    requests = []
    for case in cases:
        for version in ([3000] if case["kind"] == "axis_group" else [2000, 3000]):
            requests.append(dict(**case, source=Chem.MolToMolBlock(case["molecule"], forceV3000=version == 3000), clear=False))
            if case["kind"] == "drawn":
                requests.append(dict(**case, source=requests[-1]["source"], clear=True))
    center = next(atom.GetIdx() for atom in tetrahedron.GetAtoms() if atom.HasProp("_CIPCode"))
    for cfg in [1, 2, 3]:
        source = Chem.MolToV3KMolBlock(tetrahedron)
        source = re.sub(r" CFG=[123]", "", source)
        rows = source.splitlines()
        start = rows.index("M  V30 BEGIN ATOM") + 1
        rows[start + center] += f" CFG={cfg}"
        source = "\n".join(rows) + "\n"
        reference = Chem.MolFromMolBlock(source, removeHs=False)
        # The normal reader retains molParity; its explicit conversion helper
        # is required. This fixture has no explicit H (CTfile H-last has a
        # separately documented RDKit helper discrepancy and Rust regression).
        Chem.AssignAtomChiralTagsFromMolParity(reference)
        requests.append(dict(kind="atom_cfg", name=f"CFG={cfg}", molecule=reference, source=source, clear=False, origin=tetrahedron_origin))
    completed = subprocess.run([str(args.probe.resolve())], input="".join(json.dumps({"molfile": case["source"], "clear_double_stereo": case["clear"]}) + "\n" for case in requests), text=True, capture_output=True, check=True)
    responses = [json.loads(line) for line in completed.stdout.splitlines()]
    assert len(responses) == len(requests)
    failures = []
    evidence = []
    checks = 0
    for request_id, (case, response) in enumerate(zip(requests, responses)):
        expected = {"graph": chemical_graph(case["molecule"]),
                    "labels": expected_labels(case["molecule"], case["clear"]),
                    "groups": groups(case["molecule"])}
        record = {"request_id": request_id, "kind": case["kind"], "name": case["name"],
                  "origin": case["origin"], "clear_double_stereo": case["clear"],
                  "molfile": case["source"], "source_sha256": hashlib.sha256(case["source"].encode("utf-8")).hexdigest(),
                  "expected": expected, "probe_response": response, "outputs": {}}
        evidence.append(record)
        if response["status"] != "ok" or response["labels"] != expected["labels"]:
            failures.append(dict(request_id=request_id, name=case["name"], stage="read", response=response, expected=expected))
            continue
        for version in ["v2000", "v3000"]:
            if case["kind"] == "axis_group" and version == "v2000":
                record["outputs"][version] = {"status": "expected_rejection" if "Err" in response[version] else "unexpected_success"}
                if "Err" not in response[version]:
                    failures.append(dict(request_id=request_id, name=case["name"], stage=version, error="enhanced group was not rejected"))
                continue
            result = response[version]
            if "Err" in result:
                record["outputs"][version] = {"status": "error", "message": result["Err"]}
                failures.append(dict(request_id=request_id, name=case["name"], stage=version, error=result["Err"]))
                continue
            reread = Chem.MolFromMolBlock(result["Ok"], removeHs=False)
            checks += 1
            comparison = compare_output(expected, reread, case["clear"])
            record["outputs"][version] = comparison
            if comparison["status"] != "match":
                failures.append(dict(request_id=request_id, name=case["name"], stage=version, comparison=comparison))
    summary = dict(schema_version=2, rdkit_version=rdBase.rdkitVersion, python=platform.python_version(),
                   inputs=inputs, probe={"path": str(args.probe), "sha256": sha256(args.probe)},
                   requests=len(requests), cross_tool_outputs=checks, failures=failures, request_evidence=evidence)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as output:
        output.write(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: summary[key] for key in ("rdkit_version", "requests", "cross_tool_outputs", "failures")}))
    raise SystemExit(bool(failures))


if __name__ == "__main__":
    main()
