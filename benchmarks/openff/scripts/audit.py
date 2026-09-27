"""Optional OpenFF prerequisite audit. External tools are reference-only.

inventory: standard-library-only source verification and OFFXML inventory.
reference: executable RDKit/OpenFF/NAGL observations plus native comparisons.
No result from a missing or failed implementation is counted as agreement.
"""
from __future__ import annotations

import argparse
import ast
import gzip
import hashlib
import importlib.metadata
import inspect
import json
from enum import Enum
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone
import xml.etree.ElementTree as ET

from paths import HERE
ROOT = HERE.parents[1]
MODEL_HASH = "7981e7f5b0b1e424c9e10a40d9e7606d96dcd3dd2b095cb4eeff6829f92238ee"
# User-approved absolute tolerance, in elementary charges. Archived observations
# retain their original threshold and raw arrays; classification uses this policy.
CHARGE_ATOL = 5e-5
VERSIONS = {"rdkit": "2026.3.3", "openff-toolkit": "0.19.0", "openff-nagl": "0.6.1",
            "openff-nagl-models": "2026.9.0", "openff-interchange": "0.5.5"}


def describe_configuration(value):
    # NAGL's activation fields contain Python classes, not JSON strings.
    # Preserve their qualified identities; never silently stringify other objects.
    if isinstance(value, type):
        return dict(python_type=value.__module__ + "." + value.__qualname__)
    if isinstance(value, Enum):
        return value.value
    raise TypeError(f"Unsupported report value: {type(value).__name__}")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def generation_metadata():
    return dict(started_utc=datetime.now(timezone.utc).isoformat(),
                audit_source_sha256=hashlib.sha256(Path(__file__).read_text(encoding="utf-8").encode()).hexdigest(),
                smarts_adapter_sha256=hashlib.sha256((ROOT / "benchmarks/smarts_conformance.py")
                    .read_text(encoding="utf-8").encode()).hexdigest(),
                native_source_sha256=hashlib.sha256((ROOT / "benchmarks/src/bin/openff_prerequisites.rs")
                    .read_text(encoding="utf-8").encode()).hexdigest())


def verify_sources():
    lock = json.loads((HERE / "fixtures/sources.lock.json").read_text())
    for source in lock["sources"]:
        if digest(HERE / source["path"]) != source["sha256"]:
            raise ValueError(f"Source checksum mismatch: {source['path']}")
    # Verify reused corpora against their existing provenance, not a newly computed pin.
    smoke = ROOT / "benchmarks/corpora/smoke"
    smoke_lock = json.loads((smoke / "sources.lock.json").read_text())
    for entry in smoke_lock["entries"]:
        for source in entry["files"]:
            if source["path"].startswith("data/pubchem_smiles/"):
                if digest(smoke / source["path"]) != source["sha256"]:
                    raise ValueError(f"Source checksum mismatch: {source['path']}")
    auxiliary = ROOT / "benchmarks/smarts-fixtures/openff-smarts"
    for source in json.loads((auxiliary / "sources.lock.json").read_text())["sources"]:
        if source["path"].endswith(".smi") and digest(auxiliary / source["path"]) != source["sha256"]:
            raise ValueError(f"Source checksum mismatch: {source['path']}")
    return lock


def inventory(path):
    root = ET.parse(path).getroot()
    sections = []
    for section in root:
        if section.tag in {"Author", "Date"}:
            continue
        rows = [dict(index=i, tag=item.tag, **item.attrib) for i, item in enumerate(section)]
        sections.append(dict(name=section.tag, attributes=dict(section.attrib), parameters=rows))
    return dict(root=dict(root.attrib), sections=sections,
                pattern_count=sum("smirks" in p for s in sections for p in s["parameters"]),
                bondorder_parameters=[p for s in sections for p in s["parameters"]
                                      if any("bondorder" in k for k in p)])


def literal(source, value):
    """Require that a selected case is a literal in its pinned external source."""
    tree = ast.parse((HERE / "fixtures" / source).read_text(encoding="utf-8"))
    matches = [node for node in ast.walk(tree)
               if isinstance(node, ast.Constant) and node.value == value]
    if not matches:
        raise ValueError(f"Case is not present in {source}: {value}")
    return dict(path="fixtures/" + source, line=min(n.lineno for n in matches))


def cases():
    selections = [
        ("lookup-entry-nitro", "test_lookups.py", "[H:5][C:1]([H:6])([H:7])[N+:2](=[O:3])[O-:4]"),
        ("lookup-charge-separated", "test_lookups.py", "[H:5][C:1]([H:6])([H:7])[N+2:2](-[O-:3])[O-:4]"),
        ("lookup-entry-sulfide", "test_lookups.py", "[H:2][S:1][H:3]"),
        ("lookup-reordered-sulfide", "test_lookups.py", "[H:1][S:2][H:3]"),
        ("lookup-miss-ethane", "test_lookups.py", "CC"),
        ("carboxylate", "conftest.py", "[H:1][C:2](=[O:3])[O-:4]"),
        ("methyl-methanoate", "conftest.py", "[H:5][C:1](=[O:2])[O:3][C:4]([H:6])([H:7])[H:8]"),
        ("methane", "conftest.py", "C"),
        ("amine-alcohol", "conftest.py", "[H:6][C:1]([H:7])([H:8])[C:2]([H:9])([H:10])[N:3]([H:11])[C:4]([H:12])([H:13])[O:5][H:14]"),
        ("dimethylamine", "conftest.py", "[H:4][C:1]([H:5])([H:6])[N:2]([H:7])[C:3]([H:8])([H:9])[H:10]"),
        ("resonance-fragments", "test_resonance.py", "[O-:1][N+:2](=[O:3])[N:4](-[H:16])[c:5]1[c:6](-[H:12])[c:7](-[H:13])[c:8](-[H:14])[c:9](-[H:15])[n+:10]1[O-:11]"),
    ]
    result = [dict(id=name, smiles=s, source=literal(path, s)) for name, path, s in selections]
    # Reuse externally supplied PubChem inputs and their original provenance.
    files = sorted((ROOT / "benchmarks/corpora/smoke/data/pubchem_smiles").glob("*.txt"))
    files += sorted((ROOT / "benchmarks/smarts-fixtures/openff-smarts").glob("cid_*.smi"))
    for path in files:
        result.append(dict(id=path.stem, smiles=path.read_text().split()[0],
                           source=dict(path=path.relative_to(ROOT).as_posix(), sha256=digest(path))))
    return result


def native(binary, rows):
    result = subprocess.run([str(Path(binary).resolve())],
                            input="".join(json.dumps(row) + "\n" for row in rows),
                            capture_output=True, text=True, check=True, timeout=300)
    output = [json.loads(line) for line in result.stdout.splitlines()]
    if len(output) != len(rows):
        raise ValueError("Incomplete native record stream")
    return output


def classify(expected, actual):
    if expected.get("status") != "ok":
        return "reference_error"
    if actual.get("status") != "ok":
        return "implementation_error"
    return "equal" if expected == actual else "mismatch"


def reference_graph(molecule):
    from rdkit import Chem
    mol = Chem.Mol(molecule)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    return dict(status="ok", atoms=sorted([
        dict(map=a.GetAtomMapNum(), element=a.GetAtomicNum(), formal_charge=a.GetFormalCharge(),
             degree=a.GetDegree(), rings_3_6=[a.IsInRingSize(n) for n in range(3, 7)],
             mdl_aromatic=a.GetIsAromatic()) for a in mol.GetAtoms()], key=lambda a: a["map"]))


def mapped_rdkit(smiles):
    from rdkit import Chem
    options = Chem.SmilesParserParams()
    options.removeHs = False
    mol = Chem.MolFromSmiles(smiles, options)
    if mol is None:
        raise ValueError("RDKit rejected source molecule")
    mol = Chem.AddHs(mol)
    labels = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    if not all(labels):
        # A new complete correspondence for unmapped source inputs, before either observer.
        for i, atom in enumerate(mol.GetAtoms(), 1):
            atom.SetAtomMapNum(i)
    elif len(set(labels)) != len(labels):
        raise ValueError("Duplicate source atom maps")
    return mol


def lookup_contract():
    """Execute upstream small-table regressions, including the relaxed-map case."""
    from openff.toolkit import Molecule
    from openff.nagl.lookups import AtomPropertiesLookupTable, AtomPropertiesLookupTableEntry
    entries = [
        dict(inchi="InChI=1/CH3NO2/c1-2(3)4/h1H3",
             mapped_smiles=cases()[0]["smiles"],
             property_value=[-0.103, 0.234, -0.209, -0.209, 0.096, 0.096, 0.096]),
        dict(inchi="InChI=1/H2S/h1H2", mapped_smiles=cases()[2]["smiles"],
             property_value=[-0.441, 0.22, 0.22]),
    ]
    table = AtomPropertiesLookupTable(property_name="test", properties=[
        AtomPropertiesLookupTableEntry(**e, provenance={"description": "upstream test_lookups.py"})
        for e in entries])
    records = []
    for i, expected in [(1, entries[0]["property_value"]), (3, [0.22, -0.441, 0.22])]:
        case = cases()[i]
        mol = Molecule.from_mapped_smiles(case["smiles"], allow_undefined_stereo=True)
        value = table.lookup(mol).tolist()
        records.append(dict(case=case["id"], fixed_h_inchi=mol.to_inchi(fixed_hydrogens=True),
                            expected=expected, actual=value,
                            passed=all(abs(a-b) < 1e-7 for a, b in zip(value, expected))
                            and len(value) == len(expected)))
    try:
        table.lookup(Molecule.from_smiles(cases()[4]["smiles"]))
    except KeyError:
        missed = True
    else:
        missed = False
    records.append(dict(case="lookup-miss-ethane", passed=missed))
    left = Molecule.from_mapped_smiles(cases()[0]["smiles"])
    right = Molecule.from_mapped_smiles(cases()[1]["smiles"])
    records.append(dict(case="same-inchi-different-represented-chemistry",
                        left=left.to_inchi(fixed_hydrogens=True), right=right.to_inchi(fixed_hydrogens=True),
                        strict_isomorphic=Molecule.are_isomorphic(left, right)[0],
                        passed=left.to_inchi(fixed_hydrogens=True) == right.to_inchi(fixed_hydrogens=True)
                        and not Molecule.are_isomorphic(left, right)[0]))
    return records


def matching_observations(args, case):
    from rdkit import Chem
    sys.path.insert(0, str(ROOT / "benchmarks"))
    from smarts_conformance import query_rows, reference
    mol = mapped_rdkit(case["smiles"])
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    matching_smiles = Chem.MolToSmiles(Chem.RemoveHs(mol), canonical=False)
    rows = [dict(q, smiles=matching_smiles) for q in query_rows([HERE / "fixtures/rosemary.offxml"])]
    actual = native(args.smarts_binary, rows)
    observations = []
    for row, observed in zip(rows, actual):
        try:
            expected = reference(row, Chem)
        except Exception as error:
            expected = dict(status="error", message=str(error))
        observations.append(dict(case=case["id"], section=row["section"], row=row["row"],
            smarts=row["smarts"], expected=expected, actual=observed,
            comparison=classify(expected, observed)))
    return matching_smiles, observations


def run_primitives(args):
    """Run independently of NAGL: no missing charge engine is treated as a pass."""
    from rdkit import Chem, rdBase
    if rdBase.rdkitVersion != "2026.03.3":
        raise ValueError(f"Expected RDKit 2026.03.3, got {rdBase.rdkitVersion}")
    report = dict(schema=1, stage="primitives", complete=False, generation=generation_metadata(), rdkit=rdBase.rdkitVersion,
                  sources=verify_sources(), records=[], smarts=[],
                  implementation=dict(graph_sha256=digest(args.graph_binary),
                                      smarts_sha256=digest(args.smarts_binary)))
    for case in cases():
        record = dict(input=case)
        report["records"].append(record)
        try:
            mol = mapped_rdkit(case["smiles"])
            mapped = Chem.MolToSmiles(mol, canonical=False, allHsExplicit=True)
            expected = reference_graph(mol)
            actual = native(args.graph_binary, [dict(smiles=mapped)])[0]
            record.update(status="ok", mapped_smiles=mapped,
                          fixed_h_inchi=Chem.MolToInchi(mol, options="/FixedH"),
                          graph=dict(expected=expected, actual=actual, comparison=classify(expected, actual)))
            record["matching_smiles"], observations = matching_observations(args, case)
            report["smarts"].extend(observations)
        except Exception as error:
            record.update(status="error", message=str(error), kind=type(error).__name__)
        publish(args.output, report)
    report["complete"] = True
    report["summary"] = summarize(report)
    return report


def run_reference(args):
    import numpy as np
    from rdkit import Chem
    from openff.toolkit import Molecule, ForceField
    from openff.toolkit.utils import RDKitToolkitWrapper, ToolkitRegistry
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.toolkit.utils.toolkit_registry import toolkit_registry_manager
    from openff.nagl import GNNModel
    from openff.nagl.features.atoms import AtomAverageFormalCharge
    from openff.nagl.toolkits.openff import normalize_molecule
    from openff.nagl_models import get_model
    from openff.units import unit

    from packaging.version import Version
    versions = {name: importlib.metadata.version(name) for name in VERSIONS}
    if any(Version(versions[name]) != Version(version) for name, version in VERSIONS.items()):
        raise ValueError(f"Reference environment differs from pin: {versions}")
    versions.update({name: importlib.metadata.version(name)
                     for name in ["torch", "numpy", "openff-interchange"]})
    model_path = Path(args.model) if args.model else Path(get_model("openff-gnn-am1bcc-1.0.0.pt"))
    if digest(model_path) != MODEL_HASH:
        raise ValueError("Charge model checksum mismatch")
    model = GNNModel.load(model_path, eval_mode=True)
    forcefield = ForceField(str(HERE / "fixtures/rosemary.offxml"))
    report = dict(schema=1, stage="reference", complete=False, generation=generation_metadata(), versions=versions, sources=verify_sources(),
                  model_sha256=digest(model_path), model_config=model.config.model_dump(mode="python"),
                  chemical_domain=model.chemical_domain.model_dump(mode="python"),
                  lookup_table_sizes={name: len(table) for name, table in model.lookup_tables.items()},
                  implementation=dict(graph_sha256=digest(args.graph_binary),
                                      smarts_sha256=digest(args.smarts_binary)),
                  records=[], smarts=[])
    # Bind runtime code as well as package versions; source-review commits may be newer.
    import openff.nagl.lookups, openff.nagl.utils.resonance, openff.nagl.toolkits.openff
    report["reference_code"] = {module.__name__: digest(inspect.getfile(module)) for module in
                                [openff.nagl.lookups, openff.nagl.utils.resonance, openff.nagl.toolkits.openff]}
    registry = ToolkitRegistry([RDKitToolkitWrapper()])
    charge_wrapper = NAGLToolkitWrapper()
    system_registry = ToolkitRegistry([RDKitToolkitWrapper(), charge_wrapper])
    with toolkit_registry_manager(registry):
        report["lookup_contract"] = lookup_contract()
        for case in cases():
            record = dict(input=case)
            report["records"].append(record)
            try:
                rd = mapped_rdkit(case["smiles"])
                # Preserve maps and hydrogen vertices; NAGL orders atoms by map label.
                mapped = Chem.MolToSmiles(rd, canonical=False, allHsExplicit=True)
                off = Molecule.from_mapped_smiles(mapped, allow_undefined_stereo=True)
                record["mapped_smiles"] = mapped
                record["charge_atom_maps"] = sorted(a.GetAtomMapNum() for a in rd.GetAtoms())
                expected = reference_graph(rd)
                observed = native(args.graph_binary, [dict(smiles=mapped)])[0]
                record["graph"] = dict(expected=expected, actual=observed,
                                       comparison=classify(expected, observed))
                record["fixed_h_inchi"] = off.to_inchi(fixed_hydrogens=True)
                record["standard_inchi"] = off.to_inchi()
                normalized = normalize_molecule(off)
                record["normalized_mapped_smiles"] = normalized.to_smiles(mapped=True)
                record["average_formal_charge"] = AtomAverageFormalCharge().encode(off).flatten().tolist()
                record["domain"] = model.chemical_domain.check_molecule(off, return_error_message=True)
                # Preserve unsupported feature encodings (e.g. Na+) independently;
                # library charges may still make the force field applicable.
                try:
                    record["features"] = dict(status="ok", values=[dict(config=f.model_dump(mode="python"),
                        values=f.encode(off).tolist()) for f in model.config.atom_features])
                except Exception as error:
                    record["features"] = dict(status="error", kind=type(error).__name__, message=str(error))
                try:
                    lookup = model._check_property_lookup_table(off, "am1bcc_charges")
                    record["lookup"] = dict(status="hit", charges=lookup.flatten().tolist())
                except KeyError as error:
                    record["lookup"] = dict(status="miss", message=str(error))
                # Unsupported chemistry remains a reported charge error, not dropped input.
                try:
                    charges = model.compute_property(off, readout_name="am1bcc_charges", check_domains=True,
                                                     error_if_unsupported=True)
                    record["charges"] = dict(status="ok", raw=np.asarray(charges).flatten().tolist())
                    off.assign_partial_charges(str(model_path), toolkit_registry=charge_wrapper)
                    record["charges"]["toolkit_normalized"] = off.partial_charges.m_as(unit.elementary_charge).tolist()
                    record["charges"]["formal_total"] = off.total_charge.m_as(unit.elementary_charge)
                except Exception as error:
                    record.setdefault("charges", {}).update(status="error", kind=type(error).__name__, message=str(error))
                labels = forcefield.label_molecules(off.to_topology())[0]
                record["assigned_parameters"] = {name: [dict(atoms=list(key), id=value.id, smirks=value.smirks)
                                                          for key, value in terms.items()]
                                                  for name, terms in labels.items()}
                try:
                    system = forcefield.create_interchange(off.to_topology(), toolkit_registry=system_registry)
                    record["system"] = dict(status="ok", collections={
                        name: dict(assignments=len(collection.key_map), potentials=len(collection.potentials))
                        for name, collection in system.collections.items()},
                        charges=[dict(atoms=list(key.atom_indices), charge=value.m_as(unit.elementary_charge))
                                 for key, value in system.collections["Electrostatics"].charges.items()])
                except Exception as error:
                    record["system"] = dict(status="error", kind=type(error).__name__, message=str(error))
                # Atom permutation check exercises lookup remapping and GNN equivariance.
                permutation = {i: off.n_atoms - 1 - i for i in range(off.n_atoms)}
                reverse = off.remap(permutation)
                if record["charges"]["status"] == "ok":
                    reverse.assign_partial_charges(str(model_path), toolkit_registry=charge_wrapper)
                    remapped = reverse.partial_charges.m_as(unit.elementary_charge)[::-1]
                    original = np.array(record["charges"]["toolkit_normalized"])
                    record["charge_permutation"] = dict(actual=remapped.tolist(),
                        max_abs_error=float(np.max(np.abs(remapped-original))),
                        absolute_tolerance_e=CHARGE_ATOL,
                        passed=bool(np.allclose(remapped, original, rtol=0, atol=CHARGE_ATOL)))
                record["status"] = "ok"
            except Exception as error:
                record.update(status="error", kind=type(error).__name__, message=str(error))
            # The legacy SMARTS protocol compares traversal indices. Its RDKit adapter
            # removes input graph H and then appends H, unlike Kekule's preserving parser.
            # Supply the SAME unmapped hydrogen-suppressed serialization to both; retain
            # the original input and this explicit representation transform in the report.
            record["matching_smiles"], observations = matching_observations(args, case)
            report["smarts"].extend(observations)
            publish(args.output, report)
    report["complete"] = True
    report["summary"] = summarize(report)
    return report


def summarize(report):
    comparisons = [r["comparison"] for r in report["smarts"]]
    graph = [r.get("graph", {}).get("comparison", "not_observed") for r in report["records"]]
    return dict(cases=len(report["records"]), smarts={v: comparisons.count(v) for v in sorted(set(comparisons))},
                graph={v: graph.count(v) for v in sorted(set(graph))},
                reference_errors=sum(r["status"] != "ok" for r in report["records"]),
                charge_errors=sum(r["charges"]["status"] != "ok" for r in report["records"] if "charges" in r),
                charge_cases=sum("charges" in r for r in report["records"]),
                system_errors=sum(r["system"]["status"] != "ok" for r in report["records"] if "system" in r),
                system_cases=sum("system" in r for r in report["records"]),
                charge_permutation_failures=[r["input"]["id"] for r in report["records"]
                                             if permutation_failed(r.get("charge_permutation"))],
                lookup_contract_passed=all(r["passed"] for r in report["lookup_contract"])
                    if "lookup_contract" in report else None)


def permutation_failed(observation):
    if observation is None:
        return False
    error = observation.get("max_abs_error")
    if error is None:
        return observation.get("passed") is False
    return not (0 <= error <= CHARGE_ATOL)


def audit_failed(report):
    summary = summarize(report)
    failed = (not report["complete"] or summary["reference_errors"] > 0
              or summary["lookup_contract_passed"] is False
              or any(v != "equal" for v in summary["graph"])
              or any(v != "equal" for v in summary["smarts"])
              or bool(summary["charge_permutation_failures"]))
    if report.get("stage") == "reference":
        # A domain rejection remains visible, but library assignment can still succeed.
        # Unexpected failures and absent observations may never make the audit pass.
        failed |= summary["system_cases"] != summary["cases"] or summary["system_errors"] > 0
        failed |= summary["charge_cases"] != summary["cases"]
        failed |= any(r.get("charges", {}).get("status") != "ok" and r.get("domain", [True])[0]
                      for r in report["records"])
    return bool(failed)


def publish(path, report):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(report, indent=2, allow_nan=False, default=describe_configuration) + "\n").encode()
    if path.suffix == ".gz":
        payload = gzip.compress(payload, mtime=0)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(payload)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("mode", choices=["inventory", "primitives", "reference"])
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--graph-binary", type=Path)
    parser.add_argument("--smarts-binary", type=Path)
    parser.add_argument("--model", type=Path)
    args = parser.parse_args()
    verify_sources()
    if args.mode == "inventory":
        report = dict(schema=1, sources=verify_sources(), forcefield=inventory(HERE / "fixtures/rosemary.offxml"),
                      cases=cases())
    else:
        if not args.graph_binary or not args.smarts_binary:
            parser.error("reference requires both native observer binaries")
        report = run_reference(args) if args.mode == "reference" else run_primitives(args)
    publish(args.output, report)
    print(json.dumps(dict(cases=len(report["cases"]), patterns=report["forcefield"]["pattern_count"])
                     if args.mode == "inventory" else report["summary"]))
    if args.mode != "inventory":
        return int(audit_failed(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
