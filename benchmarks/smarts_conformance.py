"""Optional full-mapping comparison against pinned RDKit; no runtime dependency.

Accepts external query tables/OFFXML and SMILES files. Every query row and every
selected molecule is retained, including parse failures and resource errors.
The report contains inputs, complete observations, checksums, and tool versions.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def query_rows(paths):
    for path in paths:
        if path.suffix == ".offxml":
            root = ET.parse(path).getroot()
            for section in root:
                for index, item in enumerate(section):
                    if "smirks" in item.attrib:
                        yield dict(smarts=item.attrib["smirks"], source=str(path),
                                   section=section.tag, row=index, mdl=True,
                                   explicit_hydrogens=True)
        else:
            for index, line in enumerate(path.read_text().splitlines(), 1):
                if line.strip() and not line.lstrip().startswith("#"):
                    yield dict(smarts=line.split()[0], source=str(path), row=index,
                               mdl=False, explicit_hydrogens=False)


def reference(row, Chem):
    query = Chem.MolFromSmarts(row["smarts"])
    if query is None:
        return {"status": "parse_error"}
    mol = Chem.MolFromSmiles(row["smiles"])
    if mol is None:
        return {"status": "target_error"}
    if row["explicit_hydrogens"]:
        mol = Chem.AddHs(mol)
    if row["mdl"]:
        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    matches = sorted(mol.GetSubstructMatches(query, uniquify=False,
                                           useChirality=True, maxMatches=1_000_001))
    if len(matches) > 1_000_000:
        return {"status": "resource_limit", "resource": "matches"}
    tags = sorted((a.GetAtomMapNum(), a.GetIdx()) for a in query.GetAtoms()
                  if a.GetAtomMapNum())
    return dict(status="ok", atom_count=query.GetNumAtoms(), bond_count=query.GetNumBonds(),
                matches=[list(m) for m in matches], tags=[t for t, _ in tags],
                tagged_matches=[[m[i] for _, i in tags] for m in matches],
                target_aromatic_atoms=[a.GetIsAromatic() for a in mol.GetAtoms()])


def classify(expected, actual):
    if expected == actual:
        return "equal"
    if expected.get("status") == "parse_error" and actual.get("status") == "parse_error":
        return "both_rejected"
    if expected.get("status") != "ok":
        return "reference_error"
    if actual.get("status") == "parse_error" and actual.get("kind") == "Unsupported":
        return "dialect_exclusion"
    if actual.get("status") != "ok":
        return "implementation_error"
    return "mismatch"


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--queries", nargs="+", type=Path, required=True)
    parser.add_argument("--molecules", nargs="+", type=Path, required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    from rdkit import Chem, rdBase, RDLogger
    RDLogger.DisableLog("rdApp.*")
    if rdBase.rdkitVersion != "2026.03.3":
        raise SystemExit(f"Expected pinned RDKit 2026.03.3, got {rdBase.rdkitVersion}")
    molecules = []
    for path in args.molecules:
        for index, line in enumerate(path.read_text().splitlines(), 1):
            if line.strip() and not line.lstrip().startswith("#"):
                molecules.append((line.split()[0], str(path), index))
    queries = list(query_rows(args.queries))
    rows = [dict(q, smiles=s, molecule_source=p, molecule_row=i)
            for s, p, i in molecules for q in queries]
    result = subprocess.run([str(args.binary.resolve())],
                            input="".join(json.dumps(r)+"\n" for r in rows),
                            capture_output=True, text=True, check=True)
    actual = [json.loads(line) for line in result.stdout.splitlines()]
    if len(actual) != len(rows):
        raise RuntimeError("Rust observer returned an incomplete record stream")
    counts = {}
    records = []
    for row, observed in zip(rows, actual):
        try:
            expected = reference(row, Chem)
        except Exception as exc:
            expected = {"status": "reference_error", "message": str(exc)}
        status = classify(expected, observed)
        counts[status] = counts.get(status, 0)+1
        records.append(dict(input=row, comparison=status, expected=expected, actual=observed))
    sources = [dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest())
               for p in args.queries+args.molecules]
    report = dict(schema=1, rdkit=rdBase.rdkitVersion, python=sys.version,
                  rust_binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
                  sources=sources, query_rows=len(queries), molecule_rows=len(molecules),
                  summary=counts, records=records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({k: report[k] for k in ("rdkit", "query_rows", "molecule_rows", "summary")}))
    return int(any(counts.get(k) for k in ("mismatch", "implementation_error", "reference_error")))


if __name__ == "__main__":
    raise SystemExit(main())
