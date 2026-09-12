# /// script
# requires-python = ">=3.11"
# dependencies = ["rdkit==2026.3.6"]
# ///
"""Check both writers against external graphs, atom reorderings, and fixed points."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import random

from rdkit import Chem, RDLogger, rdBase

from probe_support import RDKIT_VERSION, run_probe, sha256, write_report


def parse(source):
    options = Chem.SmilesParserParams()
    options.removeHs = False
    options.parseName = False
    return Chem.MolFromSmiles(source, options)


def graph_key(molecule):
    """Normalize only removable H vertices, retaining isotope/map/stereo identity."""
    options = Chem.RemoveHsParameters()
    options.removeMapped = False
    options.removeIsotopes = False
    normalized = Chem.RemoveHs(molecule, options)
    return Chem.MolToSmiles(normalized, canonical=True, isomericSmiles=True)


def source_assertions(source):
    """Count supplied assertions before RDKit removes nonstereogenic tags."""
    if not any(marker in source for marker in ("@", "/", "\\")):
        return {"tetrahedral": 0, "double_bond": 0}
    options = Chem.SmilesParserParams()
    options.removeHs = False
    options.sanitize = False
    molecule = Chem.MolFromSmiles(source, options)
    Chem.SanitizeMol(molecule, sanitizeOps=(Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY))
    Chem.AssignStereochemistry(molecule, cleanIt=False, force=True)
    return {
        "tetrahedral": sum(atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                           for atom in molecule.GetAtoms()),
        "double_bond": sum(bond.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY)
                           for bond in molecule.GetBonds()),
    }


def compare_output(expected, result):
    if "Ok" not in result:
        return {"status": "write_error", "result": result}
    actual = parse(result["Ok"])
    if actual is None:
        return {"status": "reparse_error", "smiles": result["Ok"]}
    key = graph_key(actual)
    return {"status": "match" if key == expected else "difference",
            "smiles": result["Ok"], "graph": key}


def probe_batches(probe, requests, jobs):
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        batches = [requests[start:start + 32] for start in range(0, len(requests), 32)]
        futures = {pool.submit(run_probe, probe, [text for _, text in batch]): batch
                   for batch in batches}
        for future in as_completed(futures):
            yield futures[future], future.result()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", default=[])
    parser.add_argument("--corpus")
    parser.add_argument("--limit", type=int, help="First N distinct external inputs; omit for all")
    parser.add_argument("--variants", type=int, default=3)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--probe", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != RDKIT_VERSION:
        parser.error(f"expected RDKit {RDKIT_VERSION}, found {rdBase.rdkitVersion}")
    if args.jobs < 1 or args.variants < 0 or (args.limit is not None and args.limit <= 0):
        parser.error("jobs and limit must be positive; variants must be nonnegative")
    if args.output.exists():
        parser.error("choose a fresh output path")
    probe = args.probe.resolve()
    probe_hash = sha256(probe)
    root = Path(__file__).resolve().parents[3]
    paths = args.input + (sorted((root / "benchmarks/corpora" / args.corpus / "data/packs").glob("*.smi")) if args.corpus else [])
    sources = list(dict.fromkeys(line.split()[0] for path in paths
                                for line in path.read_text().splitlines()
                                if line.strip() and not line.lstrip().startswith("#")))
    if args.limit is not None:
        sources = sources[:args.limit]
    if not sources:
        parser.error("no external input records")
    Chem.SetUseLegacyStereoPerception(False)
    RDLogger.DisableLog("rdApp.*")
    rng = random.Random(0)
    cases = []
    requests = []
    for source in sources:
        molecule = parse(source)
        case = {"source": source, "variants": []}
        cases.append(case)
        if molecule is None:
            case["status"] = "reference_error"
            continue
        expected = graph_key(molecule)
        case["expected"] = expected
        assertions = source_assertions(source)
        case["source_assertions"] = assertions
        reference_roundtrip = graph_key(parse(expected))
        if reference_roundtrip != expected:
            case["status"] = "reference_error"
            case.setdefault("reference_errors", []).append({
                "input": expected, "graph": reference_roundtrip,
                "reason": "RDKit's own canonical output is not a fixed point",
            })
        variants = [source]
        for _ in range(args.variants):
            order = list(range(molecule.GetNumAtoms()))
            rng.shuffle(order)
            variant = Chem.MolToSmiles(Chem.RenumberAtoms(molecule, order), canonical=False, isomericSmiles=True)
            variant_molecule = parse(variant)
            variant_key = graph_key(variant_molecule) if variant_molecule is not None else None
            if variant_key != expected:
                case["status"] = "reference_error"
                case.setdefault("reference_errors", []).append({
                    "input": variant, "graph": variant_key,
                    "reason": "RDKit reordering changed its own graph identity",
                })
                continue
            variant_assertions = source_assertions(variant)
            if variant_assertions != assertions:
                case["status"] = "reference_error"
                case.setdefault("reference_errors", []).append({
                    "input": variant, "source_assertions": assertions,
                    "variant_assertions": variant_assertions,
                    "reason": "RDKit reordering removed supplied stereo assertions",
                })
                continue
            if variant not in variants:
                variants.append(variant)
        for variant in variants:
            record = {"input": variant}
            case["variants"].append(record)
            requests.append((case, record))
    probe_requests = [((case, record), record["input"]) for case, record in requests]
    for batch, outputs in probe_batches(probe, probe_requests, args.jobs):
        for ((case, record), _), output in zip(batch, outputs):
            record["response"] = output
            if output["status"] != "ok":
                record["status"] = "implementation_error"
                continue
            record["checks"] = {mode: compare_output(case["expected"], output[mode])
                                for mode in ("isomeric", "canonical")}
            record["status"] = "match" if all(value["status"] == "match" for value in record["checks"].values()) else "difference"
    fixed_points = []
    for case in cases:
        if case.get("status") == "reference_error":
            continue
        canonical = {record["checks"]["canonical"]["smiles"] for record in case["variants"]
                     if record.get("checks", {}).get("canonical", {}).get("status") == "match"}
        case["canonical_invariant"] = len(canonical) == 1
        case["status"] = "match" if case["canonical_invariant"] and all(record["status"] == "match" for record in case["variants"]) else "difference"
        if len(canonical) == 1:
            fixed_points.append((case, next(iter(canonical))))
    for batch, outputs in probe_batches(probe, fixed_points, args.jobs):
        for (case, written), output in zip(batch, outputs):
            case["fixed_point"] = output
            if output.get("canonical") != {"Ok": written}:
                case["status"] = "difference"
    counts = {status: sum(case["status"] == status for case in cases)
              for status in ("match", "difference", "reference_error")}
    if sha256(probe) != probe_hash:
        raise RuntimeError("probe executable changed during comparison; rerun with an immutable copy")
    report = {"rdkit": rdBase.rdkitVersion, "variants": args.variants, "limit": args.limit,
              "inputs": [{"path": str(path), "sha256": sha256(path)} for path in paths],
              "probe_sha256": probe_hash, "counts": counts, "cases": cases}
    write_report(args.output, report)
    print(json.dumps({"inputs": len(cases), "encodings": len(requests), **counts}))
    return int(counts["difference"] != 0 or counts["reference_error"] != 0)


if __name__ == "__main__":
    raise SystemExit(main())
