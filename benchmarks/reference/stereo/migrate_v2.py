#!/usr/bin/env python3
"""Migrate reviewed stereo goldens across the source/graph boundary.

This program never calls kekule or reads implementation output. The archived
schema-1 assertions are its authority for existing configurations. Source text,
RDKit's independently parsed graph, permutation parity, and signed geometry
provide the explicitly documented schema/behavior mappings.
"""

from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path

VERSION = "2026.03.6"
FEATURES = ("stereo.representation", "stereo.perception", "io.smiles.isomeric")
ROOT = Path(__file__).resolve().parents[3]
ARCHIVE = Path(__file__).parent / "schema-v1"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_golden(path):
    return json.loads(gzip.decompress(path.read_bytes()))


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def merged_provenance(previous, changes):
    if previous["schema_version"] != 2 or previous["rdkit"] != VERSION:
        raise ValueError("incompatible migration provenance already exists")
    merged = {entry["path"]: entry for entry in previous["fixtures"]}
    merged.update((entry["path"], entry) for entry in changes)
    return [merged[key] for key in sorted(merged)]


def changed_paths(old, new, path="$"):
    if type(old) is not type(new):
        return [path]
    if isinstance(old, dict):
        missing = old.keys() - new.keys()
        if missing:
            raise ValueError(f"meaningful reference fields removed at {path}: {sorted(missing)}")
        return [difference for key in old for difference in changed_paths(old[key], new[key], f"{path}.{key}")]
    if isinstance(old, list):
        if len(old) != len(new):
            return [path]
        return [difference for i, (left, right) in enumerate(zip(old, new))
                for difference in changed_paths(left, right, f"{path}[{i}]")]
    return [] if old == new else [path]


def archive(path, corpus, feature, write=False):
    destination = ARCHIVE / corpus / feature / path.name
    if destination.exists():
        return destination
    if read_golden(path)["schema_version"] != 1:
        raise ValueError(f"missing schema-1 archive for {path}")
    if write:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(path.read_bytes())
        return destination
    return path


def carrier_key(carrier):
    if "atom_index" in carrier:
        return (0, carrier["atom_index"])
    return (1 if "implicit_hydrogen" in carrier else 2, 0)


def invert(orientation):
    return {None: None, "clockwise": "counter_clockwise", "counter_clockwise": "clockwise",
            "together": "opposite", "opposite": "together"}[orientation]


def canonical_element(element, mol):
    result = copy.deepcopy(element)
    result.pop("source", None)
    result.pop("specifiedness", None)
    if result["type"] == "tetrahedral":
        old = result["carriers"]
        order = sorted(range(len(old)), key=lambda i: carrier_key(old[i]))
        odd = sum(a > b for i, a in enumerate(order) for b in order[i + 1:]) % 2
        result["carriers"] = [old[i] for i in order]
        if odd:
            result["orientation"] = invert(result["orientation"])
    elif result["type"] == "double_bond":
        left, right = result["left_atom_index"], result["right_atom_index"]
        if left > right:
            left, right = right, left
            result["left_atom_index"], result["right_atom_index"] = left, right
            result["left_carrier"], result["right_carrier"] = result["right_carrier"], result["left_carrier"]
        for side, endpoint, other in [("left", left, right), ("right", right, left)]:
            neighbors = sorted(a.GetIdx() for a in mol.GetAtomWithIdx(endpoint).GetNeighbors() if a.GetIdx() != other)
            carrier = {"atom_index": neighbors[0]} if neighbors else {"implicit_hydrogen": True}
            if result[f"{side}_carrier"] != carrier:
                result["orientation"] = invert(result["orientation"])
                result[f"{side}_carrier"] = carrier
    return result


def source_document(source, format_, marks, elements, origins):
    return {"format": format_, "source": source, "stereo_bond_marks": copy.deepcopy(marks),
            "stereo_sources": [{"element_index": e["index"], "source": origin,
                                "specifiedness": "unknown" if e["orientation"] is None else "specified"}
                               for e, origin in zip(elements, origins)]}


def source_records(path):
    text = path.read_text(encoding="utf-8")
    if path.suffix in {".smi", ".smiles", ".txt"}:
        return {i: line.split()[0] for i, line in enumerate(text.splitlines())
                if line.strip() and not line.lstrip().startswith("#")}
    records = {}
    lines = []
    for line in text.splitlines():
        if line.strip() == "$$$$":
            end = next(i for i, value in enumerate(lines) if value.strip() == "M  END")
            records[len(records)] = "\n".join(lines[:end + 1]) + "\n"
            lines = []
        else:
            lines.append(line)
    if any(line.strip() for line in lines):
        end = next(i for i, value in enumerate(lines) if value.strip() == "M  END")
        records[len(records)] = "\n".join(lines[:end + 1]) + "\n"
    return records


def graph(source, format_, Chem):
    if format_ == "smiles":
        options = Chem.SmilesParserParams()
        options.removeHs = False
        options.sanitize = False
        mol = Chem.MolFromSmiles(source, options)
    else:
        mol = Chem.MolFromMolBlock(source, sanitize=False, removeHs=False, strictParsing=True)
    if mol is None:
        raise ValueError("independent source parser rejected an archived successful record")
    mol.UpdatePropertyCache(strict=False)
    return mol


def vector(a, b):
    return tuple(x - y for x, y in zip(a, b))


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def tetra_from_geometry(candidate, mol):
    conformer = mol.GetConformer()
    def point(i):
        return tuple(conformer.GetAtomPosition(i))
    points = [point(c.get("atom_index", candidate["center_atom_index"])) for c in candidate["carriers"]]
    if len(points) != 4:
        raise ValueError("tetrahedron must have four spatial carriers")
    vectors = [vector(p, points[3]) for p in points[:3]]
    product = dot(vectors[0], cross(vectors[1], vectors[2]))
    scale = 1.0
    for v in vectors:
        scale *= dot(v, v) ** 0.5
    if scale == 0 or abs(product / scale) <= 1e-10:
        return None
    return {**copy.deepcopy(candidate), "orientation": "clockwise" if product > 0 else "counter_clockwise"}


def double_from_geometry(candidate, mol):
    sides = [[c for c in candidate[f"{side}_carriers"] if "atom_index" in c] for side in ("left", "right")]
    if not all(sides):
        return None
    conformer = mol.GetConformer()
    point = lambda i: tuple(conformer.GetAtomPosition(i))
    left, right = point(candidate["left_atom_index"]), point(candidate["right_atom_index"])
    axis = vector(right, left)
    normals = [cross(axis, vector(point(sides[0][0]["atom_index"]), left)),
               cross(axis, vector(point(sides[1][0]["atom_index"]), right))]
    scale = (dot(normals[0], normals[0]) * dot(normals[1], normals[1])) ** 0.5
    if scale == 0 or abs(dot(*normals) / scale) <= 1e-10:
        return None
    result = {k: copy.deepcopy(v) for k, v in candidate.items() if k not in {"left_carriers", "right_carriers"}}
    result.update(left_carrier=sides[0][0], right_carrier=sides[1][0],
                  orientation="together" if dot(*normals) > 0 else "opposite")
    return result


def add_sulfur_candidates(candidates, mol):
    present = {c["center_atom_index"] for c in candidates if c["type"] == "tetrahedral"}
    for atom in mol.GetAtoms():
        orders = [b.GetBondTypeAsDouble() for b in atom.GetBonds()]
        if (atom.GetAtomicNum() in (16, 34) and atom.GetDegree() == 3
                and atom.GetTotalNumHs() == 0 and all(o in (1, 2) for o in orders)
                and orders.count(2) <= 1 and atom.GetIdx() not in present):
            candidates.append({"type": "tetrahedral", "center_atom_index": atom.GetIdx(),
                               "carriers": [{"atom_index": i} for i in sorted(n.GetIdx() for n in atom.GetNeighbors())]
                                           + [{"implicit_lone_pair": True}]})
    candidates.sort(key=lambda c: (0, c["center_atom_index"]) if c["type"] == "tetrahedral"
                    else (1, c["center_bond_index"]))


def geometric_candidates(mol, Chem):
    """Local spatial candidates, independently enumerated on RDKit's graph."""
    result = []
    for atom in mol.GetAtoms():
        hydrogens = atom.GetTotalNumHs()
        if (atom.GetAtomicNum() != 1 and hydrogens <= 1 and atom.GetDegree() + hydrogens == 4
                and all(b.GetBondType() == Chem.BondType.SINGLE for b in atom.GetBonds())):
            carriers = [{"atom_index": i} for i in sorted(n.GetIdx() for n in atom.GetNeighbors())]
            if hydrogens:
                carriers.append({"implicit_hydrogen": True})
            result.append({"type": "tetrahedral", "center_atom_index": atom.GetIdx(), "carriers": carriers})
    rings = mol.GetRingInfo().BondRings()
    for bond in mol.GetBonds():
        left, right = sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        if (bond.GetBondType() != Chem.BondType.DOUBLE or bond.GetIsAromatic()
                or (mol.GetAtomWithIdx(left).GetIsAromatic() and mol.GetAtomWithIdx(right).GetIsAromatic())
                or any(bond.GetIdx() in ring and len(ring) < 8 for ring in rings)
                or (bond.IsInRing() and any(mol.GetAtomWithIdx(i).GetAtomicNum() != 6 for i in (left, right)))):
            continue
        def carriers(endpoint, other):
            atom = mol.GetAtomWithIdx(endpoint)
            values = [{"atom_index": n.GetIdx()} for n in atom.GetNeighbors()
                      if n.GetIdx() != other and mol.GetBondBetweenAtoms(endpoint, n.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
            values.sort(key=carrier_key)
            if atom.GetTotalNumHs() == 1:
                values.append({"implicit_hydrogen": True})
            return values
        a, b = carriers(left, right), carriers(right, left)
        if a and b:
            result.append({"type": "double_bond", "center_bond_index": bond.GetIdx(),
                           "left_atom_index": left, "right_atom_index": right,
                           "left_carriers": a, "right_carriers": b})
    add_sulfur_candidates(result, mol)
    return result


def migrate_record(representation, perception, source, format_, feature, Chem):
    mol = graph(source, format_, Chem)
    if perception.get("report", {}).get("issues") or perception.get("report", {}).get("is_ok") is False:
        raise ValueError("archived perception issues require explicit migration review")
    candidates = copy.deepcopy(perception.get("report", {}).get("candidates", []))
    recovered_perception = False
    if perception["status"] == "normalization_or_perception_error":
        sanitized = Chem.Mol(mol)
        if Chem.SanitizeMol(sanitized, catchErrors=True) == Chem.SanitizeFlags.SANITIZE_NONE:
            if sanitized.GetNumAtoms() != mol.GetNumAtoms() or sanitized.GetNumBonds() != mol.GetNumBonds():
                raise ValueError("reference cleanup changed graph identity")
            candidates = geometric_candidates(sanitized, Chem)
            recovered_perception = True
    add_sulfur_candidates(candidates, mol)
    prior = representation.get("stereo_elements", []) + [
        e for e in perception.get("report", {}).get("assembled_elements", [])
        if e.get("source") != "coordinates_3d"]
    elements, origins, focuses = [], [], set()
    def add(element, origin):
        focus = (element["type"], element.get("center_atom_index", element.get("center_bond_index")))
        normalized = canonical_element(element, mol)
        if focus in focuses:
            previous = next(e for e in elements if (e["type"], e.get("center_atom_index", e.get("center_bond_index"))) == focus)
            normalized["index"] = previous["index"]
            if normalized != previous:
                raise ValueError(f"conflicting archived assertions at {focus}")
            return
        normalized["index"] = len(elements)
        elements.append(normalized)
        origins.append(origin)
        focuses.add(focus)
    for element in prior:
        add(element, element["source"])
    if format_ != "smiles":
        for candidate in candidates:
            if candidate["type"] == "double_bond" and ("double_bond", candidate["center_bond_index"]) not in focuses:
                element = double_from_geometry(candidate, mol)
                if element is not None:
                    add(element, format_)
    assembled = copy.deepcopy(elements)
    created = []
    if feature == "stereo.perception" and format_ != "smiles":
        for candidate in candidates:
            if candidate["type"] == "tetrahedral" and ("tetrahedral", candidate["center_atom_index"]) not in focuses:
                element = tetra_from_geometry(candidate, mol)
                if element is not None:
                    before = len(elements)
                    add(element, "coordinates_3d")
                    if len(elements) > before:
                        created.append(before)
    legacy = representation if feature == "stereo.representation" else perception
    result = {key: copy.deepcopy(legacy[key]) for key in ("record_index", "status", "title", "atom_count", "bond_count")}
    result["document"] = source_document(source, format_, representation.get("stereo_bond_marks", []), elements, origins)
    if feature == "stereo.perception" and recovered_perception:
        result["status"] = "ok"
    if result["status"] != "ok":
        if result["status"] == "normalization_or_perception_error":
            result["status"] = "perception_error"
        return result
    result.update(stereo_elements=elements, stereo_groups=copy.deepcopy(legacy.get("stereo_groups", representation["stereo_groups"])))
    if feature == "stereo.perception":
        result["report"] = {"is_ok": True, "candidates": candidates, "issues": [],
                            "assembled_elements": assembled, "created_element_indices": created}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", action="append")
    parser.add_argument("--fixture", help="Republish one pinned fixture within one selected corpus.")
    parser.add_argument("--write", action="store_true", help="Publish independently migrated goldens and metadata.")
    args = parser.parse_args()
    if args.fixture and (not args.corpus or len(args.corpus) != 1):
        parser.error("--fixture requires exactly one --corpus")
    spec = importlib.util.spec_from_file_location("rdkit_reference", ROOT / "benchmarks/reference/rdkit/run_feature.py")
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    rdkit = reference.import_rdkit()
    if rdkit["version"] != VERSION:
        raise SystemExit(f"requires RDKit {VERSION}; found {rdkit['version']}")
    changes = []
    for corpus_dir in sorted((ROOT / "benchmarks/corpora").iterdir()):
        if args.corpus and corpus_dir.name not in args.corpus:
            continue
        for feature in FEATURES:
            manifest_path = corpus_dir / "features" / f"{feature}.toml"
            if not manifest_path.exists():
                continue
            manifest = reference.read_manifest(manifest_path)
            for fixture in manifest["fixtures"]:
                if args.fixture and fixture != args.fixture:
                    continue
                path = corpus_dir / "golden" / feature / f"{reference.slugify_fixture(fixture)}.json.gz"
                old_path = archive(path, corpus_dir.name, feature, args.write)
                old = read_golden(old_path)
                fixture_path = corpus_dir / fixture
                if digest(fixture_path.read_bytes()) != old["input_sha256"]:
                    raise ValueError(f"source drift: {fixture_path}")
                if feature == "io.smiles.isomeric":
                    new = reference.generate_document(feature, corpus_dir.name, fixture, fixture_path, rdkit)
                else:
                    pairs = []
                    for paired_feature in FEATURES[:2]:
                        paired_path = corpus_dir / "golden" / paired_feature / path.name
                        if paired_path.exists():
                            pairs.append(read_golden(archive(paired_path, corpus_dir.name, paired_feature, args.write)))
                        elif feature == "stereo.perception" and not old["expected"]["records"][0].get("stereo_elements"):
                            pairs.append(copy.deepcopy(old))
                        else:
                            raise ValueError(f"missing paired reviewed golden: {paired_path}")
                    rep = {r["record_index"]: r for r in pairs[0]["expected"]["records"]}
                    perc = {r["record_index"]: r for r in pairs[1]["expected"]["records"]}
                    sources = source_records(fixture_path)
                    format_ = "smiles" if fixture_path.suffix in {".smi", ".smiles", ".txt"} else "molfile_v2000"
                    new = copy.deepcopy(old)
                    new["expected"]["records"] = [migrate_record(rep[i], perc[i], sources[i], format_, feature, rdkit["Chem"])
                                                       for i in rep]
                    new["reference"] = {"tool": "reviewed-stereo-schema", "version": f"v2 + RDKit {VERSION}", "runtime_dependency": False}
                new["schema_version"] = 2
                archived_path = ARCHIVE / corpus_dir.name / feature / path.name
                new["migration"] = {"from_schema": 1, "archived_golden": str(archived_path.relative_to(ROOT)).replace("\\", "/"),
                                    "archived_sha256": digest(old_path.read_bytes()),
                                    "method": "source-graph-boundary-v2", "reference_version": VERSION}
                compressed = gzip.compress(json_bytes(new), mtime=0)
                changes.append({"path": str(path.relative_to(ROOT)).replace("\\", "/"),
                                "old_sha256": digest(old_path.read_bytes()), "new_sha256": digest(compressed),
                                "old_records": len(old["expected"]["records"]), "new_records": len(new["expected"]["records"])})
                if feature == "io.smiles.isomeric":
                    by_index = {record["record_index"]: record for record in new["expected"]["records"]}
                    overlaps = []
                    for record in old["expected"]["records"]:
                        differences = changed_paths(record, by_index[record["record_index"]])
                        if differences:
                            overlaps.append({"record_index": record["record_index"], "changed_paths": differences})
                    changes[-1]["overlap_changed_records"] = overlaps
                if args.write:
                    path.write_bytes(compressed)
                print(f"{corpus_dir.name} {feature} {fixture}: {changes[-1]['old_records']} -> {changes[-1]['new_records']} records", flush=True)
            if args.write and (not args.fixture or args.fixture in manifest["fixtures"]):
                text = manifest_path.read_text()
                import re
                version = f"RDKit {VERSION}" if feature == "io.smiles.isomeric" else f"v2 + RDKit {VERSION}"
                text = re.sub(r'^reference_version = .*$', f'reference_version = "{version}"', text, flags=re.M)
                if feature != "io.smiles.isomeric":
                    text = re.sub(r'^reference_tool = .*$', 'reference_tool = "reviewed-stereo-schema"', text, flags=re.M)
                text = text.replace("Goldens compare molecules' semantic stereo elements and source marks, not bytewise SMILES spelling.",
                                    "Schema 2 asserts source-document evidence separately from canonical molecular stereo elements; schema 1 reference bytes and hashes are archived.")
                text = text.replace("Goldens compare semantic stereo elements, stereo groups, and source bond marks, not bytewise SMILES spelling.",
                                    "Schema 2 asserts canonical stereo elements and groups plus separately retained source-document evidence; schema 1 bytes and hashes are archived.")
                text = text.replace("Goldens compare semantic stereo elements, groups, and source bond marks, not bytewise SDF spelling or parser-local IDs.",
                                    "Schema 2 asserts canonical stereo elements and groups plus separately retained source-document evidence; schema 1 bytes and hashes are archived.")
                text = text.replace("Goldens compare semantic stereo perception reports and resulting stereo elements, not bytewise SMILES spelling.",
                                    "Schema 2 uses retained Model coordinates for explicit stereo inference and asserts source elements, created elements, candidates and groups.")
                text = text.replace("Goldens compare semantic stereo reports and resulting stereo elements after normalization and default perception, not bytewise SDF spelling.",
                                    "Schema 2 uses retained Model coordinates for explicit stereo inference and asserts source elements, created elements, candidates and groups.")
                text = text.replace("normalized/perceived records with source stereo syntax", "all source records, including non-stereo controls")
                text = text.replace("Reference RDKit goldens compare semantic output after noncanonical isomeric SMILES write and reparse.",
                                    "The decoded isomeric roundtrip must preserve independently parsed source chemistry and hydrogen declarations; RDKit emission is separately retained and checked for whole-graph stereo identity.")
                if feature == "io.smiles.isomeric":
                    text = text.replace("Goldens compare semantic output after write and reparse, including normalized/perceived graph semantics plus CIP descriptor-bearing stereo semantics.",
                                        "Decoded output is compared with source graph, hydrogen declarations and CIP stereo; RDKit raw emission and its whole-graph equivalence check remain separate reference evidence.")
                    text = text.replace("source chemistry and hydrogen declarations;", "source chemistry and hydrogen declarations subject to mandatory charge/bracket syntax;")
                    text = text.replace("source graph, hydrogen declarations and CIP stereo;", "source graph, hydrogen declarations subject to mandatory charge/bracket syntax, and CIP stereo;")
                manifest_path.write_text(text)
    report = {"schema_version": 2, "rdkit": VERSION, "method": "source-graph-boundary-v2", "fixtures": changes}
    if args.fixture and not changes:
        raise ValueError("selected fixture is not pinned by an affected benchmark manifest")
    if args.write:
        report_path = Path(__file__).parent / "migration-v2.json"
        if report_path.exists():
            previous = json.loads(report_path.read_text())
            report["fixtures"] = merged_provenance(previous, changes)
        report_path.write_bytes(json_bytes(report))


if __name__ == "__main__":
    main()
