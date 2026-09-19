#!/usr/bin/env python3
"""Generate RDKit-backed golden data for small-molecule benchmark features."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any
from reference.rdkit import molecule
from reference.rdkit.molecule import atom_json, graph as molecular_graph, records as molecular_records


RDKIT_STRICT_ROTATABLE_BOND_SMARTS = (
    "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)"
    "&!$(C([CH3])([CH3])[CH3])&!$([CH3])"
    "&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])"
    "&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])"
    "&!$([CD3](=[N+])-!@[#7!D1])"
    "&!$([#7!D1]-!@[CD3]=[N+])]"
    "-,:;!@"
    "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)"
    "&!$(C([CH3])([CH3])[CH3])&!$([CH3])]"
)

SUBSTRUCTURE_QUERIES = (Path(__file__).parents[2] / 'queries.smarts').read_text().splitlines()


def import_rdkit() -> dict[str, Any]:
    try:
        from rdkit import Chem, RDLogger, rdBase
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError as error:
        raise SystemExit(
            "RDKit is not importable. Create the environment from "
            "benchmarks/reference/rdkit/environment.yml before generating goldens."
        ) from error
    RDLogger.DisableLog("rdApp.*")
    return {
        "Chem": Chem,
        "Descriptors": Descriptors,
        "rdMolDescriptors": rdMolDescriptors,
        "version": rdBase.rdkitVersion,
    }


def evaluate(feature_id, fixture_path, rdkit):
    # All available components are evaluated. Reference errors stay errors.
    Chem = rdkit['Chem']
    if feature_id == 'query.smarts':
        return {'records': smarts_query_records(fixture_path, Chem)}
    if feature_id == 'stereo.cip':
        records = []
        for index, (title, mol, fields) in enumerate(molecular_records(fixture_path)) :
            records.append({'record_index':index,'title':title,'status':'ok','mol':mol})
        return {'records':[stereo_cip_record(record, Chem) for record in records]}
    records = []
    for title, mol, fields in molecular_records(fixture_path):
        fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        for fragment in fragments:
            fragment.RemoveAllConformers()
            records.append({'record_index':len(records),'title':title,'status':'ok','mol':fragment,
                'properties':fields})
    functions = {
        'descriptor.molecular': lambda r: molecular_descriptor_record(r, Chem, rdkit['Descriptors']),
        'descriptor.rotatable-bonds.rdkit-strict': lambda r: rotatable_bond_record(r, Chem, rdkit['rdMolDescriptors'], Chem.MolFromSmarts(RDKIT_STRICT_ROTATABLE_BOND_SMARTS)),
        'algo.substructure.vf2': lambda r: substructure_record(r, Chem),
        'algo.rings.fast': ring_record,
        'algo.rings.sssr': ring_set_record,
        'algo.valence.rdkit-like': valence_record,
        'algo.aromaticity.rdkit-like': aromaticity_record,
        'algo.canonical-ranking': canonical_ranking_record,
        'chem.perception.default': perceived_atom_record,
        'chem.hydrogen-transforms': hydrogen_transform_record,
    }
    if feature_id not in functions:
        raise ValueError(f'no independent reference for {feature_id}')
    return {'records':[functions[feature_id](record) for record in records]}


def smarts_query_records(fixture_path: Path, Chem: Any) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index, raw_line in enumerate(fixture_path.read_text(encoding="utf-8").splitlines()):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        smarts, title = parts[0], parts[1] if len(parts) == 2 else ''
        query = Chem.MolFromSmarts(smarts)
        records.append(
            {
                "record_index": index,
                "status": "ok" if query is not None else "parse_error",
                "smarts": smarts,
                "title": title,
                "atom_count": query.GetNumAtoms() if query is not None else None,
                "bond_count": query.GetNumBonds() if query is not None else None,
            }
        )
    return records


def substructure_record(record: dict[str, Any], Chem: Any) -> dict[str, Any]:
    if record["status"] != "ok" or record["mol"] is None:
        return {
            "record_index": record["record_index"],
            "status": "parse_error",
            "title": record["title"],
            "queries": [],
        }
    mol = Chem.Mol(record["mol"])
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return {
            "record_index": record["record_index"],
            "status": "normalization_or_perception_error",
            "title": record["title"],
            "queries": [],
        }

    queries: list[dict[str, Any]] = []
    for smarts in SUBSTRUCTURE_QUERIES:
        query = Chem.MolFromSmarts(smarts)
        if query is None:
            raise RuntimeError(f"benchmark SMARTS did not parse in RDKit: {smarts}")
        matches = mol.GetSubstructMatches(query, uniquify=False, maxMatches=0)
        queries.append(
            {
                "smarts": smarts,
                "matches": [list(match) for match in sorted(matches)],
            }
        )
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "queries": queries,
    }


def molecular_descriptor_record(
    record: dict[str, Any], Chem: Any, Descriptors: Any
) -> dict[str, Any]:
    if record["status"] != "ok" or record["mol"] is None:
        return {
            "record_index": record["record_index"],
            "status": "parse_error",
            "title": record["title"],
        }
    mol = Chem.Mol(record["mol"])
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return {
            "record_index": record["record_index"],
            "status": "normalization_or_perception_error",
            "title": record["title"],
        }

    counts: dict[tuple[str, int | None], int] = {}
    charge = 0
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope() or None
        key = (atom.GetSymbol(), isotope)
        counts[key] = counts.get(key, 0) + 1
        hydrogens = atom.GetNumExplicitHs() + atom.GetNumImplicitHs()
        if hydrogens:
            hydrogen_key = ("H", None)
            counts[hydrogen_key] = counts.get(hydrogen_key, 0) + hydrogens
        charge += atom.GetFormalCharge()

    contains_carbon = any(symbol == "C" for symbol, _ in counts)

    def hill_key(term: tuple[str, int | None]) -> tuple[int, str, int]:
        symbol, isotope = term
        if contains_carbon and symbol == "C":
            element_key = (0, "")
        elif contains_carbon and symbol == "H":
            element_key = (1, "")
        else:
            element_key = (2 if contains_carbon else 0, symbol)
        return (*element_key, isotope or 0)

    terms = [
        {"element": symbol, "isotope": isotope, "count": count}
        for (symbol, isotope), count in sorted(counts.items(), key=lambda item: hill_key(item[0]))
    ]
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "formula": {"terms": terms, "formal_charge": charge},
        "average_mass_da": Descriptors.MolWt(mol),
        # RDKit ExactMolWt already applies its electron-mass correction.
        "monoisotopic_mass_da": Descriptors.ExactMolWt(mol),
    }


def rotatable_bond_record(
    record: dict[str, Any], Chem: Any, rdMolDescriptors: Any, query: Any
) -> dict[str, Any]:
    if record["status"] != "ok" or record["mol"] is None:
        return {
            "record_index": record["record_index"],
            "status": "parse_error",
            "title": record["title"],
        }

    mol = Chem.Mol(record["mol"])
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return {
            "record_index": record["record_index"],
            "status": "normalization_or_perception_error",
            "title": record["title"],
        }

    bond_endpoints = sorted(
        {
            tuple(sorted(match))
            for match in mol.GetSubstructMatches(query, uniquify=True, maxMatches=0)
        }
    )
    descriptor_count = rdMolDescriptors.CalcNumRotatableBonds(
        mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
    )
    if len(bond_endpoints) != descriptor_count:
        raise RuntimeError(
            "RDKit strict SMARTS/count disagreement for record "
            f"{record['record_index']}: {len(bond_endpoints)} endpoints versus "
            f"descriptor count {descriptor_count}"
        )

    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "count": descriptor_count,
        "bonds": [
            {"begin_atom_index": begin, "end_atom_index": end}
            for begin, end in bond_endpoints
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
    Chem.FastFindRings(mol)
    rings = mol.GetRingInfo()
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_in_ring": [atom.IsInRing() for atom in mol.GetAtoms()],
        "bond_in_ring": bond_values(mol, lambda bond: rings.NumBondRings(bond.GetIdx()) > 0),
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


def perceived_atom_record(record: dict[str, Any]) -> dict[str, Any]:
    sanitized = clone_and_sanitize(record["mol"]) if record["mol"] is not None else None
    if sanitized is None:
        return {
            "record_index": record["record_index"],
            "status": "normalization_or_perception_error",
            "title": record["title"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atoms": [basic_atom_json(atom) for atom in sanitized.GetAtoms()],
        "graph": molecular_graph(sanitized),
        "valence": [valence_atom_json(atom) for atom in sanitized.GetAtoms()],
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
        # Compare valence on normalized represented chemistry, as published by
        # Kekule. Cleanup does not run aromaticity or radical perception.
        Chem.Cleanup(prepared)
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


def hydrogen_transform_record(record: dict[str, Any]) -> dict[str, Any]:
    from rdkit import Chem

    mol = record["mol"]
    prepared = clone_and_sanitize(mol) if mol is not None else None
    if prepared is None:
        return {
            "record_index": record["record_index"],
            "status": "normalization_or_perception_error",
            "title": record["title"],
        }

    original_atom_count = prepared.GetNumAtoms()
    try:
        expanded = Chem.AddHs(prepared, addCoords=False)
    except Exception:
        return {
            "record_index": record["record_index"],
            "status": "add_error",
            "title": record["title"],
        }

    added_by_parent: dict[int, int] = {}
    for atom in list(expanded.GetAtoms())[original_atom_count:]:
        if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
            return {
                "record_index": record["record_index"],
                "status": "add_error",
                "title": record["title"],
            }
        parent = atom.GetNeighbors()[0].GetIdx()
        added_by_parent[parent] = added_by_parent.get(parent, 0) + 1

    try:
        collapsed = Chem.RemoveHs(expanded, sanitize=True)
    except Exception:
        return {
            "record_index": record["record_index"],
            "status": "remove_error",
            "title": record["title"],
        }

    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_count_after_add": expanded.GetNumAtoms(),
        "added_hydrogens_by_parent": [
            {"parent_atom_index": parent, "count": count}
            for parent, count in sorted(added_by_parent.items())
        ],
        "round_trip": molecular_graph(collapsed),
        "added_graph": molecular_graph(expanded),
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
            "status": "normalization_or_perception_error",
            "title": record["title"],
        }
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_aromatic": [atom.GetIsAromatic() for atom in sanitized.GetAtoms()],
        "bond_aromatic": bond_values(sanitized, lambda bond: bond.GetIsAromatic()),
    }


def canonical_ranking_record(record: dict[str, Any]) -> dict[str, Any]:
    from rdkit import Chem

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
            "status": "normalization_or_perception_error",
            "title": record["title"],
        }
    ranks = Chem.CanonicalRankAtoms(
        sanitized,
        breakTies=False,
        includeChirality=False,
        includeIsotopes=True,
        includeAtomMaps=True,
        includeChiralPresence=False,
    )
    classes: dict[int, list[int]] = {}
    for atom_index, rank in enumerate(ranks):
        classes.setdefault(int(rank), []).append(atom_index)
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "classes": sorted(classes.values()),
    }


def stereo_cip_record(record: dict[str, Any], Chem: Any) -> dict[str, Any]:
    failure = {"record_index": record["record_index"], "title": record["title"]}
    mol = record["mol"]
    if mol is None:
        return {**failure, "status": "parse_error"}
    prepared = clone_and_sanitize(mol)
    if prepared is None:
        return {**failure, "status": "sanitize_error"}
    try:
        Chem.AssignCIPLabels(prepared, maxRecursiveIterations=1_000_000)
    except Exception:
        return {**failure, "status": "cip_error"}
    atom_descriptors = cip_atom_descriptors(prepared)
    bond_descriptors = cip_bond_descriptors(prepared)
    # Kekule's topology groups connected components. Express source atom IDs in
    # that same component order, retaining every explicit hydrogen vertex.
    atom_index = {source: dense for dense, source in enumerate(
        atom for component in Chem.GetMolFrags(prepared) for atom in component)}
    for descriptor in atom_descriptors:
        descriptor['atom_index'] = atom_index[descriptor['atom_index']]
    for descriptor in bond_descriptors:
        ends = sorted([atom_index[descriptor['begin_atom_index']], atom_index[descriptor['end_atom_index']]])
        descriptor['begin_atom_index'], descriptor['end_atom_index'] = ends
    atom_descriptors.sort(key=lambda d: d['atom_index'])
    bond_descriptors.sort(key=lambda d: (d['begin_atom_index'], d['end_atom_index']))
    return {
        "record_index": record["record_index"],
        "status": "ok",
        "title": record["title"],
        "atom_count": prepared.GetNumAtoms(),
        "bond_count": prepared.GetNumBonds(),
        "atom_descriptors": atom_descriptors,
        "bond_descriptors": bond_descriptors,
    }


def cip_atom_descriptors(mol: Any) -> list[dict[str, Any]]:
    descriptors = []
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            descriptors.append(
                {
                    "atom_index": atom.GetIdx(),
                    "descriptor": atom.GetProp("_CIPCode"),
                }
            )
    descriptors.sort(key=lambda item: item["atom_index"])
    return descriptors


def cip_bond_descriptors(mol: Any) -> list[dict[str, Any]]:
    descriptors = []
    for bond in mol.GetBonds():
        if bond.HasProp("_CIPCode"):
            descriptors.append(
                {
                    "begin_atom_index": bond.GetBeginAtomIdx(),
                    "end_atom_index": bond.GetEndAtomIdx(),
                    "descriptor": bond.GetProp("_CIPCode"),
                }
            )
    descriptors.sort(
        key=lambda item: (
            item["begin_atom_index"],
            item["end_atom_index"],
            item["descriptor"],
        )
    )
    return descriptors


def clone_and_sanitize(mol: Any) -> Any | None:
    from rdkit import Chem

    cloned = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(cloned)
        Chem.AssignStereochemistry(cloned, cleanIt=True, force=True)
    except Exception:
        return None
    return cloned


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



def bond_values(mol, value):
    bonds=[]
    for bond in mol.GetBonds():
        ends=[bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()]
        if str(bond.GetBondType()) != 'DATIVE':
            ends.sort()
        bonds.append({'begin_atom_index':ends[0],'end_atom_index':ends[1],'value':value(bond)})
    return sorted(bonds,key=lambda b:(b['begin_atom_index'],b['end_atom_index']))
