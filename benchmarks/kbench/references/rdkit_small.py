"""RDKit reference observer for small-molecule tasks.

Emits the same facts as `benchmarks/observer/src/small.rs`, keyed by source
atom order, which RDKit keeps when hydrogens are not removed. Every value is
computed by RDKit; nothing here reimplements chemistry.

Reading: parse without sanitization (no normalization, explicit hydrogens
kept), then sanitize a copy. Represented atoms, bonds and hydrogen counts
come from the unsanitized molecule; perceived results come from the
sanitized copy, with stereo assigned from coordinates or bond directions and
non-stereogenic assertions cleaned, as Kekule does before comparison.
"""
from __future__ import annotations

import io
from functools import lru_cache

from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import rdCIPLabeler

RDLogger.DisableLog("rdApp.*")

TOOL = f"RDKit {rdBase.rdkitVersion}"
MAX_MATCHES = 100_000
# Resonance enumeration steps before RDKit gives up on a fragment; see
# WorkLimit. On the small dataset the median fragment takes 7 steps and the
# largest that finishes 2,770; a few porphyrins and polycations never finish.
RESONANCE_WORK_LIMIT = 5_000
# RDKit's strict rotatable-bond definition (rdMolDescriptors, NumRotatableBonds Strict).
STRICT_ROTATABLE = Chem.MolFromSmarts(
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
ORDERS = {
    Chem.BondType.ZERO: 0,
    Chem.BondType.SINGLE: 1,
    Chem.BondType.DOUBLE: 2,
    Chem.BondType.TRIPLE: 3,
    Chem.BondType.QUADRUPLE: 4,
    Chem.BondType.AROMATIC: "aromatic",
    Chem.BondType.DATIVE: "dative",
}


class Failure(Exception):
    """A reference failure for one input, with a stable kind."""

    def __init__(self, kind: str, message: str):
        super().__init__(message)
        self.kind = kind


def read(fmt: str, text: str) -> Chem.Mol:
    """Unsanitized molecule with explicit hydrogens kept and properties cached."""
    if fmt == "smiles":
        params = Chem.SmilesParserParams()
        params.sanitize = False
        params.removeHs = False
        mol = Chem.MolFromSmiles((text.splitlines() or [""])[0].strip(), params)
    elif fmt == "mol":
        mol = Chem.MolFromMolBlock(text, sanitize=False, removeHs=False, strictParsing=True)
    elif fmt == "sdf":
        supplier = Chem.ForwardSDMolSupplier(
            io.BytesIO(text.encode("utf-8")), sanitize=False, removeHs=False, strictParsing=True
        )
        records = list(supplier)
        if len(records) != 1:
            raise Failure("parse", f"expected one SDF record, found {len(records)}")
        mol = records[0]
    else:
        raise Failure("request", f"unknown format {fmt!r}")
    if mol is None:
        raise Failure("parse", "RDKit could not read the input")
    mol.UpdatePropertyCache(strict=False)
    if fmt == "smiles":
        Chem.SetBondStereoFromDirections(mol)
    return mol


@lru_cache(maxsize=2)
def prepared(fmt: str, text: str) -> tuple[Chem.Mol, Chem.Mol]:
    """The unsanitized and sanitized molecule of one input, read once per
    worker for all of its tasks. Callers must not modify either."""
    source = read(fmt, text)
    return source, sanitized(source, fmt)


def is_3d(mol: Chem.Mol) -> bool:
    return mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()


def sanitized(source: Chem.Mol, fmt: str) -> Chem.Mol:
    mol = Chem.Mol(source)
    try:
        Chem.SanitizeMol(mol)
    except Exception as error:  # RDKit raises several sanitization exception types
        raise Failure("sanitize", str(error)) from None
    if is_3d(mol):
        Chem.AssignStereochemistryFrom3D(mol)
    else:
        if fmt in ("mol", "sdf"):
            Chem.AssignChiralTypesFromBondDirs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol


def ends(bond: Chem.Bond) -> tuple[int, int]:
    a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    if bond.GetBondType() == Chem.BondType.DATIVE:
        return a, b
    return (a, b) if a < b else (b, a)


def atom_context(atom: Chem.Atom) -> str:
    parts = [atom.GetSymbol()]
    if atom.GetIsAromatic():
        parts.append("ar")
    if atom.GetFormalCharge():
        parts.append(f"{atom.GetFormalCharge():+d}")
    return "|".join(parts)


def bond_context(mol: Chem.Mol, bond: Chem.Bond) -> str:
    a, b = ends(bond)
    order = ORDERS.get(bond.GetBondType(), str(bond.GetBondType()))
    order = "ar" if bond.GetIsAromatic() else order
    return f"{mol.GetAtomWithIdx(a).GetSymbol()}-{mol.GetAtomWithIdx(b).GetSymbol()}|{order}"


def normalized(source: Chem.Mol, mol: Chem.Mol) -> bool:
    """Whether sanitization's cleanup rewrote charges or localized bonds."""
    for a, b in zip(source.GetAtoms(), mol.GetAtoms()):
        if a.GetFormalCharge() != b.GetFormalCharge():
            return True
    for a, b in zip(source.GetBonds(), mol.GetBonds()):
        if not (a.GetIsAromatic() or b.GetIsAromatic()) and a.GetBondType() != b.GetBondType():
            return True
    return False


def parity(values) -> int:
    return sum(a > b for i, a in enumerate(values) for b in values[i + 1 :]) % 2


def reference_carrier(mol: Chem.Mol, center: int, other: int) -> int:
    neighbors = (a.GetIdx() for a in mol.GetAtomWithIdx(center).GetNeighbors() if a.GetIdx() != other)
    return min(neighbors, default=-1)


def stereo_facts(mol: Chem.Mol, facts: list) -> None:
    for atom in mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            raise Failure("stereo", f"no common representation for {tag}")
        carriers = [a.GetIdx() for a in atom.GetNeighbors()]
        if len(carriers) == 3:
            carriers.append(-1 if atom.GetTotalNumHs() else -2)
        if len(carriers) != 4:
            raise Failure("stereo", "tetrahedral centre without four carriers")
        value = int(tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW) ^ parity(carriers)
        focus = [atom.GetIdx()]
        facts.append(["stereo", "tetrahedral", focus, "carriers", sorted(carriers)])
        facts.append(["stereo", "tetrahedral", focus, "parity", value])
    for bond in mol.GetBonds():
        tag = bond.GetStereo()
        if tag == Chem.BondStereo.STEREONONE:
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        reference = [reference_carrier(mol, a, b), reference_carrier(mol, b, a)]
        if tag in (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW):
            # RDKit defines atropisomer orientation from the lowest-numbered
            # neighbour at each end, independent of bond direction.
            kind, value = "axis", int(tag == Chem.BondStereo.STEREOATROPCCW)
        elif tag == Chem.BondStereo.STEREOANY:
            kind, value = "double_bond", None
        elif tag in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOE,
                     Chem.BondStereo.STEREOZ):
            selected = list(bond.GetStereoAtoms())
            if len(selected) != 2:
                raise Failure("stereo", "double-bond stereo without reference atoms")
            opposite = tag in (Chem.BondStereo.STEREOTRANS, Chem.BondStereo.STEREOE)
            kind = "double_bond"
            value = int(opposite) ^ int(selected[0] != reference[0]) ^ int(selected[1] != reference[1])
        else:
            raise Failure("stereo", f"no common representation for {tag}")
        focus = [a, b]
        if a > b:
            focus.reverse()
            reference.reverse()
        facts.append(["stereo", kind, focus, "carriers", reference])
        facts.append(["stereo", kind, focus, "parity", value])
    kinds = {
        Chem.StereoGroupType.STEREO_ABSOLUTE: "absolute",
        Chem.StereoGroupType.STEREO_AND: "and",
        Chem.StereoGroupType.STEREO_OR: "or",
    }
    for group in mol.GetStereoGroups():
        members = [["tetrahedral", [a.GetIdx()]] for a in group.GetAtoms()]
        members += [["axis", sorted([b.GetBeginAtomIdx(), b.GetEndAtomIdx()])] for b in group.GetBonds()]
        members.sort(key=str)
        facts.append(["group", members, "kind", kinds[group.GetGroupType()]])


def parse(source: Chem.Mol, mol: Chem.Mol) -> tuple[list, dict]:
    facts = []
    for atom in source.GetAtoms():
        i = atom.GetIdx()
        hydrogens = None if atom.GetAtomicNum() == 0 else atom.GetTotalNumHs(includeNeighbors=True)
        # Radicals implied by an explicit valence are assigned during sanitization.
        radicals = mol.GetAtomWithIdx(i).GetNumRadicalElectrons()
        facts += [
            ["atom", i, "element", atom.GetSymbol()],
            ["atom", i, "isotope", atom.GetIsotope() or None],
            ["atom", i, "charge", atom.GetFormalCharge()],
            ["atom", i, "radicals", radicals],
            ["atom", i, "hydrogens", hydrogens],
            ["atom", i, "map", atom.GetAtomMapNum() or None],
        ]
    aromatic = []
    for bond in source.GetBonds():
        a, b = ends(bond)
        order = ORDERS.get(bond.GetBondType(), str(bond.GetBondType()))
        if order == "aromatic":
            aromatic.append([a, b])
        facts.append(["bond", a, b, "order", order])
    for atoms in Chem.GetMolFrags(source):
        facts.append(["component", min(atoms), "atoms", sorted(atoms)])
    stereo_facts(mol, facts)
    context = {
        "atoms": {str(a.GetIdx()): atom_context(a) for a in mol.GetAtoms()},
        "bonds": {"-".join(map(str, ends(b))): bond_context(mol, b) for b in mol.GetBonds()},
        "source_aromatic": aromatic,
        "normalized": normalized(source, mol),
    }
    return facts, context


def fragments(mol: Chem.Mol):
    """Each connected component with the map from fragment to source atom index."""
    mapping: list = []
    parts = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False, fragsMolAtomMapping=mapping)
    return list(zip(parts, mapping))


def rings(mol: Chem.Mol) -> list:
    facts = []
    info = mol.GetRingInfo()
    for atom in mol.GetAtoms():
        facts.append(["atom", atom.GetIdx(), "in_ring", info.NumAtomRings(atom.GetIdx()) > 0])
    for bond in mol.GetBonds():
        facts.append(["bond", *ends(bond), "in_ring", info.NumBondRings(bond.GetIdx()) > 0])
    sizes: dict[int, int] = {}
    for ring in Chem.GetSymmSSSR(mol):
        sizes[len(ring)] = sizes.get(len(ring), 0) + 1
    facts += [["rings", size, "count", count] for size, count in sorted(sizes.items())]
    return facts


def aromaticity(mol: Chem.Mol, model: str) -> list:
    if model == "mdl":
        mol = Chem.Mol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MDL)
    elif model != "rdkit":
        raise Failure("request", f"unknown model {model!r}")
    facts = [["atom", a.GetIdx(), "aromatic", a.GetIsAromatic()] for a in mol.GetAtoms()]
    facts += [["bond", *ends(b), "aromatic", b.GetIsAromatic()] for b in mol.GetBonds()]
    return facts


def conjugation(mol: Chem.Mol) -> list:
    return [["bond", *ends(b), "conjugated", b.GetIsConjugated()] for b in mol.GetBonds()]


class WorkLimit(Chem.ResonanceMolSupplierCallback):
    """Cancels enumeration after a fixed number of progress steps. Single-
    threaded enumeration reports the same steps on every run, so the limit is
    deterministic, unlike a timeout."""

    def __init__(self):
        super().__init__()
        self.steps = 0

    def __call__(self) -> bool:
        self.steps += 1
        return self.steps < RESONANCE_WORK_LIMIT


def resonance(mol: Chem.Mol) -> list:
    facts = []
    for fragment, mapping in fragments(mol):
        supplier = Chem.ResonanceMolSupplier(fragment, 0, 1000)
        supplier.SetProgressCallback(WorkLimit())
        len(supplier)  # enumerates
        if supplier.WasCanceled():
            raise Failure("resonance", f"enumeration work limit exceeded ({RESONANCE_WORK_LIMIT} progress steps)")
        for group in range(supplier.GetNumConjGrps()):
            atoms = sorted(mapping[a.GetIdx()] for a in fragment.GetAtoms()
                           if supplier.GetAtomConjGrpIdx(a.GetIdx()) == group)
            bonds = sorted(
                sorted([mapping[b.GetBeginAtomIdx()], mapping[b.GetEndAtomIdx()]])
                for b in fragment.GetBonds() if supplier.GetBondConjGrpIdx(b.GetIdx()) == group
            )
            facts.append(["resonance", "group", atoms[0], "atoms", atoms])
            facts.append(["resonance", "group", atoms[0], "bonds", bonds])
        contributors = set()
        count = 0
        for form in supplier:
            count += 1
            charges = sorted((mapping[a.GetIdx()], a.GetFormalCharge()) for a in form.GetAtoms()
                             if a.GetFormalCharge())
            bonds = []
            for bond in form.GetBonds():
                order = ORDERS.get(bond.GetBondType(), str(bond.GetBondType()))
                if order != 1:
                    a, b = mapping[bond.GetBeginAtomIdx()], mapping[bond.GetEndAtomIdx()]
                    if order != "dative" and a > b:
                        a, b = b, a
                    bonds.append(((a, b), json_order(order)))
            bonds.sort()
            contributors.add(
                ",".join(f"{i}:{q}" for i, q in charges) + "|" + ",".join(f"{a}-{b}:{o}" for (a, b), o in bonds)
            )
        first = min(mapping)
        facts.append(["resonance", "contributors", first, "count", count])
        facts.append(["resonance", "contributors", first, "set", sorted(contributors)])
    return facts


def json_order(order) -> str:
    """Bond order as the Rust observer prints it inside contributor strings."""
    return f'"{order}"' if isinstance(order, str) else str(order)


def candidates(mol: Chem.Mol) -> list:
    atom_kinds = {
        "Atom_Tetrahedral": "tetrahedral",
        "Atom_SquarePlanar": "square_planar",
        "Atom_TrigonalBipyramidal": "trigonal_bipyramidal",
        "Atom_Octahedral": "octahedral",
    }
    bond_kinds = {"Bond_Double": "double_bond", "Bond_Atropisomer": "axis", "Bond_Cumulene_Even": "cumulene"}
    mol = Chem.Mol(mol)
    facts = []
    for info in Chem.FindPotentialStereo(mol):
        kind = str(info.type)
        if kind in atom_kinds:
            facts.append(["candidate", atom_kinds[kind], [info.centeredOn], "present", True])
        elif kind in bond_kinds:
            bond = mol.GetBondWithIdx(info.centeredOn)
            facts.append(["candidate", bond_kinds[kind], sorted(ends(bond)), "present", True])
        else:
            facts.append(["candidate", kind, [info.centeredOn], "present", True])
    return facts


def cip(mol: Chem.Mol) -> list:
    mol = Chem.Mol(mol)
    try:
        rdCIPLabeler.AssignCIPLabels(mol, maxRecursiveIterations=1_000_000)
    except Exception as error:
        raise Failure("cip", str(error)) from None
    facts = []
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            facts.append(["cip", "atom", [atom.GetIdx()], "label", atom.GetProp("_CIPCode")])
    for bond in mol.GetBonds():
        if bond.HasProp("_CIPCode"):
            facts.append(["cip", "bond", sorted(ends(bond)), "label", bond.GetProp("_CIPCode")])
    return facts


def symmetry(mol: Chem.Mol) -> list:
    facts = []
    for fragment, mapping in fragments(mol):
        ranks = Chem.CanonicalRankAtoms(
            fragment, breakTies=False, includeChirality=False, includeIsotopes=True, includeAtomMaps=True
        )
        classes: dict[int, list] = {}
        for index, rank in enumerate(ranks):
            classes.setdefault(rank, []).append(mapping[index])
        for atoms in classes.values():
            atoms.sort()
            facts.append(["class", atoms[0], "atoms", atoms])
    return facts


def formula(mol: Chem.Mol) -> list:
    facts = []
    for fragment, mapping in fragments(mol):
        counts: dict[tuple, int] = {}
        charge = 0
        for atom in fragment.GetAtoms():
            key = (atom.GetSymbol(), atom.GetIsotope() or None)
            counts[key] = counts.get(key, 0) + 1
            hydrogens = atom.GetTotalNumHs(includeNeighbors=False)
            if hydrogens:
                counts[("H", None)] = counts.get(("H", None), 0) + hydrogens
            charge += atom.GetFormalCharge()
        first = min(mapping)
        for (symbol, isotope), count in sorted(counts.items(), key=lambda item: (item[0][0], item[0][1] or 0)):
            facts.append(["formula", first, "count", symbol, isotope, count])
        facts.append(["formula", first, "charge", charge])
    return facts


def rotatable(mol: Chem.Mol) -> list:
    """Strict rotatable bonds after removing explicit hydrogens, like for like
    with Kekule's heavy-atom definition."""
    mol = Chem.Mol(mol)
    for atom in mol.GetAtoms():
        atom.SetIntProp("source", atom.GetIdx())
    heavy = Chem.RemoveHs(mol)
    source = [atom.GetIntProp("source") for atom in heavy.GetAtoms()]
    bonds = set()
    for a, b in heavy.GetSubstructMatches(STRICT_ROTATABLE, uniquify=True, maxMatches=MAX_MATCHES):
        bonds.add(tuple(sorted((source[a], source[b]))))
    return [["bond", a, b, "rotatable", True] for a, b in sorted(bonds)]


def smarts(mol: Chem.Mol, queries: list[str]) -> list:
    facts = []
    for index, text in enumerate(queries):
        query = Chem.MolFromSmarts(text)
        if query is None:
            facts.append(["query", index, "error", "parse"])
            continue
        matches = sorted(
            list(match)
            for match in mol.GetSubstructMatches(query, uniquify=False, useChirality=True, maxMatches=MAX_MATCHES)
        )
        facts.append(["query", index, "count", len(matches)])
        facts.append(["query", index, "matches", matches])
    return facts


def observe(task: str, fmt: str, text: str, options: dict) -> dict:
    """One reference observation: `{"facts": [...]}`, plus `context` for parse.
    `parse` with `coordinates` also reports each atom's source position."""
    source, mol = prepared(fmt, text)
    if task == "parse":
        facts, context = parse(source, mol)
        if options.get("coordinates") and source.GetNumConformers():
            conformer = source.GetConformer()
            for atom in source.GetAtoms():
                point = conformer.GetAtomPosition(atom.GetIdx())
                facts.append(["atom", atom.GetIdx(), "xyz", [point.x, point.y, point.z]])
        return {"facts": facts, "context": context}
    functions = {
        "rings": lambda: rings(mol),
        "aromaticity": lambda: aromaticity(mol, options.get("model", "rdkit")),
        "conjugation": lambda: conjugation(mol),
        "resonance": lambda: resonance(mol),
        "stereo.candidates": lambda: candidates(mol),
        "stereo.cip": lambda: cip(mol),
        "symmetry": lambda: symmetry(mol),
        "formula": lambda: formula(mol),
        "rotatable": lambda: rotatable(mol),
        "smarts": lambda: smarts(mol, options["queries"]),
    }
    if task not in functions:
        raise Failure("request", f"unknown task small.{task}")
    return {"facts": functions[task]()}
