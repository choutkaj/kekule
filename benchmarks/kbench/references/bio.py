"""gemmi, Biotite and mkdssp reference observer for biomolecule tasks.

Emits the facts of `benchmarks/observer/src/bio.rs`, keyed by `_atom_site.id`.
- CIF syntax, atom sites and superposition: gemmi.
- Covalent connectivity and residue classes: Biotite, which reads intra-residue
  bonds from the file's `_chem_comp_bond` and links from `_struct_conn` and
  standard polymer linkage.
- Secondary structure: mkdssp on the same single-model, single-conformer file
  Kekule reads, written by gemmi, so altloc policies cannot confound DSSP.
"""
from __future__ import annotations

import io
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import biotite
import biotite.structure as struc
import biotite.structure.info as info
import biotite.structure.io.pdbx as pdbx
import gemmi
import numpy as np


def mkdssp() -> str:
    found = shutil.which("mkdssp")
    if found is None:
        raise RuntimeError("mkdssp is not on PATH")
    return found


def mkdssp_version() -> str:
    try:
        return subprocess.run([mkdssp(), "--version"], capture_output=True, text=True).stdout.split("\n")[0].strip()
    except (RuntimeError, OSError):
        return "mkdssp unavailable"


TOOL = f"gemmi {gemmi.__version__}, Biotite {biotite.__version__}, {mkdssp_version()}"


class Failure(Exception):
    def __init__(self, kind: str, message: str):
        super().__init__(message)
        self.kind = kind


def fingerprint(values) -> str:
    """FNV-1a over values separated by 0x1f, as the Rust observer computes it."""
    digest = 0xCBF29CE484222325
    for index, value in enumerate(values):
        data = value.encode("utf-8")
        if index:
            data = b"\x1f" + data
        for byte in data:
            digest ^= byte
            digest = (digest * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return f"{digest:016x}"


def token(raw: str) -> str:
    return "\0" + raw if raw in ("?", ".") else gemmi.cif.as_string(raw)


def site_order(identifier: str):
    return (0, int(identifier), "") if identifier.isdigit() else (1, 0, identifier)


def cif_syntax(text: str) -> list:
    document = gemmi.cif.read_string(text)
    facts = []
    for index, block in enumerate(document):
        facts.append(["block", index, "name", block.name])
        for item in block:
            if item.pair is not None:
                tag, raw = item.pair
                facts.append(["item", index, tag.lower(), "value", token(raw)])
            elif item.loop is not None:
                loop = item.loop
                width, rows = loop.width(), loop.length()
                flat = loop.values
                for column, tag in enumerate(loop.tags):
                    values = [token(flat[row * width + column]) for row in range(rows)]
                    facts.append(["column", index, tag.lower(), "rows", rows])
                    facts.append(["column", index, tag.lower(), "fingerprint", fingerprint(values)])
    return facts


def optional(raw: str):
    return None if raw in ("?", ".") else gemmi.cif.as_string(raw)


def number(raw: str, kind=float):
    return None if raw in ("?", ".") else kind(raw)


def hierarchy(text: str) -> list:
    """Every atom-site row of the first model; the comparison keeps the rows Kekule kept."""
    block = gemmi.cif.read_string(text).sole_block()
    columns = ["id", "label_asym_id", "auth_asym_id", "label_seq_id", "auth_seq_id", "pdbx_PDB_ins_code",
               "label_comp_id", "label_atom_id", "type_symbol", "label_alt_id", "occupancy", "B_iso_or_equiv",
               "Cartn_x", "Cartn_y", "Cartn_z", "pdbx_PDB_model_num"]
    table = block.find("_atom_site.", columns)
    facts = []
    first_model = None
    for row in table:
        model = row[15]
        first_model = first_model or model
        if model != first_model:
            continue
        facts.append(["site", gemmi.cif.as_string(row[0]), [
            optional(row[1]), optional(row[2]), number(row[3], int), optional(row[4]), optional(row[5]),
            optional(row[6]), optional(row[7]), (optional(row[8]) or "").upper(), optional(row[9]),
            number(row[10]), number(row[11]), [float(row[12]), float(row[13]), float(row[14])],
        ]])
    return facts


def biotite_structure(text: str):
    file = pdbx.CIFFile.read(io.StringIO(text))
    try:
        return pdbx.get_structure(file, model=1, altloc="occupancy", extra_fields=["atom_id"], include_bonds=True)
    except Exception as error:
        raise Failure("biotite", str(error)) from None


BOND_ORDERS = {
    struc.BondType.ANY: "any",
    struc.BondType.SINGLE: 1,
    struc.BondType.DOUBLE: 2,
    struc.BondType.TRIPLE: 3,
    struc.BondType.QUADRUPLE: 4,
    struc.BondType.AROMATIC_SINGLE: 1,
    struc.BondType.AROMATIC_DOUBLE: 2,
    struc.BondType.AROMATIC_TRIPLE: 3,
    struc.BondType.AROMATIC: "aromatic",
    struc.BondType.COORDINATION: "dative",
}


def connectivity(text: str) -> list:
    atoms = biotite_structure(text)
    ids = [str(value) for value in atoms.atom_id]
    facts = []
    parent = list(range(len(ids)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j, kind in atoms.bonds.as_array():
        a, b = sorted((ids[i], ids[j]), key=site_order)
        facts.append(["bond", a, b, "order", BOND_ORDERS.get(struc.BondType(kind), str(kind))])
        parent[root(i)] = root(j)
    molecules: dict[int, list[str]] = {}
    for index, identifier in enumerate(ids):
        molecules.setdefault(root(index), []).append(identifier)
    for sites in molecules.values():
        sites.sort(key=site_order)
        facts.append(["molecule", sites[0], "sites", sites])
    return facts


def classification(text: str) -> list:
    atoms = biotite_structure(text)
    masks = [
        ("amino-acid", struc.filter_amino_acids(atoms)),
        ("nucleotide", struc.filter_nucleotides(atoms)),
        ("carbohydrate", struc.filter_carbohydrates(atoms)),
        ("water", struc.filter_solvent(atoms)),
        ("ion", struc.filter_monoatomic_ions(atoms)),
    ]
    facts = []
    starts = struc.get_residue_starts(atoms, add_exclusive_stop=True)
    for start, stop in zip(starts[:-1], starts[1:]):
        name = atoms.res_name[start]
        label = "other"
        for candidate, mask in masks:
            if mask[start]:
                label = candidate
                break
        if label == "nucleotide":
            link = (info.link_type(name) or "").upper()
            label = "dna" if "DNA" in link else "rna" if "RNA" in link else "other"
        insertion = atoms.ins_code[start] or ""
        key = f"{atoms.chain_id[start]}/{atoms.res_id[start]}{insertion}/{name}"
        facts.append(["residue", key, "class", label])
    return sorted(facts, key=str)


def single_conformer(text: str) -> str:
    """The same document reduced to its first model and, per residue, the first
    alternate location in file order. Every other category is kept, including
    `_chem_comp_bond`; anisotropic rows are dropped with their atom sites."""
    document = gemmi.cif.read_string(text)
    block = document.sole_block()
    table = block.find("_atom_site.", ["pdbx_PDB_model_num", "label_alt_id", "label_asym_id", "label_seq_id",
                                        "auth_seq_id", "pdbx_PDB_ins_code"])
    first_model = table[0][0] if len(table) else None
    chosen: dict[tuple, str] = {}
    remove = []
    for index, row in enumerate(table):
        if row[0] != first_model:
            remove.append(index)
            continue
        if row[1] in (".", "?"):
            continue
        if chosen.setdefault((row[2], row[3], row[4], row[5]), row[1]) != row[1]:
            remove.append(index)
    for index in reversed(remove):
        table.remove_row(index)
    block.find_mmcif_category("_atom_site_anisotrop.").erase()
    return document.as_string()


def dssp(text: str) -> list:
    directory = Path(tempfile.mkdtemp(prefix="kekule-dssp-"))
    try:
        source, output = directory / "input.cif", directory / "output.cif"
        source.write_text(text, encoding="utf-8")
        run = subprocess.run([mkdssp(), "--output-format", "mmcif", str(source), str(output)],
                             capture_output=True, text=True, timeout=300)
        if run.returncode != 0 or not output.exists():
            message = run.stderr.strip().splitlines()[-1] if run.stderr.strip() else f"exit {run.returncode}"
            raise Failure("mkdssp", message)
        block = gemmi.cif.read(str(output)).sole_block()
    finally:
        shutil.rmtree(directory, ignore_errors=True)
    summary = block.find("_dssp_struct_summary.", ["label_asym_id", "label_seq_id", "secondary_structure",
                                                     "phi", "psi", "kappa", "alpha", "TCO"])
    pairs = {}
    names = ["label_asym_id", "label_seq_id"]
    for side in ("acceptor_1", "acceptor_2", "donor_1", "donor_2"):
        names += [f"{side}_label_asym_id", f"{side}_label_seq_id", f"{side}_energy"]
    for row in block.find("_dssp_struct_bridge_pairs.", names):
        values = []
        for offset in range(2, len(names), 3):
            asym, seq, energy = row[offset], row[offset + 1], row[offset + 2]
            partner = None if asym in ("?", ".") or seq in ("?", ".") else f"{asym}:{seq}"
            values += [partner, None if partner is None else number(energy)]
        pairs[(row[0], row[1])] = values
    facts = []
    for row in summary:
        ss = row[2] if row[2] not in ("?",) else "."
        angles = [number(row[i]) for i in range(3, 8)]
        facts.append(["dssp", row[0], int(row[1]), [ss, *angles, *pairs.get((row[0], row[1]), [None] * 8)]])
    return facts


def superposition(text: str) -> list:
    structure = gemmi.make_structure_from_block(gemmi.cif.read_string(text).sole_block())
    keys = []
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.name == "CA" and atom.element.name == "C":
                    keys.append((chain.name, str(residue.seqid), atom.altloc))

    def positions(model):
        found = {}
        for chain in model:
            for residue in chain:
                for atom in residue:
                    key = (chain.name, str(residue.seqid), atom.altloc)
                    if atom.name == "CA" and atom.element.name == "C":
                        found[key] = atom.pos
        return [found[key] for key in keys]

    fixed = positions(structure[0])
    facts = [["alpha", "count", len(keys)]]
    for model in structure:
        result = gemmi.superpose_positions(fixed, positions(model))
        facts.append(["model", str(model.num), "rmsd", result.rmsd])
    return facts


def observe(task: str, fmt: str, text: str, options: dict) -> dict:
    if fmt != "mmcif":
        raise Failure("request", f"bio tasks read mmCIF, not {fmt!r}")
    functions = {
        "cif.syntax": cif_syntax,
        "hierarchy": hierarchy,
        "connectivity": connectivity,
        "classification": classification,
        "dssp": dssp,
        "superposition": superposition,
    }
    if task not in functions:
        raise Failure("request", f"unknown task bio.{task}")
    try:
        return {"facts": functions[task](text)}
    except Failure:
        raise
    except Exception as error:
        raise Failure(task, f"{type(error).__name__}: {error}") from None
