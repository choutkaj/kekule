"""Separate RDKit's VS132 SMILES interpretation from its CIP ranking.

Run from the repository root with:
    uv run --with rdkit==2026.3.6 --python 3.13 python \
        benchmarks/reference/rdkit/vs132_reproducer.py

RDKit and NumPy are optional external references, never Rust dependencies.
The source fixture is pinned in tests/fixtures/cip/provenance.toml.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase


SMILES = "CC=1C=CC=2[N@@]3CC=4C=C(C=CC4[N@](CC2C1)C3)C"
FIXTURE_SHA256 = "00964798d1f40a9ae73e6a44bcb249c61cc5c2fa5cd5290e8f8399ecfdda15fd"


def nitrogen_state(molecule):
    return [
        {
            "atom": atom.GetIdx(),
            "tag": str(atom.GetChiralTag()),
            "neighbors": [neighbor.GetIdx() for neighbor in atom.GetNeighbors()],
            "label": atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None,
            "priority": list(
                atom.GetPropsAsDict(includePrivate=True).get("_CIPNeighborOrder", [])
            ),
        }
        for atom in molecule.GetAtoms()
        if atom.GetSymbol() == "N"
    ]


def assign(molecule):
    for obj in list(molecule.GetAtoms()) + list(molecule.GetBonds()):
        if obj.HasProp("_CIPCode"):
            obj.ClearProp("_CIPCode")
    Chem.AssignCIPLabels(molecule, maxRecursiveIterations=1_000_000)
    return nitrogen_state(molecule)


def main():
    if rdBase.rdkitVersion != "2026.03.6":
        raise SystemExit("This reproduction is pinned to RDKit 2026.03.6")
    Chem.SetUseLegacyStereoPerception(False)
    fixture = (
        Path(__file__).resolve().parents[3]
        / "crates/kekule/tests/fixtures/cip/VS132.sdf"
    )
    fixture_bytes = fixture.read_bytes()
    assert hashlib.sha256(fixture_bytes).hexdigest() == FIXTURE_SHA256
    keep_chirality = (
        Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
    )

    original = Chem.MolFromSmiles(SMILES, sanitize=False)
    parsed = nitrogen_state(original)
    cleaned = Chem.Mol(original)
    Chem.SanitizeMol(cleaned)
    cleaned_labels = assign(cleaned)
    Chem.SanitizeMol(original, sanitizeOps=keep_chirality)
    original_labels = assign(original)

    geometry = Chem.MolFromMolBlock(
        fixture_bytes.decode("utf-8"), sanitize=False, removeHs=False
    )
    Chem.SanitizeMol(geometry, sanitizeOps=keep_chirality)
    xyz = np.array(geometry.GetConformer().GetPositions())
    centers = []
    for atom in geometry.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        center = atom.GetIdx()
        neighbors = list(atom.GetNeighbors())
        # Rule 1a at both nitrogens is decided in the first ligand sphere:
        # bridge [N,H,H] > aryl [C,C,C] > benzyl [C,H,H] > lone pair.
        bridge = next(
            neighbor
            for neighbor in neighbors
            if any(
                other.GetSymbol() == "N" and other.GetIdx() != center
                for other in neighbor.GetNeighbors()
            )
        ).GetIdx()
        aryl = next(n.GetIdx() for n in neighbors if n.GetIsAromatic())
        benzyl = next(n.GetIdx() for n in neighbors if n.GetIdx() not in (bridge, aryl))
        priority = [bridge, aryl, benzyl]
        determinant = float(np.linalg.det(xyz[priority] - xyz[center]))
        # The center can stand in for the lone-pair vertex: moving that
        # vertex outwards, opposite the three bonds, preserves the sign.
        centers.append(
            {
                "atom": center,
                "priority": priority,
                "oriented_determinant": determinant,
                "independent_label": "S" if determinant > 0 else "R",
            }
        )

    # Transport the published coordinates onto the ORIGINAL SMILES atom
    # ordering, so no change of graph or ligand ranking can explain a label
    # difference. Explicit SDF hydrogens follow all nineteen heavy atoms.
    heavy = Chem.RemoveHs(geometry, sanitize=False)
    mapping = heavy.GetSubstructMatch(original, useChirality=False)
    assert len(mapping) == original.GetNumAtoms() == 19
    from_geometry = Chem.Mol(original)
    for atom in from_geometry.GetAtoms():
        if atom.GetSymbol() != "N":
            continue
        center = mapping[atom.GetIdx()]
        order = [mapping[n.GetIdx()] for n in atom.GetNeighbors()]
        volume = float(np.linalg.det(xyz[order] - xyz[center]))
        assert abs(volume) > 1.0
        atom.SetChiralTag(
            Chem.ChiralType.CHI_TETRAHEDRAL_CCW
            if volume > 0
            else Chem.ChiralType.CHI_TETRAHEDRAL_CW
        )
    geometry_labels = assign(from_geometry)
    assert [center["independent_label"] for center in centers] == ["S", "S"]
    assert [atom["label"] for atom in geometry_labels] == ["S", "S"]
    assert [atom["priority"] for atom in original_labels] == [
        atom["priority"] for atom in geometry_labels
    ]
    assert [atom["label"] for atom in original_labels] == ["R", "S"]

    print(
        json.dumps(
            {
                "rdkit": rdBase.rdkitVersion,
                "smiles": SMILES,
                "parsed": parsed,
                "default_sanitization": cleaned_labels,
                "preserved_smiles_tags": original_labels,
                "published_3d_independent_calculation": centers,
                "smiles_to_sdf_heavy_atom_mapping": mapping,
                "same_graph_with_published_3d_tags": geometry_labels,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
