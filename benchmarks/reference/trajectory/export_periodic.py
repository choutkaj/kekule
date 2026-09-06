"""Export optional scientific references from externally supplied PDB/XTC files.

Development references only: MDAnalysis==2.9.0 and mdtraj==1.11.1.post1.
No fixtures are synthesized. The input must be one connected molecule, with
MDTraj's PDB atom order matching the XTC coordinate order. Bond connectivity is
exported for numerical reconstruction; chemical bond orders are not inferred.

Run this script with --topology SYSTEM.pdb --trajectory INPUT.xtc --output DIR,
then run the trajectory_periodic_reference Rust example with
DIR/topology.txt INPUT.xtc DIR. No Python libraries are Rust runtime dependencies.
"""

import argparse
import hashlib
import json
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.transformations.nojump import NoJump
import mdtraj as md
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if mda.__version__ != "2.9.0" or md.__version__ != "1.11.1.post1":
        raise ValueError("Use the pinned reference versions documented above")
    source = md.load_xtc(str(args.trajectory), top=str(args.topology))
    molecules = source.topology.find_molecules()
    if len(molecules) != 1:
        raise ValueError("This reference profile requires one connected molecule")
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "topology.txt").open("w", encoding="utf-8") as output:
        output.write(f"{source.n_atoms} {source.topology.n_bonds}\n")
        output.write(" ".join(atom.element.symbol for atom in source.topology.atoms) + "\n")
        for bond in source.topology.bonds:
            output.write(f"{bond.atom1.index} {bond.atom2.index}\n")
    whole = source.make_molecules_whole(inplace=False)
    image = source.image_molecules(inplace=False, anchor_molecules=molecules)
    universe = mda.Universe(str(args.topology), str(args.trajectory))
    nojump = NoJump()
    unwrapped = np.asarray([nojump(ts).positions.copy().astype(np.float64) / 10.0 for ts in universe.trajectory])
    if unwrapped.shape != source.xyz.shape:
        raise ValueError("Independent readers disagree about trajectory dimensions")
    for name, values in {"raw": source.xyz, "whole": whole.xyz, "image": image.xyz, "unwrap": unwrapped}.items():
        np.savetxt(args.output / f"{name}.txt", values.reshape(-1, 3), fmt="%.17g")
    metadata = {
        "schema": 1,
        "frames": source.n_frames,
        "atoms": source.n_atoms,
        "bonds": source.topology.n_bonds,
        "coordinate_unit": "nm",
        "tolerance_nm": 0.0001,
        "references": {"mdtraj": md.__version__, "MDAnalysis": mda.__version__, "numpy": np.__version__},
        "operations": {"whole": "MDTraj make_molecules_whole", "image": "MDTraj image_molecules, entire molecule as explicit anchor", "unwrap": "MDAnalysis NoJump, sequential fractional-coordinate continuity"},
        "inputs": [{"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in (args.topology, args.trajectory)],
    }
    (args.output / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
