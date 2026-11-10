"""Generate the committed TRR interoperability fixture with MDAnalysis."""

from pathlib import Path

import MDAnalysis as mda
import numpy as np


OUTPUT = Path(__file__).with_name("mdanalysis-2.9.0-three-atoms.trr")

universe = mda.Universe.empty(
    3,
    trajectory=True,
    velocities=True,
    forces=True,
)
universe.add_TopologyAttr("names", ["C", "H", "O"])
universe.add_TopologyAttr("types", ["C", "H", "O"])

frames = [
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]], dtype=np.float32),
    np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float32),
]

with mda.coordinates.TRR.TRRWriter(str(OUTPUT), n_atoms=3) as writer:
    for step, positions in enumerate(frames):
        universe.atoms.positions = positions
        universe.atoms.velocities = positions + 0.5
        universe.atoms.forces = positions + 1.0
        universe.dimensions = np.array([10, 11, 12, 90, 90, 90], dtype=np.float32)
        universe.trajectory.ts.time = step * 0.25
        universe.trajectory.ts.data["step"] = step
        universe.trajectory.ts.data["lambda"] = 0.125 + step * 0.125
        writer.write(universe.atoms)
