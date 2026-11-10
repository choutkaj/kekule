"""Generate the committed DCD interoperability fixture with MDAnalysis."""

from pathlib import Path

import MDAnalysis as mda
import numpy as np


OUTPUT = Path(__file__).with_name("mdanalysis-2.9.0-three-atoms.dcd")

universe = mda.Universe.empty(3, trajectory=True)
universe.add_TopologyAttr("names", ["C", "H", "O"])
universe.add_TopologyAttr("types", ["C", "H", "O"])

frames = [
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]], dtype=np.float32),
    np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype=np.float32),
]

with mda.coordinates.DCD.DCDWriter(str(OUTPUT), n_atoms=3, dt=0.5) as writer:
    for positions in frames:
        universe.atoms.positions = positions
        universe.dimensions = np.array([10, 11, 12, 90, 90, 90], dtype=np.float32)
        writer.write(universe.atoms)
