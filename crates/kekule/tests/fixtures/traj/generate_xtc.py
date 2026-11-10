"""Generate the committed compressed XTC fixture with MDAnalysis."""

from pathlib import Path

import MDAnalysis as mda
import numpy as np


OUTPUT = Path(__file__).with_name("mdanalysis-2.9.0-twelve-atoms.xtc")

universe = mda.Universe.empty(12, trajectory=True)
universe.add_TopologyAttr("names", [f"C{i + 1}" for i in range(12)])
universe.add_TopologyAttr("types", ["C"] * 12)

base = np.array(
    [[0.1 * i, 0.2 * i, 0.3 * i] for i in range(12)],
    dtype=np.float32,
)

with mda.coordinates.XTC.XTCWriter(str(OUTPUT), n_atoms=12, precision=3) as writer:
    for step, shift in enumerate([0.0, 0.01]):
        universe.atoms.positions = base + shift
        universe.dimensions = np.array([20, 21, 22, 80, 85, 75], dtype=np.float32)
        universe.trajectory.ts.time = step * 0.25
        universe.trajectory.ts.frame = step
        universe.trajectory.ts.data["step"] = step
        writer.write(universe.atoms)
