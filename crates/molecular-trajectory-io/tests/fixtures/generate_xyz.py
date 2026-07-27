"""Generate the independently written strict XYZ interoperability fixture."""

from pathlib import Path

from ase import Atoms
from ase.io import write


OUTPUT = Path(__file__).with_name("ase-3.26.0-water.xyz")

frames = [
    Atoms(
        symbols=["O", "H", "H"],
        positions=[
            [0.000000, 0.000000, 0.000000],
            [0.957200, 0.000000, 0.000000],
            [-0.239987, 0.927297, 0.000000],
        ],
    ),
    Atoms(
        symbols=["O", "H", "H"],
        positions=[
            [0.100000, 0.200000, 0.300000],
            [1.057200, 0.200000, 0.300000],
            [-0.139987, 1.127297, 0.300000],
        ],
    ),
]

write(OUTPUT, frames, format="xyz")
