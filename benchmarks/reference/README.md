# Independent benchmark references

Use the commands and scientific contract in [GUIDE.md](../GUIDE.md). The supported
entry point is `cargo benchmark generate`; `run.py` provides its JSON protocol.
`rdkit/molecule.py` observes indexed molecular structures and validates writer
output; the engine `run_feature.py` modules implement the remaining observations.
No module reads Kekule results to construct an expected value.

The combined pinned environment is `environment.yml`. Engine-specific RDKit and
Biopython environments support narrower installations. RDKit, Biopython, mkdssp,
MDTraj and MDAnalysis are development references, never Rust runtime dependencies.

`rdkit/vs132_reproducer.py` is a focused external-reference investigation with its
own documented version and fixture provenance. The optional trajectory workflow
is documented in [trajectory/VALIDATION.md](trajectory/VALIDATION.md).
