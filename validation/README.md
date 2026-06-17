# Validation

Reference validation is separated from the pure Rust library.

- `validation/reference/rdkit/` contains small-molecule reference generators.
- `validation/reference/biopython/` contains macromolecular reference generators.
- `validation/golden/` contains normalized JSON expectations checked by Rust tests or `xtask` commands.

Golden files should record the input, expected normalized behavior, reference implementation, and reference version.
