# Validation

## Validation Level

`core.graph` is a data-structure feature. Initial validation is unit-test based and does not require RDKit or Biopython golden data.

Reference tools are not runtime dependencies and should not be used to define hidden behavior for this feature. Later chemistry features may generate normalized JSON golden data from RDKit or Biopython, but this feature only needs deterministic graph behavior tests.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
```

## Unit Test Matrix

The implementation PR should include tests for:

- Creating an empty molecule.
- Adding one atom and fetching it by `AtomId`.
- Adding multiple atoms and verifying monotonic, typed IDs.
- Iterating atoms skips deleted atoms.
- Deleting an atom returns the atom and makes later access fail.
- Deleting an atom does not change IDs for remaining atoms.
- Deleting an atom removes all incident bonds.
- Adding a bond between valid atoms returns a stable `BondId`.
- Fetching and iterating bonds skips deleted bonds.
- Deleting a bond returns the bond and makes later access fail.
- Deleting a bond does not change IDs for remaining bonds.
- Rejecting self-bonds.
- Rejecting duplicate undirected bonds in both endpoint orders.
- Rejecting bonds with missing or deleted atom endpoints.
- Neighbor iteration for isolated, terminal, middle, and branched atoms.
- Incident-bond iteration after bond deletion and atom deletion.
- `bond_between` returns the expected bond ID or `None`.
- Topology mutations turn `Fresh` perception state into `Stale`.
- Topology mutations leave `Absent` perception state as `Absent`.
- Molecule, atom, and bond property maps can be read and mutated without changing topology.
- `SmallMolecule` and `MacroMolecule` still share the same core `Molecule` topology type.

## Manual Validation

Manual validation is acceptable only for documenting public API review before implementation. It cannot mark the feature as validated.

## Golden Data

No golden data is required for this feature. If future validation fixtures are added, they should be normalized JSON and include:

- Tool name and version, if generated from a reference tool.
- Input molecule identifier.
- Atom and bond counts.
- Edge list with stable endpoint IDs.
- Expected error or mutation behavior when applicable.

Do not mark `validated = true` until either the unit-test matrix is implemented with review evidence or a documented manual validation record exists.
