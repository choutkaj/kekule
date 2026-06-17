# Implementation Plan

## Feature ID

`core.graph`

## Goal

Implement the first real shared molecular graph for the crate without adding chemistry perception algorithms or file parsers. The graph should be suitable for both `SmallMolecule` and `MacroMolecule` wrappers and should preserve the architecture boundary between raw topology, biomolecular hierarchy, and computed chemistry state.

## Public API

Add or complete:

- Typed `AtomId` and `BondId` newtypes with `raw()` and `index()` helpers.
- `Molecule::add_atom`, `delete_atom`, `atom`, `atom_mut`, `atoms`, and `atom_ids`.
- `Molecule::add_bond`, `delete_bond`, `bond`, `bond_mut`, `bonds`, and `bond_ids`.
- `Molecule::neighbors`, `incident_bonds`, and `bond_between`.
- `Molecule::atom_count` and `bond_count`.
- Property-map accessors for molecule, atom, and bond records.
- Read-only and internal mutation access for perception state.

Avoid adding public chemistry APIs for validation, sanitization, rings, aromaticity, stereochemistry, SDF, or mmCIF in this feature.

## Internal Modules Touched

Expected initial scope:

- `crates/molecules/src/lib.rs` for the scaffold implementation, unless the file becomes large enough to justify a small `graph` module split.
- Unit tests colocated with the graph implementation.
- `features/core.graph/*` docs.

Do not add RDKit or Biopython dependencies to the Rust workspace.

## Data Model Changes

Use:

- `Vec<Option<Atom>>` for atom slots.
- `Vec<Option<Bond>>` for bond slots.
- `Vec<Vec<BondId>>` adjacency lists keyed by atom slot.
- `PropMap` on molecule, atom, and bond.
- `PerceptionState` on molecule.

IDs are stable handles to slot positions. Deleted slots are tombstones. The first implementation should not reuse deleted IDs.

## Behavior Details

Atom behavior:

- Add atom appends a live atom slot and an empty adjacency list.
- Delete atom rejects invalid or deleted IDs, deletes incident bonds, tombstones the atom slot, and keeps all other IDs stable.
- Atom iterators skip tombstones.

Bond behavior:

- Add bond verifies both endpoints exist and are live.
- Add bond rejects self-bonds and duplicate unordered edges.
- Add bond appends a live bond slot and records the bond in both endpoint adjacency lists.
- Delete bond rejects invalid or deleted IDs, removes adjacency references, tombstones the bond slot, and keeps all other IDs stable.
- Bond iterators skip tombstones.

Neighbor and bond iteration:

- `neighbors(atom)` yields live adjacent atom IDs.
- `incident_bonds(atom)` yields live incident `(BondId, &Bond)` pairs.
- `bond_between(a, b)` treats the graph as undirected and returns `Ok(Some(id))`, `Ok(None)`, or an endpoint error.

Perception invalidation:

- Add/delete atom and add/delete bond invalidate computed chemistry state.
- Mutating molecule, atom, or bond properties does not invalidate topology-derived state unless future docs define a chemically interpreted property.
- Parsing and perception remain separate; no parser should call sanitization implicitly because of this feature.

## Error Cases

Represent and test:

- Invalid atom ID.
- Deleted atom ID.
- Invalid bond ID.
- Deleted bond ID.
- Self-bond.
- Duplicate undirected bond.
- Bond creation with missing or deleted endpoint.
- Optional internal invariant errors only if adjacency consistency cannot be guaranteed by construction.

It is acceptable for invalid and deleted IDs to share the existing `InvalidAtomId` or `InvalidBondId` errors if the public docs state that deleted IDs are no longer live.

## Unit Tests

Implement the matrix from `validation.md`, with special attention to:

- Stable IDs after deletion.
- Cascading bond deletion when deleting atoms.
- Neighbor and incident-bond iteration before and after deletion.
- Perception invalidation for every topology mutation.
- Property map behavior.
- Wrapper sharing through `SmallMolecule` and `MacroMolecule`.

## Validation

Run the standard checks:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
```

No RDKit or Biopython validation is required for the first graph implementation. The feature should remain `validated = false` until tests or documented manual validation are complete and reviewed.

## Risks

- Borrowing ergonomics for iterator APIs may require named iterator structs or collected IDs before yielding references.
- Cascading atom deletion must avoid mutating adjacency while iterating borrowed adjacency lists.
- Tombstones keep IDs stable but can make count and capacity semantics easy to confuse.
- Future tombstone reuse would require generation semantics or stronger lifetime documentation.

## Explicitly Out of Scope

- SDF, Molfile, PDB, and mmCIF parsing.
- Ring detection.
- Aromaticity perception.
- Stereochemistry assignment.
- Valence or sanitization algorithms.
- Validation golden-data generation.
- Biomolecular model, chain, residue, atom-site, alternate-location, occupancy, B-factor, author, or label fields on core `Atom`.
- Runtime RDKit or Biopython integration.
