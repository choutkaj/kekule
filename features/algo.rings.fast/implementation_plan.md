# Implementation Plan

## Feature ID

`algo.rings.fast`

## Goal

Detect whether atoms and bonds are members of any graph cycle without computing or exposing a canonical ring basis. This feature provides fast ring membership state for later aromaticity and chemistry algorithms.

## Public API

Add an API equivalent to:

- `RingMembership`
- `Molecule::ring_membership() -> Option<&RingMembership>` or a perception-state accessor.
- `perceive_ring_membership(mol: &mut Molecule) -> RingMembership`
- `RingMembership::atom_in_ring(atom: AtomId) -> bool`
- `RingMembership::bond_in_ring(bond: BondId) -> bool`
- `RingMembership::ring_atom_ids() -> impl Iterator<Item = AtomId>`
- `RingMembership::ring_bond_ids() -> impl Iterator<Item = BondId>`

Exact placement can be adjusted to match the emerging module structure, but the API should make clear that this is membership only, not SSSR, cycle basis, ring families, or aromaticity.

## Internal Modules Touched

Expected scope:

- Core graph algorithms module or `crates/molecules/src/lib.rs` while the crate is small.
- `PerceptionState` or a computed-state sidecar for ring membership.
- Unit tests for graph shapes.
- Feature docs under `features/algo.rings.fast/`.

Do not touch SDF or mmCIF parsers except to ensure they do not implicitly run this algorithm.

## Data Model

Store computed ring membership as:

- Set or vector of ring atom flags keyed by `AtomId` slot.
- Set or vector of ring bond flags keyed by `BondId` slot.
- Computed state linked to `PerceptionState::rings`.

The result should tolerate stable IDs and tombstones from `core.graph`. Deleted atoms and bonds should not appear in results.

## Algorithm Outline

Use bridge detection on the undirected molecular graph:

1. Treat live atoms as graph vertices and live bonds as undirected edges.
2. Run DFS over each connected component.
3. Compute discovery indices and low-link values.
4. Mark a bond as not ring-member if it is a bridge.
5. Mark all non-bridge bonds as ring bonds.
6. Mark atoms incident to at least one ring bond as ring atoms.

In an undirected graph, an edge is in at least one cycle exactly when it is not a bridge. This matches the feature goal without constructing a cycle basis.

## Tests

Add tests for:

- Empty molecule has no ring atoms or bonds.
- Single atom has no rings.
- Single bond has no rings.
- Linear chain has no rings.
- Triangle marks all atoms and bonds.
- Square marks all atoms and bonds.
- Ring with a tail marks only ring atoms and ring bonds.
- Fused rings mark shared atoms and all cyclic bonds.
- Disconnected components are handled independently.
- Deleted atoms and bonds are ignored.
- Running perception sets rings computed state to fresh.
- Topology mutation after perception invalidates rings to stale.

## Reference Validation

RDKit can be used to generate golden data for ring membership only:

- Use fixtures with simple chains, monocyclic rings, fused rings, bridged systems, and disconnected molecules.
- Golden JSON should include RDKit version, atom ring membership flags, and bond ring membership flags.
- Do not compare ring basis ordering or SSSR membership in this feature.

If RDKit disagrees because of query bonds or unsupported graph constructs, document fixture exclusions.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature algo.rings.fast
```

## Risks

- Accidentally exposing a ring basis would overcommit this feature.
- Deleted graph slots can cause false positives if not skipped.
- Parallel bonds are currently rejected by `core.graph`, so the algorithm can assume a simple undirected graph for now.
- Future query bonds or zero-order bonds may need explicit inclusion rules.
- Storing results must not become stale after mutation.

## Edge Cases

- Disconnected graph.
- Isolated atoms.
- Ring plus acyclic substituent.
- Fused and bridged cyclic systems.
- Deleted atom with old incident bond tombstones.
- Molecule with no bonds.
- Self-bonds are impossible through `core.graph` but should not be assumed if tests construct internals later.

## Explicitly Out of Scope

- SSSR.
- Minimum cycle basis.
- Ring enumeration.
- Aromaticity.
- Valence perception.
- Stereochemistry.
- Parser behavior.
- RDKit runtime dependency.
