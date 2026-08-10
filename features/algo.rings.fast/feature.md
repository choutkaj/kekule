# Fast Ring Membership Detection

## Summary

Detect whether atoms and bonds in a connected molecular graph are members of
any cycle without computing a canonical ring basis.

## Behavior/API

- Exposes `perception::rings::{RingMembership, perceive_ring_membership}`.
- Reports ring membership for live atoms and bonds.
- Ignores deleted graph slots.
- Sets ring perception state to fresh after successful perception.
- Cached membership is accessible only while ring perception remains fresh.
- `RingMembership` exposes and reconstructs complete stable atom/bond slot
  flags, including false tombstone positions, for checked whole-state restore.

## Implementation Notes

- Uses bridge detection on the undirected molecular graph.
- Uses an explicit DFS frame stack, so graph depth does not consume the Rust call stack.
- A bond is a ring bond exactly when it is not a bridge.
- Ring atoms are atoms incident to at least one ring bond.
- Traversal remains component-safe internally, but the public molecule boundary
  supplies one connected graph.
- Topology or chemistry mutation clears cached membership rather than exposing stale results.

## Tests

- Unit tests cover empty/single-atom boundary cases, linear graphs, rings, fused
  rings with acyclic tails, tombstones, and core graph-cycle membership.

## Benchmarks

- RDKit-generated goldens compare ring membership for external PubChem fixtures.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- SSSR, minimum cycle basis, ring enumeration, aromaticity, valence perception, stereochemistry, and parser behavior.
- Runtime RDKit dependency.

## Revision Notes

- v1: Graph-cycle membership perception.
- v2: Hide and clear cached membership after invalidating mutations.
- v3: Replace recursive bridge traversal with an explicit stack for very large graphs.
- v4: Move the public expert API under the `perception::rings` facade.
- v5: Add PubChem-100k as required broad-corpus external-parity evidence.
- v6: Mark current broad-corpus external-parity evidence as matching in feature metadata.
- v8: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v9: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v10: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v11: Add exact stable-slot flag inspection/construction for canonical
  perception reconstruction.
- v12: Align the public contract and regressions with the connected
  `Molecule` invariant while retaining component-safe internal traversal.
