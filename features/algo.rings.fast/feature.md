# Fast Ring Membership Detection

## Summary

Detect whether atoms and bonds are members of any graph cycle without computing a canonical ring basis.

## Behavior/API

- Exposes `perception::rings::{RingMembership, perceive_ring_membership}`.
- Reports ring membership for live atoms and bonds.
- Ignores deleted graph slots.
- Sets ring perception state to fresh after successful perception.
- Cached membership is accessible only while ring perception remains fresh.

## Implementation Notes

- Uses bridge detection on the undirected molecular graph.
- Uses an explicit DFS frame stack, so graph depth does not consume the Rust call stack.
- A bond is a ring bond exactly when it is not a bridge.
- Ring atoms are atoms incident to at least one ring bond.
- Handles disconnected components.
- Topology or chemistry mutation clears cached membership rather than exposing stale results.

## Validation

- Unit tests cover core graph-cycle membership behavior.
- RDKit-generated goldens compare ring membership for external PubChem fixtures.

## Out Of Scope

- SSSR, minimum cycle basis, ring enumeration, aromaticity, valence perception, stereochemistry, and parser behavior.
- Runtime RDKit dependency.

## Revision Notes

- v1: Graph-cycle membership perception.
- v2: Hide and clear cached membership after invalidating mutations.
- v3: Replace recursive bridge traversal with an explicit stack for very large graphs.
- v4: Move the public expert API under the `perception::rings` facade.
