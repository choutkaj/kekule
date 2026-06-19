# Fast Ring Membership Detection

## Summary

Detect whether atoms and bonds are members of any graph cycle without computing a canonical ring basis.

## Behavior/API

- Exposes `RingMembership` and `perceive_ring_membership`.
- Reports ring membership for live atoms and bonds.
- Ignores deleted graph slots.
- Sets ring perception state to fresh after successful perception.

## Implementation Notes

- Uses bridge detection on the undirected molecular graph.
- A bond is a ring bond exactly when it is not a bridge.
- Ring atoms are atoms incident to at least one ring bond.
- Handles disconnected components.

## Validation

- Current coverage is unit-test based.
- RDKit golden validation is planned through `validation.harness`.
- Fixtures live under `validation/features/algo.rings.fast/`.

## Out Of Scope

- SSSR, minimum cycle basis, ring enumeration, aromaticity, valence perception, stereochemistry, and parser behavior.
- Runtime RDKit dependency.

## Revision Notes

- v1: Graph-cycle membership perception.
