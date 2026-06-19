# Smallest-Set Ring Basis

## Summary

Compute a compact ring basis for downstream small-molecule perception.

## Behavior/API

- Exposes `RingSet`, `Ring`, and `perceive_ring_set`.
- Reports ring atom and bond IDs for a deterministic cycle basis.
- Sets ring perception state through the existing ring membership machinery.

## Implementation Notes

- Uses graph traversal to derive a deterministic fundamental cycle basis.
- The basis is intended for common small molecules and fixture validation, not full RDKit SymmSSSR parity.

## Validation

- Unit tests cover monocyclic, fused, and disconnected cases.
- RDKit reference generator support is included for fixture/golden generation.

## Out Of Scope

- Exact SymmSSSR parity, ring families, ring aromaticity classification, and exhaustive cycle enumeration.

## Revision Notes

- v1: Deterministic ring basis.
