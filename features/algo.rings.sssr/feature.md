# Smallest-Set Ring Basis

## Summary

Compute a compact ring basis for downstream small-molecule perception.

## Behavior/API

- Exposes `RingSet`, `Ring`, and `perceive_ring_set`.
- Reports ring atom and bond IDs for a deterministic cycle basis.
- Sets ring perception state through the existing ring membership machinery.
- Cached ring sets are accessible only while ring perception remains fresh.

## Implementation Notes

- Enumerates deterministic shortest-cycle candidates and greedily keeps linearly independent cycles.
- Adds RDKit-like symmetric extra rings when a ring component has exactly one more shortest-cycle candidate than its cycle rank.
- The basis is intended for common small molecules and fixture validation, not full RDKit SymmSSSR parity for all graph families.

## Validation

- Unit tests cover monocyclic, fused, and disconnected cases.
- RDKit-generated goldens compare ring atom sets for external PubChem fixtures.

## Out Of Scope

- Exact SymmSSSR parity, ring families, ring aromaticity classification, and exhaustive cycle enumeration.

## Revision Notes

- v1: Deterministic ring basis.
- v2: Shortest-cycle basis passes the RDKit-backed `tiny` corpus; broader required corpora remain pending.
- v3: Fixed bridged and symmetric ring selection exposed by external PubChem validation.
- v4: Hide and clear cached ring sets after invalidating mutations.
