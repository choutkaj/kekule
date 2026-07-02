# Smallest-Set Ring Basis

## Summary

Compute a compact ring basis for downstream small-molecule perception.

## Behavior/API

- Exposes `perception::rings::{RingSet, Ring, RingWork, RingPerceptionOptions, RingPerceptionError, perceive_ring_set, perceive_ring_set_with_options}`.
- Reports ring atom and bond IDs for a deterministic cycle basis.
- Reports graph size, candidate cycles, equivalent shortest paths, path expansions, queue/stack peaks, and total work.
- Returns a structured resource-limit error without caching a partial ring set.
- Sets ring perception state through the existing ring membership machinery.
- Cached ring sets are accessible only while ring perception remains fresh.

## Implementation Notes

- Enumerates deterministic shortest-cycle candidates with bounded iterative path reconstruction and greedily keeps linearly independent cycles.
- Adds RDKit-like symmetric extra rings when a ring component has exactly one more shortest-cycle candidate than its cycle rank.
- Defaults allow 1,000,000 atoms, 2,000,000 bonds, 100,000 candidates, 2,000,000 path expansions, 100,000 equivalent shortest paths, cycles up to 4,096 atoms, and 5,000,000 total work units.
- The graph limits accommodate large sparse molecular inputs; candidate/path limits bound symmetric-cycle growth well above observed required corpora.

## Validation

- Unit tests cover monocyclic, fused, and disconnected cases.
- Adversarial tests cover long chains, ladders, theta graphs, fused/bridged systems, symmetric cages, and disconnected mixtures using work counters rather than timing.
- RDKit-generated goldens compare ring atom sets for external PubChem fixtures.

## Out Of Scope

- Exact SymmSSSR parity, ring families, ring aromaticity classification, and exhaustive cycle enumeration.

## Revision Notes

- v1: Deterministic ring basis.
- v2: Shortest-cycle basis passes the RDKit-backed `tiny` corpus; broader required corpora remain pending.
- v3: Fixed bridged and symmetric ring selection exposed by external PubChem validation.
- v4: Hide and clear cached ring sets after invalidating mutations.
- v5: Add bounded work instrumentation, structured resource errors, configurable limits, and iterative shortest-path reconstruction.
- v6: Move the public expert API under the `perception::rings` facade.
- v7: Add PubChem-100k as required broad-corpus validation evidence.
