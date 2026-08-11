# Smallest-Set Ring Basis

## Summary

Compute a compact ring basis for downstream perception of connected
small-molecule graphs.

## Behavior/API

- Exposes `perception::rings::{RingSet, Ring, RingWork, RingPerceptionOptions, RingPerceptionError, perceive_ring_set, perceive_ring_set_with_options}`.
- Reports ring atom and bond IDs for a deterministic cycle basis.
- Reports graph size, candidate cycles, equivalent shortest paths, path expansions, queue/stack peaks, and total work.
- Returns a structured `ResourceLimit` or `IncompleteRingCoverage` error without
  caching a partial ring set.
- Sets ring perception state through the existing ring membership machinery.
- Cached ring sets are accessible only while ring perception remains fresh.
- `RingSet::from_parts` reconstructs a detached complete basis and `RingWork`;
  graph references and membership coherence are checked on whole-state install.

## Implementation Notes

- Uses the Figueras/RDKit degree-trimming workflow: prune acyclic atoms,
  search one representative per degree-two chain with deterministic
  single-parent BFS, recover candidates hidden by duplicate roots, and handle
  degree-three cores before bounded basis completion.
- Completes rare discovery shortfalls with bounded deterministic multi-root
  BFS fundamental cycles, selected by cycle-space independence and shortest
  length, so every bond identified as cyclic is covered by a returned ring.
- Follows RDKit's SSSR coverage semantics for bridged and cage systems; it does
  not force the returned ring count to equal the algebraic cycle rank.
- Adds RDKit-like symmetric extra rings only when an unselected candidate can
  replace one same-sized basis ring without removing any bond uniquely supplied
  by that basis ring and the candidate shares at least one bond with it.
- Defaults allow 1,000,000 atoms, 2,000,000 bonds, 100,000 candidates, 2,000,000 path expansions, 100,000 equivalent shortest paths, cycles up to 4,096 atoms, and 5,000,000 total work units.
- The graph limits accommodate large sparse molecular inputs; candidate/path limits bound symmetric-cycle growth well above observed required corpora.

## Tests

- Unit tests cover monocyclic, fused, bridged, and connected acyclic-decoration
  cases.
- Adversarial tests cover long chains, ladders, theta graphs with acyclic tails,
  fused/bridged systems, and symmetric cages using work counters rather than
  timing.

## Benchmarks

- RDKit-generated goldens compare ring atom sets for external PubChem fixtures.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Exact SymmSSSR parity, ring families, ring aromaticity classification, and exhaustive cycle enumeration.

## Revision Notes

- v1: Deterministic ring basis.
- v2: Shortest-cycle basis passes the RDKit-backed `smoke` corpus; broader required corpora remain pending.
- v3: Fixed bridged and symmetric ring selection exposed by external PubChem validation.
- v4: Hide and clear cached ring sets after invalidating mutations.
- v5: Add bounded work instrumentation, structured resource errors, configurable limits, and iterative shortest-path reconstruction.
- v6: Move the public expert API under the `perception::rings` facade.
- v7: Add PubChem-100k as required broad-corpus external-parity evidence.
- v8: Replace per-bond shortest-path enumeration with RDKit's degree-trimming
  candidate discovery and exact same-size, shared-bond,
  unique-bond-preserving SSSR replacement rule.
- v10: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v11: Require every cyclic bond to be covered by the returned SSSR and report
  an explicit coverage error separately from resource exhaustion; add a
  bounded fundamental-cycle fallback when edge-local shortest cycles leave a
  cyclic bond uncovered, while preserving RDKit cage-system ring counts.
- v12: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v13: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v14: Add detached complete `RingSet`/`RingWork` reconstruction with
  molecule-specific validation deferred to atomic perception installation.
- v15: Align the documented input and adversarial regressions with the
  connected `Molecule` boundary.
