# Core Molecular Graph

## Summary

Store one asserted chemical entity as a possibly disconnected graph with stable
typed IDs, graph-adjacent stereo, properties, conformers, and private perception.

## Behavior/API

- Provides one shared `Molecule` graph used by both `SmallMolecule` and `MacroMolecule`.
- Permits disconnected topology; connectedness is queried from the graph and is
  not an asserted-entity invariant.
- Supports adding and deleting atoms and bonds.
- Supports first-class stereo elements, stereo groups, and source bond marks attached to stable graph IDs.
- Replaces stereo elements through a validating transactional operation; direct
  mutable access cannot bypass graph-reference or stereo-group invariants.
- Rejects empty stereo groups and duplicate group members.
- Rejects invalid atom IDs, invalid bond IDs, self-bonds, and duplicate bonds.
- Iterates live atoms, live bonds, neighbors, and incident bonds.
- Reports the `i64` sum of asserted formal charges across live atoms without
  requiring sanitization or perception.
- Preserves stable `AtomId` and `BondId` values after deletion.
- Returns scoped `AtomMut` and `BondMut` guards from mutable graph access.
- Owns one internally consistent `PerceptionState` with read-only valence,
  implicit-H, ring, aromaticity/provenance, and CIP queries.
- Owns the stored perception-result vocabulary (`ValenceModel`,
  `AromaticityModel`, `RingMembership`, `Ring`, `RingSet`, and `RingWork`);
  algorithm implementations depend on core, never the reverse.

## Implementation Notes

- Uses slot storage with tombstones so IDs remain stable.
- Maintains adjacency for neighbor and incident-bond iteration.
- Deleting an atom removes its incident bonds.
- Deleting atoms or bonds prunes stereo elements and source bond marks that reference removed topology and drops pruned elements from stereo groups.
- Topology and chemistry-relevant changes immediately clear affected perception;
  stale/fresh flags and mutable cache setters do not exist.
- Property-only and coordinate-only edits do not invalidate chemistry state.
- Mutation guards compare chemistry-relevant fields when released, so obtaining mutable access alone does not stale perception.
- Wrapper `graph_mut()` access is likewise state-neutral; concrete `Molecule`
  mutators remain solely responsible for targeted invalidation.
- Molecule, atom, and bond property maps are stored on the core data structures.
- Local stereo state is graph-adjacent storage on `Molecule`, separate from atom and bond payloads and from derived CIP descriptors.

## Validation

- Current coverage is unit-test based.
- Tests cover empty molecules, insertion, deletion, invalid IDs, self-bonds, duplicates, iteration, stable IDs, counts, chemistry invalidation, state-neutral property/coordinate edits, stereo CRUD, and stereo pruning.
- Reference-tool golden data is not required for this data-structure feature.

## Out Of Scope

- SDF, PDB, or mmCIF parsing.
- Ring detection, aromaticity, valence perception, stereochemistry perception, canonicalization, and validation generation.
- Runtime RDKit or Biopython dependency.

## Revision Notes

- v1: Stable-ID molecular graph and wrapper integration.
- v2: Centralize chemistry invalidation in scoped mutation guards, remove mutable perception-state access, clear stale ring caches, and preserve state across property/coordinate edits.
- v3: Hide perception freshness/cache state from public core API while retaining internal invalidation checks.
- v4: Add graph-adjacent stereo elements, stereo groups, source bond marks, typed stereo IDs, mutation invalidation, and topology-aware stereo pruning.
- v5: Make asserted entity boundaries independent of graph connectedness and
  consolidate all derived chemistry in one private optional `PerceptionState`.
- v6: Keep wrapper mutable access state-neutral so chained perception
  operations retain their prerequisite state; concrete graph mutations still
  invalidate immediately.
- v7: Replace unchecked mutable stereo-element access with transactional
  replacement and enforce nonempty, duplicate-free stereo groups.
- v8: Move all kernel-stored perception vocabulary into core so the graph has
  no physical dependency on algorithm implementations.
- v9: Add `Molecule::formal_charge` as an overflow-safe aggregate over live
  asserted atom payloads.
