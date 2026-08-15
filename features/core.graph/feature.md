# Core Molecular Graph

## Summary

Store one asserted chemical entity as a connected graph with stable typed IDs,
graph-adjacent stereo, properties, conformers, and private perception. Empty and
single-atom values are valid connected boundary cases.

## Behavior/API

- Provides one shared `Molecule` graph used by both `SmallMolecule` and `MacroMolecule`.
- Requires every completed nontrivial `Molecule` to contain exactly one graph
  component.
- Constructs topology through `MoleculeBuilder` and changes topology through a
  transactional `MoleculeEditor`. Their working copies may be temporarily
  disconnected, but the staging `Molecule` cannot be publicly borrowed, cloned,
  replaced, or otherwise extracted. `build` and `commit` are the only public
  publication routes and reject disconnected results with a structured error.
- Failed builds and edit commits do not expose or install an invalid graph; a
  failed edit leaves the original molecule unchanged.
- Keeps raw atom/bond insertion and deletion plus builder/editor staging access
  crate-private so public callers cannot bypass the connectedness boundary.
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
  implicit-H, ring, model-perceived aromaticity, and CIP queries.
- Exposes immutable exact perception-section views, detached checked
  construction, and one whole-state validated atomic installation operation;
  incremental perception mutators remain crate-private.
- Owns the stored perception-result vocabulary (`ValenceModel`,
  `AromaticityModel`, `RingMembership`, `Ring`, and `RingSet`);
  algorithm implementations depend on core, never the reverse.
- Keeps `PerceptionState` limited to semantic installed chemical perception;
  algorithm diagnostics and work counters are ephemeral or sidecar state and
  do not participate in molecule identity, equality, or reconstruction.

## Implementation Notes

- Uses slot storage with tombstones so IDs remain stable.
- Builder/editor staging owns the only temporarily disconnected public workflow,
  but exposes only focused topology operations rather than a staging-graph
  reference. Ordinary conformer, property, and chemistry operations occur after
  `build`; all completed values exposed through the core API satisfy
  connectedness.
- Checks atom, bond, conformer, stereo-element, and stereo-group collection
  slots before insertion and returns `MoleculeError::IdentifierCapacityExceeded`
  without changing graph state.
- Live-ID iterators advance in the native `u32` identifier domain and never
  reconstruct an ID through narrowing from `usize`.
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
- Whole-state perception installation validates stable-slot dimensions, live
  atom/bond/stereo references, duplicate assignments, and ring coherence before
  replacing prior state. It does not validate or retain algorithm work
  diagnostics. Failed installation is state-preserving.

## Tests

- Current coverage is unit-test based.
- Tests cover empty and single-atom boundary cases, connected builder success,
  disconnected builder rejection, editor rollback on a disconnecting change,
  invalid IDs, self-bonds, duplicates, iteration, stable IDs, counts, checked
  synthetic capacity boundaries, transactional capacity rejection, chemistry
  invalidation, state-neutral property/coordinate edits, stereo CRUD, and
  stereo pruning.
- Downstream compile-fail rustdoc regressions prove neither immutable cloning
  nor mutable `mem::take` can extract builder/editor staging as a `Molecule`.
- Reference-tool golden data is not required for this data-structure feature.
- Downstream regressions cover exact installed-perception export/install,
  malformed transactional rejection, and normal post-install invalidation.

## Out Of Scope

- SDF, PDB, or mmCIF parsing.
- Ring detection, aromaticity, valence perception, stereochemistry perception, canonicalization, and benchmark generation.
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
- v10: Make atom insertion fallible and enforce one checked fixed-width ID
  strategy across atoms, bonds, conformers, stereo elements, stereo groups,
  and live-ID iterators.
- v11: Add immutable exact perception-section views, detached construction,
  and checked atomic whole-state installation for canonical persistence.
- v12: Make connectedness a completed-`Molecule` invariant, route public
  topology construction and mutation through checked builder/editor staging,
  and guarantee rollback when an edit would disconnect the graph.
- v13: Remove public builder/editor staging-graph access, make `build` and
  `commit` the only publication routes, migrate ordinary mutation after
  finalization, and add downstream compile-fail extraction regressions.
- v14: Restrict installed ring perception to semantic membership and basis
  state, excluding algorithm diagnostics from molecule identity and restore.
- v15: Restrict installed aromaticity to semantic model-perceived membership,
  remove imported representation provenance, and expose the installed
  `AromaticityModel` directly.
