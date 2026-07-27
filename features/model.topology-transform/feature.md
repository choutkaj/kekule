# Topology Instance-Subset Transformations

## Summary

Create immutable subsets of complete molecule instances and explicitly remap
compatible topology-bound structure state through checked lineage.

## Behavior/API

- `topology::transform` retains or removes complete molecule instances without
  mutating the source topology.
- Duplicate requests are normalized; source definition and instance order is
  filtered deterministically; empty results are rejected; and exact identity is
  preserved only for no-op edits.
- Transform results contain one exact target topology and a complete checked
  `TopologyMapping` across definitions, instances, qualified atoms and bonds,
  and dense atom and bond indices.
- Topology-bound positions, configurations, observations, models, atom
  selections, ensembles, trajectory frames, in-memory trajectories, and
  reusable frame buffers can be remapped explicitly.
- Every remap validates exact source, mapping-source, mapping-target, and target
  topology identity.
- Removed selected atoms require an explicit strict or drop policy.
- Periodic cells, observation provenance and properties, ensemble weights and
  properties, velocities, forces, time, step, frame properties, and source
  order are preserved.
- Target atoms lacking mapped source state are rejected. Prepared potentials
  remain bound to their original exact topology and are never remapped.

## Implementation Notes

- Subset membership is normalized before construction. Source definition and
  instance order determine target order, and retained reusable definitions are
  cloned once. The implementation has no quadratic builder behavior, while
  topology index construction and checked lineage currently use ordered maps
  and sets, giving practical O(n log n) scaling. No graph correspondence is
  inferred.
- Dense or hash-backed topology indices and lineage maps remain a future
  optimization if strict or expected linear construction becomes necessary.
- Local atom and bond identifiers are preserved; new topology-level semantic
  and dense identifiers are recorded by the returned mapping.
- Remapping is complete and transactional. Missing target atom state is an
  error rather than a sentinel coordinate or vector.

## Tests

- Focused unit and downstream public-API tests cover deterministic subsetting,
  complete lineage, identity checks, state preservation, transactional
  failures, strict/drop selection policy, member/frame error context, and
  reusable buffer allocation stability. Buffer regressions preserve every
  pre-existing destination field after a later unmapped-target-atom failure and
  clear stale optional state when a positions-only frame follows a full frame.

## Benchmarks

- A synthetic regression subsets 20,000 solvent instances sharing one
  definition, verifies definition reuse and complete mapping cardinality, and
  guards against quadratic builder cloning. It is a practical large-input
  regression, not proof of linear asymptotics.

## Out Of Scope

- Atom-level splitting, instance merging, append or replace edits, hydrogen
  transforms, inferred isomorphism or correspondence, reactive trajectories,
  prepared-potential remapping, production trajectory codecs, and coordinate
  generation.

## Revision Notes

- v1: Track the complete-instance subset and explicit state-remapping milestone.
- v2: Correct the ordered-map complexity contract and strengthen reusable-buffer
  transactionality and stale-state regressions.
