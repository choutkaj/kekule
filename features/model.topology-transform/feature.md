# Topology Instance-Subset Transformations

## Summary

Create immutable subsets of complete molecule instances and explicitly remap
compatible topology-bound structure state through checked lineage.

## Behavior/API

- `topology::transform` retains or removes complete molecule instances without
  mutating the source topology.
- Transform results contain one exact target topology and a complete checked
  `TopologyMapping`.
- Topology-bound positions, configurations, observations, models, atom
  selections, ensembles, trajectory frames, in-memory trajectories, and
  reusable frame buffers can be remapped explicitly.
- Every remap validates exact source, mapping-source, mapping-target, and target
  topology identity.
- Removed selected atoms require an explicit strict or drop policy.

## Implementation Notes

- Subset membership is normalized before construction. Source definition and
  instance order determine target order, and retained reusable definitions are
  cloned once.
- Local atom and bond identifiers are preserved; new topology-level semantic
  and dense identifiers are recorded by the returned mapping.
- Remapping is complete and transactional. Missing target atom state is an
  error rather than a sentinel coordinate or vector.

## Tests

- Focused unit and downstream public-API tests cover deterministic subsetting,
  complete lineage, identity checks, state preservation, transactional
  failures, reusable buffers, and solvent-rich linear scaling.

## Benchmarks

- A synthetic solvent-rich regression exercises many instances sharing one
  definition and records practical linear-scaling evidence without requiring
  external benchmark data.

## Out Of Scope

- Atom-level splitting, instance merging, append or replace edits, hydrogen
  transforms, inferred isomorphism or correspondence, reactive trajectories,
  prepared-potential remapping, production trajectory codecs, and coordinate
  generation.

## Revision Notes

- v1: Track the complete-instance subset and explicit state-remapping milestone.
