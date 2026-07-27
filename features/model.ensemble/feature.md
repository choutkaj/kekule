# Shared-Topology Structural Ensembles

## Summary

Represent a finite ordered collection of non-temporal structural realizations
that share one exact topology.

## Behavior/API

- `Ensemble` owns one topology and stable-order
  `EnsembleMember` values.
- Every member contains one complete compatible configuration, optional typed
  observation state, optional finite non-negative weight, and member metadata.
- Weight normalization is explicit and rejects invalid or zero-total weights.
- Members expose borrowed structural views suitable for analyses and prepared
  potentials without creating owned models or copying coordinates.
- `Ensemble::remap_to` stages every member against one exact target topology,
  preserves member order, cells, observations, weights, and properties without
  renormalization, and reports the failing member index.

## Implementation Notes

- Ensemble order has no implicit temporal meaning.
- Independently constructed topologies, including complete layouts for which
  `Topology::same_layout` is true, are not merged without an explicit validated
  mapping.

## Tests

- Tests cover topology mismatch, stable order, weight validation and
  normalization, conformer conversion, observation state, and borrowed-view
  analysis.
- Transformation regressions cover shared target identity, variable cells,
  value-preserved weights and properties, member error context, and source
  immutability.

## Out Of Scope

- Sparse members, automatic topology reconciliation, consensus structures,
  temporal semantics, and trajectory I/O.

## Revision Notes

- v1: Track the finite shared-topology ensemble contract.
- v2: Implement exact-topology ensemble members, observations, metadata,
  explicit weight normalization, conformer conversion, and borrowed views.
- v3: Add transactional finite-ensemble remapping through exact topology
  lineage with member-index failure context.
