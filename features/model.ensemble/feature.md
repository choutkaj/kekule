# Shared-Topology Structural Ensembles

## Summary

Represent a finite ordered collection of non-temporal structural realizations
that share one exact topology.

## Behavior/API

- The planned `Ensemble` owns one topology and stable-order
  `EnsembleMember` values.
- Every member contains one complete compatible configuration, optional typed
  observation state, optional finite non-negative weight, and member metadata.
- Weight normalization is explicit and rejects invalid or zero-total weights.
- Members expose borrowed structural views suitable for analyses and prepared
  potentials without creating owned models or copying coordinates.

## Implementation Notes

- Ensemble order has no implicit temporal meaning.
- Independently constructed structurally equivalent topologies are not merged
  without an explicit validated mapping.

## Tests

- Planned tests cover topology mismatch, stable order, weight validation and
  normalization, conformer conversion, observation state, and borrowed-view
  analysis.

## Out Of Scope

- Sparse members, automatic topology reconciliation, consensus structures,
  temporal semantics, and trajectory I/O.

## Revision Notes

- v1: Track the finite shared-topology ensemble contract.
