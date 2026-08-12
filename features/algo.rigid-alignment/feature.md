# Rigid Molecular Alignment

## Summary

Fit one molecular coordinate snapshot onto another with a proper weighted rigid
transform over an exact-topology atom selection.

## Behavior/API

- Exposes `alignment::kabsch` and `alignment::kabsch_with_options` over two
  borrowed `ModelView` values and one topology-bound `AtomSelection`.
- Requires moving, reference, and selection to retain the same `Arc<Topology>`
  allocation. Correspondence follows the selection's sorted dense-index order.
- Returns a `RigidAlignment` containing the existing `geometry::RigidTransform`
  mapping moving coordinates into reference coordinates, post-fit weighted RMSD
  in `MODEL_LENGTH_UNIT`, and the selected atom count.
- Uses uniform weights by default. Explicit weights correspond to selection
  order, must match its length, and must be finite and strictly positive.
- Rejects periodic configurations by default.
  `PeriodicAlignmentPolicy::UseStoredCoordinates` ignores cells and performs no
  imaging, wrapping, unwrapping, or minimum-image correction.
- Requires at least three selected atoms and rank-two geometry on both sides.
  Planar non-collinear selections are valid; coincident, collinear, and
  scale-relatively near-collinear selections are rejected structurally.
- Leaves every input topology, coordinate container, cell, observation,
  property, and selection unchanged.

## Implementation Notes

- The weighted kernel uses `f64` weighted online centroid, centered scatter, and
  cross-covariance accumulation with compensated scalar sums.
- A private fixed-size one-sided Jacobi SVD solves the 3-by-3 orthogonal
  Procrustes problem. Orthonormal basis completion selects the optimal proper
  right-handed rotation, and `RigidTransform::new` validates the result.
- A centered point set is rank deficient when its second-largest weighted
  scatter eigenvalue is not greater than `1e-12` times its largest eigenvalue.
  This variance-relative threshold is independent of coordinate unit and
  overall geometric scale. The cross-covariance must likewise determine two
  singular directions above a `1e-12` relative ratio.
- Explicit weights are normalized by their maximum before accumulation, making
  the result invariant to a common positive weight scale within floating-point
  precision.
- Runtime is O(n) in selected atoms. Numerical workspace is fixed size; the
  implementation does not clone complete coordinate arrays or construct owned
  models.

## Tests

- Focused tests cover identity, translation, asymmetric non-axis rotation and
  translation, transform direction, subset fitting, noisy RMSD, weighted
  fitting, weight-scale invariance, reflected coordinates, planar geometry,
  point-count and rank failures, threshold-adjacent degeneracy, topology and
  stale-selection mismatches, invalid weights, periodic policy, large coordinate
  offsets, units, and input immutability.
- An external integration test compiles the focused public module without
  expanding the crate prelude.
- Rustdoc examples exercise the moving-to-reference contract.

## Out Of Scope

- Topology reconciliation, inferred or symmetry-aware correspondence, sequence
  or structural alignment, atomic-mass weighting, periodic imaging, coordinate
  mutation, batch or consensus fitting, non-rigid fitting, and GPU execution.

## Revision Notes

- v1: Add same-topology selection-based uniform or explicitly weighted proper
  Kabsch alignment with explicit periodic policy and structured failures.
- v2: Use exact shared-allocation checks
  across both views and the compiled selection.
