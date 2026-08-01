# Trajectory RMSD Series

## Summary

Calculate unit-bearing RMSD series from every frame of one finite trajectory to
one selected reference frame, with explicit separation between direct
measurement and rigidly aligned measurement.

## Behavior/API

- `Trajectory::rmsd_to_frame` and `rmsd_to_frame_with_options` calculate direct
  RMSD from coordinates exactly as stored. They never center, fit, transform,
  or mutate trajectory state.
- `Trajectory::aligned_rmsd_to_frame` and its options overload fit each frame
  without materializing a transformed trajectory, then measure RMSD. Fit and
  measurement selections are independent so callers may fit a backbone and
  measure a ligand or domain.
- Results are `Quantity<Vec<f64>>` values in `MODEL_LENGTH_UNIT`, ordered one to
  one with input frames.
- Direct RMSD supports uniform or explicit positive finite selection-order
  weights. Weights are normalized by their maximum and results divide by total
  normalized weight, making them invariant to a common positive scale.
- Direct RMSD requires at least one selected atom. Aligned RMSD inherits the
  rigid-alignment point-count, rank, weighting, topology, and periodic policies.
- Direct periodic measurement is rejected by default. Explicit
  `UseStoredCoordinates` ignores cells without imaging or unwrapping.
- Bounds, selection identity, weight validation, periodic rejection, alignment
  failure, and numerical failure are structured and retain frame context where
  applicable.

## Implementation Notes

- Direct RMSD uses compensated f64 accumulation. Runtime is O(frames × measured
  atoms) with O(frames) result storage.
- Fused aligned RMSD delegates fitting to Kekule's Kabsch kernel and applies
  each returned transform only while measuring selected points. Runtime is
  O(frames × (fit atoms + measured atoms)); it does not clone frames or mutate
  coordinates.

## Tests

- Focused regressions distinguish raw from aligned RMSD; cover uniform and
  explicit weights, weight-scale semantics, distinct fit/measurement
  selections, split/fused agreement, units, reference bounds, empty and stale
  selections, invalid weights, periodic policy, frame-indexed alignment
  failures, and input immutability.
- A downstream integration test compiles direct, split, and fused public usage.

## Out Of Scope

- Atomic-mass convenience weighting, periodic imaging, symmetry-aware or
  cross-topology correspondence, RMSD matrices, RMSF, parallel execution, and
  streaming analysis adapters.

## Revision Notes

- v1: Add explicit direct RMSD plus fused zero-copy aligned RMSD with separate
  fit and measurement selections.
