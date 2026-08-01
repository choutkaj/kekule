# Trajectory Rigid Superposition

## Summary

Superpose every frame of one finite fixed-topology trajectory onto a selected
reference frame through Kekule's proper rigid-alignment kernel.

## Behavior/API

- `Trajectory::superpose_to_frame` performs uniform-weight fitting;
  `superpose_to_frame_with_options` accepts the existing Kabsch weighting and
  periodic-cell policy through `analysis::SuperpositionOptions`.
- The fit selection is topology-bound and determines one proper rigid transform
  per frame. The transform is applied to every atom, not only selected atoms.
- Positions are rotated and translated. Velocities, forces, and periodic-cell
  basis vectors are rotated without translation. Topology, time, step,
  observation provenance, and frame properties are preserved.
- The operation is transactional across the complete trajectory: all fits and
  replacement frames validate before the original trajectory changes.
- `SuperpositionReport` retains the reference index and one applied
  `RigidAlignment`, including transform and fit RMSD, per frame.
- Periodic configurations are rejected by default. Explicit stored-coordinate
  handling performs no imaging, wrapping, unwrapping, minimum-image correction,
  or molecule reconstruction; retained cells rotate with their frames.
- Reference bounds and per-frame fit/transformation failures retain structured
  frame context.

## Implementation Notes

- Fitting delegates to `kekule::alignment::kabsch_with_options` over zero-copy
  frame `ModelView` values.
- The implementation first computes every alignment, then builds and validates
  a complete replacement trajectory before one final publication. Runtime is
  O(frames × (fit atoms + all atoms)); transactional replacement temporarily
  owns one additional complete set of frame state.
- Source raw-coordinate observation strings remain preserved provenance rather
  than being rewritten as generated text.

## Tests

- Focused regressions cover rigid translation/rotation recovery, complete
  position/vector/cell transformation, metadata preservation, periodic opt-in,
  reference bounds, late-frame rank failure, and complete source immutability
  after failure.
- A downstream integration test compiles the public options, report, method,
  and split superpose-then-RMSD workflow.

## Out Of Scope

- Periodic imaging or molecule reconstruction, progressive fitting, non-rigid
  fitting, cross-topology correspondence, parallel execution, and streaming
  transformation adapters.

## Revision Notes

- v1: Add transactional finite-trajectory superposition with complete dynamic
  geometric-state transformation and per-frame alignment reporting.
