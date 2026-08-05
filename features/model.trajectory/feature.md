# Shared-Topology Trajectory Streaming

## Summary

Provide the `kekule-traj` companion as Kekule's coherent trajectory layer:
ordered frames over one immutable topology, reusable storage, streaming file
I/O, and focused trajectory-oriented analysis workflows.

## Behavior/API

- The `kekule-traj` package and `kekule_traj` Rust import root depend one-way
  on the foundational `kekule` package. Trajectory types are exported from the
  `kekule_traj` root; codecs and path factories live under `kekule_traj::io`.
- The crate provides topology-bound frames with required
  positions and optional cell, velocities, forces, time, step, observation
  state, and frame metadata.
- In-memory `Trajectory` stores ordered frames directly rather than as
  `Vec<Model>`.
- `TrajectoryReader` fills a caller-owned `FrameBuffer`;
  `SeekableTrajectoryReader` adds explicit random access; `TrajectoryWriter`
  validates topology and supported state.
- `FrameBuffer::replace_from_data` transactionally publishes complete borrowed
  decoder state after validating every field, converts units once, reuses all
  dense-array allocations, clears absent optional state, and leaves the
  destination unchanged on failure.
- `FrameBuffer::reset_dynamic_state` clears cell, velocities, forces, time,
  step, observation, and properties without replacing its positions allocation;
  coordinate-only readers use this reset after filling positions.
- Frame and buffer views can be consumed by structural analyses and prepared
  potentials without owned-model construction or coordinate copying.
- Owned frames and finite in-memory trajectories remap through exact checked
  topology lineage while preserving order, positions, cells, velocities,
  forces, time, step, observations, and properties.
- `FrameBuffer::copy_remapped_from` transactionally copies borrowed source
  frame state into an exact target-bound buffer without constructing a model
  and reuses position, velocity, and force allocations.
- `AtomOrderAssertion` distinguishes an explicit semantic-order proof from the
  caller's explicit assertion that a topology-free file uses authoritative
  topology order, and remains bound to exact topology identity.
- Non-exhaustive `TrajectoryFormat`, `TrajectoryIoOperation`, and
  `TrajectoryCodecErrorKind` values plus cloneable I/O/codec contexts report
  typed format, operation, source, frame, byte offset, count, and underlying
  I/O information.
- General single-configuration geometry, selections, alignment, and potential
  kernels remain in `kekule`. `kekule-traj` supplies zero-copy frame views and
  owns trajectory-scale orchestration. Initial superposition and direct/fused
  RMSD workflows are tracked by `algo.trajectory-superposition` and
  `algo.trajectory-rmsd`; slicing, distance/contact series, RMSF, and related
  MDTraj-like workflows remain incremental companion features.

## Implementation Notes

- Ordinary trajectories retain one fixed exact topology. Reactive trajectories
  remain a future segmented-topology concept.
- A minimal in-memory/reference reader and writer validate the streaming
  contracts. Production codecs are part of the same `kekule-traj` crate and
  are tracked independently as `io.trajectory.file`, `.xyz`, `.dcd`, `.trr`,
  and `.xtc`.

## Tests

- Tests cover complete optional arrays, variable cells, topology and
  atom-count mismatches, buffer allocation reuse, end-of-stream, random-access
  separation, writer rejection, reference round trips, and complete dynamic
  state clearing by coordinate-only readers while preserving the positions
  allocation.
- Transformation regressions cover complete owned and borrowed frame transfer,
  frame-index error context, target-buffer identity rejection, unchanged
  positions, cell, vectors, time, step, observation, and properties after a
  later validation failure, stale optional-state clearing after positions-only
  remaps, and stable dense-array pointers and capacities over repeated remaps.
- Downstream tests prove external applications can implement sequential,
  seekable, and writer traits and publish complete data without private access.
- Cross-crate tests prove `kekule` potentials and the
  `kekule_potentials::dreiding` adapter consume `kekule-traj` frame and buffer
  views without copying coordinates.
- Publication regressions cover late validation failure, complete destination
  transactionality, stale property/optional-field clearing, exact order-token
  identity, typed error context, and stable position/vector pointers and
  capacities after warm-up.

## Out Of Scope

- Reactive trajectories, dynamics integration, neighbor lists, and higher-level
  MDTraj-like slicing, distance/contact, and RMSF operations not tracked by
  focused implemented features.

## Revision Notes

- v1: Track the fixed-topology frame and streaming-I/O contract.
- v2: Implement owned frames, reusable frame buffers, fixed-topology in-memory
  trajectories, sequential/seekable/writer traits, memory adapters, and an
  atom-order-asserted coordinate-only reference reader.
- v3: Centralize reusable-buffer dynamic-state reset so coordinate-only reads
  clear stale properties and every optional field without reallocating
  positions.
- v4: Add exact-lineage owned frame and finite trajectory remapping plus
  transactional allocation-reusing borrowed-frame copies into target buffers.
- v5: Strengthen reusable-buffer tests for complete destination transactionality
  and stale optional-state clearing after positions-only remaps.
- v6: Add complete transactional borrowed-data publication, explicit exact
  atom-order binding helpers, shared format identity, typed file/codec error
  context, and downstream codec-trait implementation coverage.
- v7: Record the supported companion-codec integration contract and its use of
  the same exact-topology, transactional reusable-buffer API.
- v8: Move the full trajectory model and production codecs into the one-way
  `kekule-traj` companion, remove `kekule::trajectory`, and reserve this crate
  as the home for future trajectory-oriented workflows over Kekule kernels.
- v9: Establish the companion analysis namespace and delegate implemented
  trajectory superposition and RMSD behavior to focused feature contracts.
