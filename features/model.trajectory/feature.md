# Shared-Topology Trajectory Streaming

## Summary

Represent ordered frames over one immutable topology and process large
trajectories through reusable caller-owned buffers.

## Behavior/API

- The `trajectory` module provides topology-bound frames with required
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
  `TrajectoryCodecErrorKind` values plus cloneable I/O/codec contexts let
  companion codecs report typed format, operation, source, frame, byte offset,
  count, and underlying I/O information without a dependency cycle.

## Implementation Notes

- Ordinary trajectories retain one fixed exact topology. Reactive trajectories
  remain a future segmented-topology concept.
- A minimal in-memory/reference reader and writer validate the streaming
  contracts. Production codecs live in the one-way
  `molecular-trajectory-io` companion crate and are tracked separately.

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
- Downstream tests prove an external companion crate can implement sequential,
  seekable, and writer traits and publish complete data without private access.
- Publication regressions cover late validation failure, complete destination
  transactionality, stale property/optional-field clearing, exact order-token
  identity, typed error context, and stable position/vector pointers and
  capacities after warm-up.

## Out Of Scope

- File codec implementations, reactive trajectories, dynamics integration, and
  neighbor lists.

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
