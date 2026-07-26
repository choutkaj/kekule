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
- `FrameBuffer::reset_dynamic_state` clears cell, velocities, forces, time,
  step, observation, and properties without replacing its positions allocation;
  coordinate-only readers use this reset after filling positions.
- Frame and buffer views can be consumed by structural analyses and prepared
  potentials without owned-model construction or coordinate copying.

## Implementation Notes

- Ordinary trajectories retain one fixed exact topology. Reactive trajectories
  remain a future segmented-topology concept.
- A minimal in-memory/reference reader and writer validate the streaming
  contracts; production binary codecs are not part of this feature.

## Tests

- Tests cover complete optional arrays, variable cells, topology and
  atom-count mismatches, buffer allocation reuse, end-of-stream, random-access
  separation, writer rejection, reference round trips, and complete dynamic
  state clearing by coordinate-only readers while preserving the positions
  allocation.

## Out Of Scope

- XTC, DCD, TRR, NetCDF, reactive trajectories, dynamics integration, and
  neighbor lists.

## Revision Notes

- v1: Track the fixed-topology frame and streaming-I/O contract.
- v2: Implement owned frames, reusable frame buffers, fixed-topology in-memory
  trajectories, sequential/seekable/writer traits, memory adapters, and an
  atom-order-asserted coordinate-only reference reader.
- v3: Centralize reusable-buffer dynamic-state reset so coordinate-only reads
  clear stale properties and every optional field without reallocating
  positions.
