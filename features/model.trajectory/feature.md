# Shared-Topology Trajectory Streaming

## Summary

Represent ordered frames over one immutable topology and process large
trajectories through reusable caller-owned buffers.

## Behavior/API

- The planned trajectory module provides topology-bound frames with required
  positions and optional cell, velocities, forces, time, step, observation
  state, and frame metadata.
- In-memory `Trajectory` stores ordered frames directly rather than as
  `Vec<Model>`.
- `TrajectoryReader` fills a caller-owned `FrameBuffer`;
  `SeekableTrajectoryReader` adds explicit random access; `TrajectoryWriter`
  validates topology and supported state.
- Frame and buffer views can be consumed by structural analyses and prepared
  potentials without owned-model construction or coordinate copying.

## Implementation Notes

- Ordinary trajectories retain one fixed exact topology. Reactive trajectories
  remain a future segmented-topology concept.
- A minimal in-memory/reference reader and writer validate the streaming
  contracts; production binary codecs are not part of this feature.

## Tests

- Planned tests cover complete optional arrays, variable cells, topology and
  atom-count mismatches, buffer allocation reuse, end-of-stream, random-access
  separation, writer rejection, and reference round trips.

## Out Of Scope

- XTC, DCD, TRR, NetCDF, reactive trajectories, dynamics integration, and
  neighbor lists.

## Revision Notes

- v1: Track the fixed-topology frame and streaming-I/O contract.
