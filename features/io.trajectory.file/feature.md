# File-Backed Trajectory I/O

## Summary

Provide bounded format-agnostic file and stream I/O for production
fixed-topology trajectories without adding codec dependencies to the
foundational `molecular` crate.

## Behavior/API

- `molecular-trajectory-io` depends one-way on `molecular` and implements its
  sequential, seekable, and writer traits without defining another topology,
  frame, coordinate, or unit model.
- Opening requires an exact `Topology` plus an explicit
  `AtomOrderAssertion`; equal atom count is never sufficient evidence.
- Auto-detection combines a bounded signature inspection with extension
  evidence. Explicit format selection still validates the selected header.
- Sequential opening retains one handle and does not scan the full file.
  Indexed opening verifies every frame and stores checked compact offsets.
- Metadata distinguishes declared and verified counts, native precision,
  lossy encoding, available fields, format variants, and random-access
  capability.
- Default resource limits bound detection bytes, atoms, frames, record and
  scratch sizes, text lines, and index memory before allocation.
- Path writers use a temporary sibling and publish only through a successful
  consuming `finish()`. Unsupported state is rejected rather than discarded.
- XYZ, DCD, TRR, and XTC are tracked independently. Amber NetCDF, multi-model
  PDB, GRO/G96, Amber ASCII, LAMMPS dump, and compressed wrappers remain
  structured unsupported formats until their own contracts land.
- The current format-agnostic factory dispatches the complete XYZ vertical
  slice. Recognized DCD, TRR, and XTC signatures remain typed unsupported
  results until their independently tracked implementations land.

## Implementation Notes

- Molecular-owned source forbids unsafe code and has no native Chemfiles,
  C/C++, or CMake runtime.
- The companion crate is a publishable workspace member depending one-way on
  `molecular`; format detection and XYZ have no production dependency beyond
  the semantic core.
- Codecs decode complete frames into reusable private scratch and publish
  transactionally into caller-owned `FrameBuffer` storage.
- Clean EOF is recognized only between frames. Partial headers, records, or
  payloads are truncation errors.

## Tests

- Public downstream tests prove the companion crate implements
  `TrajectoryReader`, `SeekableTrajectoryReader`, and `TrajectoryWriter`.
- Detection, limits, topology/order binding, metadata, transactionality,
  allocation reuse, indexing, atomic finish, and malformed/truncated input
  receive focused regressions.
- Public bounded detection is covered by the `trajectory_detection` fuzz
  target, which also asserts that the inspected stream position is restored.
- Each supported format supplies provenance-pinned interoperability fixtures,
  sequential-versus-indexed equality, and strict writer round trips.
- Bounded fuzz targets cover detection and every format parser or adapter.

## Benchmarks

- Record sequential frames/s and MB/s, index time and memory, random-frame
  latency, allocation reuse after warm-up, and writer throughput for
  representative small and large atom counts.
- External-tool comparisons are informational and never imply broad format or
  corpus parity.

## Out Of Scope

- Reactive trajectories, topology inference, atom correspondence from equal
  count, atom subsets, concatenation, persisted indices, compression wrappers,
  async/network I/O, playback, caching, renderer packing, and project storage.

## Revision Notes

- v1: Register the planned file-backed trajectory I/O boundary and milestone
  acceptance contract.
- v2: Add the experimental companion crate, bounded signature/extension
  detection, metadata/reports/limits, one-handle sequential and indexed
  wrappers, strict atomic path writing, the XYZ dispatch vertical slice, and
  bounded detection and XYZ fuzz targets.
