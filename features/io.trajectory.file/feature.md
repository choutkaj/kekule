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
  consuming `finish()`. A failed frame write poisons the writer, so a later
  `finish()` cannot publish a valid prefix. Unsupported state is rejected
  rather than discarded.
- The supported initial factory profiles are strict multi-frame XYZ, common
  CHARMM/NAMD/OpenMM DCD, GROMACS TRR with f32/f64 XDR payloads, and GROMACS
  XTC with 1995/2023 magic values.
- Amber NetCDF, multi-model PDB, GRO/G96, Amber ASCII, LAMMPS dump, reactive
  trajectories, and compressed wrappers return structured unsupported-format
  or unsupported-variant errors.

### Supported compatibility matrix

| Format | Supported read profile | Supported write profile | Native units and precision |
|---|---|---|---|
| XYZ | strict constant-count multi-frame element/x/y/z text | deterministic strict multi-frame text | configured length unit, angstrom by default; decimal text |
| DCD | 32-bit little/big-endian `CORD`, common cell records, all-atom and fixed-atom frames | canonical all-atom `CORD`, little or big endian, optional cell | angstrom positions/cell; f32 coordinates; step from `ISTART`/`NSAVC`; time only by explicit `DELTA` policy |
| TRR | GROMACS `GMX_trn_file`, f32/f64 XDR, optional box/velocity/force blocks | one explicit f32 or f64 precision with per-frame optional blocks | nm, ps, nm/ps, and kJ mol-1 nm-1 converted to Molecular units; lambda by explicit policy |
| XTC | magic 1995 and 2023, small uncompressed and ordinary compressed coordinates | magic 1995 or 2023 at explicit positive inverse-nm precision | nm and ps converted to Molecular units; lossy resolution nominally `1 / precision` nm |

The matrix is deliberately narrower than all historical files using these
extensions. A matching extension or atom count never expands the profile.

## Implementation Notes

- Molecular-owned source forbids unsafe code and has no native Chemfiles,
  C/C++, or CMake runtime.
- The companion crate is a publishable workspace member depending one-way on
  `molecular`; XYZ, DCD, and TRR have no production dependency beyond the
  semantic core, while XTC adds only the audited pure-Rust `molly` adapter.
- Codecs decode complete frames into reusable private scratch and publish
  transactionally into caller-owned `FrameBuffer` storage.
- Clean EOF is recognized only between frames. Partial headers, records, or
  payloads are truncation errors.
- Sequential readers retain one file handle. Indexed readers retain one handle,
  fully verify every frame, and store only bounded checked offsets; random
  access is therefore O(one frame decode) after an O(file size) index build.
- Defaults cap atoms at 10,000,000; frames and index entries at 100,000,000;
  frame, record, and scratch payloads at 4 GiB; index storage at 800,000,000
  bytes; text lines/comments at 1 MiB; and format detection at 4096 bytes.
  Callers may lower these limits.

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
- Regression tests prove exact-limit EOF behavior, late-failure destination
  transactionality, stale optional-state clearing, frame/offset error context,
  unsupported writer-state rejection, poisoned atomic output, and stable
  position/vector allocations after warm-up.

## Benchmarks

- The lightweight `codec_throughput` bench records sequential frames/s and
  MiB/s, index time and offset memory, random-frame latency, allocation reuse,
  and writer MiB/s without a benchmark framework dependency.
- An informational Windows x86-64 release run on 2026-07-27 (12 logical
  workers) measured, for 64 frames x 256 atoms, sequential/index/random/writer:
  XYZ 68,505 frames/s, 4.814 ms, 19.88 us, 116.5 MiB/s; DCD 847,177 frames/s,
  0.116 ms, 0.59 us, 918.2 MiB/s; TRR-f32 829,929 frames/s, 0.343 ms, 1.40 us,
  232.9 MiB/s; and XTC 162,914 frames/s, 1.116 ms, 3.37 us, 24.3 MiB/s. These
  local measurements are evidence of the exercised paths, not portable
  performance guarantees.
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
- v3: Support the bounded initial XYZ/DCD/TRR/XTC compatibility matrix; harden
  checked lengths, offsets, EOF boundaries, late-failure transactionality,
  reusable scratch, error context, poisoned atomic publication, checked XTC
  decoding, fuzz corpora, interoperability fixtures, and lightweight
  performance evidence.
