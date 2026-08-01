# File-Backed Trajectory I/O

## Summary

Provide bounded format-agnostic file and stream I/O for production
fixed-topology trajectories without adding codec dependencies to the
foundational `kekule` crate.

## Behavior/API

- The `kekule-traj` package and `kekule_traj` Rust import root depend one-way
  on `kekule`. Trajectory types and streaming traits live at the package root;
  codecs and format-agnostic factories live under `kekule_traj::io`, without
  defining another topology, coordinate, or unit model.
- Opening requires an exact `Topology` plus an explicit
  `AtomOrderAssertion`; equal atom count is never sufficient evidence.
- Auto-detection combines a bounded signature inspection with extension
  evidence. Explicit format selection still validates the selected header.
- Sequential opening retains one handle and does not scan the full file.
  Indexed opening verifies every frame and stores checked compact offsets.
- Metadata distinguishes declared and verified counts, native precision,
  lossy encoding, available fields, format variants, and random-access
  capability. Sequential TRR metadata reports the precision verified so far
  and changes to mixed when a later frame uses the other scalar width.
- Default resource limits bound detection bytes, atoms, frames, record and
  scratch sizes, text lines, and index memory before allocation.
- Path writers use a temporary sibling and publish only through a successful
  consuming `finish()`. Production files must contain at least one frame. An
  empty finish or failed frame write cannot publish a destination, and
  unsupported state is rejected rather than discarded.
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
| TRR | GROMACS `GMX_trn_file`, f32/f64 XDR, optional box/velocity/force blocks | one explicit f32 or f64 precision with per-frame optional blocks | nm, ps, nm/ps, and kJ mol-1 nm-1 converted to Kekule units; lambda by explicit policy |
| XTC | magic 1995 and 2023, small uncompressed and ordinary compressed coordinates | magic 1995 or 2023 at explicit positive inverse-nm precision | nm and ps converted to Kekule units; lossy resolution nominally `1 / precision` nm |

The matrix is deliberately narrower than all historical files using these
extensions. A matching extension or atom count never expands the profile.

## Implementation Notes

- Kekule-owned source forbids unsafe code and has no native Chemfiles,
  C/C++, or CMake runtime.
- The trajectory crate is a publishable workspace member depending one-way on
  `kekule`; XYZ, DCD, and TRR have no production dependency beyond the
  semantic core, while XTC adds only the audited pure-Rust `molly` adapter.
- Codecs decode complete frames into reusable private scratch and publish
  transactionally into caller-owned `FrameBuffer` storage.
- Clean EOF is recognized only between frames. Partial headers, records, or
  payloads are truncation errors.
- At an exact frame or index limit, readers perform only a bounded
  frame-start/EOF probe. They neither decode the next frame nor grow the index
  before reporting a configured-limit failure.
- Sequential readers retain one file handle. Indexed readers retain one handle,
  fully verify every frame, and store only bounded checked offsets. The offset
  vector grows geometrically up to the smallest configured frame, entry, or
  byte capacity, so index construction does not copy the accumulated index for
  each frame. Random access is therefore O(one frame decode) after an O(file
  size) index build.
- Defaults cap atoms at 10,000,000; frames and index entries at 100,000,000;
  frame, record, and scratch payloads at 4 GiB; index storage at 800,000,000
  bytes; text lines/comments at 1 MiB; and format detection at 4096 bytes.
  Callers may lower these limits.

## Tests

- Public downstream tests prove codecs and external implementations satisfy
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
  position/vector allocations after warm-up. Instrumented N/N+1 streams prove
  limit probes do not decode or consume frame N+1, and restoration-seek fault
  streams prove indexed reads publish only after reader state is restored.
- Shared helper regressions build at least 100,000 offsets and bound capacity
  growth to logarithmically many events. Direct and atomic path writers reject
  empty finishes, and path tests verify that neither destinations nor temporary
  siblings remain.

## Benchmarks

- The lightweight `codec_throughput` bench records sequential frames/s and
  MiB/s, index time and offset memory, random-frame latency, allocation reuse,
  and writer MiB/s without a benchmark framework dependency.
- An informational Windows x86-64 release run on 2026-07-28 (12 logical
  workers) measured sequential frames/s, index ms, random us/frame, and writer
  MiB/s. For 256 frames x 64 atoms: XYZ 82,274 / 2.808 / 13.37 / 150.0; DCD
  1,730,431 / 0.081 / 0.49 / 901.2; TRR-f32 659,930 / 0.349 / 1.59 / 277.1;
  XTC 237,865 / 1.425 / 6.85 / 30.2. For 32 frames x 4096 atoms: XYZ 814 /
  37.988 / 954.32 / 85.6; DCD 26,125 / 0.700 / 25.24 / 515.2; TRR-f32 37,219 /
  0.675 / 22.10 / 220.1; XTC 3,785 / 12.136 / 386.18 / 38.5. A high-frame
  profile of 100,000 frames x 1 atom measured index construction at 49.174 ms
  for XYZ, 18.999 ms for DCD, 111.519 ms for TRR, and 13.376 ms for XTC, with
  781.25 KiB of logical offsets for each. These local measurements are evidence
  of the exercised paths, not portable performance guarantees.
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
- v4: Enforce projected index limits before parsing or growth, use bounded
  exact-limit EOF probes, and make DCD/TRR/XTC random reads restore all stream
  and codec state before transactional destination publication.
- v5: Use shared bounded geometric index reservation after successful parsing,
  require nonempty concrete and atomic writer finishes, and add a 100,000-frame
  small-frame index benchmark.
- v6: Hard-rename the pre-release codec companion to the Kekule namespace.
- v7: Fold the codecs into the broader `kekule-traj` package and move the
  public I/O surface under `kekule_traj::io` before the first crates.io
  release.
