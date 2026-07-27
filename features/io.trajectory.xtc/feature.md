# XTC Trajectory Codec

## Summary

Read and write compact lossy GROMACS XTC through a private adapter over an
audited pure-Rust codec while retaining Molecular validation, limits, units,
metadata, and errors.

## Behavior/API

- Validate supported 1995 and 2023 magic variants, constant atom count,
  triclinic box, nonnegative representable step, finite time, and compressed
  frame structure.
- Convert nanometre positions/boxes and picosecond time once to Molecular
  units, while exposing native lossy precision in metadata and reports.
- Sequential access does not scan the full file. Indexed access validates
  checked offsets and decodes only the requested frame.
- Writers require a positive finite physical coordinate resolution and verify
  round-trip error against its declared quantization bound.
- Malformed small- and large-atom frames return structured corruption or
  truncation errors without panicking.

## Implementation Notes

- The selected dependency is `molly` 0.6.1, crates.io source revision
  `600e27a5ce17ad3a1579d54a943e1cf6f1ccd485`: MIT licensed, default-feature
  dependency-free, pure Rust, and MSRV 1.74.1 versus this workspace's 1.89.
- Audit findings: its public convenience path can allocate from an untrusted
  compressed byte count and contains assertions/expectations reachable from
  malformed data. Its optional file-buffered fast path contains one
  `get_unchecked` access. The private Molecular adapter will not call that
  buffered path.
- The adapter preflights bounded frame lengths and structural fields, uses only
  the allocation-reusing unbuffered low-level path, and contains dependency
  panics as typed codec errors. `molly` types do not appear in public APIs.
- Molecular-owned source remains `#![forbid(unsafe_code)]`; there is no native
  library, C/C++, or CMake runtime.

## Tests

- Cover both magic values, small and large atom paths, triclinic cells,
  step/time, precision, corrupt compressed blocks, count/limit failures,
  transactionality, allocation reuse, indexed equality, and precision-aware
  round trips.
- Gate supported status on GROMACS and independent-reader fixtures plus a
  bounded fuzz target around adapter preflight and decode containment.

## Benchmarks

- Record sequential decode, index construction, random access, scratch reuse,
  and writer throughput for small and representative compressed frames.

## Out Of Scope

- Leaking `molly` APIs publicly, atom-subset semantics, unbounded convenience
  methods, or claiming variants beyond the tested magic/precision matrix.

## Revision Notes

- v1: Register the planned XTC contract and pin the initial `molly` 0.6.1
  dependency audit.
