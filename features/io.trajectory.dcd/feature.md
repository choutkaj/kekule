# DCD Trajectory Codec

## Summary

Read common production CHARMM/NAMD/OpenMM-style DCD files and write one
deterministic little-endian compatibility profile.

## Behavior/API

- Support little- and big-endian 32-bit Fortran record markers, `CORD`
  headers, bounded title/atom-count records, all-atom frames, common fixed-atom
  trajectories, and documented unit-cell variants.
- Validate record sizes and trailing markers, constant atom count, declared
  versus indexed frame count, fixed/free atom lists, and every checked offset.
- Preserve defensible step metadata. Time remains absent unless an explicit
  `DcdTimePolicy` establishes its units and conversion.
- Writers emit positions and optional cells with explicit start step, save
  interval, and time convention; they do not emit fixed-atom optimization.
- Unknown record or cell dialects return `UnsupportedVariant`.

## Implementation Notes

- Coordinate values use the documented angstrom profile and convert once into
  Molecular model units.
- Sequential readers retain one handle and do not build an index; indexed
  readers structurally scan records without materializing frames.

## Tests

- Cover both endiannesses, unit cells, fixed atoms, producer variants, count
  disagreement, marker corruption, truncation, overflow/limits,
  transactionality, allocation reuse, indexing equality, and strict round trips.
- The supported producer/version matrix is gated by provenance-pinned fixtures
  and independent writer-output reads.

## Benchmarks

- Record sequential decode, structural indexing, random reads, allocation reuse,
  and canonical writer throughput for all-atom and fixed-atom inputs.

## Out Of Scope

- Undocumented historical dialects, 64-bit record markers, automatic time-unit
  invention, and fixed-atom writer optimization.

## Revision Notes

- v1: Register the planned common DCD compatibility profile.
