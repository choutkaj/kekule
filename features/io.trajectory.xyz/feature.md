# XYZ Trajectory Codec

## Summary

Read and write strict multi-frame XYZ as the transparent text vertical slice of
file-backed trajectory I/O.

## Behavior/API

- Every open requires an exact topology and explicit atom-order assertion.
- Each frame has a checked nonzero constant atom count and an element sequence
  matching authoritative topology order.
- Coordinates are finite and use an explicit length-unit policy; the documented
  default is angstrom because ordinary XYZ carries no reliable unit tag.
- Comments, ordinary whitespace, final-newline absence, and line endings are
  handled within configured text limits. Extended XYZ schemas are not inferred.
- Writers use deterministic locale-independent precision and reject cell,
  velocity, force, time, step, observation, or properties they cannot preserve.

## Implementation Notes

- Parsing is line-bounded and distinguishes clean EOF before a frame from
  truncation after a frame begins.
- Later coordinate-only frames clear all stale optional destination state.

## Tests

- Cover comments, whitespace, CRLF/LF, symbols, units, malformed counts,
  changing counts, element mismatch, non-finite values, truncation, limits,
  transactionality, allocation reuse, indexing, and multi-frame round trips.
- Read committed independently generated fixtures and validate writer output
  with an independent reader.

## Benchmarks

- Record sequential parsing, index construction, random access, allocation
  reuse, and deterministic writer throughput for small and large frames.

## Out Of Scope

- Extended XYZ lattice/property schemas, inferred time or step, and topology
  construction from element labels.

## Revision Notes

- v1: Register the planned strict multi-frame XYZ read/write contract.
