# TRR Trajectory Codec

## Summary

Read and write the tested GROMACS TRR profile through a private bounded
pure-Rust XDR layer.

## Behavior/API

- Validate XDR framing, magic/version, single- or double-precision payloads,
  constant atom count, and per-frame block sizes.
- Preserve positions, triclinic box, optional velocities and forces, time,
  nonnegative step, and explicitly handled lambda metadata.
- Convert native GROMACS units once to Molecular `f64` model units.
- Writers require explicit scalar precision and stable field presence, and
  reject every unsupported block or frame field instead of dropping it.
- Indexed access scans checked frame offsets without decoding full arrays.

## Implementation Notes

- XDR arithmetic and padding are checked before allocation or seek.
- Recognized optional blocks not representable by the shipped frame contract
  are either preserved through a documented policy or rejected explicitly.

## Tests

- Cover f32/f64, every supported optional-vector combination, triclinic cells,
  time, step, lambda, block consistency, non-finite values, oversized/truncated
  XDR, transactionality, allocation reuse, indexed equality, and strict round
  trips.
- Gate the supported profile on GROMACS and independent-reader fixtures.

## Benchmarks

- Record f32/f64 sequential decode, index construction, random access,
  allocation reuse, and writer throughput.

## Out Of Scope

- Arbitrary unmodeled TRR blocks, negative steps, changing atom counts, and
  opaque metadata loss.

## Revision Notes

- v1: Register the planned pure-Rust TRR/XDR read/write contract.
