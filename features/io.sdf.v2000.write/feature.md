# SDF V2000 Writer

## Summary

Write small molecules as SDF V2000 records while preserving molecule data fields.

## Behavior/API

- Exposes `write_sdf_v2000`.
- Writes Molfile V2000 blocks followed by `sdf.field.*` properties as SDF data fields.
- Does not run sanitization, canonicalization, or perception.

## Implementation Notes

- The writer uses existing molecule properties for title/program/comment and SDF fields.
- Records are emitted in input slice order.

## Validation

- Unit tests cover multi-record round trips and data fields.
- RDKit reference generator support is included for fixture/golden generation.

## Out Of Scope

- Compression, streaming, V3000, and canonical output ordering beyond current graph order.

## Revision Notes

- v1: Initial SDF V2000 writer.
