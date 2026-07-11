# SDF V2000 Writer

## Summary

Write small molecules as SDF V2000 records while preserving molecule data fields.

## Behavior/API

- Exposes `sdf::write_v2000`.
- Writes Molfile V2000 blocks followed by `sdf.field.*` properties as SDF data fields.
- Inherits exact radical and supported source bond-stereo mark handling from the Molfile V2000 writer.
- Inherits V2000 atom-block valence/no-implicit output and structured rejection
  of quartet/quintet radical multiplicity.
- Does not run sanitization, canonicalization, or perception.

## Implementation Notes

- The writer uses existing molecule properties for title/program/comment and SDF fields.
- Records are emitted in input slice order.
- Unsupported Molfile representations in any record return a structured error and no SDF text is returned.

## Validation

- Unit tests cover multi-record round trips, multiline data fields, and Molfile metadata symmetry.
- RDKit-generated goldens compare SDF writer records for external PubChem fixtures.

## Out Of Scope

- Compression, streaming, V3000, and canonical output ordering beyond current graph order.

## Revision Notes

- v1: SDF V2000 writer.
- v2: Preserve exact Molfile radical and supported bond-stereo semantics in SDF records.
- v3: Move the public writer API under the `sdf` facade.
- v4: Add PubChem-100k as required broad-corpus validation evidence.
- v5: Inherit first-class source bond stereo marks from Molfile V2000 writing.
- v6: Inherit atom-block valence/no-implicit output and lossless high-spin
  radical rejection from Molfile V2000 writing.
