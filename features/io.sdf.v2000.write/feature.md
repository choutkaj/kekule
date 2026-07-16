# SDF V2000 Writer

## Summary

Write canonical `SdfRecord` values as ordered SDF V2000 records.

## Behavior/API

- Exposes `sdf::write_v2000`.
- Accepts `&[SdfRecord]`, writes each record title and ordered `SdfDataField`
  values, and never reads SDF metadata from molecule properties.
- Inherits exact radical and supported source bond-stereo mark handling from the Molfile V2000 writer.
- Inherits V2000 atom-block valence/no-implicit output and structured rejection
  of quartet/quintet radical multiplicity.
- Rejects titles with line breaks, malformed field names, field values with
  blank lines or carriage returns, and field values containing a `$$$$`
  delimiter line instead of emitting SDF text that cannot round-trip safely.
- Does not run sanitization, canonicalization, or perception.

## Implementation Notes

- Molfile headers and SDF fields remain record/document concerns.
- Records are emitted in input slice order.
- Unsupported Molfile representations in any record return a structured error and no SDF text is returned.

## Validation

- Unit tests cover multi-record round trips, multiline data fields, and Molfile metadata symmetry.
- RDKit-generated goldens compare SDF writer records for external PubChem fixtures.
- PubChem-1k is required baseline evidence; manifest-backed broader corpora
  remain available for deliberate local parity checks.

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
- v7: Hard break to canonical `SdfRecord` input and record-owned metadata.
- v8: Make the committed smoke corpus the CI-reproducible required evidence
  tier while retaining every ignored corpus on demand.
- v9: Reject record metadata that is not representable by the line-oriented SDF
  grammar and preserve value lines beginning with `>` on round trip.
- v10: Use PubChem-1k as the required baseline validation corpus after retiring the former smoke corpus from public validation.
