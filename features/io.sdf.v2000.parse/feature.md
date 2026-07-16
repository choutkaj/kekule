# SDF Document Parsing and Interpretation

## Summary

Parse ordered SDF records into a loss-preserving `SdfDocument`, then interpret
their Molfile documents into canonical `SdfRecord` values.

## Behavior/API

- Exposes `SdfDocument`, `SdfRecordDocument`, `SdfDataField`, canonical
  `SdfRecord`, `sdf::parse_str`, and `sdf::interpret`.
- Each raw record owns a `MolfileDocument` and ordered raw data fields with source
  lines. Molfile versions are auto-detected per record.
- Interpretation returns one `SmallMolecule` plus record title/data fields and
  keeps syntax errors distinct from chemical interpretation errors.
- Headers and SDF fields are record metadata and are never injected into
  `Molecule::props`.
- Parsing and interpretation never sanitize or run chemical perception.

## Validation

- Tests cover ordered records, raw headers/data fields, record round trips,
  V2000 metadata/stereo inheritance, malformed input, and absence of implicit
  perception. Existing external corpora remain the reference evidence.
- All ignored non-smoke corpora remain available for explicit local-only
  parity checks but are not selected as required routine evidence.

## Out Of Scope

- Sanitization, canonicalization, and automatic promotion of data fields to
  molecule properties.

## Revision Notes

- v1-v9: Direct V2000 SDF-to-molecule reader and expanded metadata coverage.
- v10: Hard break to `SdfDocument` parse/interpret and canonical record metadata.
- v11: Inherit fixed-width V2000 counts parsing for adjacent three-digit atom
  and bond counts in raw SDF record documents.
- v12: Make the committed smoke corpus the CI-reproducible required evidence
  tier while retaining every ignored corpus on demand.
