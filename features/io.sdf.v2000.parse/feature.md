# SDF Document Parsing and Interpretation

## Summary

Parse ordered SDF records into a loss-preserving `SdfDocument`, then interpret
their Molfile documents into canonical `SdfRecord` values plus qualified
per-record reports.

## Behavior/API

- Exposes `SdfDocument`, `SdfRecordDocument`, `SdfDataField`, canonical
  `SdfRecord`, `SdfInterpretation`, `sdf::parse_str`, and `sdf::interpret`.
- Each raw record owns a `MolfileDocument` and ordered raw data fields with source
  lines. Molfile versions are auto-detected per record.
- Interpretation returns ordered records plus an `SdfInterpretationReport`.
  Each `SdfRecordInterpretationReport` qualifies the underlying Molfile
  atom/bond mappings with its record index and source start line.
- Headers and SDF fields are record metadata and are never injected into
  `Molecule::props`.
- The final record requires its own `$$$$` delimiter unless
  `SdfParseOptions::allow_missing_final_delimiter` is enabled; an earlier
  record delimiter never makes a later unterminated record valid.
- `SdfParseOptions` also bounds total input bytes, record count, and bytes per
  record. Defaults support the repository's broad corpora while preventing
  unbounded document and per-record allocation.
- Data values continue until a blank line or the record's `$$$$` delimiter, so
  a value line beginning with `>` is preserved as data rather than misread as a
  new field header. Bare end-of-input does not implicitly terminate a field.
- Nonblank text after the Molfile block must be a data-field header, and every
  field must end at a blank line or record delimiter; malformed post-CTAB
  content is rejected rather than discarded.
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
- v13: Return a first-class SDF interpretation result with qualified per-record
  Molfile reports and source-to-canonical mappings.
- v14: Enforce the final delimiter per record and preserve data-value lines
  beginning with `>` until the field's blank-line terminator.
- v15: Add configurable whole-input, record-count, and per-record resource
  limits with broad-corpus-safe defaults.
- v16: Reject unstructured post-CTAB text and data fields without their required
  blank-line terminator instead of silently dropping or ambiguously accepting
  malformed record content.
- v17: Accept the common, unambiguous SDF form where `$$$$` directly terminates
  the final data field while still rejecting a field that reaches bare
  end-of-input without a blank line.
