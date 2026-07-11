# SDF V2000 Parser

## Summary

Parse multi-record SDF V2000 input into small molecules using the Molfile V2000 parser for each record.

## Behavior/API

- Exposes `sdf::{SdfParseOptions, SdfRecord, SdfParseError, read_v2000_str, read_v2000_records}`.
- Parses V2000 Molfile blocks, coordinates, common atom metadata, supported source bond stereo marks, and SDF data fields.
- Preserves SDF title, program, comment, and data fields as molecule properties.
- Rejects V3000 input and malformed graph endpoints.
- Returns located `SdfParseError` values for truncated, non-ASCII, overflowing, or inconsistent V2000 structures and M records.
- Does not run sanitization or perception.

## Implementation Notes

- Delegates Molfile block parsing to `io.mol.v2000.parse`.
- Preserves common V2000 `M  CHG`, `M  ISO`, and exact `M  RAD` metadata.
- Inherits exact V2000 single-/double-bond source mark parsing from the Molfile parser.
- Inherits V2000 atom-block valence/no-implicit semantics and valence-implied
  tetrahedral hydrogen-carrier preservation.
- Inherits fixed-width V2000 count and bond parsing from the Molfile parser.
- Uses checked record/block offsets and validates declared M-record pair counts before mutation.
- Treats raw parsing as separate from chemistry interpretation.

## Validation

- Unit tests cover multi-record parsing, data fields, radical/source-stereo metadata inherited from Molfile blocks, malformed blocks and M records, and explicit parse-without-perceive behavior.
- The standalone `sdf_v2000` fuzz target checks panic safety and successful parse/write/parse paths.
- RDKit-generated goldens compare SDF records and preserved properties for external PubChem fixtures.

## Out Of Scope

- SDF V3000 parsing, Molfile writing, sanitization, valence perception, ring detection, aromaticity, stereochemistry, and canonicalization.
- Runtime RDKit dependency.

## Revision Notes

- v1: Raw multi-record SDF V2000 parser.
- v2: Delegate Molfile parsing and preserve coordinates plus common atom metadata.
- v3: Handle fixed-width three-digit V2000 count and bond fields in larger external PubChem records.
- v4: Inherit panic-free checked V2000 parsing and add SDF fuzz coverage.
- v5: Inherit exact V2000 radical multiplicity and supported bond stereo parsing.
- v6: Move the public parser API under the `sdf` facade.
- v7: Add PubChem-100k as required broad-corpus validation evidence.
- v8: Inherit first-class source bond stereo marks from Molfile V2000 parsing.
- v9: Inherit atom-block valence/no-implicit semantics and valence-implied
  tetrahedral hydrogen-carrier preservation from Molfile V2000 parsing.
