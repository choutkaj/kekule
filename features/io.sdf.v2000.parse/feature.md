# SDF V2000 Parser

## Summary

Parse multi-record SDF V2000 input into small molecules using the Molfile V2000 parser for each record.

## Behavior/API

- Exposes `SdfParseOptions`, `SdfRecord`, `SdfParseError`, `read_sdf_v2000_str`, and `read_sdf_v2000_records`.
- Parses V2000 Molfile blocks, coordinates, common atom metadata, bond blocks, and SDF data fields.
- Preserves SDF title, program, comment, and data fields as molecule properties.
- Rejects V3000 input and malformed graph endpoints.
- Does not run sanitization or perception.

## Implementation Notes

- Delegates Molfile block parsing to `io.mol.v2000.parse`.
- Preserves common V2000 `M  CHG`, `M  ISO`, and `M  RAD` metadata.
- Treats raw parsing as separate from chemistry interpretation.

## Validation

- Current coverage is unit-test based.
- RDKit golden validation is planned through `validation.harness`.
- Fixtures live under `validation/features/io.sdf.v2000.parse/`.

## Out Of Scope

- SDF V3000 parsing, Molfile writing, sanitization, valence perception, ring detection, aromaticity, stereochemistry, and canonicalization.
- Runtime RDKit dependency.

## Revision Notes

- v1: Raw multi-record SDF V2000 parser.
- v2: Delegate Molfile parsing and preserve coordinates plus common atom metadata.
