# SDF V2000 Parser

## Summary

Parse multi-record SDF V2000 input into small molecules using the shared core graph.

## Behavior/API

- Exposes `SdfParseOptions`, `SdfRecord`, `SdfParseError`, `read_sdf_v2000_str`, and `read_sdf_v2000_records`.
- Parses V2000 counts, atom blocks, bond blocks, and SDF data fields.
- Preserves SDF title, program, comment, and data fields as molecule properties.
- Rejects V3000 input and malformed graph endpoints.
- Does not run sanitization or perception.

## Implementation Notes

- Converts SDF one-based atom indices to stable `AtomId`s.
- Supports the initial documented bond order mappings.
- Treats raw parsing as separate from chemistry interpretation.
- Some V2000 metadata records such as `M  CHG` and `M  ISO` are fixture pressure cases for future parser expansion.

## Validation

- Current coverage is unit-test based.
- RDKit golden validation is planned through `validation.harness`.
- Fixtures live under `validation/features/io.sdf.v2000.parse/`.

## Out Of Scope

- SDF V3000 parsing, Molfile writing, sanitization, valence perception, ring detection, aromaticity, stereochemistry, and canonicalization.
- Runtime RDKit dependency.

## Revision Notes

- v1: Initial raw multi-record SDF V2000 parser.
