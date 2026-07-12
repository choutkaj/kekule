# mmCIF Document Parser

## Summary

Parse a structural mmCIF file into a format-level document before chemical interpretation.

## Behavior/API

- Exposes `mmcif::parse_str` returning `MmcifDocument`.
- Preserves every data block, scalar item, loop table, unknown category, missing-value marker, and value source line.
- Offers case-insensitive lookup while preserving original data names and values.
- Distinguishes quoted values from bare syntax controls and supports explicit `stop_` loop terminators.
- Rejects content outside a data block, unnamed blocks, missing scalar values, duplicate data names, duplicate loop tags, and ragged loop rows with structured parse errors.
- Applies the existing mmCIF tokenizer resource limits.

## Implementation Notes

- The document layer owns strings and has no live dependency on the input buffer.
- It does not infer entities, molecule boundaries, bonds, coordinate models, or chemistry.

## Validation

- Unit tests cover unknown scalar/loop content, quoted controls, multiline values, multiple blocks, lookup, and malformed structure.
- No external full-document golden evidence currently exists, so the feature remains unvalidated.

## Out Of Scope

- Dictionary save frames, DDL validation, category schemas, binary CIF, writing, and chemical interpretation.

## Revision Notes

- v1: Add the format-level multi-block mmCIF document representation and parser.
