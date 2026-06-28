# Molfile V3000 Writer

## Summary

Write deterministic Molfile V3000 CTAB output for the supported raw graph subset.

## Behavior/API

- Exposes `write_mol_v3000`.
- Emits three-line Molfile headers, V3000 `CTAB`, `COUNTS`, `ATOM`, and `BOND` sections in graph order.
- Preserves title/program/comment properties, coordinates from the first conformer, bond orders, atom map numbers, formal charges, isotopes via `MASS`, radical multiplicities, and supported bond `CFG` stereo.
- Successful output is accepted by `read_mol_v3000_str`.
- Rejects atom stereochemistry, perceived E/Z bond stereo, bond `CFG` values incompatible with the bond order, and quadruple bonds with structured writer errors.
- Does not canonicalize, sanitize, or perceive chemistry before writing.

## Implementation Notes

- Atom and bond output follows current graph insertion order.
- Missing coordinates write as zero coordinates.
- Supported bond orders match the parser subset: zero, single, double, triple, aromatic, and dative.
- Unsupported chemistry is rejected instead of silently being coerced into a different representation.

## Validation

- Unit tests cover writer parse-back preservation for metadata, coordinates, atom maps, charges, isotopes, radicals, and supported bond stereo.
- Future golden validation should add RDKit-compatible V3000 fixtures.

## Out Of Scope

SDF V3000 writing, canonical atom ordering, query atom/bond semantics, atom stereochemistry, enhanced stereochemistry collections, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Feature contract reserved.
- v2: Molfile V3000 writer for the parser-supported raw graph subset.
