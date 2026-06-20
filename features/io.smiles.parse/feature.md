# SMILES Parser

## Summary

Parse a practical subset of SMILES into `SmallMolecule`.

## Behavior/API

- Exposes `SmilesParseOptions`, `SmilesParseError`, and `read_smiles_str`.
- Supports organic subset atoms, bracket atoms, branches, ring closures, dot fragments, common bond symbols, charges, isotopes, explicit hydrogens, aromatic lowercase atoms, and atom maps.
- Supports one-digit ring labels and `%10` through `%99`.
- Rejects empty or malformed brackets, overflowed isotope/charge/map values, incomplete bonds/branches/rings, cross-component closures, conflicting ring bond symbols, and unsupported stereo/query syntax.
- Does not run sanitization or perception.

## Implementation Notes

- The parser is intentionally non-query and non-reaction.
- Cursor offsets remain UTF-8 character boundaries; bracket grammar is consumed strictly as ASCII rather than skipping unknown bytes.
- Aromatic lowercase atoms set aromatic flags but do not replace explicit aromaticity perception.

## Validation

- Unit tests cover branches, ring closures, bracket atoms, dot fragments, malformed percent labels, non-ASCII brackets, pending bonds, and unsupported syntax.
- The standalone `smiles` fuzz target checks panic safety and successful parse/write/parse paths.
- RDKit-generated goldens compare parsed graph records for external PubChem SMILES fixtures.

## Out Of Scope

- SMARTS, reactions, full stereochemistry, wildcard/query atoms, and full grammar parity.

## Revision Notes

- v1: SMILES parser.
- v2: Strictly consume bracket and structural syntax, add `%10`-`%99` ring labels, and add fuzz coverage.
