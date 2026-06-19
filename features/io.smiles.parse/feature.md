# Basic SMILES Parser

## Summary

Parse a practical initial subset of SMILES into `SmallMolecule`.

## Behavior/API

- Exposes `SmilesParseOptions`, `SmilesParseError`, and `read_smiles_str`.
- Supports organic subset atoms, bracket atoms, branches, ring closures, dot fragments, common bond symbols, charges, isotopes, explicit hydrogens, aromatic lowercase atoms, and atom maps.
- Does not run sanitization or perception.

## Implementation Notes

- The parser is intentionally non-query and non-reaction.
- Aromatic lowercase atoms set initial aromatic flags but do not replace explicit aromaticity perception.

## Validation

- Unit tests cover branches, ring closures, bracket atoms, and dot fragments.
- RDKit reference generator support is included for fixture/golden generation.

## Out Of Scope

- SMARTS, reactions, full stereochemistry, wildcard/query atoms, and full grammar parity.

## Revision Notes

- v1: Initial basic SMILES parser.
