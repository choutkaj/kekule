# SMILES Document Parsing and Interpretation

## Summary

Parse SMILES syntax into a loss-preserving `SmilesDocument`, then explicitly
interpret one record as one asserted `SmallMolecule`.

## Behavior/API

- Exposes `smiles::parse_str(&str) -> SmilesDocument` and
  `smiles::interpret(&SmilesDocument) -> SmallMolecule` with distinct errors.
- Preserves source text, tokens and UTF-8 spans, branch/ring/stereo marks, and
  dot-component token ranges.
- Interpretation supports the established atom, bond, branch, ring, charge,
  isotope, map, radical, aromatic-token, and local-stereo subset.
- Dots create disconnected graph components but do not split the asserted
  molecule; `[Na+].[Cl-]` is one `SmallMolecule`.
- Interpretation may install imported aromatic annotations with input
  provenance. It never sanitizes or runs full perception.
- `SmallMolecule::from_smiles` and `from_smiles_sanitized` remain deliberate
  convenience orchestrators. Empty parse options and direct reader aliases are
  removed.

## Validation

- Unit/fuzz tests cover raw document spans and components, malformed syntax,
  parse/interpret separation, disconnected molecules, local stereo, and
  parse-write-parse behavior. Existing external validation corpora remain the
  semantic reference.

## Out Of Scope

- SMARTS, reactions, sanitization during parsing, and full grammar parity.

## Revision Notes

- v1-v12: Practical raw SMILES reader and expanded chemistry/stereo coverage.
- v13: Hard break to `SmilesDocument` parse/interpret and remove direct readers
  and `SmilesParseOptions`.
- v14: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
