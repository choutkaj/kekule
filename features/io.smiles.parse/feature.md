# SMILES Parser

## Summary

Parse a practical subset of SMILES into `SmallMolecule`.

## Behavior/API

- Exposes `SmilesParseOptions`, `SmilesParseError`, and `read_smiles_str`.
- Supports organic subset atoms, bracket atoms, branches, ring closures, dot fragments, common bond symbols, directional single-bond markers, charges, isotopes, explicit hydrogens, bracket tetrahedral markers, aromatic lowercase atoms including bracketed Se/Te chalcogens, and atom maps.
- Supports one-digit ring labels and `%10` through `%99`.
- Resolves omitted bonds between aromatic atoms, including ring closures, as aromatic bonds in the supported subset.
- Rejects empty or malformed brackets, overflowed isotope/charge/map values, incomplete bonds/branches/rings, cross-component closures, conflicting ring bond symbols, and unsupported query syntax.
- Does not run sanitization or perception.

## Implementation Notes

- The parser is intentionally non-query and non-reaction.
- Bracket `@` and `@@` tetrahedral markers are stored as atom stereo metadata but are not assigned, normalized, or validated against neighbor order.
- Directional `/` and `\` bond markers are stored as bond stereo metadata but are not assigned, normalized, or validated across double-bond stereo systems.
- Metal-bound organic-subset atoms are marked as no-implicit-hydrogen atoms during parsing when their parsed bond valence is already filled, preserving RDKit-like valence semantics after canonical reparse without suppressing hydrogens on lower-valence organometallic shorthand.
- Cursor offsets remain UTF-8 character boundaries; bracket grammar is consumed strictly as ASCII rather than skipping unknown bytes.
- Aromatic lowercase atoms set aromatic flags and aromatic omitted bonds, but they do not replace explicit sanitization/perception. Bracketed aromatic element support follows the current RDKit-like aromatic donor set used by perception.

## Validation

- Unit tests cover branches, ring closures, bracket atoms, dot fragments, malformed percent labels, non-ASCII brackets, pending bonds, and unsupported syntax.
- The standalone `smiles` fuzz target checks panic safety and successful parse/write/parse paths.
- RDKit-generated goldens compare raw parsed graph records plus parse-sanitize-write-reparse atom identity, labeled-neighbor topology, bond order/aromaticity, charge, isotope, hydrogen, map, and valence records for external PubChem SMILES fixtures.

## Out Of Scope

- SMARTS, reactions, full stereochemistry assignment, wildcard/query atoms, and full grammar parity.

## Revision Notes

- v1: SMILES parser.
- v2: Strictly consume bracket and structural syntax, add `%10`-`%99` ring labels, and add fuzz coverage.
- v3: Resolve omitted aromatic bonds, preserve writer-compatible ring labels through validation, and compare RDKit-backed sanitize/write/reparse semantics for supported SMILES while recording unsupported stereo/query inputs explicitly.
- v4: Parse bracket tetrahedral `@`/`@@` markers as atom stereo metadata without assigning stereochemistry.
- v5: Parse directional `/` and `\` bond markers as bond stereo metadata, and preserve no-implicit hydrogen semantics for metal-bound organic halogens.
- v6: Parse bracketed aromatic Se/Te atoms and generalize metal-bound no-implicit preservation to valence-filled organic-subset atoms so canonical organometallic output can be reparsed and sanitized.
