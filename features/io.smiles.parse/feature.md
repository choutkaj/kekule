# SMILES Parser

## Summary

Parse a practical subset of SMILES into `SmallMolecule`.

## Behavior/API

- Exposes `smiles::{SmilesParseOptions, SmilesParseError, read_str, read_str_with_options, read_sanitized_str}`.
- Supports organic subset atoms, bracket atoms, branches, ring closures, dot
  fragments, common bond symbols, directional single-bond markers, charges,
  isotopes, explicit hydrogens, bracket tetrahedral markers, aromatic lowercase
  atoms including bracketed Se/Te chalcogens, atom maps, and bracket-atom
  singlet-through-quintet radical states implied by charge-adjusted default
  valence.
- Supports one-digit ring labels and `%10` through `%99`.
- Resolves omitted bonds between aromatic atoms, including ring closures, as aromatic bonds in the supported subset.
- Rejects empty or malformed brackets, overflowed isotope/charge/map values, incomplete bonds/branches/rings, cross-component closures, conflicting ring bond symbols, and unsupported query syntax.
- Does not run sanitization or perception.

## Implementation Notes

- The parser is intentionally non-query and non-reaction.
- Bracket `@` and `@@` tetrahedral markers are stored as `Molecule` stereo elements with SMILES local orientation, carrier order, source, and specifiedness; supported three-carrier no-H N/P/As/Sb/O/S/Se/Te centers receive an implicit-lone-pair carrier. Parsed markers are not assigned CIP descriptors, normalized, or validated against stereogenicity.
- Directional `/` and `\` bond markers are stored as source bond marks on `Molecule`; they are not assembled, normalized, or validated across double-bond stereo systems.
- Cursor offsets remain UTF-8 character boundaries; bracket grammar is consumed strictly as ASCII rather than skipping unknown bytes.
- Aromatic lowercase atoms set aromatic flags and aromatic omitted bonds, but they do not replace explicit sanitization/perception. Bracketed aromatic element support follows the current RDKit-like aromatic donor set used by perception.
- Radical inference is graph-wide after parsing, so a bracket atom's explicit
  bonds, hydrogens, formal charge, and aromatic bond valence determine its
  unpaired-electron count without element-specific SMILES motifs.

## Validation

- Unit tests cover branches, ring closures, bracket atoms, dot fragments, malformed percent labels, non-ASCII brackets, pending bonds, supported local stereo preservation, and unsupported syntax.
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
- v7: Move the public parser API under the `smiles` facade and keep sanitizing reads explicit by name.
- v8: Add PubChem-100k as required broad-corpus validation evidence.
- v9: Remove parser-side metal-bound organic no-implicit preservation so sanitized SMILES semantics keep RDKit-like no-implicit flags while valence perception still assigns zero implicit hydrogens.
- v10: Preserve supported local SMILES stereo syntax in the first-class stereo representation instead of atom/bond payload fields.
- v11: Preserve supported three-carrier heteroatom tetrahedral markers with an implicit-lone-pair carrier rather than dropping the local stereo element.
- v12: Infer bracket-atom doublet, triplet, quartet, and quintet multiplicity
  from charge-adjusted default valence after graph construction.
