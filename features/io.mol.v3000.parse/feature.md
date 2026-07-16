# Molfile V3000 Parser

## Summary

Interpret V3000 records from the version-autodetected loss-preserving
`MolfileDocument` as one `SmallMolecule`.

## Behavior/API

- Shares `molfile::parse_str` and `molfile::interpret` with V2000; no
  version-specific reader remains.
- Parses three-line Molfile headers, V3000 `CTAB`, `COUNTS`, `ATOM`, and `BOND` sections, and line continuations.
- Preserves title/program/comment properties, atom coordinates, bond orders, atom map numbers, formal charges, isotopes via `MASS`, radical multiplicities, and supported V3000 bond `CFG` stereo as source bond marks.
- Preserves valence-implied hydrogen carriers on degree-three tetrahedral `CFG`
  centers as explicit atom hydrogens using the same rule as V2000 parsing.
- Rejects malformed sections, count mismatches, duplicate atom indices, out-of-range bond endpoints, unknown elements, non-finite coordinates, unsupported bond orders, atom stereochemistry, and enhanced stereo collections with structured parse errors.
- Does not run sanitization, valence perception, ring perception, aromaticity, or stereochemistry perception.

## Implementation Notes

- V3000 atom indices are mapped to stable `AtomId`s.
- Supported bond orders use the existing core `BondOrder` representation: zero, single, double, triple, aromatic, and dative.
- Supported bond `CFG` mappings are stored as source bond marks for wedge/either cases.
- Coordinates are interpreted as angstrom quantities and stored in the first
  conformer with that explicit unit.

## Validation

- Unit tests cover successful raw parsing, line continuations, metadata fields, no-perception behavior, malformed counts, count mismatches, non-finite coordinates, bad endpoints, supported source bond stereo marks, unsupported atom stereo, and unsupported bond types.
- RDKit-generated goldens compare Molfile-preserved atom, bond, metadata, and coordinate records for the same external PubChem fixtures used by the V2000 parser tier.
- All ignored non-smoke corpora remain available for explicit local-only
  validation but do not determine repository-wide validation state.

## Out Of Scope

SDF V3000 parsing, V3000 writing, query atom/bond semantics, atom stereochemistry, enhanced stereochemistry collections, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Feature contract reserved.
- v2: Raw Molfile V3000 parser for CTAB atoms, bonds, coordinates, common atom metadata, and supported bond stereo.
- v3: Declare the same required small-molecule validation corpora as V2000 Molfile parsing.
- v4: Move the public parser API under the `molfile` facade.
- v5: Add PubChem-100k as required broad-corpus validation evidence.
- v6: Store supported V3000 bond `CFG` values as first-class source bond marks.
- v7: Preserve valence-implied tetrahedral hydrogen carriers from V3000 `CFG`
  syntax using the shared RDKit-like allowed-valence table.
- v8: Migrate V3000 input to `MolfileDocument` parse/interpret with distinct
  syntax and chemical errors.
- v9: Make the committed smoke corpus the CI-reproducible required evidence
  tier while retaining every ignored corpus on demand.
- v10: Record the CTfile coordinate convention as an explicit angstrom
  conformer unit.
