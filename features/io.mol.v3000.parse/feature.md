# Molfile V3000 Parser

## Summary

Interpret V3000 records from the version-autodetected loss-preserving
`MolfileDocument` as one `SmallMolecule` plus source mappings.

## Behavior/API

- Shares `molfile::parse_str` and `molfile::interpret` with V2000; no
  version-specific reader remains.
- Shares `molfile::parse_str_with_options` and `MolfileParseOptions`; default
  limits bound input bytes, V3000 atoms, V3000 bonds, and continued logical-line
  size before allocation.
- Parses three-line Molfile headers, V3000 `CTAB`, `COUNTS`, `ATOM`, and `BOND`
  sections, and line continuations into a validated private typed syntax
  representation.
- Interpretation returns `MolfileInterpretation`; its report maps source atom
  and bond records to stable canonical `AtomId` and `BondId` values.
- Preserves title/program/comment document metadata, atom coordinates, bond
  orders, atom map numbers, formal charges, isotopes via `MASS`, radical
  multiplicities, and supported V3000 bond `CFG` stereo as source bond marks.
- Preserves valence-implied hydrogen carriers on degree-three tetrahedral `CFG`
  centers as explicit atom hydrogens using the same rule as V2000 parsing.
- Rejects malformed sections, incomplete or malformed count records, count
  mismatches, duplicate count or section-control records, misplaced count
  records, structural records outside the CTAB,
  zero or duplicate atom/bond indices, out-of-range bond endpoints, unknown
  elements, non-finite coordinates, malformed, duplicate, or unsupported
  embedded atom/bond options, unsupported bond orders, and atom
  stereochemistry with structured parse errors.
- Loss-preserves unsupported nonstructural CTAB records such as enhanced stereo
  collections in the document and reports their source lines as ignored during
  interpretation; required CTAB control records are not misreported as ignored
  chemistry.
- Does not run sanitization, valence perception, ring perception, aromaticity, or stereochemistry perception.

## Implementation Notes

- V3000 atom indices are mapped to stable `AtomId`s.
- Supported bond orders use the existing core `BondOrder` representation: zero, single, double, triple, aromatic, and dative.
- Supported bond `CFG` mappings are stored as source bond marks for wedge/either cases.
- Coordinates are interpreted as angstrom quantities and stored in the first
  conformer with that explicit unit.

## Validation

- Unit tests cover successful raw parsing, line continuations, metadata fields, no-perception behavior, malformed counts, count mismatches, non-finite coordinates, bad endpoints, supported source bond stereo marks, unsupported atom stereo, and unsupported bond types.
- A dedicated bounded fuzz target exercises V3000 parse, interpretation, write,
  and reparse in CI smoke tests and scheduled campaigns.
- RDKit-generated goldens compare Molfile-preserved atom, bond, metadata, and coordinate records for the same external PubChem fixtures used by the V2000 parser tier.
- PubChem-1k is required baseline evidence; manifest-backed broader corpora
  remain available for deliberate local parity checks.

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
- v11: Make parsing own the typed V3000 atom/bond representation and validation;
  interpretation consumes it without reparsing source lines and returns
  canonical mappings.
- v12: Validate complete count records and positive unique atom/bond indices,
  preserve unsupported nonstructural CTAB records explicitly, and keep required
  V3000 control lines out of the ignored-record report.
- v13: Add configurable V3000 atom, bond, logical-line, and whole-document
  resource bounds through the shared Molfile parser options.
- v14: Add a dedicated bounded V3000 parse/interpret/write round-trip fuzz
  target to CI and scheduled campaigns.
- v15: Reject duplicate structural records and malformed, duplicate, or
  unsupported embedded atom/bond options instead of silently discarding them.
- v16: Require exactly one CTAB, ATOM, and BOND section-control pair and place
  the sole COUNTS record before ATOM, closing remaining structural
  record-discard paths.
- v17: Use PubChem-1k as the required baseline validation corpus after retiring the former smoke corpus from public validation.
