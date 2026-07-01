# Molfile V3000 Parser

## Summary

Parse a single Molfile V3000 CTAB into `SmallMolecule` using raw parsing semantics.

## Behavior/API

- Exposes `molfile::read_v3000_str`.
- Parses three-line Molfile headers, V3000 `CTAB`, `COUNTS`, `ATOM`, and `BOND` sections, and line continuations.
- Preserves title/program/comment properties, atom coordinates, bond orders, atom map numbers, formal charges, isotopes via `MASS`, radical multiplicities, and supported V3000 bond `CFG` stereo.
- Rejects malformed sections, count mismatches, duplicate atom indices, out-of-range bond endpoints, unknown elements, non-finite coordinates, unsupported bond orders, and atom stereochemistry with structured parse errors.
- Does not run sanitization, valence perception, ring perception, aromaticity, or stereochemistry perception.

## Implementation Notes

- V3000 atom indices are mapped to stable `AtomId`s.
- Supported bond orders use the existing core `BondOrder` representation: zero, single, double, triple, aromatic, and dative.
- Supported bond `CFG` mappings are stored in the existing `BondStereo` representation for wedge/either cases.
- Coordinates are stored in the first conformer.

## Validation

- Unit tests cover successful raw parsing, line continuations, metadata fields, no-perception behavior, malformed counts, count mismatches, non-finite coordinates, bad endpoints, unsupported atom stereo, and unsupported bond types.
- RDKit-generated goldens compare Molfile-preserved atom, bond, metadata, and coordinate records for the same external PubChem fixtures used by the V2000 parser tier.

## Out Of Scope

SDF V3000 parsing, V3000 writing, query atom/bond semantics, atom stereochemistry, enhanced stereochemistry collections, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Feature contract reserved.
- v2: Raw Molfile V3000 parser for CTAB atoms, bonds, coordinates, common atom metadata, and supported bond stereo.
- v3: Declare the same required small-molecule validation corpora as V2000 Molfile parsing.
- v4: Move the public parser API under the `molfile` facade.
