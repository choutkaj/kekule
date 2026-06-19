# Molfile V2000 Parser

## Summary

Parse a single Molfile V2000 block into `SmallMolecule` using raw parsing semantics.

## Behavior/API

- Exposes `read_mol_v2000_str`.
- Preserves title/program/comment properties, atom coordinates, bond orders, atom map numbers, formal charges, isotopes, and radical electron counts where represented by common V2000 fields.
- Rejects V3000 and malformed graph endpoints.
- Does not run sanitization, valence perception, ring perception, aromaticity, or stereochemistry perception.

## Implementation Notes

- SDF V2000 parsing delegates Molfile-block parsing to this feature.
- V2000 one-based atom indices are mapped to stable `AtomId`s.
- Coordinates are stored in the first conformer.

## Validation

- Unit tests cover coordinates, `M  CHG`, `M  ISO`, radicals, atom maps, and malformed blocks.
- RDKit-generated goldens compare raw Molfile-preserved atom, bond, metadata, and coordinate records for external PubChem fixtures.

## Out Of Scope

- V3000, query atom/bond semantics, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Raw Molfile V2000 parser.
