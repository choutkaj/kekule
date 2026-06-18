# mmCIF Parser

## Summary

Parse mmCIF atom-site data into `MacroMolecule` while preserving label and author identifiers in `BioHierarchy`.

## Behavior/API

- Exposes `MmcifParseOptions`, `MmcifParseError`, and `read_mmcif_str`.
- Parses `_atom_site` loops into core atoms plus SMCRA hierarchy records.
- Preserves model, chain, residue, atom-site labels, alternate locations, occupancy, B-factor, and insertion code where present.
- Does not infer bonds or run chemical perception.

## Implementation Notes

- Tokenizes mmCIF with support for quotes, comments, missing values, and semicolon text blocks.
- Keeps biomolecular labels in `BioHierarchy`, not core `Atom`.
- Unknown categories are ignored by the initial parser.

## Validation

- Current coverage is unit-test based.
- Biopython golden validation is planned through `validation.harness`.
- Fixtures live under `validation/features/io.mmcif.parse/`.

## Out Of Scope

- Full mmCIF category coverage, PDB parsing, distance-based bond perception, polymer template bonds, rings, aromaticity, valence, and stereochemistry.
- Runtime Biopython dependency.

## Revision Notes

- v1: Initial raw `_atom_site` parser into `MacroMolecule`.
