# mmCIF Parser

## Summary

Parse mmCIF atom-site data into `MacroMolecule` while preserving label and author identifiers in `BioHierarchy`.

## Behavior/API

- Exposes `MmcifParseOptions`, `MmcifParseError`, and `read_mmcif_str`.
- Parses `_atom_site` loops into core atoms plus SMCRA hierarchy records.
- Preserves model, chain, residue, atom-site labels, alternate locations, occupancy, B-factor, and insertion code where present.
- Preserves `_atom_site.group_PDB`, `_atom_site.id`, label and author identifiers, and complete Cartesian coordinate triplets.
- Stores parsed coordinates in the first conformer; each atom-site row maps to the corresponding atom position in that conformer.
- In strict mode, rejects rows with partial coordinate triplets and residues lacking both label and author sequence identifiers.
- Does not infer bonds or run chemical perception.

## Implementation Notes

- Tokenizes mmCIF with support for quotes, comments, missing values, and semicolon text blocks.
- Keeps biomolecular labels in `BioHierarchy`, not core `Atom`.
- Unknown categories are ignored by this parser version.

## Validation

- Current coverage combines unit tests with Biopython golden validation through `validation.harness`.
- Fixtures, goldens, and evidence live under each applicable corpus directory.

## Out Of Scope

- Full mmCIF category coverage, PDB parsing, distance-based bond perception, polymer template bonds, rings, aromaticity, valence, and stereochemistry.
- Runtime Biopython dependency.

## Revision Notes

- v1: Raw `_atom_site` parser into `MacroMolecule`.
- v2: Preserve atom-site row metadata, author-keyed residues without label sequence IDs, and Cartesian coordinates.
