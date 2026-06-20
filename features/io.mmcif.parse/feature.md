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
- Rejects unterminated quoted/semicolon values, ragged loops, integer overflow, and non-finite floating-point values with source-line errors.
- `MmcifParseOptions` bounds input bytes, token count, token bytes, and atom-site row count. Defaults are 256 MiB, 10,000,000 tokens, 16 MiB per token, and 5,000,000 atom-site rows.
- Does not infer bonds or run chemical perception.

## Implementation Notes

- Tokenizes mmCIF with support for quotes, comments, missing values, and semicolon text blocks.
- Enforces resource limits before unbounded token or atom-site allocation.
- Keeps biomolecular labels in `BioHierarchy`, not core `Atom`.
- Unknown categories are ignored by this parser version.

## Validation

- Current coverage combines unit tests with Biopython golden validation through `validation.harness`.
- Fixtures, goldens, and evidence live under each applicable corpus directory.
- Focused malformed-input tests and the standalone `mmcif` fuzz target cover panic safety and successful hierarchy access.

## Out Of Scope

- Full mmCIF category coverage, PDB parsing, distance-based bond perception, polymer template bonds, rings, aromaticity, valence, and stereochemistry.
- Runtime Biopython dependency.

## Revision Notes

- v1: Raw `_atom_site` parser into `MacroMolecule`.
- v2: Preserve atom-site row metadata, author-keyed residues without label sequence IDs, and Cartesian coordinates.
- v3: Add structured malformed-input handling, finite numeric checks, configurable resource limits, and fuzz coverage.
