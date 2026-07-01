# mmCIF Parser

## Summary

Parse mmCIF atom-site data into `MacroMolecule` while preserving label and author identifiers in `BioHierarchy`.

## Behavior/API

- Exposes `bio::{MmcifParseOptions, MmcifParseError, read_mmcif_str}`.
- Parses `_atom_site` loops into core atoms plus SMCRA hierarchy records.
- Preserves model, chain, residue, atom-site labels, alternate locations, occupancy, B-factor, and insertion code where present.
- Preserves `_atom_site.group_PDB`, `_atom_site.id`, label and author identifiers, and complete Cartesian coordinate triplets.
- Stores parsed coordinates in one shared conformer; atom-site rows from every model are distinct graph atoms whose model identity remains in `BioHierarchy`.
- In strict mode, rejects rows with partial coordinate triplets and residues lacking both label and author sequence identifiers.
- In lenient mode, groups contiguous sequence-less rows until an atom identity repeats; alternate-location variants remain in one residue, while repeated waters and ligands begin a new occurrence.
- Rejects unterminated quoted/semicolon values, ragged loops, integer overflow, and non-finite floating-point values with source-line errors.
- Rejects duplicate `_atom_site` loop tags instead of silently choosing one column.
- `MmcifParseOptions` bounds input bytes, token count, token bytes, and atom-site row count. Defaults are 256 MiB, 10,000,000 tokens, 16 MiB per token, and 5,000,000 atom-site rows.
- Does not infer bonds or run chemical perception.

## Implementation Notes

- Tokenizes mmCIF with support for quotes, comments, missing values, and semicolon text blocks.
- Enforces resource limits before unbounded token or atom-site allocation.
- Keeps biomolecular labels in `BioHierarchy`, not core `Atom`.
- Preserves label and author chain, component, sequence, and atom identifiers separately.
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
- v4: Preserve label/auth component IDs separately and group ambiguous lenient residues by conservative contiguous occurrences.
- v5: Reject duplicate `_atom_site` loop tags as malformed input.
- v6: Clarify the public parser API under the `bio` facade.
