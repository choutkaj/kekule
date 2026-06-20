# SMCRA-like Biomolecular Hierarchy

## Summary

Represent model, chain, residue, and atom-site hierarchy as a sidecar over the shared core molecule graph.

## Behavior/API

- Exposes typed hierarchy IDs for models, chains, residues, and atom sites.
- Stores biomolecular labels and atom-site metadata in `BioHierarchy`, not core `Atom`.
- `MacroMolecule` owns one core `Molecule` plus one `BioHierarchy`.
- Atom-site insertion validates that referenced core atoms exist.
- Atom-site metadata preserves `_atom_site.group_PDB`, `_atom_site.id`, label/auth atom IDs, alternate location, occupancy, and B-factor.

## Implementation Notes

- Preserves insertion order for hierarchy iteration.
- Tracks label and author identifiers separately.
- Supports alternate locations, occupancy, B-factor, insertion code, and model identifiers.
- Stores label and author component IDs separately on residues.
- mmCIF residue grouping uses label sequence identity when available, author sequence identity otherwise, and strict mode rejects ambiguous sequence-less residues.
- Lenient sequence-less grouping is occurrence-based and keeps alternate locations together without merging repeated waters or ligands.

## Validation

- Current coverage combines unit tests with Biopython golden validation through `validation.harness`.
- Fixtures, goldens, and evidence live under each applicable corpus directory.

## Out Of Scope

- PDB parsing, full mmCIF category coverage, polymer connectivity, sequence extraction, and chemical perception.
- Runtime Biopython dependency.

## Revision Notes

- v1: SMCRA sidecar hierarchy for macromolecular parsing.
- v2: Preserve atom-site row metadata and distinguish author-keyed residues when label sequence IDs are absent.
- v3: Preserve label/auth component IDs separately and support conservative lenient occurrence grouping.
