# SMCRA-like Biomolecular Hierarchy

## Summary

Represent model, chain, residue, and atom-site hierarchy as a sidecar over the shared core molecule graph.

## Behavior/API

- Exposes typed hierarchy IDs for models, chains, residues, and atom sites.
- Stores biomolecular labels and atom-site metadata in `BioHierarchy`, not core `Atom`.
- `MacroMolecule` owns one core `Molecule` plus one `BioHierarchy`.
- Atom-site insertion validates that referenced core atoms exist.

## Implementation Notes

- Preserves insertion order for hierarchy iteration.
- Tracks label and author identifiers separately.
- Supports alternate locations, occupancy, B-factor, insertion code, and model identifiers.

## Validation

- Current coverage is unit-test based.
- Biopython golden validation is planned through `validation.harness`.
- Fixtures live under `validation/features/bio.hierarchy.smcra/`.

## Out Of Scope

- PDB parsing, full mmCIF category coverage, polymer connectivity, sequence extraction, and chemical perception.
- Runtime Biopython dependency.

## Revision Notes

- v1: Initial SMCRA sidecar hierarchy for macromolecular parsing.
