# Interpreted Molecular Contents

## Summary

Represent the chemically meaningful output of a multi-entity structure import without treating the entire file as one molecule.

## Behavior/API

- Exposes `bio::MolecularContents` with separate iterators over `SmallMolecule` and `MacroMolecule` values.
- Exposes `bio::Solvent` as a collection of individual solvent `SmallMolecule` instances.
- Single-atom ions remain ordinary `SmallMolecule` values.
- The containers own their molecules and support consuming decomposition through `into_parts` and `into_molecules`.
- The types do not retain a live dependency on a source file document.

## Implementation Notes

- A molecule remains a chemical object rather than a structural-file record or disconnected mixture.
- `Solvent` is a semantic collection, not a disconnected mega-molecule, so each water can later be selected, removed, transformed, or modeled independently.
- Population is crate-controlled so import interpreters establish category boundaries consistently.

## Validation

- Unit tests cover mixed polymer, ligand, ion, and repeated-water contents plus consuming and read-only access.
- No external corpus evidence currently validates decomposition boundaries, so the feature remains unvalidated.

## Out Of Scope

- Generic solvent species, solvent boxes, periodic cells, bulk-solvent models, molecule roles, preparation, serialization, and direct `MolecularModel` conversion.

## Revision Notes

- v1: Add clean molecular-content and individual-solvent containers for staged structure import.
