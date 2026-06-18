# Core Molecular Graph

## Summary

Store atoms, bonds, properties, and computed-state invalidation behind stable typed atom and bond IDs.

## Behavior/API

- Provides one shared `Molecule` graph used by both `SmallMolecule` and `MacroMolecule`.
- Supports adding and deleting atoms and bonds.
- Rejects invalid atom IDs, invalid bond IDs, self-bonds, and duplicate bonds.
- Iterates live atoms, live bonds, neighbors, and incident bonds.
- Preserves stable `AtomId` and `BondId` values after deletion.

## Implementation Notes

- Uses slot storage with tombstones so IDs remain stable.
- Maintains adjacency for neighbor and incident-bond iteration.
- Deleting an atom removes its incident bonds.
- Every topology mutation invalidates computed perception state.
- Molecule, atom, and bond property maps are stored on the core data structures.

## Validation

- Current coverage is unit-test based.
- Tests cover empty molecules, insertion, deletion, invalid IDs, self-bonds, duplicates, iteration, stable IDs, counts, and perception invalidation.
- Reference-tool golden data is not required for this data-structure feature.

## Out Of Scope

- SDF, PDB, or mmCIF parsing.
- Ring detection, aromaticity, valence perception, stereochemistry, canonicalization, and validation generation.
- Runtime RDKit or Biopython dependency.

## Revision Notes

- v1: Initial stable-ID molecular graph and wrapper integration.
