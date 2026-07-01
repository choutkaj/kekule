# Core Molecular Graph

## Summary

Store atoms, bonds, properties, and computed-state invalidation behind stable typed atom and bond IDs.

## Behavior/API

- Provides one shared `Molecule` graph used by both `SmallMolecule` and `MacroMolecule`.
- Supports adding and deleting atoms and bonds.
- Rejects invalid atom IDs, invalid bond IDs, self-bonds, and duplicate bonds.
- Iterates live atoms, live bonds, neighbors, and incident bonds.
- Preserves stable `AtomId` and `BondId` values after deletion.
- Returns scoped `AtomMut` and `BondMut` guards from mutable graph access.
- Tracks computed perception state internally without exposing cache freshness as public API.

## Implementation Notes

- Uses slot storage with tombstones so IDs remain stable.
- Maintains adjacency for neighbor and incident-bond iteration.
- Deleting an atom removes its incident bonds.
- Topology and chemistry-relevant atom or bond changes invalidate computed perception state and clear cached ring objects.
- Property-only and coordinate-only edits do not invalidate chemistry state.
- Mutation guards compare chemistry-relevant fields when released, so obtaining mutable access alone does not stale perception.
- Molecule, atom, and bond property maps are stored on the core data structures.

## Validation

- Current coverage is unit-test based.
- Tests cover empty molecules, insertion, deletion, invalid IDs, self-bonds, duplicates, iteration, stable IDs, counts, chemistry invalidation, and state-neutral property/coordinate edits.
- Reference-tool golden data is not required for this data-structure feature.

## Out Of Scope

- SDF, PDB, or mmCIF parsing.
- Ring detection, aromaticity, valence perception, stereochemistry, canonicalization, and validation generation.
- Runtime RDKit or Biopython dependency.

## Revision Notes

- v1: Stable-ID molecular graph and wrapper integration.
- v2: Centralize chemistry invalidation in scoped mutation guards, remove mutable perception-state access, clear stale ring caches, and preserve state across property/coordinate edits.
- v3: Hide perception freshness/cache state from public core API while retaining internal invalidation checks.
