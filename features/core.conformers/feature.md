# Core Conformer Coordinate Storage

## Summary

Store 2D or 3D atom coordinates as conformers on the shared core `Molecule` graph.

## Behavior/API

- Exposes `Point3`, `Conformer`, `ConformerId`, and conformer accessors on `Molecule`.
- Conformers store optional coordinates keyed by `AtomId`.
- Adding or deleting topology invalidates coordinate-bearing conformers only when the topology operation removes atoms.

## Implementation Notes

- Coordinates are chemically general and live in core, not the small-molecule wrapper.
- Stable conformer IDs use slot storage, matching atom and bond ID behavior.
- Parsers may attach a conformer without running sanitization or perception.
- The mmCIF parser stores coordinates from all models in one conformer because each atom-site row is a distinct graph atom; model identity remains in `BioHierarchy`.

## Validation

- Unit tests cover insertion, lookup, and SDF/Molfile coordinate preservation.
- RDKit-generated goldens compare conformer coordinate preservation for external PubChem fixtures.

## Out Of Scope

- Coordinate generation, alignment, RMSD, force-field minimization, and conformer ensembles from external tools.

## Revision Notes

- v1: Shared conformer storage.
