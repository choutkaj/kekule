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

## Validation

- Unit tests cover insertion, lookup, and SDF/Molfile coordinate preservation.
- No external reference golden is required for the initial storage feature.

## Out Of Scope

- Coordinate generation, alignment, RMSD, force-field minimization, and conformer ensembles from external tools.

## Revision Notes

- v1: Initial shared conformer storage.
