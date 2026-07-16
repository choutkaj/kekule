# Core Conformer Coordinate Storage

## Summary

Store 2D or 3D atom coordinates as conformers on the shared core `Molecule` graph.

## Behavior/API

- Exposes `Point3`, `Conformer`, `ConformerId`, and conformer accessors on `Molecule`.
- Conformers store optional coordinates keyed by `AtomId` plus one explicit
  compatible length unit for the complete coordinate array.
- Position setters accept `Quantity<Point3>`, convert to the conformer's unit,
  and position accessors return quantities retaining that unit.
- `Molecule::add_conformer` is fallible and transactionally rejects coordinates
  assigned to invalid, deleted, or otherwise non-live atom IDs.
- Adding or deleting topology invalidates coordinate-bearing conformers only when the topology operation removes atoms.

## Implementation Notes

- Coordinates are chemically general and live in core, not the small-molecule wrapper.
- A conformer selects its storage unit at construction; it does not assume a
  hidden coordinate convention.
- Stable conformer IDs use slot storage, matching atom and bond ID behavior.
- Parsers may attach a conformer without running sanitization or perception.
- The mmCIF parser stores coordinates from all models in one conformer because each atom-site row is a distinct graph atom; model identity remains in `SmcraHierarchy`.

## Validation

- Unit tests cover insertion, lookup, and SDF/Molfile coordinate preservation.
- RDKit-generated goldens compare conformer coordinate preservation for external PubChem fixtures.

## Out Of Scope

- Coordinate generation, alignment, RMSD, force-field minimization, and conformer ensembles from external tools.

## Revision Notes

- v1: Shared conformer storage.
- v2: Add PubChem-100k as required broad-corpus validation evidence.
- v3: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v4: Make conformer attachment fallible and reject coordinates for non-live
  graph atoms without inserting a partial conformer.
- v5: Require explicit length units for conformer construction and coordinate
  access through `Quantity<Point3>`.
- v6: Use PubChem-1k as the required baseline validation corpus after retiring the former smoke corpus from public validation.
