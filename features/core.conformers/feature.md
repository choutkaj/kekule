# Core Conformer Coordinate Storage

## Summary

Store 2D or 3D atom coordinates as conformers on the shared core `Molecule` graph.

## Behavior/API

- Exposes `Conformer`, `ConformerId`, and conformer accessors on `Molecule`;
  the shared `Point3` type lives at `geometry::Point3`.
- Conformers store optional coordinates keyed by `AtomId` plus one explicit
  compatible length unit for the complete coordinate array.
- Position setters accept `Quantity<Point3>`, convert to the conformer's unit,
  and position accessors return quantities retaining that unit.
- `Molecule::add_conformer` is fallible and transactionally rejects coordinates
  assigned to invalid, deleted, or otherwise non-live atom IDs.
- `Conformer::new`, `with_atom_capacity`, and `set_position` return
  `ConformerError`; position arrays outside the `AtomId` addressable range fail
  with `PositionCapacityExceeded` before allocation or mutation.
- Adding or deleting topology invalidates coordinate-bearing conformers only when the topology operation removes atoms.

## Implementation Notes

- Coordinates are chemically general and live in core, not the small-molecule wrapper.
- A conformer selects its storage unit at construction; it does not assume a
  hidden coordinate convention.
- Stable conformer IDs use slot storage, matching atom and bond ID behavior.
- Parsers may attach a conformer without running sanitization or perception.
- Molecule-local conformers remain convenience coordinate storage. System-wide
  models, ensembles, and trajectories use the topology-bound structure layer.

## Tests

- Unit tests cover insertion, lookup, synthetic position-capacity boundaries,
  and SDF/Molfile coordinate preservation.

## Benchmarks

- RDKit-generated goldens compare conformer coordinate preservation for external PubChem fixtures.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Coordinate generation, alignment, RMSD, force-field minimization, and conformer ensembles from external tools.

## Revision Notes

- v1: Shared conformer storage.
- v2: Add PubChem-100k as required broad-corpus external-parity evidence.
- v3: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v4: Make conformer attachment fallible and reject coordinates for non-live
  graph atoms without inserting a partial conformer.
- v5: Require explicit length units for conformer construction and coordinate
  access through `Quantity<Point3>`.
- v6: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v7: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v8: Move `Point3` to the shared `geometry` module and distinguish local
  molecule conformers from topology-bound configurations and trajectories.
- v9: Add structured conformer position-capacity errors and remove unchecked
  slot reconstruction from conformer coordinate iteration.
