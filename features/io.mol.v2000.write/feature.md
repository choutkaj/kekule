# Molfile V2000 Writer

## Summary

Write `SmallMolecule` values to Molfile V2000 text for round-trip oriented workflows.

## Behavior/API

- Exposes `write_mol_v2000`.
- Emits atom coordinates from the first conformer when present.
- Emits common bond orders plus `M  CHG`, `M  ISO`, and `M  RAD` records.
- Does not sanitize or canonicalize before writing.

## Implementation Notes

- Writer preserves current graph iteration order.
- Unsupported bond-order details are rejected rather than silently downgraded.

## Validation

- Unit tests cover Molfile parse/write/parse round trips.
- RDKit-generated goldens compare Molfile-preservable atoms, bonds, coordinates, charges, isotopes, atom maps, and headers for external PubChem fixtures.

## Out Of Scope

- V3000 writing, canonical atom ordering, query features, stereo parity, and runtime RDKit.

## Revision Notes

- v1: V2000 writer.
- v2: Validation contract excludes SDF data fields and passes the RDKit-backed `tiny` corpus.
