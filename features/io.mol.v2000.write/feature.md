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
- RDKit reference generator support is included for fixture/golden generation.

## Out Of Scope

- V3000 writing, canonical atom ordering, query features, stereo parity, and runtime RDKit.

## Revision Notes

- v1: Initial V2000 writer.
