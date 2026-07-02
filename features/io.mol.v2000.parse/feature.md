# Molfile V2000 Parser

## Summary

Parse a single Molfile V2000 block into `SmallMolecule` using raw parsing semantics.

## Behavior/API

- Exposes `molfile::read_v2000_str`.
- Preserves title/program/comment properties, atom coordinates, bond orders, atom map numbers, formal charges, isotopes, radical multiplicities, and supported V2000 bond stereo codes.
- Rejects V3000, zero or out-of-range graph endpoints, non-ASCII structural fields, truncated records, malformed M records, non-finite coordinates, and counts above the V2000 limit.
- Does not run sanitization, valence perception, ring perception, aromaticity, or stereochemistry perception.

## Implementation Notes

- SDF V2000 parsing delegates Molfile-block parsing to this feature.
- V2000 one-based atom indices are mapped to stable `AtomId`s.
- The code tables are pinned to BIOVIA CTfile Formats 2020 V2000 CTAB bond-block and properties-block definitions.
- `M  RAD` codes map exactly to radical multiplicities: code 1 is singlet, 2 is doublet, and 3 is triplet.
- Supported bond stereo codes are code 1/4/6 on single bonds for wedge-up/either/wedge-down and code 3 on double bonds for double-bond either; unsupported code/order combinations are rejected.
- Counts, atom fields, and bond lines use ASCII byte-field helpers with checked block arithmetic and whitespace fallback for short permissive inputs.
- Allocation is bounded by the V2000 limit of 999 atoms and 999 bonds.
- Coordinates are stored in the first conformer.

## Validation

- Unit tests cover coordinates, `M  CHG`, `M  ISO`, radical multiplicities, supported and unsupported bond stereo, atom maps, zero endpoints, non-ASCII/truncated fields, extreme counts, and malformed blocks.
- The standalone `mol_v2000` fuzz target checks panic safety and successful parse/write/parse paths.
- RDKit-generated goldens compare raw Molfile-preserved atom, bond, metadata, and coordinate records for external PubChem fixtures.

## Out Of Scope

- V3000, query atom/bond semantics, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Raw Molfile V2000 parser.
- v2: Preserve RDKit-compatible V2000 double-bond `STEREOANY` markers and parse fixed-width three-digit counts/endpoints.
- v3: Make structural parsing byte-safe and checked, reject malformed counts/endpoints/M records, and add fuzz coverage.
- v4: Preserve exact V2000 radical multiplicity and supported single-/double-bond stereo mappings.
- v5: Move the public parser API under the `molfile` facade.
- v6: Add PubChem-100k as required broad-corpus validation evidence.
