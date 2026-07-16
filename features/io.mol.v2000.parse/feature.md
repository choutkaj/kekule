# Molfile V2000 Parser

## Summary

Parse a Molfile into `MolfileDocument`, then explicitly interpret V2000 chemistry
as one `SmallMolecule`.

## Behavior/API

- Exposes version-autodetecting `molfile::parse_str` and `molfile::interpret`
  with distinct `MolfileParseError` and `MolfileInterpretError`.
- The document preserves headers, atom/bond/property records, unsupported
  records, source lines, and original text without constructing a molecule.
- Interpretation preserves atom coordinates, bond orders, atom maps, charges,
  isotopes, radicals, and supported source stereo marks. Headers remain document
  metadata rather than molecule properties.
- Preserves the valence-implied hydrogen carrier on degree-three tetrahedral
  wedge centers as an explicit atom hydrogen, matching RDKit's CTfile parsing
  semantics without running general valence perception.
- Preserves a populated V2000 atom-block valence field as
  `no_implicit_hydrogens`, including code 15 for zero valence.
- Rejects V3000, zero or out-of-range graph endpoints, non-ASCII structural fields, truncated records, malformed M records, non-finite coordinates, and counts above the V2000 limit.
- Does not run sanitization, valence perception, ring perception, aromaticity, or stereochemistry perception.

## Implementation Notes

- SDF V2000 parsing delegates Molfile-block parsing to this feature.
- V2000 one-based atom indices are mapped to stable `AtomId`s.
- The code tables are pinned to BIOVIA CTfile Formats 2020 V2000 CTAB bond-block and properties-block definitions.
- `M  RAD` codes map exactly to radical multiplicities: code 1 is singlet, 2 is doublet, and 3 is triplet.
- Supported bond stereo codes are code 1/4/6 on single bonds for wedge-up/either/wedge-down and code 3 on double bonds for double-bond either; they are stored as source bond marks, and unsupported code/order combinations are rejected.
- Counts, atom fields, and bond lines use ASCII byte-field helpers with checked block arithmetic and whitespace fallback for short permissive inputs.
- Allocation is bounded by the V2000 limit of 999 atoms and 999 bonds.
- Coordinates are interpreted as angstrom quantities and stored in the first
  conformer with that explicit unit.

## Validation

- Unit tests cover coordinates, `M  CHG`, `M  ISO`, radical multiplicities, supported and unsupported source bond stereo marks, atom maps, zero endpoints, non-ASCII/truncated fields, extreme counts, and malformed blocks.
- The standalone `mol_v2000` fuzz target checks panic safety and successful parse/write/parse paths.
- RDKit-generated goldens compare raw Molfile-preserved atom, bond, metadata, and coordinate records for external PubChem fixtures.
- All ignored non-smoke corpora remain available for explicit local-only
  parity checks but are not selected as required routine evidence.

## Out Of Scope

- V3000, query atom/bond semantics, atom parity, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Raw Molfile V2000 parser.
- v2: Preserve RDKit-compatible V2000 double-bond `STEREOANY` markers and parse fixed-width three-digit counts/endpoints.
- v3: Make structural parsing byte-safe and checked, reject malformed counts/endpoints/M records, and add fuzz coverage.
- v4: Preserve exact V2000 radical multiplicity and supported single-/double-bond stereo mappings.
- v5: Move the public parser API under the `molfile` facade.
- v6: Add PubChem-100k as required broad-corpus validation evidence.
- v7: Store supported V2000 bond stereo codes as first-class source bond marks.
- v8: Preserve populated atom-block valence/no-implicit semantics and
  valence-implied tetrahedral hydrogen carriers from V2000 wedge syntax using
  the shared RDKit-like allowed-valence table.
- v9: Hard break to loss-preserving `MolfileDocument` parsing, explicit
  interpretation, and separate syntax/chemistry errors.
- v10: Parse `MolfileDocument` V2000 counts from fixed-width fields so adjacent
  three-digit atom and bond counts remain unambiguous.
- v11: Make the committed smoke corpus the CI-reproducible required evidence
  tier while retaining every ignored corpus on demand.
- v12: Record the CTfile coordinate convention as an explicit angstrom
  conformer unit.
