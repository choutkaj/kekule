# Molfile V2000 Parser

## Summary

Parse a Molfile into `MolfileDocument`, then explicitly interpret V2000 chemistry
as one `SmallMolecule` plus source mappings.

## Behavior/API

- Exposes version-autodetecting `molfile::parse_str` and `molfile::interpret`
  with distinct `MolfileParseError` and `MolfileInterpretError`.
- `molfile::parse_str_with_options` accepts `MolfileParseOptions` for input-byte
  and V3000 allocation bounds; `parse_str` uses the documented defaults.
- The document preserves headers, atom/bond/property records, unsupported
  records, source lines, original text, and a validated private typed CTAB
  representation without constructing a molecule.
- Nonempty content after `M  END` remains visible as unsupported document
  records and is included in the interpretation report instead of being
  silently omitted from the structured view.
- Interpretation returns `MolfileInterpretation`; its report maps source atom
  and bond records to stable canonical `AtomId` and `BondId` values.
- Interpretation preserves atom coordinates, bond orders, atom maps, charges,
  isotopes, radicals, and supported source stereo marks. Headers remain document
  metadata rather than molecule properties.
- Atom-block charge code 4 is preserved as a doublet radical; `M  RAD` records
  continue to preserve explicit singlet, doublet, and triplet multiplicities.
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
- v13: Parse and validate atom, bond, and supported property records into a
  typed V2000 syntax representation; interpretation consumes it without
  reparsing source lines and returns canonical mappings.
- v14: Preserve atom-block charge code 4 as a doublet radical and reject
  malformed or unsupported atom-block charge codes instead of silently
  coercing them to neutral atoms.
- v15: Add configurable Molfile input bounds and shared V3000 atom, bond, and
  logical-line limits while retaining the existing default-bounded parser.
- v16: Preserve and report nonempty records after `M  END` as unsupported
  document content instead of retaining them only in the raw source string.
