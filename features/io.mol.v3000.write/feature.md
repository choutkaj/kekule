# Molfile V3000 Writer

## Summary

Write deterministic Molfile V3000 CTAB output for the supported raw graph subset.

## Behavior/API

- Exposes `molfile::write_v3000`.
- Emits three-line Molfile headers, V3000 `CTAB`, `COUNTS`, `ATOM`, and `BOND` sections in graph order.
- Emits neutral generated headers plus coordinates from the first conformer,
  bond orders, maps, charges, isotopes, radicals, and supported source `CFG`
  marks. Format headers are not read from molecule properties.
- Successful output is accepted by `molfile::parse_str` then
  `molfile::interpret`.
- Rejects stored stereo elements, perceived E/Z bond stereo, bond `CFG` source marks incompatible with the bond order, enhanced stereo, and quadruple bonds with structured writer errors.
- Does not canonicalize, sanitize, or perceive chemistry before writing.

## Implementation Notes

- Atom and bond output follows current graph insertion order.
- Compatible conformer length quantities are converted to angstroms; missing
  coordinates write as zero coordinates.
- Supported bond orders match the parser subset: zero, single, double, triple, aromatic, and dative.
- Supported V3000 bond `CFG` output is read from source bond marks, not from atom or bond payload fields.
- Unsupported chemistry is rejected instead of silently being coerced into a different representation.
- `RAD=1`, `RAD=2`, and `RAD=3` preserve singlet, doublet, and triplet;
  quartet/quintet radicals return a structured error because V3000 defines no
  lossless code for them.

## Validation

- Unit tests cover writer parse-back preservation for metadata, coordinates, atom maps, charges, isotopes, radicals, and supported source bond stereo marks.
- A dedicated bounded fuzz target exercises V3000 parse, interpretation, write,
  and reparse in CI smoke tests and scheduled campaigns.
- RDKit-generated goldens compare Molfile-preservable content for the same external PubChem fixtures used by the V2000 writer tier.
- PubChem-1k is required baseline evidence; manifest-backed broader corpora
  remain available for deliberate local parity checks.

## Out Of Scope

SDF V3000 writing, canonical atom ordering, query atom/bond semantics, atom stereochemistry, enhanced stereochemistry collections, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Feature contract reserved.
- v2: Molfile V3000 writer for the parser-supported raw graph subset.
- v3: Declare the same required small-molecule validation corpora as V2000 Molfile writing.
- v4: Move the public writer API under the `molfile` facade.
- v5: Add PubChem-100k as required broad-corpus validation evidence.
- v6: Read supported V3000 bond `CFG` output from source bond marks and reject stored stereo elements.
- v7: Reject quartet/quintet radical multiplicity explicitly instead of
  silently mapping an unrepresentable high-spin state.
- v8: Migrate parse-back validation to `MolfileDocument` and remove
  molecule-property header coupling.
- v9: Make the committed smoke corpus the CI-reproducible required evidence
  tier while retaining every ignored corpus on demand.
- v10: Convert explicit conformer length units to the CTfile angstrom convention.
- v11: Add a dedicated bounded V3000 parse/interpret/write round-trip fuzz
  target to CI and scheduled campaigns.
- v12: Use PubChem-1k as the required baseline validation corpus after retiring the former smoke corpus from public validation.
