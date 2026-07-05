# Molfile V3000 Writer

## Summary

Write deterministic Molfile V3000 CTAB output for the supported raw graph subset.

## Behavior/API

- Exposes `molfile::write_v3000`.
- Emits three-line Molfile headers, V3000 `CTAB`, `COUNTS`, `ATOM`, and `BOND` sections in graph order.
- Preserves title/program/comment properties, coordinates from the first conformer, bond orders, atom map numbers, formal charges, isotopes via `MASS`, radical multiplicities, and supported source bond `CFG` stereo marks.
- Successful output is accepted by `molfile::read_v3000_str`.
- Rejects stored stereo elements, perceived E/Z bond stereo, bond `CFG` source marks incompatible with the bond order, enhanced stereo, and quadruple bonds with structured writer errors.
- Does not canonicalize, sanitize, or perceive chemistry before writing.

## Implementation Notes

- Atom and bond output follows current graph insertion order.
- Missing coordinates write as zero coordinates.
- Supported bond orders match the parser subset: zero, single, double, triple, aromatic, and dative.
- Supported V3000 bond `CFG` output is read from source bond marks, not from atom or bond payload fields.
- Unsupported chemistry is rejected instead of silently being coerced into a different representation.

## Validation

- Unit tests cover writer parse-back preservation for metadata, coordinates, atom maps, charges, isotopes, radicals, and supported source bond stereo marks.
- RDKit-generated goldens compare Molfile-preservable content for the same external PubChem fixtures used by the V2000 writer tier.

## Out Of Scope

SDF V3000 writing, canonical atom ordering, query atom/bond semantics, atom stereochemistry, enhanced stereochemistry collections, full MDL property coverage, sanitization, and runtime RDKit.

## Revision Notes

- v1: Feature contract reserved.
- v2: Molfile V3000 writer for the parser-supported raw graph subset.
- v3: Declare the same required small-molecule validation corpora as V2000 Molfile writing.
- v4: Move the public writer API under the `molfile` facade.
- v5: Add PubChem-100k as required broad-corpus validation evidence.
- v6: Read supported V3000 bond `CFG` output from source bond marks and reject stored stereo elements.
