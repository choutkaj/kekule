# Molfile V2000 Writer

## Summary

Write `SmallMolecule` values to Molfile V2000 text for round-trip oriented workflows.

## Behavior/API

- Exposes `molfile::write_v2000`.
- Emits atom coordinates from the first conformer when present.
- Emits common bond orders plus `M  CHG`, `M  ISO`, and exact `M  RAD` records.
- Emits the atom-block valence field when `no_implicit_hydrogens` is set,
  including code 15 for explicit zero valence, and rejects values beyond the
  V2000 field range.
- Emits supported V2000 bond stereo codes from source bond marks without conflating wedge direction and double-bond either stereo.
- Does not sanitize or canonicalize before writing.

## Implementation Notes

- Writer preserves current graph iteration order.
- Unsupported bond-order and source-mark combinations are rejected rather than silently downgraded.
- Stored stereo elements are rejected until atom stereo and enhanced stereo writing are explicitly implemented.
- Radical and stereo code tables are pinned to BIOVIA CTfile Formats 2020 V2000 CTAB bond-block and properties-block definitions.
- Radical multiplicities write as the inverse of the parser mapping: singlet
  code 1, doublet code 2, and triplet code 3. Quartet and quintet are rejected
  because V2000 has no lossless code for them.

## Validation

- Unit tests cover Molfile parse/write/parse round trips for radical multiplicity, supported source bond stereo marks, charge codes, isotope/map records, coordinates, and unsupported representations.
- RDKit-generated goldens compare Molfile-preservable atoms, bonds, coordinates, charges, isotopes, atom maps, and headers for external PubChem fixtures.

## Out Of Scope

- V3000 writing, canonical atom ordering, query features beyond supported V2000 stereo fields, atom stereo elements, and runtime RDKit.

## Revision Notes

- v1: V2000 writer.
- v2: Validation contract excludes SDF data fields and passes the RDKit-backed `smoke` corpus.
- v3: Write exact radical multiplicities and supported V2000 bond stereo codes; reject unsupported stereo/order combinations.
- v4: Move the public writer API under the `molfile` facade.
- v5: Add PubChem-100k as required broad-corpus validation evidence.
- v6: Read supported V2000 stereo output from source bond marks and reject stored stereo elements.
- v7: Preserve explicit no-implicit valence through the atom-block valence
  field and return a structured error for quartet/quintet radicals rather than
  coercing them to a CTfile radical code.
