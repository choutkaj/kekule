# mmCIF Molecular Interpretation

## Summary

Interpret one selected coordinate model from a loss-preserving `MmcifDocument`
into a `MolecularModel` with distinct molecule instances and a report.

## Behavior/API

- `mmcif::interpret` returns `MmcifInterpretation { model, report }` and never
  sanitizes or prepares chemistry.
- `MmcifModelSelection::RequireSingle` is the default and rejects multiple model
  IDs; `Select(String)` and `First` are explicit alternatives.
- Requires one complete finite position per interpreted atom after deterministic
  alternate-location selection.
- Uses entity, structural-instance, atom-site, and declared-connection metadata
  to build Small/Macro instances. Only declared covalent links merge boundaries.
- Assigns conservative evidence-backed roles and exposes exact source
  classifications through graph properties/report data.
- Reports selected and ignored models, altloc omissions, inferred entity kinds,
  applied/ignored/unresolved connections, and pending template connectivity.

## Implementation Notes

- `BioHierarchy` maps labels to local `AtomId`; model insertion provides the
  instance-qualified view.
- Polymer/branched instances establish conservative macro boundaries; nonpolymer
  and water occurrences remain distinct unless a declared covalent link joins
  them.
- The document remains the loss-preserving source representation.

## Validation

- Tests cover mixed typed instances and roles, complete positions, default
  multi-model rejection, explicit selection, altloc policy/reporting, missing
  coordinates, covalent merging, and noncovalent separation.
- Successful bounded fuzz parses traverse the loss-preserving document and then
  exercise explicit selected-model interpretation and qualified model lookup.

## Out Of Scope

- CCD/template lookup, inferred polymer bonds, assembly generation, sanitization,
  force-field preparation, and serialization.

## Revision Notes

- v1: Staged interpretation into molecular-content containers.
- v2: Remove direct whole-file molecule reader.
- v3: Hard break to selected-model `MolecularModel` output and remove
  `MolecularContents`/`Solvent`.
