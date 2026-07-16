# mmCIF Molecular Interpretation

## Summary

Interpret one selected coordinate model from a loss-preserving `MmcifDocument`
into a `Model` with distinct molecule instances and a report.

## Behavior/API

- `mmcif::interpret` returns `MmcifInterpretation { model, report }` and never
  sanitizes or prepares chemistry.
- `MmcifModelSelection::RequireSingle` is the default and rejects multiple model
  IDs; `Select(String)` and `First` are explicit alternatives.
- Requires one complete finite position per interpreted atom after deterministic
  alternate-location selection.
- Treats `_atom_site.Cartn_*` coordinates as explicit angstrom quantities before
  model construction.
- Uses entity, structural-instance, atom-site, and declared-connection metadata
  to build Small/Macro instances. Only declared covalent links merge boundaries.
- Preserves `sing`, `doub`, `trip`, and `quad` values from
  `_struct_conn.pdbx_value_order`; a missing value defaults to single and an
  unsupported explicit value returns a structured interpretation error.
- Assigns conservative evidence-backed roles and exposes exact source
  classifications through report/provenance data.
- Reports selected and ignored models, altloc omissions, inferred entity kinds,
  applied/ignored/unresolved connections, and pending template connectivity.
- Reports every interpreted atom through `MmcifAtomProvenance`, qualified by
  `MoleculeInstanceId` and `InstanceAtomId`, with source line, atom-site,
  component, asymmetry, entity, and coordinate-model identifiers.
- Never writes mmCIF-specific labels into generic atom, molecule, or conformer
  property maps.

## Implementation Notes

- `SmcraHierarchy` maps labels to local `AtomId`; model insertion provides the
  instance-qualified view.
- Polymer/branched instances establish conservative macro boundaries; nonpolymer
  and water occurrences remain distinct unless a declared covalent link joins
  them.
- The document remains the loss-preserving source representation.

## Validation

- Tests cover mixed typed instances and roles, complete positions, default
  multi-model rejection, explicit selection, altloc policy/reporting, missing
  coordinates, covalent merging, noncovalent separation, supported connection
  order interpretation, and unknown-order rejection.
- Successful bounded fuzz parses traverse the loss-preserving document and then
  exercise explicit selected-model interpretation and qualified model lookup.

## Out Of Scope

- CCD/template lookup, inferred polymer bonds, assembly generation, sanitization,
  force-field preparation, and serialization.

## Revision Notes

- v1: Staged interpretation into molecular-content containers.
- v2: Remove direct whole-file molecule reader.
- v3: Hard break to selected-model `Model` output and remove
  `MolecularContents`/`Solvent`.
- v4: Preserve the four PDBx/mmCIF covalent bond orders carried by
  `_struct_conn.pdbx_value_order` instead of coercing every connection to single.
- v5: Return the renamed canonical `Model` and populate the fully prefixed
  `SmcraHierarchy` API without changing interpretation semantics.
- v6: Carry the mmCIF Cartesian angstrom convention through explicit conformer
  and model quantities.
- v7: Move all mmCIF labels and source identity into structured,
  instance-qualified interpretation provenance and keep core property maps
  format-neutral.
