# mmCIF Molecular Interpretation

## Summary

Interpret one explicitly selected coordinate model from a loss-preserving
`MmcifDocument` into separated topology, configuration, observation state, and
report, or use a separate path for a verified shared-topology ensemble.

## Behavior/API

- `mmcif::interpret` returns `MmcifInterpretation { model, report }` and never
  sanitizes or prepares chemistry.
- `MmcifInterpretation::into_model` consumes the interpretation and returns its
  canonical `Model` when the report is not needed.
- `mmcif::interpret_ensemble` separately interprets selected or all coordinate
  models, proves molecule partition, chemistry, connectivity, semantic atom
  identity, and dense order consistency, and returns one shared-topology
  `Ensemble`; inconsistent atom sets or topology fail structurally.
- `MmcifModelSelection::RequireSingle` is the default and rejects multiple model
  IDs; `Select(String)` and `First` are explicit alternatives.
- Requires one complete finite position per interpreted atom after deterministic
  alternate-location selection.
- Treats `_atom_site.Cartn_*` coordinates as explicit angstrom quantities before
  model construction.
- Uses entity, structural-instance, atom-site, polymer-sequence, and
  declared-connection metadata to build Small/Macro instances. Only declared
  local covalent links merge boundaries; symmetry-mate links remain unresolved.
- Preserves `sing`, `doub`, `trip`, and `quad` values from
  `_struct_conn.pdbx_value_order`; a missing value defaults to single and an
  unsupported explicit value returns a structured interpretation error.
- Assigns conservative evidence-backed roles and exposes exact source
  classifications through report/provenance data.
- Reports selected and ignored models, altloc omissions, inferred entity kinds,
  applied/ignored/unresolved connections, pending template connectivity, and
  distance-based connectivity candidates.
- Reports every interpreted atom through `MmcifAtomProvenance`, qualified by
  `MoleculeInstanceId` and `InstanceAtomId`, with source line, atom-site,
  component, asymmetry, entity, and coordinate-model identifiers.
- Never writes mmCIF-specific labels into generic atom, molecule, or conformer
  property maps.
- Preserves source coordinate-model ID, selected alternate location, occupancy,
  B-factor, source atom-site ID, and raw Cartesian fields in topology-bound
  `StructureObservation` records rather than static hierarchy.
- Distance heuristics report connectivity candidates but do not assert
  authoritative single bonds without evidence-backed bond order.

## Implementation Notes

- `SmcraHierarchy` maps labels to local `AtomId`; model insertion provides the
  instance-qualified view.
- Polymer/branched instances establish conservative macro boundaries; nonpolymer
  and water occurrences remain distinct unless a declared covalent link joins
  them.
- Within merged macro instances, hierarchy chains follow their first
  `_pdbx_poly_seq_scheme.asym_id` occurrence rather than incidental atom-row or
  map order.
- The document remains the loss-preserving source representation.

## Tests

- Tests cover mixed typed instances and roles, complete positions, default
  multi-model rejection, explicit selection, altloc policy/reporting, missing
  coordinates, covalent merging, noncovalent separation, symmetry-mate
  rejection, deposited polymer-chain ordering, supported connection order
  interpretation, and unknown-order rejection.
- Multi-model tests cover shared-topology ensemble construction, distinct
  per-member source IDs/occupancy/B-factors, and structured inconsistent atom
  set rejection.
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
- v8: Preserve deposited polymer-chain order within merged instances and report
  symmetry-mate connections as unresolved instead of creating local self-bonds.
- v9: Add `MmcifInterpretation::into_model` for direct consuming access to the
  canonical model while retaining the report-bearing interpretation contract.
- v10: Separate topology, configuration, and observation state; report
  distance-based connectivity only as candidates; and add a distinct,
  consistency-proving multi-model ensemble interpretation path.
