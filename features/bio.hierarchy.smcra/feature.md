# SMCRA-like Biomolecular Hierarchy

## Summary

Represent coordinate-independent chain, residue, and atom-site hierarchy as a
sidecar over one connected shared core molecule graph.

## Behavior/API

- Exposes `SmcraChain`, `SmcraResidue`, and `SmcraAtomSite` nodes plus
  correspondingly prefixed typed IDs.
- Stores biomolecular labels and atom-site metadata in `SmcraHierarchy`, not core `Atom`.
- `MacroMolecule` owns one connected core `Molecule` plus one
  `SmcraHierarchy`.
- `MacroMolecule` exposes chain, residue, and atom-site iterators plus
  `atom_site_for_atom`.
- `MacroMolecule::validate` checks graph/hierarchy consistency and all retained
  conformer coordinates without mutation.
- Static validation rejects a graph with multiple components through
  `MacroValidateError::DisconnectedGraph`.
- `MacroMolecule::validate_with_options` can restrict work to static
  graph/hierarchy validation for coordinate-independent consumers.
- `MacroMoleculeBuilder` and `MacroMolecule::try_from_parts` are the only raw
  assembly paths and reject invalid graph/hierarchy pairs, including a
  disconnected graph.
- Public macro graph growth starts with `add_atom` and then uses the atomic
  `add_atom_bonded_to` operation; a second unattached atom returns
  `MacroGraphEditError::DisconnectedAtomInsertion` without mutation.
- `MacroMolecule::edit` provides coordinated transactional graph/hierarchy
  mutation; its working copy may be temporarily disconnected, but commit
  validates connectedness before atomically replacing the original.
- Every live graph atom must have exactly one atom site, and every atom site
  must reference a live graph atom.
- Atom-site insertion validates that referenced core atoms exist.
- Chain, residue, and atom-site insertion checks fixed-width capacity before
  mutating parent/lookup state and reports
  `SmcraHierarchyError::IdentifierCapacityExceeded` with the affected ID kind.
- Atom-site metadata preserves static type, label/auth chain, and label/auth
  atom identity. Occupancy and B-factor live in topology-bound `AtomData`;
  format-specific coordinate-model identity remains interpretation provenance.
- Checked targeted restoration sets distinct residue label/author component
  IDs and returns mutable chain, residue, or atom-site property maps without
  exposing parentage, child arrays, complete mutable records, or atom lookup.

## Implementation Notes

- Preserves insertion order for hierarchy iteration.
- Tracks label and author identifiers separately.
- Supports insertion codes and distinct label/author identifiers without
  structural coordinate-model parent nodes.
- Stores label and author component IDs separately on residues.
- Child enrichment validates the requested typed ID before mutation and is
  usable through both builder and transactional editor hierarchy access.
- mmCIF interpretation populates hierarchy only after molecular boundaries and alternate locations have been resolved.

## Tests

- Unit tests cover hierarchy construction, checked assembly, disconnected
  builder/parts rejection, transactional mutation, disconnecting-edit rollback,
  lookup, full versus static-only validation reports, unused conformer
  skipping, synthetic ID-capacity boundaries, capacity rollback, and
  failed-commit rollback.
- The former Biopython evidence exercised the removed whole-file reader rather
  than the format-neutral hierarchy contract, so no current hierarchy parity
  evidence is recorded pending a replacement comparison.
- Canonical reconstruction tests preserve arbitrary child properties, distinct
  residue name/label/author component IDs, all atom-site metadata, complete
  hierarchy equality, and invalid-ID rollback.

## Out Of Scope

- PDB parsing, full mmCIF category coverage, polymer connectivity, sequence extraction, and chemical perception.
- Runtime Biopython dependency.

## Revision Notes

- v1: SMCRA sidecar hierarchy for macromolecular parsing.
- v2: Preserve atom-site row metadata and distinguish author-keyed residues when label sequence IDs are absent.
- v3: Preserve label/auth component IDs separately and support conservative lenient occurrence grouping.
- v4: Add direct macro hierarchy accessors plus conservative macro validation and sanitization APIs.
- v5: Make macro sanitization defaults honest by enabling only implemented validation behavior and rejecting requested unimplemented stages.
- v6: Remove validation coupling to the deleted direct mmCIF reader and keep `SmcraHierarchy` format-neutral.
- v7: Hard-break the complete hierarchy vocabulary to `Smcra*` names so
  structural hierarchy nodes cannot be confused with concrete configurations.
- v8: Enforce a valid-state `MacroMolecule` boundary with checked builders,
  complete graph-to-atom-site coverage, and transactional coordinated editing;
  remove the placeholder macromolecule sanitization surface.
- v9: Remove coordinate-model nodes from the structural hierarchy and move
  occupancy, B-factor, alternate-location, source-row, and raw-coordinate state
  into topology-bound structure observations.
- v10: Document static-only validation for coordinate-independent topology
  consumers while retaining all-conformer validation as the standalone
  `MacroMolecule::validate` default.
- v11: Add structured chain, residue, and atom-site capacity errors and verify
  hierarchy insertion remains transactional at the fixed-width boundary.
- v12: Add checked targeted restoration for residue component IDs and
  chain/residue/atom-site property maps.
- v13: Require a connected graph at every completed `MacroMolecule` boundary,
  expose structured disconnected validation, and preserve the original value
  when a staged edit would disconnect it.
