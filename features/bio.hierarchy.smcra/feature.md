# SMCRA-like Biomolecular Hierarchy

## Summary

Represent model, chain, residue, and atom-site hierarchy as a sidecar over the shared core molecule graph.

## Behavior/API

- Exposes `SmcraModel`, `SmcraChain`, `SmcraResidue`, and `SmcraAtomSite`
  nodes plus correspondingly prefixed typed IDs.
- Stores biomolecular labels and atom-site metadata in `SmcraHierarchy`, not core `Atom`.
- `MacroMolecule` owns one core `Molecule` plus one `SmcraHierarchy`.
- `MacroMolecule` exposes model, chain, residue, and atom-site iterators plus `atom_site_for_atom`.
- `MacroMolecule::validate` checks graph/hierarchy/coordinate consistency without mutation.
- `MacroMolecule::sanitize` uses separate macro options, reports, and errors; defaults only validate current graph/hierarchy/coordinate state.
- Unsupported normalization, residue recognition, connectivity, alternate-location selection, and ligand sanitization requests fail explicitly.
- Atom-site insertion validates that referenced core atoms exist.
- Atom-site metadata preserves `_atom_site.group_PDB`, `_atom_site.id`, label/auth atom IDs, alternate location, occupancy, and B-factor.

## Implementation Notes

- Preserves insertion order for hierarchy iteration.
- Tracks label and author identifiers separately.
- Supports alternate locations, occupancy, B-factor, insertion code, and model identifiers.
- Stores label and author component IDs separately on residues.
- mmCIF interpretation populates hierarchy only after molecular boundaries and alternate locations have been resolved.

## Validation

- Unit tests cover hierarchy construction, mutation, lookup, validation, and sanitization behavior.
- The former Biopython evidence exercised the removed whole-file reader rather
  than the format-neutral hierarchy contract, so no current hierarchy parity
  evidence is recorded pending a replacement comparison.

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
- v7: Hard-break the complete hierarchy vocabulary to `Smcra*` names so its
  structural model cannot be confused with `modeling::Model`.
