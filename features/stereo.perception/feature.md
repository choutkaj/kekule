# Stereochemistry Perception

## Summary

Validate graph-local stereo elements, detect candidate stereochemical units,
and assemble supported source marks into local stereo elements.

## Behavior/API

- Exposes `perception::stereo::{StereoPerceptionOptions,
  StereoPerceptionReport, StereoCandidate, StereoPerceptionIssue,
  validate_stereo, validate_stereo_with_options, perceive_stereo,
  perceive_stereo_with_options}`.
- `validate_stereo` is read-only. It reports invalid local stereo elements,
  potential tetrahedral atom candidates, potential double-bond candidates, and
  source-mark assembly diagnostics without mutating the graph.
- `perceive_stereo` is mutating. It runs the same checks, assembles supported
  source marks into first-class stereo elements, records created stereo element
  IDs, and marks stereo perception fresh.
- Candidate detection uses current graph and hydrogen state. Sanitization or
  valence perception should run first when implicit hydrogens matter.
- Tetrahedral candidates are local geometry candidates only; no exact ligand
  ranking, symmetry pruning, or CIP descriptor assignment is performed.
- Double-bond candidates include atom and implicit-hydrogen carriers. Paired
  directional `/` and `\` source marks can create specified double-bond stereo
  elements when one marked single bond is available on each end.
- Molfile wedge up/down marks can create specified tetrahedral stereo elements.
  Molfile wedge/either marks can create explicit unknown tetrahedral stereo
  elements. In both cases the marked bond's first endpoint is treated as the
  local stereo center and the marked carrier is placed first in carrier order.
- Coordinate-derived assignment remains out of scope.

## Implementation Notes

This feature identifies candidate tetrahedral atoms and double bonds, validates
existing local stereo elements against current topology and hydrogen semantics,
assembles SMILES-style paired directional bond marks, and assembles supported
Molfile tetrahedral wedge/either source marks. It does not assign exact CIP
descriptors directly; that belongs to `stereo.cip`.

Small-molecule perception should run as an explicit staged workflow and may be
integrated into sanitization later. It should not run over whole
`MacroMolecule` structures by default.

## Validation

- Unit tests cover read-only validation, candidate detection after sanitization,
  directional double-bond assembly, Molfile wedge/either assembly, unsupported
  source-mark diagnostics, and preservation of explicit unknown versus absent
  stereo.
- Smoke validation records semantic perception JSON for externally pinned
  PubChem fixtures covering absent stereo, stored tetrahedral stereo, and
  directional double-bond source-mark assembly.

## Out Of Scope

Exact CIP descriptors, isomeric SMILES writing, enhanced stereo serialization,
coordinate-derived stereo assignment, stereo enumeration, and reaction stereo
transfer.

## Revision Notes

- v1: Feature contract reserved for stereo candidate detection and local
  validation.
- v2: Add public stereo perception API, local element validation, conservative
  tetrahedral/double-bond candidate detection, paired directional source-mark
  assembly, unit coverage, and smoke semantic validation.
- v3: Assemble supported Molfile wedge up/down source marks into specified
  tetrahedral elements and wedge/either source marks into explicit unknown
  tetrahedral elements.
