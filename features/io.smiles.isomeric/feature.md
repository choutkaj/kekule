# Isomeric SMILES Writer

## Summary

Write noncanonical isomeric SMILES from the first-class stereo representation.

## Behavior/API

- Exposes `smiles::{IsomericSmilesWriteOptions, write_isomeric, write_isomeric_with_options}`.
- Exposes `SmallMolecule::to_isomeric_smiles()`.
- Emits tetrahedral `@`/`@@` markers for specified tetrahedral stereo elements.
- Computes the marker from the writer's actual emitted carrier order, so branch,
  continuation, hydrogen, and ring-closure ordering preserve the stored local
  tetrahedral parity instead of assuming the original parse order.
- Emits `/` and `\` markers for specified double-bond stereo elements when both
  stored carriers are explicit atoms connected to the double-bond endpoints by
  printable single bonds.
- Allows source directional bond marks only when they are covered by stored
  double-bond stereo elements; unassembled source marks are rejected rather than
  treated as authoritative stereo.
- Rejects implicit-carrier double-bond elements, axial elements, enhanced
  stereo groups, non-directional source bond marks, and
  unknown/unspecified/invalid-cleared stereo until those layers have explicit
  writer support.
- Does not sanitize, perceive stereo, or assign CIP descriptors before writing.

## Implementation Notes

This feature is an opt-in writer surface layered over `io.smiles.write`. The
plain noncanonical writer continues to reject stereo, and canonical output
continues to be explicitly non-isomeric.

The tetrahedral writer builds a small per-atom stereo context from
`StereoElementKind::Tetrahedral`. For each emitted chiral atom, it reconstructs
the carrier order that a parser would see in the generated SMILES: parent atom,
explicit bracket hydrogen when needed, ring closures, branch children, main
continuation child, and implicit lone pair when present. If that order is an
odd permutation of the stored carrier order, the emitted `@`/`@@` marker is
flipped.

The double-bond writer builds endpoint-local directional constraints from
`StereoElementKind::DoubleBond`. It chooses a local direction for the left
endpoint carrier, chooses the same or opposite local direction on the right
endpoint from the stored `Together`/`Opposite` relation, and resolves the
concrete slash or backslash at the moment the adjacent single bond is emitted.
This mirrors the perception layer's endpoint-relative normalization and lets a
shared directional bond report a conflict instead of silently dropping one
constraint.

## Validation

Unit tests cover tetrahedral marker emission, marker flipping under odd writer
carrier-order permutations, round-trip parsing of emitted tetrahedral SMILES,
directional double-bond slash/backslash emission for `Together` and `Opposite`
elements, semantic round-trip perception of emitted double-bond stereo, and
rejection of unperceived source bond marks plus explicit unknown stereo.

Broad RDKit parity validation is not enabled for this writer feature yet. The
next validation step is to add an `io.smiles.isomeric` implementation
comparison mode plus externally generated RDKit semantic goldens for smoke and
the existing stereo SMILES corpora.

## Out Of Scope

Canonical isomeric SMILES, implicit-hydrogen double-bond carrier output, axial
SMILES extensions, enhanced stereo groups, stereo enumeration, query
stereochemistry, and broad RDKit isomeric SMILES parity validation.

## Revision Notes

- v1: Add opt-in tetrahedral isomeric writer API over first-class stereo
  elements, with parity-aware marker emission and structured rejection for
  unsupported stereo layers.
- v2: Emit endpoint-normalized slash/backslash marks for explicit-carrier
  double-bond stereo elements while rejecting unassembled source marks and
  unsupported implicit-carrier double-bond output.
