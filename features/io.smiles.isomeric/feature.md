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
- Rejects source bond stereo marks, double-bond elements, axial elements,
  enhanced stereo groups, and unknown/unspecified/invalid-cleared stereo until
  those layers have explicit writer support.
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

## Validation

Unit tests cover tetrahedral marker emission, marker flipping under odd writer
carrier-order permutations, round-trip parsing of emitted tetrahedral SMILES,
and rejection of source bond marks plus explicit unknown stereo.

Broad RDKit parity validation is not enabled for this first writer slice. The
next validation step is to add externally generated RDKit isomeric SMILES
semantic goldens once double-bond slash/backslash output is implemented, so the
corpus does not mostly exercise expected writer rejections.

## Out Of Scope

Canonical isomeric SMILES, double-bond slash/backslash output, axial SMILES
extensions, enhanced stereo groups, stereo enumeration, query stereochemistry,
and broad RDKit isomeric SMILES parity validation.

## Revision Notes

- v1: Add opt-in tetrahedral isomeric writer API over first-class stereo
  elements, with parity-aware marker emission and structured rejection for
  unsupported stereo layers.
