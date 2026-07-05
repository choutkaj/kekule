# Stereochemistry Perception

## Summary

Plan the explicit perception layer that turns graph, coordinate, and source
bond-mark information into validated local stereo elements.

## Behavior/API

No public API is implemented yet.

## Implementation Notes

This feature should identify candidate tetrahedral atoms, double bonds, and
reserved axial units; validate existing local stereo elements against current
topology and hydrogen semantics; assign local orientation from 2D wedges or 3D
coordinates; and repair or clear invalid stereo with explicit diagnostics. It
should not assign exact CIP descriptors directly; that belongs to
`stereo.cip`.

Small-molecule perception should run as an explicit staged workflow and may be
integrated into sanitization later. It should not run over whole
`MacroMolecule` structures by default.

## Validation

Future validation should include manually reviewed fixtures for tetrahedral
centers, double-bond systems, partial/unknown stereo, coordinate-derived wedge
assignment, invalid-cleared stereo, and stereo-preserving topology edits.

## Out Of Scope

Exact CIP descriptors, isomeric SMILES writing, enhanced stereo serialization,
stereo enumeration, and reaction stereo transfer.

## Revision Notes

- v1: Feature contract reserved for stereo candidate detection and local
  validation.
