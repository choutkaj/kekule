# Noncanonical SMILES Writer

## Summary

Write small molecules as deterministic noncanonical SMILES for round-trip workflows.

## Behavior/API

- Exposes `SmilesWriteOptions` and `write_smiles`.
- Emits graph-order-based noncanonical SMILES with branches, ring closures, dot fragments, common bond symbols, and bracket atoms when needed.
- Does not canonicalize or sanitize before writing.

## Implementation Notes

- The writer targets readability and deterministic output, not canonical ranking.
- A deterministic DFS tree is rendered with preassigned ring closures at both endpoints and branch children before the selected continuation path.
- Unsupported advanced stereo/query details are omitted until later feature work.

## Validation

- Unit tests cover parse/write/parse round trips for branches, rings, brackets, and fragments.
- RDKit-generated goldens compare noncanonical output SMILES for external PubChem SMILES fixtures.

## Out Of Scope

- Canonical SMILES, isomeric SMILES parity, SMARTS, reactions, and full stereochemical output.

## Revision Notes

- v1: Noncanonical writer.
- v2: Deterministic ring-closure and branch emission validated against RDKit goldens for current external fixtures.
