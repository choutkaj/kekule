# Basic RDKit-like Aromaticity Perception

## Summary

Assign aromatic atom and bond flags for common organic ring systems using an initial RDKit-like model. This is a perception step and must remain separate from parsing.

## Behavior/API

- Exposes `AromaticityModel::RdkitLikeBasic` and `perceive_aromaticity`.
- Requires or computes ring membership before assigning aromaticity.
- Marks supported aromatic atoms and bonds and sets aromaticity perception state to fresh.
- Clears prior aromatic flags deterministically before assignment.

## Implementation Notes

- Operates on the shared core `Molecule` graph.
- Uses ring membership from `algo.rings.fast`.
- Applies a basic 4n+2 electron-count model for common C, N, O, S, and P rings.
- Leaves unsupported or ambiguous systems non-aromatic rather than claiming full RDKit parity.

## Validation

- Current coverage is unit-test based.
- RDKit golden validation is planned through `validation.harness`.
- Fixtures live under `validation/features/algo.aromaticity.rdkit-like-basic/`.

## Out Of Scope

- Full RDKit aromaticity parity.
- Valence perception, sanitization, kekulization, stereochemistry, and parser behavior.
- Runtime RDKit dependency.

## Revision Notes

- v1: Initial basic aromaticity perception for common organic rings.
