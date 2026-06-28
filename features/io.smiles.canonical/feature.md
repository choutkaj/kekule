# Canonical SMILES Writer

## Summary

Write deterministic non-stereo canonical SMILES for supported small-molecule graphs.

## Behavior/API

- Exposes `CanonicalSmilesWriteOptions` and `write_canonical_smiles`.
- Reuses the noncanonical writer's supported chemistry subset and structured write errors.
- Chooses a deterministic representation by ranking atoms, trying every atom in each connected component as a root, rendering rank-ordered branches and ring closures, and selecting the lexicographically smallest component string.
- Sorts disconnected component strings before joining with `.`.
- Does not sanitize or perceive chemistry before writing.

## Implementation Notes

- Builds on `canonical_atom_ranking` for atom symmetry classes.
- Symmetric ties are handled by candidate string selection, with `AtomId` only as a final deterministic fallback inside rank-equivalent traversal choices.
- The implementation is intentionally non-isomeric until stereochemistry perception and canonical stereo policy are available.

## Validation

- Unit tests cover atom-order-independent tree output, component sorting, branch/ring round trips, and inherited unsupported-chemistry errors through the noncanonical writer contract.
- Future golden validation should compare compact RDKit canonical SMILES fixtures for the supported non-stereo subset.

## Out Of Scope

Isomeric SMILES, SMARTS, reactions, query atoms/bonds, radicals, unsupported bond orders, and full RDKit canonicalization parity for every symmetry edge case.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement deterministic non-stereo canonical SMILES for the existing writer subset.
