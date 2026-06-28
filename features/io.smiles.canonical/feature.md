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
- The implementation is intentionally non-isomeric until stereochemistry perception and canonical stereo policy are available; stored atom and bond stereo metadata is ignored when writing canonical output.

## Validation

- Unit tests cover atom-order-independent tree output, component sorting, branch/ring round trips, and inherited unsupported-chemistry errors through the noncanonical writer contract.
- RDKit-generated tiny goldens compare exact non-isomeric canonical SMILES plus sanitized reparse semantics for external PubChem SMILES fixtures in the current non-fused-ring subset.
- RDKit-generated PubChem-100 and PubChem-1000 goldens compare sanitized reparse semantics for canonical output across all declared records. Validation sanitizes parsed fixtures before canonical writing to match RDKit's canonicalization input model. It does not apply a feature-specific unsupported-chemistry filter; parser, sanitizer, or writer gaps surface as validation failures.

## Out Of Scope

Isomeric SMILES, fused-ring canonical traversal parity, SMARTS, reactions, query atoms/bonds, radicals, unsupported bond orders, and full RDKit canonicalization parity for every symmetry edge case.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement deterministic non-stereo canonical SMILES for the existing writer subset.
- v3: Declare tiny RDKit canonical SMILES validation.
- v4: Declare PubChem-100 and PubChem-1000 semantic canonical-output validation.
- v5: Remove canonical-specific unsupported filtering so broad-corpus gaps are reported as validation failures.
- v6: Allow non-isomeric canonical output from stereo-bearing graphs by ignoring stored stereo metadata.
- v7: Sanitize canonical validation fixtures before writing so Kekule/aromatic normalization matches RDKit.
- v8: Preserve bracket hydrogens on metal-bound organic atoms so canonical output reparses with RDKit-like valence semantics.
