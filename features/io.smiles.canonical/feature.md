# Canonical SMILES Writer

## Summary

Write deterministic non-stereo canonical SMILES for supported small-molecule graphs.

## Behavior/API

- Exposes `smiles::{CanonicalSmilesWriteOptions, write_canonical, write_canonical_with_options}`.
- Reuses the noncanonical writer's supported chemistry subset and structured write errors.
- Chooses a deterministic representation by ranking atoms, trying every atom in each connected component as a root, rendering rank-ordered branches and ring closures, and selecting the smallest component string that preserves sanitized local atom-neighbor semantics when one is available.
- Sorts disconnected component strings before joining with `.`.
- Does not sanitize or perceive chemistry before writing.

## Implementation Notes

- Builds on `canon::atom_ranking` for atom symmetry classes.
- Symmetric ties are handled by candidate string selection, with `AtomId` only as a final deterministic fallback inside rank-equivalent traversal choices.
- Ranks normal aromatic SMILES candidates first, then lazily tries a stored-Kekule candidate family for aromatic components without aromatic heteroatoms when aromatic spelling would alter sanitized reparse topology.
- Preserves explicit hydrogens on organic atoms bound to a broad metal-like element set when organic shorthand would alter sanitized valence semantics after reparse.
- Preserves no-implicit organic atoms bound to broad metal-like/main-group neighbors when organic shorthand would alter sanitized atom semantics after reparse.
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
- v9: Prefer aromatic continuations as canonical main paths so fused heteroaromatic branches reparse with stable aromaticity.
- v10: Rank canonical candidates by sanitized semantic preservation before string shape, improving lactone, fused aromatic, and topology-sensitive round trips.
- v11: Broaden metal-like neighbor detection for bracketed organic hydrogens, preserving organothallium and related main-group/actinide valence semantics.
- v12: Incorporate directional-bond parse support, metal-bound halogen reparse semantics, and neutral imide aromaticity cleanup exposed by PubChem-1000 canonical validation.
- v13: Advance PubChem-1000 semantic validation through additional fused saturated-ring and cyclic amidine aromaticity cases.
- v14: Advance PubChem-1000 validation through aromatic Se/Te bracket reparse support and valence-filled metal-bound organic no-implicit preservation shared with the SMILES parser.
- v15: Advance PubChem-1000 validation through saturated sulfonamide fused-ring aromaticity cleanup.
- v16: Advance PubChem-1000 validation through exocyclic alkene deactivation in fused nitrogen/chalcogen aromatic systems.
- v17: Advance PubChem-1000 validation through canonical reparse support for thione-rich imported nitrogen/chalcogen aromatic-order rings.
- v18: Advance PubChem-1000 validation through RDKit-like fused lactam/enone and saturated oxygen bridge aromaticity cleanup for canonical reparse semantics.
- v19: Advance PubChem-1000 validation through RDKit-like saturated fused nitrogen carbonyl aromaticity cleanup for benzodiazepinone lactam canonical reparse semantics.
- v20: Preserve no-implicit aromatic organic atoms bound to germanium-like main-group centers in canonical output, advancing PubChem-1000 validation past aryl germanium trichloride.
- v21: Advance PubChem-1000 validation through imported aromatic-order support for five-member nitrogen/chalcogen rings with exocyclic cationic imine pi bonds.
- v22: Rank canonical candidates with local atom-neighbor semantics and add a lazy stored-Kekule fallback for aromatic carbocyclic components without aromatic heteroatoms, advancing PubChem-1000 past oxygen-rich multicomponent lactone/carboxyl mixtures while keeping heteroaromatic nitrogen output in aromatic form.
- v23: Move the public canonical writer API under the `smiles` facade and the ranking dependency under `canon`.
