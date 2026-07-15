# Canonical SMILES Writer

## Summary

Write deterministic non-stereo canonical SMILES for supported small-molecule graphs.

## Behavior/API

- Exposes `smiles::{CanonicalSmilesWriteOptions, write_canonical, write_canonical_with_options}`.
- Reuses the noncanonical writer's supported chemistry subset and structured write errors.
- Chooses a deterministic representation by ranking atoms, trying every atom in each connected component as a root, rendering rank-ordered branches and ring closures, and selecting the smallest rank-guided component string using SMILES syntax tie-breakers.
- Sorts disconnected component strings before joining with `.`.
- Does not sanitize or perceive chemistry before writing.
- Applies non-isomeric normalization on a clone: isotope labels are suppressed,
  ordinary explicit hydrogen atoms are collapsed into parent hydrogen counts
  when safe, isotopic hydrogen atoms remain graph atoms but lose their isotope
  label, and the caller's molecule and perception caches are unchanged.

## Implementation Notes

- Builds on `canon::atom_ranking` for atom symmetry classes.
- Follows the RDKit-inspired split between canonical atom ranking and canonical traversal/output: graph-derived atom and bond invariants drive traversal, branch order, ring closure order, and disconnected-component ordering.
- Symmetric ties are handled by candidate string selection across roots and bond-order traversal preferences, with `AtomId` only as a final deterministic fallback inside rank-equivalent traversal choices.
- Canonical candidate ranking is derived from the graph and emitted SMILES syntax only; it does not reparse candidates, run sanitization, or switch to motif-specific stored-Kekule fallback spellings.
- Uses stored Kekule atom/bond spelling when a mixed aromatic/aliphatic pi component has non-aromatic multiple-bonded framework atoms that aromatic shorthand cannot represent without seeding a different aromaticity partition on reparse, and only when concrete stored single/double bond orders are available.
- Preserves explicit hydrogens and no-implicit organic atoms when organic shorthand cannot represent the stored atom state, including metal/main-group neighbor cases handled by the shared SMILES parse/write fidelity rules.
- Uses the shared charge-adjusted valence rules to choose organic shorthand and
  bracket hydrogens, and retains stored triple/quadruple orders even on bonds
  that belong to an aromaticized framework.
- The implementation is intentionally non-isomeric until stereochemistry perception and canonical stereo policy are available; first-class stereo elements and source bond marks are ignored only in this explicit canonical non-isomeric output mode.

## Validation

- Unit tests cover atom-order-independent tree output, component sorting, branch/ring round trips, and inherited unsupported-chemistry errors through the noncanonical writer contract.
- RDKit-generated smoke goldens compare exact non-isomeric canonical SMILES plus sanitized reparse semantics for external PubChem SMILES fixtures in the current non-fused-ring subset.
- RDKit-generated PubChem-100 and pubchem-1k goldens compare sanitized reparse semantics for canonical output across all declared records. Validation sanitizes parsed fixtures before canonical writing to match RDKit's canonicalization input model. It does not apply a feature-specific unsupported-chemistry filter; parser, sanitizer, or writer gaps surface as validation failures.

## Out Of Scope

Isomeric SMILES, fused-ring canonical traversal parity, SMARTS, reactions,
query atoms/bonds, radical states that are not implied by bracket valence,
unsupported bond orders, and full RDKit canonicalization parity for every
symmetry edge case.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement deterministic non-stereo canonical SMILES for the existing writer subset.
- v3: Declare smoke RDKit canonical SMILES validation.
- v4: Declare PubChem-100 and pubchem-1k semantic canonical-output validation.
- v5: Remove canonical-specific unsupported filtering so broad-corpus gaps are reported as validation failures.
- v6: Allow non-isomeric canonical output from stereo-bearing graphs by ignoring stored stereo metadata.
- v7: Sanitize canonical validation fixtures before writing so Kekule/aromatic normalization matches RDKit.
- v8: Preserve bracket hydrogens on metal-bound organic atoms so canonical output reparses with RDKit-like valence semantics.
- v9: Prefer aromatic continuations as canonical main paths so fused heteroaromatic branches reparse with stable aromaticity.
- v10: Rank canonical candidates by sanitized semantic preservation before string shape, improving lactone, fused aromatic, and topology-sensitive round trips.
- v11: Broaden metal-like neighbor detection for bracketed organic hydrogens, preserving organothallium and related main-group/actinide valence semantics.
- v12: Incorporate directional-bond parse support, metal-bound halogen reparse semantics, and neutral imide aromaticity cleanup exposed by pubchem-1k canonical validation.
- v13: Advance pubchem-1k semantic validation through additional fused saturated-ring and cyclic amidine aromaticity cases.
- v14: Advance pubchem-1k validation through aromatic Se/Te bracket reparse support and valence-filled metal-bound organic no-implicit preservation shared with the SMILES parser.
- v15: Advance pubchem-1k validation through saturated sulfonamide fused-ring aromaticity cleanup.
- v16: Advance pubchem-1k validation through exocyclic alkene deactivation in fused nitrogen/chalcogen aromatic systems.
- v17: Advance pubchem-1k validation through canonical reparse support for thione-rich imported nitrogen/chalcogen aromatic-order rings.
- v18: Advance pubchem-1k validation through RDKit-like fused lactam/enone and saturated oxygen bridge aromaticity cleanup for canonical reparse semantics.
- v19: Advance pubchem-1k validation through RDKit-like saturated fused nitrogen carbonyl aromaticity cleanup for benzodiazepinone lactam canonical reparse semantics.
- v20: Preserve no-implicit aromatic organic atoms bound to germanium-like main-group centers in canonical output, advancing pubchem-1k validation past aryl germanium trichloride.
- v21: Advance pubchem-1k validation through imported aromatic-order support for five-member nitrogen/chalcogen rings with exocyclic cationic imine pi bonds.
- v22: Rank canonical candidates with local atom-neighbor semantics and add a lazy stored-Kekule fallback for aromatic carbocyclic components without aromatic heteroatoms, advancing pubchem-1k past oxygen-rich multicomponent lactone/carboxyl mixtures while keeping heteroaromatic nitrogen output in aromatic form.
- v23: Move the public canonical writer API under the `smiles` facade and the ranking dependency under `canon`.
- v24: Remove candidate reparse/sanitize ranking and motif-specific stored-Kekule cleanup fallback; canonical output is selected by RDKit-inspired rank-guided graph traversal plus syntax tie-breakers, with graph-derived stored-Kekule emission for mixed aromatic/aliphatic pi components that cannot be faithfully represented by aromatic shorthand.
- v25: Complete pubchem-1k canonical validation by broadening graph-derived representability rules for aromatic carbonyls, charged aromatic carbon components, and zero-hydrogen metal/main-group-bound organic atoms, plus aligning validation semantics for anionic aromatic nitrogen and charged aromatic carbon hydrogens.
- v26: Add PubChem-100k as required broad-corpus validation evidence.
- v27: Document that canonical non-isomeric output ignores first-class stereo elements and source bond marks by policy, rather than reading removed atom/bond stereo metadata.
- v28: Complete PubChem-100k non-isomeric semantics by normalizing isotope and
  explicit-hydrogen representation on a clone, using periodic-table valence for
  bracket/organic spelling, preserving high-order aromatic-framework bonds, and
  writing valence-implied radical states.
- v29: Preserve collapsed explicit-hydrogen counts while canonicalizing an
  unsanitized molecule without falsely marking full valence perception as
  installed.
- v30: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
