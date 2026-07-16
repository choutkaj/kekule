# RDKit-like Aromaticity Perception

## Summary

Assign aromatic atom and bond membership in canonical `PerceptionState` for
supported organic ring systems using an RDKit-like graph aromaticity model.

## Behavior/API

- Exposes `perception::aromaticity::{AromaticityModel::RdkitLike, perceive_aromaticity, perceive_aromaticity_with_ring_options}`.
- Requires or computes ring perception before assigning aromaticity.
- Replaces existing imported/perceived aromatic membership transactionally and
  records `AromaticityModel::RdkitLike` provenance.
- Can be run directly or through the explicit sanitization pipeline.
- Treats unsupported ring elements as non-candidates, allowing a supported
  aromatic subring to remain aromatic when fused or attached to a nonaromatic
  unsupported-element ring.
- Converts an imported aromatic-order component that has no valid aromatic
  atoms to deterministic localized single/double bonds when a bounded
  valence-demand matching exists; otherwise returns
  `InvalidAromaticRepresentation` transactionally.
- Reports imported-aromatic matching budget exhaustion separately as
  `ImportedAromaticKekulizationLimit`; it is never presented as invalid
  chemistry.
- Propagates bounded ring-perception failures as `AromaticityError::RingPerception`.

## Implementation Notes

- Operates on the shared core `Molecule` graph and uses ring data from the `algo.rings.sssr` stack.
- Localizes imported aromatic-order components through one general
  valence-demand matching step, then sends imported and already-localized input
  through the same donor/candidate/Huckel engine.
- Uses one RDKit-style donor classifier for aromatic candidate checks and
  simple- and fused-component electron counting.
- Supports common RDKit organic aromatic elements: C, N, O, P, S, Se, and Te.
- Applies RDKit-like candidate gates for atom degree, explicit pi-bond count, triple bonds, exocyclic multiple-bond options, charge-adjusted default valence, and radical eligibility.
- Counts localized saturated, vacant, lone-pair, anionic, and pi-bond donors with a `countAtomElec`-style helper using default valence, outer-shell electrons, charge, radical electrons, effective hydrogens, and exocyclic electronegativity.
- Evaluates localized simple rings of arbitrary size through the same Huckel donor-count path, including two-electron rings such as cyclopropenyl cation.
- Treats exocyclic pi bonds through electronegativity-aware donor logic rather than raw hetero-atom symbol checks.
- Applies a bounded RDKit-style fused-system pass: fused candidate rings are grouped by shared bonds, connected subsets are evaluated from small to large, subset atom sets use fused-ring multiplicity, and accepted subsets mark perimeter bonds.
- Uses RDKit-like fused-system atom multiplicity, selected-subsystem perimeter bonds, additive accepted subsets, and the 24-atom fused-ring candidate cap.
- Does not run molecule- or functional-group-specific cleanup. Carbonyl,
  heteroatom, radical, charge, and fused-ring behavior emerges from the shared
  valence, donor, candidate, and graph-topology rules.
- Keeps parsing separate from aromaticity perception. Canonical SMILES normalization issues exposed by these flags belong to `io.smiles.canonical`, not hidden aromaticity cleanup.
- Treats unsupported or ambiguous systems conservatively rather than claiming full RDKit parity.
- The direct public API stages the complete operation, including ring
  perception and imported-order localization, and commits only on success.

## Validation

- Unit tests cover localized donor analysis, candidate gates, radical and valence eligibility, imported aromatic-order handling, fused-subsystem search, and SMILES sanitize/reparse smoke cases.
- RDKit-generated goldens compare aromatic atom and bond flags for external PubChem fixtures.

## Out Of Scope

- Full RDKit aromaticity parity.
- Runtime RDKit dependency.
- Valence perception, sanitization policy, general-purpose Kekule-form
  enumeration, stereochemistry, and parser behavior.
- Canonical SMILES normalization for every valid aromaticity assignment.

## Revision Notes

- v1-v83: Built the RDKit-like donor classifier, fused-subsystem search, validation workflow, and public expert facade.
- v84: Removed the post-Huckel motif cleanup passes and their direct tests so aromaticity perception is driven by the shared RDKit-like donor/candidate/fused rules.
- v85: Reworked fused-system perception around a single RDKit-style connected-subset Huckel evaluator and removed the separate exocyclic fused fallback marking pass.
- v86: Add PubChem-100k as required broad-corpus validation evidence.
- v87: Narrow fused-neighbor nonaromatic bond suppression so accepted simple aromatic rings are not vetoed by adjacent nonaromatic rings.
- v88: Localize fused-system bond suppression to accepted simple rings and fused subsets, and admit exocyclic-pi chalcogen fused candidates into the subset Huckel evaluator.
- v89: Add fused-topology handling for ring-local exocyclic pi links: veto lone-pair-rescued six-member rings that RDKit keeps aliphatic, admit lone-pair five-member macrocycle partners that RDKit keeps aromatic, and allow accepted fused subsets to mark shared bonds through candidate-compatible four-electron dione partners.
- v90: Split fused support from perimeter assignment so hetero five-electron support rings fused to a large accepted member can contribute to accepted fused systems without marking their non-shared outer perimeter.
- v91: Replace localized motif gates with global RDKit-style donor
  classification and connected fused-subset marking, keep unsupported ring
  atoms as non-candidates, and add bounded valence-demand localization for
  imported aromatic components that are valid chemistry but not aromatic.
- v92: Store derived membership and model provenance only in `PerceptionState`.
- v93: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v94: Remove the parallel imported-aromatic perception engine and all
  motif-specific fused exceptions, localize imported components before one
  shared donor/candidate/Huckel pass, make the direct API transactional, widen
  electron counts, and distinguish matching limits from invalid chemistry.
- v95: Use PubChem-1k as the required baseline validation corpus after retiring the former smoke corpus from public validation.
