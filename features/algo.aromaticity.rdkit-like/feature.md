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
- Propagates bounded ring-perception failures as `AromaticityError::RingPerception`.

## Implementation Notes

- Operates on the shared core `Molecule` graph and uses ring data from the `algo.rings.sssr` stack.
- Uses one RDKit-style donor classifier for localized rings, imported aromatic-order rings, aromatic candidate checks, and fused-component electron counting.
- Supports common RDKit organic aromatic elements: C, N, O, P, S, Se, and Te.
- Applies RDKit-like candidate gates for atom degree, explicit pi-bond count, triple bonds, exocyclic multiple-bond options, charge-adjusted default valence, and radical eligibility.
- Counts localized saturated, vacant, lone-pair, anionic, and pi-bond donors with a `countAtomElec`-style helper using default valence, outer-shell electrons, charge, radical electrons, effective hydrogens, and exocyclic electronegativity.
- Evaluates localized simple rings of arbitrary size through the same Huckel donor-count path, including two-electron rings such as cyclopropenyl cation.
- Handles imported aromatic-bond rings with variable donor ranges and accepts them when the donor set contains a valid 4n+2 count.
- Treats exocyclic pi bonds through electronegativity-aware donor logic rather than raw hetero-atom symbol checks.
- Uses cached per-ring donor analysis so initial ring gates, Huckel counting, fused candidate admission, and fused single-bond protection share the same candidate state.
- Applies a bounded RDKit-style fused-system pass: fused candidate rings are grouped by shared bonds, connected subsets are evaluated from small to large, subset atom sets use fused-ring multiplicity, and accepted subsets mark perimeter bonds.
- Uses RDKit-like fused-system atom multiplicity, selected-subsystem perimeter bonds, additive accepted subsets, and the 24-atom fused-ring candidate cap.
- Keeps simple-ring nonaromatic fused-bond suppression local to simple-ring assignment, then lets accepted fused subsets decide perimeter and internal bond flags from the accepted subset topology.
- Keeps multi-protected fused-subset perimeter singles and explicitly rejected internal shared singles aliphatic, while allowing accepted fused subsets to aromatize shared bonds when member rings remain candidate-compatible and have individual Huckel or one-electron-deficient fused support.
- Distinguishes fused support rings from perimeter-marking rings, so candidate-compatible hetero five-electron rings fused to a large accepted member can support internal fused aromaticity without aromatizing their non-shared perimeter.
- Allows accepted fused subsets to aromatize internal shared bonds through candidate-compatible four-electron dione partners when the fused-system Huckel evaluator accepts the larger system.
- Allows low-unsaturation chalcogen-containing fused candidates with exocyclic pi links into the surrounding fused system to reach the fused-subset Huckel evaluator.
- Allows candidate-compatible five-member rings with a nitrogen lone-pair Huckel count and two fused-system-local exocyclic pi links to remain aromatic inside accepted macrocyclic/fused systems.
- Distinguishes terminal exocyclic pi bonds from exocyclic pi bonds that are ring-local in an adjacent fused system, so six-member rings whose Huckel count depends on a nitrogen lone pair are not admitted when the fused topology keeps the neighboring pi bond aliphatic.
- Does not run a post-Huckel molecule-specific cleanup pass. Carbonyl, imide, lactam, lactone, amidine, chalcogen-oxo, terminal-imine, and orphan-atom corrections are expected to emerge from the general donor/candidate/fused rules rather than separate motif clearing.
- Keeps parsing separate from aromaticity perception. Canonical SMILES normalization issues exposed by these flags belong to `io.smiles.canonical`, not hidden aromaticity cleanup.
- Treats unsupported or ambiguous systems conservatively rather than claiming full RDKit parity.
- Uses the source-faithful localized donor path directly when the input has no
  imported aromatic-order bonds, avoiding legacy ring-local motif gates.

## Validation

- Unit tests cover localized donor analysis, candidate gates, radical and valence eligibility, imported aromatic-order handling, fused-subsystem search, and SMILES sanitize/reparse smoke cases.
- RDKit-generated goldens compare aromatic atom and bond flags for external PubChem fixtures.

## Out Of Scope

- Full RDKit aromaticity parity.
- Runtime RDKit dependency.
- Valence perception, sanitization policy, kekulization, stereochemistry, and parser behavior.
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
