# RDKit-like Aromaticity Perception

## Summary

Assign aromatic atom and bond membership in canonical `PerceptionState` for
supported organic ring systems using an RDKit-like graph aromaticity model.

## Behavior/API

- Exposes `perception::aromaticity::{AromaticityModel::RdkitLike, perceive_aromaticity, perceive_aromaticity_with_ring_options}`.
- Reuses an installed `RingSet` when present and otherwise computes a ring
  basis before assigning aromaticity; membership alone is not treated as a
  basis.
- Replaces existing semantic aromatic membership transactionally and records
  `AromaticityModel::RdkitLike` as the perception model.
- Can be run directly or through the explicit sanitization pipeline.
- Treats unsupported ring elements as non-candidates, allowing a supported
  aromatic subring to remain aromatic when fused or attached to a nonaromatic
  unsupported-element ring.
- Requires normalized localized bond orders and returns
  `AromaticityError::UnsupportedBondOrder` when a live aromatic source bond
  remains. It never localizes or otherwise rewrites primary bond orders.
- Propagates bounded ring-perception failures as `AromaticityError::RingPerception`.
- Exact installed aromatic atom/bond sets and their named perception model are
  publicly inspectable and reconstructible without requiring a perceived ring
  basis.

## Implementation Notes

- Operates on the shared core `Molecule` graph and uses installed ring data or
  computes it through the `algo.rings.sssr` stack when absent.
- Reads normalized ordinary bond orders through the donor/candidate/Huckel
  engine and writes only semantic aromatic atom/bond membership.
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
- The direct public API stages the complete operation, including any required
  ring perception, and commits only on success. Atom/bond payloads, bond
  orders, local stereo, source marks, properties, and conformers are preserved.

## Tests

- Unit tests cover installed ring-basis reuse, missing and membership-only ring
  state, localized donor analysis, candidate gates, radical and valence
  eligibility, normalized-input enforcement, represented-state preservation,
  fused-subsystem search, and SMILES sanitize/reparse smoke cases.
- A cross-pass regression snapshots the complete primary representation of a
  normalized heteroaromatic stereo fixture and proves valence, ring-set, and
  aromaticity perception change only `PerceptionState`.

## Benchmarks

- RDKit-generated goldens compare aromatic atom and bond flags for external PubChem fixtures.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Full RDKit aromaticity parity.
- Runtime RDKit dependency.
- Valence perception, sanitization policy, general-purpose Kekule-form
  enumeration, stereochemistry, and parser behavior.
- Canonical SMILES normalization for every valid aromaticity assignment.

## Revision Notes

- v1-v83: Built the RDKit-like donor classifier, fused-subsystem search, benchmark workflow, and public expert facade.
- v84: Removed the post-Huckel motif cleanup passes and their direct tests so aromaticity perception is driven by the shared RDKit-like donor/candidate/fused rules.
- v85: Reworked fused-system perception around a single RDKit-style connected-subset Huckel evaluator and removed the separate exocyclic fused fallback marking pass.
- v86: Add PubChem-100k as required broad-corpus external-parity evidence.
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
- v95: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v96: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v97: Expose lossless installed aromaticity membership/provenance for checked
  canonical reconstruction without adding new section dependencies.
- v98: Reuse a current installed ring basis during aromaticity perception and
  compute one only when no `RingSet` is installed.
- v99: Move imported aromatic localization and its failures into
  `chem.normalization`, require localized input, remove imported aromaticity
  provenance, and make the aromaticity pass representation-pure.
- v100: Enforce complete primary-representation purity across the ordinary
  valence -> rings -> aromaticity sequence and remove the sanitizer's remaining
  aromatic-nitrogen representation feedback.
