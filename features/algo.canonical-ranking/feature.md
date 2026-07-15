# Canonical Atom Ranking

## Summary

Assign deterministic non-stereo atom symmetry ranks for later canonical SMILES and deterministic graph workflows.

## Behavior/API

- Exposes `canon::{CanonicalAtomRanking, atom_ranking(&Molecule)}`.
- Returns one rank per live atom. Atoms in the same unresolved symmetry class share a rank.
- Uses RDKit's non-stereo base invariants: atom map, graph degree, atomic
  number, isotope, total hydrogen count, and formal charge, followed by bond
  order/aromaticity and iterative neighbor refinement.
- Treats equivalent total hydrogen counts identically whether hydrogens entered
  through explicit atom metadata or valence perception; parser representation
  details do not create false chemical asymmetry.
- Treats a perceived aromatic bond independently of its localized Kekule
  single/double choice.
- Adds rooted cyclic-topology refinement for regular fused and bridged graphs
  that ordinary local color refinement cannot separate.
- Does not mutate the molecule or require a particular input atom order.

## Implementation Notes

- Uses deterministic Weisfeiler-Lehman-style local refinement followed by an
  RDKit-like ring-topology pass. Breadth-first level widths encode cyclic
  distances and revisit multiplicities encode path reconvergence.
- Callers should sanitize first when rank invariants should include perceived valence, hydrogens, rings, and aromaticity.
- Non-stereo ranking and its reference validation do not require the stereo
  perception stage to accept the input molecule.
- Symmetric ties are intentionally preserved as equal ranks. Later canonical SMILES work must add traversal/backtracking policy instead of treating rank plus `AtomId` as chemically canonical.

## Validation

- Unit tests cover symmetric atom grouping, atom-order-independent rank class
  counts, payload fields that break symmetry, Kekule-choice independence, and
  regular cyclic graphs that require global topology refinement.
- RDKit-generated goldens compare equivalence-class partitions rather than
  implementation-specific numeric rank labels across the external molecular
  corpora.

## Out Of Scope

Stereochemistry-sensitive ranking, total canonical atom ordering, canonical SMILES traversal, and RDKit parity for every tie-breaking edge case.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement non-stereo atom symmetry ranking.
- v3: Move the public ranking API under the `canon` facade.
- v4: Give quartet and quintet radicals distinct atom invariants rather than
  collapsing newly representable high-spin states.
- v5: Add RDKit-backed symmetry-partition validation, normalize perceived
  aromatic bond codes across Kekule forms, align base invariants with RDKit's
  total-H semantics, preserve graph-sized degrees, and add conditionally
  triggered rooted cyclic-topology refinement for WL-hard ring systems.
