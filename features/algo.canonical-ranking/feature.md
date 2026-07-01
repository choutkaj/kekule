# Canonical Atom Ranking

## Summary

Assign deterministic non-stereo atom symmetry ranks for later canonical SMILES and deterministic graph workflows.

## Behavior/API

- Exposes `canon::{CanonicalAtomRanking, atom_ranking(&Molecule)}`.
- Returns one rank per live atom. Atoms in the same unresolved symmetry class share a rank.
- Uses atom identity fields, hydrogen counts, aromaticity flags, atom-map numbers, graph degree, bond order/aromaticity, and iterative neighbor refinement.
- Does not mutate the molecule or require a particular input atom order.

## Implementation Notes

- The first implementation is a deterministic Weisfeiler-Lehman-style refinement over the current graph state.
- Callers should sanitize first when rank invariants should include perceived valence, hydrogens, rings, and aromaticity.
- Symmetric ties are intentionally preserved as equal ranks. Later canonical SMILES work must add traversal/backtracking policy instead of treating rank plus `AtomId` as chemically canonical.

## Validation

- Unit tests cover symmetric atom grouping, atom-order-independent rank class counts, and payload fields that break symmetry.
- Future golden validation should compare compact RDKit-generated canonical ranking cases.

## Out Of Scope

Stereochemistry-sensitive ranking, total canonical atom ordering, canonical SMILES traversal, and RDKit parity for every tie-breaking edge case.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement non-stereo atom symmetry ranking.
- v3: Move the public ranking API under the `canon` facade.
