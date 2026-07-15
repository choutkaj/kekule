# Hydrogen Topology Normalization

## Summary

Materialize stored and perceived hydrogens as graph atoms, or conservatively
collapse ordinary graph hydrogens back into explicit counts or perceived
implicit hydrogens.

## Behavior/API

- Exposes `hydrogens::{add_hydrogens, add_hydrogens_with_options,
  remove_hydrogens}` plus matching `SmallMolecule` convenience methods.
- Both transforms are transactional: failure leaves the input unchanged.
- `add_hydrogens` materializes bracket-style explicit counts and current
  RDKit-like implicit-hydrogen assignments as appended neutral protium atoms
  with single bonds. `explicit_only` limits materialization to stored counts.
- `remove_hydrogens` collapses neutral, ordinary, degree-one graph hydrogens
  while preserving the total hydrogen count at each parent. It uses implicit
  hydrogens where current valence rules reproduce the count and explicit counts
  where they do not, including aromatic bracket-hydrogen cases.
- Retained atom IDs remain unchanged. Reports map every added or removed
  hydrogen to its parent and explain every retained hydrogen.
- Tetrahedral and double-bond hydrogen carriers are converted between explicit
  atom and implicit-hydrogen carriers without changing local orientation.
- Added hydrogen coordinates are missing in every existing conformer; the
  transform does not invent 2D or 3D geometry.
- Both operations require current valence perception whenever implicit counts
  are consumed or reconstructed. They do not sanitize implicitly.

## Implementation Notes

- Work is staged on a clone and committed only after topology, stereo, and
  hydrogen-count postconditions succeed.
- Addition is linear in the input graph plus added topology and is guarded by a
  configurable default limit of 1,000,000 new hydrogen atoms.
- Removal is conservative and lossless. Isotopic, mapped, charged, radical,
  property-bearing, isolated, multiply connected, hydrogen-only, non-single-
  bonded, and unsupported stereo-role hydrogens are retained. Source-marked
  hydrogen bonds are removed only when an installed semantic stereo element
  supplies a lossless explicit-to-implicit carrier conversion.
- Retained topology uses stable IDs. Removed IDs become tombstones through the
  ordinary graph deletion path, which also clears their conformer positions.
- Topology changes invalidate installed perception; callers may explicitly
  sanitize or perceive the result afterward.

## Validation

- Unit and downstream API tests cover implicit and explicit-count addition,
  resource-limit and missing-perception transactionality, methane and aromatic
  bracket-hydrogen collapse, stereo carrier preservation, conformer behavior,
  and conservative retention of lossy hydrogen cases.
- RDKit-generated smoke, PubChem-100, PubChem-1k, and PL-REX goldens compare
  per-parent added hydrogen counts, retained graph topology, atom identity,
  charges, isotopes, maps, stable original-index connectivity, and total
  encoded hydrogen counts.
  Explicit-versus-implicit count storage and unrelated perception caches are
  deliberately normalized away.

## Out Of Scope

- Hydrogen coordinate generation, protonation-state changes, pH models,
  isotopic-hydrogen tracking across collapse, query atoms, macromolecule
  hierarchy updates, removal of charged or otherwise lossy hydrogen atoms, and
  an unsafe remove-all-hydrogens mode.

## Revision Notes

- v1: Add bounded transactional hydrogen materialization and conservative
  collapse with mappings, stereo preservation, and explicit coordinate policy.
