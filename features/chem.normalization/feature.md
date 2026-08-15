# Small-Molecule Representation Normalization

## Summary

Provide an explicit deterministic, idempotent layer for meaning-preserving
small-molecule representation normalization.

## Behavior/API

- Exposes `normalization::{normalize, NormalizationError}` and the thin
  `SmallMolecule::normalize()` convenience.
- Normalization stages work on a clone and publish only after every rewrite
  succeeds. Failure leaves the complete input molecule unchanged.
- Every successful call clears the complete installed `PerceptionState`, even
  when the represented chemistry was already canonical.
- Canonicalizes supported neutral hypervalent chlorine, bromine, and iodine
  oxo patterns into charge-separated single-bond representations when a
  terminal single-bond oxygen establishes that representation.
- Localizes every accepted imported `BondOrder::Aromatic` component into a
  deterministic ordinary single/double-bond representation. No live aromatic
  bond order remains after successful normalization.
- Preserves total formal charge and all chemistry outside the recognized
  representation pattern.
- Reports `NormalizationError::FormalChargeOutOfRange` rather than saturating
  when the canonical central-atom charge cannot fit the core charge type.
- Reports invalid imported aromatic representations and bounded matching
  exhaustion as focused normalization errors. Either failure preserves the
  complete original molecule.

## Implementation Notes

- Candidate atoms and bonds are visited in stable identifier order, making the
  rewrite deterministic.
- Aromatic localization derives its matching demands only from represented
  atom/bond fields and fixed normalization rules. It never calls valence or
  aromaticity perception, reads installed implicit hydrogens, selects a named
  perception model, or installs temporary perception state.
- Aromatic matching is bounded to 100,000 deterministic search states per
  imported component.
- The rewrite no longer belongs to sanitization. The compatibility sanitizer
  delegates to normalization before running its existing perception stages.
- No success report exists because Stage 1 produces no warnings or other
  useful sidecar output.

## Tests

- Focused unit tests cover the meaning-preserving oxo-halide rewrite, aromatic
  SMILES localization, fused-component idempotence, complete perception
  clearing, structured matching limits, and exact rollback on invalid source
  representation or formal-charge overflow.
- A downstream integration test compiles both the focused facade and the
  `SmallMolecule` convenience and verifies perception clearing in a production
  library build.

## Out Of Scope

- Source-stereo assembly, general parser redesign, chemical standardization,
  perception-pipeline frameworks, and removal of the compatibility sanitizer
  APIs.

## Revision Notes

- v1: Introduce first-class transactional representation normalization with
  hypervalent oxo-halide cleanup and perception clearing.
- v2: Add model-independent represented-state aromatic localization, guarantee
  ordinary bond orders after success, and move invalid-source and matching
  limit failures out of aromaticity perception.
