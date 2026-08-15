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
- Preserves total formal charge and all chemistry outside the recognized
  representation pattern.
- Reports `NormalizationError::FormalChargeOutOfRange` rather than saturating
  when the canonical central-atom charge cannot fit the core charge type.

## Implementation Notes

- Candidate atoms and bonds are visited in stable identifier order, making the
  rewrite deterministic.
- The rewrite no longer belongs to sanitization. The compatibility sanitizer
  delegates to normalization before running its existing perception stages.
- No success report exists because Stage 1 produces no warnings or other
  useful sidecar output.

## Tests

- Focused unit tests cover the meaning-preserving oxo-halide rewrite,
  idempotence, complete perception clearing, and exact rollback on formal-charge
  overflow.
- A downstream integration test compiles both the focused facade and the
  `SmallMolecule` convenience and verifies perception clearing in a production
  library build.

## Out Of Scope

- Imported aromatic localization, source-stereo assembly, parser or
  interpretation changes, chemical standardization, perception redesign, and
  removal of the compatibility sanitizer APIs.

## Revision Notes

- v1: Introduce first-class transactional representation normalization with
  hypervalent oxo-halide cleanup and perception clearing.
