# Default Discrete Perception Pipeline

## Summary

Provide one explicit, opinionated operation for installing Kekule's default
discrete small-molecule perception state after representation normalization.

## Behavior/API

- Exposes `perception::perceive(&mut Molecule) -> Result<(), PerceptionError>`
  and a thin `SmallMolecule::perceive()` convenience.
- Runs exactly RDKit-like valence and implicit-H perception, default ring-set
  perception, then RDKit-like aromaticity perception.
- Requires normalized represented chemistry. Remaining
  `BondOrder::Aromatic` source bonds fail through `PerceptionError::Valence`;
  the operation never normalizes implicitly.
- Stages the complete operation on a clone. Any valence, ring, or aromaticity
  failure preserves the exact input molecule and prior `PerceptionState`.
- On success, installs complete fresh valence, ring-set, and aromaticity state
  while leaving all primary represented chemistry unchanged.
- Does not assemble source stereo, perceive coordinate stereo, or assign CIP
  descriptors.
- Has no options or success report. Expert partial/custom workflows use the
  focused `perception::valence`, `perception::rings`, and
  `perception::aromaticity` operations.

## Implementation Notes

- `PerceptionError` is a focused enum wrapping `ValenceError`,
  `RingPerceptionError`, or `AromaticityError`.
- Normalization and perception remain separate public operations. The ordinary
  one-call parse-to-ready convenience remains intentionally deferred.

## Tests

- Focused tests prove successful default installation, exact primary-state
  purity, idempotence, and aromatic/heteroaromatic chemistry behavior.
- Failure tests prove exact rollback for unnormalized aromatic input and for a
  downstream ring resource failure after staged valence perception succeeds.
- Public API and reconstruction tests exercise the explicit
  `normalize()` then `perceive()` workflow.

## Benchmarks

- Existing RDKit-generated atom-state comparisons are retained under this
  feature identity for the optional external molecular corpora.
- Benchmark observations remain informational and never determine repository
  health or release status.

## Out Of Scope

- Normalization, configurable partial perception, source/coordinate stereo,
  CIP assignment, parser redesign, and an ordinary fast-forward constructor.

## Revision Notes

- v1: Replace the former sanitizer abstraction with one transactional default
  discrete perception operation over normalized chemistry.
