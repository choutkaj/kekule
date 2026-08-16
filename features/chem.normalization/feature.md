# Small-Molecule Representation Normalization

## Summary

Provide an explicit deterministic, idempotent layer for meaning-preserving
small-molecule representation normalization.

## Behavior/API

- Exposes `normalization::{normalize, NormalizationReport,
  NormalizationWarning, NormalizationError}` plus focused source-stereo issue
  types and the thin `SmallMolecule::normalize()` convenience.
- A standalone normalization call owns one staging clone and publishes only
  after every rewrite succeeds. Failure leaves the complete input molecule
  unchanged.
- Every successful call clears the complete installed `PerceptionState`, even
  when the represented chemistry was already canonical.
- Canonicalizes supported neutral hypervalent chlorine, bromine, and iodine
  oxo patterns into charge-separated single-bond representations when a
  terminal single-bond oxygen establishes that representation.
- Localizes every accepted imported `BondOrder::Aromatic` component into a
  deterministic ordinary single/double-bond representation. No live aromatic
  bond order remains after successful normalization.
- Resolves supported directional, Molfile wedge/hash/either, double-bond
  either, and source-declared axis marks into canonical `StereoElement` state
  after aromatic localization. Successfully represented source marks are
  consumed, so the canonical element is the single represented assertion.
- Preserves already-canonical SMILES `@`/`@@` elements without creating a
  duplicate.
- Preserves total formal charge and all chemistry outside the recognized
  representation pattern.
- Reports `NormalizationError::FormalChargeOutOfRange` rather than saturating
  when the canonical central-atom charge cannot fit the core charge type.
- Reports invalid imported aromatic representations and bounded matching
  exhaustion as focused normalization errors. Either failure preserves the
  complete original molecule.
- Reports unpaired, ambiguous, unsupported, or unassemblable source stereo as
  collected `SourceStereoNormalizationIssue` values under
  `NormalizationError::SourceStereo`. Conflicting multiple tetrahedral wedge
  marks remain the one nonfatal `NormalizationWarning`, preserving the
  established warning-and-ignore compatibility behavior.
- `NormalizationReport` returns created stereo-element IDs and concrete
  nonfatal warnings; it contains no generic pipeline output.

## Implementation Notes

- Candidate atoms and bonds are visited in stable identifier order, making the
  rewrite deterministic.
- Aromatic localization derives its matching demands only from represented
  atom/bond fields and fixed normalization rules. It never calls valence or
  aromaticity perception, reads installed implicit hydrogens, selects a named
  perception model, or installs temporary perception state.
- Aromatic matching is bounded to 100,000 deterministic search states per
  imported component.
- Source-stereo decoding visits source bonds in stable identifier order and
  derives implicit carriers only from represented source declarations such as
  `explicit_hydrogens`. It does not read installed implicit-H assignments,
  semantic aromaticity, or installed rings, call valence perception, or
  install temporary perception state.
- Source-axis and small-ring exclusions use private temporary graph-theoretical
  ring work. Coordinates are consulted only to decode the drawing-local sense
  of an explicit source wedge/hash assertion; unmarked coordinate stereo is
  never inferred by normalization.
- Normalization remains a separate public operation. Default perception
  expects its localized represented output and never invokes it implicitly.
- The combined small-molecule workflow invokes the same normalization stages
  in its one outer staging molecule instead of opening a nested transaction.

## Tests

- Focused unit tests cover the meaning-preserving oxo-halide rewrite, aromatic
  SMILES localization, fused-component idempotence, source wedge/either,
  directional, double-bond-either, and axis assembly, mark consumption,
  direct-SMILES preservation, complete perception clearing, independence from
  arbitrary installed perception, structured matching limits, and exact
  rollback on invalid or ambiguous source representation and charge overflow.
- A downstream integration test compiles both the focused facade and the
  `SmallMolecule` convenience and verifies perception clearing in a production
  library build.

## Out Of Scope

- General parser redesign, coordinate-only stereo perception publication,
  chemical standardization, perception-pipeline frameworks, and an ordinary
  parse-to-ready convenience API.

## Revision Notes

- v1: Introduce first-class transactional representation normalization with
  hypervalent oxo-halide cleanup and perception clearing.
- v2: Add model-independent represented-state aromatic localization, guarantee
  ordinary bond orders after success, and move invalid-source and matching
  limit failures out of aromaticity perception.
- v3: Canonicalize source-declared stereo transactionally after chemistry
  normalization, consume resolved marks, expose focused source diagnostics,
  and keep decoding independent of installed chemical perception.
- v4: Centralize standalone transaction ownership around one staging clone and
  let the atomic combined workflow reuse its existing outer staging molecule.
