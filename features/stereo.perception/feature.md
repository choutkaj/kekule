# Stereo Inference and Materialization

## Summary

Validate graph-local represented stereo, detect candidate stereochemical units,
infer coordinate-derived stereo read-only, and materialize it only through an
explicit representation transform.

## Behavior/API

- Exposes focused operations under top-level `stereo` for
  `validate_stereo`, `detect_stereo_candidates`, read-only
  `infer_coordinate_stereo*`, and explicit
  `materialize_coordinate_stereo*`.
- `validate_stereo(&Molecule) -> Result<(), StereoValidationError>` validates
  represented `StereoElement` structure only: live focus references, carrier
  counts and uniqueness, adjacency, focus/end-point coherence, bond order, and
  carrier forms supported by each represented stereo kind. It does not require
  valence, hydrogen, ring, or aromaticity perception.
- Implicit hydrogen and lone-pair carriers are structurally valid tetrahedral
  forms, and implicit hydrogen is a structurally valid double-bond carrier.
  Whether such a carrier is chemically available or stereogenic is deferred to
  perception or CIP. Double-bond lone-pair and implicit axis carriers remain
  structurally unsupported.
- `detect_stereo_candidates(&Molecule) -> Vec<StereoCandidate>` remains a
  read-only exploratory query over the current graph and installed hydrogen or
  ring state where its conservative chemical heuristics require them.
- Coordinate inference and materialization never read, assemble, consume, warn
  about, or error on `StereoBondMark` state. Source-declared stereo belongs to
  normalization.
- `infer_coordinate_stereo*(&Molecule)` returns detached proposed
  `StereoElement` values and leaves the complete molecule, including
  `PerceptionState`, unchanged.
- `CoordinateStereoOptions::infer_axes` enables the default-off conservative
  3D axis subset. Tetrahedral and double-bond inference are always considered.
- `materialize_coordinate_stereo*(&mut Molecule)` is the separately named
  transactional transform. It stages the inferred elements on a clone,
  skips already represented foci, validates the complete final stereo
  representation, publishes only on success, and reports created element IDs.
- `CoordinateStereoError` preserves structural-validation and insertion
  failures. Source-mark warning and issue vocabulary remains absent.
- Coordinate-derived tetrahedral assignment requires four explicit atom
  carriers with nondegenerate 3D coordinates. Double-bond assignment requires
  one explicit atom carrier on each side with nondegenerate 2D or 3D geometry.
- Opt-in coordinate-axis assignment requires a single-bond axis with two
  SP2-like endpoints, exactly two explicit atom carriers per endpoint, no
  existing stored-axis element, and nondegenerate 3D handedness from the
  lowest-ID endpoint references.
- Conservative double-bond candidate exclusions retain the RDKit-like
  small-ring boundary and current unsupported aromatic/endocyclic hetero cases.

## Implementation Notes

Stored stereo validation is representation-structural and independent of
`PerceptionState`. Candidate detection and coordinate inference retain their
existing conservative chemistry heuristics because they are true perception
work. Inference returns a snapshot result; only explicit materialization
changes represented stereo.

Default perception does not call coordinate inference or materialization.
Source marks and their warnings or errors remain owned by `normalization`.


## Tests

- Focused unit and regression tests cover perception-independent structural
  validation, standalone candidate detection, completely read-only
  tetrahedral/double-bond/axis inference, coordinate-axis opt-in behavior,
  explicit transactional materialization, duplicate avoidance, orientation
  and source preservation, and exact rollback.
- Normalization regressions own wedge/hash/either, directional, source-axis,
  implicit source carrier, source-warning, and source-error coverage.

## Benchmarks

- PubChem, Enamine, and PL-REX semantic fixtures continue to record represented
  stereo plus candidate and explicitly materialized coordinate-stereo output
  after explicit normalization and default perception.
  Source-declared elements are now present in the normalization report and
  canonical molecule before coordinate inference runs.
- Exact descriptor parity remains the responsibility of `stereo.cip`.
- Optional external-reference manifests are available for `pubchem-1k`,
  `pubchem-100k`, `pl-rex`, and `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's
  release status or repository health.


## Out Of Scope
Exact CIP descriptors, default broad axis candidate inference, non-opt-in
coordinate-only axis inference, 2D coordinate-only axis inference,
source-mark normalization, CXSMILES atropisomeric syntax, isomeric SMILES
writing, enhanced stereo serialization, implicit-hydrogen coordinate
reconstruction, stereo enumeration, and reaction stereo transfer.

## Revision Notes

- v1: Feature contract reserved for stereo candidate detection and local
  validation.
- v2: Add public stereo perception API, local element validation, conservative
  tetrahedral/double-bond candidate detection, paired directional source-mark
  assembly, unit coverage, and smoke semantic validation.
- v3: Assemble supported Molfile wedge up/down source marks into specified
  tetrahedral elements and wedge/either source marks into explicit unknown
  tetrahedral elements.
- v4: Integrate stereo perception into the explicit small-molecule sanitization
  pipeline with opt-out options, reporting, freshness-state handling, and
  transactional rollback on stereo issues.
- v5: Add conservative coordinate-derived local assignment for explicit-atom
  tetrahedral centers and double bonds using the first conformer.
- v6: Normalize SMILES directional bond marks relative to alkene endpoints and
  accept redundant two-mark substituted endpoints when the marks cover both atom
  carriers with opposite normalized directions.
- v7: Validate supported implicit lone-pair tetrahedral carriers and skip
  unsupported aromatic or endocyclic hetero double-bond stereo candidates before
  source-mark assembly.
- v8: Add PubChem 100 and PubChem 1k semantic regression requirements for
  stereo perception over externally supplied isomeric SMILES.
- v9: Add PubChem 100k and Enamine diversity semantic regression requirements
  and preserve unsupported-sanitization records as per-record validation output
  instead of aborting broad perception fixtures.
- v10: Exclude double-bond stereo candidates in rings smaller than eight atoms
  using the RDKit-like stereogenic-bond boundary while preserving cyclooctene
  and larger ring alkene candidates.
- v11: Derive Molfile wedge up/down tetrahedral orientation from conformer
  coordinates when present so coordinate-bearing V2000 records preserve
  RDKit-like local stereo sense.
- v12: Validate stored axis elements structurally instead of treating every
  axis element as unsupported, enabling the CIP layer to assign descriptors for
  explicitly stored axes.
- v13: Assemble a conservative Molfile atropisomeric wedge subset into stored
  axis elements using ring-atom axis eligibility and 2D opposite-side carrier
  selection.
- v14: Restrict Molfile wedge-derived atrop axis candidates to non-ring axis
  bonds, matching the exocyclic subset covered by official RDKit RP-6306
  variants and avoiding ambiguous ring-internal candidates.
- v15: Use virtual implicit-H geometry for coordinate-bearing Molfile
  tetrahedral wedges and store Molfile atrop axes with RDKit-style
  lowest-neighbor endpoint references plus coordinate-derived handedness.
- v16: Broaden Molfile atrop axis assembly to RDKit-like SP2 endpoint
  eligibility, consume redundant same-axis wedge marks, and cover official
  BMS/ZM atropisomer regressions.
- v17: Allow ring-internal macrocyclic Molfile atrop axes when no non-ring
  candidate is available from the same source mark, with official RDKit
  macrocycle regressions.
- v18: Add PL-REX ligand SDF packs to the perception benchmark contract for
  coordinate-bearing Molfile stereo and source-mark assembly regression
  coverage.
- v19: Add default-off conservative 3D coordinate-derived axis assignment for
  explicit SP2-like single-bond axes with two atom carriers per endpoint and
  lowest-neighbor endpoint references.
- v20: Assemble Molfile double-bond either marks into explicit unknown
  double-bond stereo elements when both endpoints have valid carriers, and
  treat consumed conflicting multi-wedge input as a non-fatal sanitization
  ambiguity rather than rejecting valid chemistry. Preserve explicit
  hypervalent and pyramidal heteroatom wedge centers without broadening
  unmarked candidate detection.
- v22: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v23: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v24: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v25: Split stored stereo validation, candidate detection, and mutating
  perception into focused operations; make direct perception transactional;
  separate successful warnings from fatal structured errors; and narrow
  perception options and successful report output. Use carrier-accurate public
  validation issues for unavailable implicit and unsupported axis carriers.
- v26: Move all source-mark assembly and diagnostics into normalization,
  remove source warnings from perception output, and narrow validation to
  represented structural invariants independent of installed perception.
- v27: Replace mutating coordinate-stereo perception with detached read-only
  inference results and a separately named transactional materialization
  transform.
- v28: Move the complete stereo-domain public surface to the top-level
  `stereo` facade and remove the mixed-responsibility `perception::stereo`
  namespace.
