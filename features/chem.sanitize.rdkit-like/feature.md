# RDKit-like Sanitization Pipeline

## Summary

Provide an explicit opt-in transactional normalization + perception workflow
for common small molecules.

## Behavior/API

- Exposes `perception::{SanitizeOptions, SanitizeReport, SanitizeError, sanitize, sanitize_with_options, sanitize_with_ring_options}`.
- Normalizes represented chemistry, then runs valence, ring set, aromaticity,
  and stereo perception according to options.
- Normalization assembles supported source marks into canonical local stereo
  and consumes them before valence perception. The stereo-perception stage no
  longer decodes source marks; coordinate-only stereo assignment remains an
  explicit `stereo.perception` operation.
- Commits changes only after every requested pass succeeds; any error leaves the input exactly unchanged.
- Propagates valence, ring, aromaticity, and stereo failures through
  `SanitizeError` without committing staged mutations.
- Propagates canonical representation failures through
  `SanitizeError::Normalization` without committing staged mutations.
- `SanitizeError::Stereo` carries coordinate-perception validation or insertion
  failures. Source-mark diagnostics propagate through
  `SanitizeError::Normalization`.
- `SanitizeError::Valence` carries `ValenceError`. Successful
  `SanitizeReport` values retain the useful `NormalizationReport` plus optional
  stereo-stage output; installed valence and ring state are
  inspected through the molecule rather than duplicated in the report.
- Installs requested successful perception state and clears skipped state.
  Aromaticity reuses an already-installed ring basis. It may compute rings
  internally when no basis is installed, but an unrequested ring result is not
  retained or exposed.
- Performs no representation cleanup after perception. Represented explicit
  hydrogens remain primary state, while valence-derived implicit hydrogens stay
  in `PerceptionState`.
- Does not run automatically from file parsers.

## Implementation Notes

- The workflow stages work on a clone and atomically replace the caller's
  molecule after success.
- The workflow returns normalization sidecar information because created
  source-stereo elements and concrete warnings remain useful to callers.
- It operates on `SmallMolecule` while using shared core graph algorithms internally.
- The public facade is `perception`; lower-level sanitizer internals are not root-level API.
- Delegates canonical representation cleanup to the focused normalization
  layer before valence perception. The sanitizer owns neither the hypervalent
  oxyhalogen rewrite nor imported aromatic localization, and it does not
  restore perception cleared by normalization.
- Does not translate perceived aromatic-nitrogen hydrogens into represented
  explicit-H declarations. Source `[nH]` remains represented from
  interpretation/normalization; other inferred donor hydrogens remain derived.
- Runs stereo perception after aromaticity for compatibility, but source-mark
  assembly has already completed without installed perception during
  normalization.
- Preserves representable Molfile double-bond either marks as explicit unknown
  stereo elements instead of rejecting the whole molecule.
- Retains conflicting multi-wedge input as a normalization warning while
  allowing otherwise valid chemistry to sanitize; lone unassemblable marks and
  structural stereo errors remain fatal and transactional normalization
  failures.
- Retains valence-implied explicit hydrogen carriers established by Molfile
  wedge parsing so tetrahedral local stereo survives the full pipeline.
- Its valence, ring, aromaticity, and stereo passes are compared together against each required corpus.
- Inherits the current valence and aromaticity improvements, including radical implicit-hydrogen handling, imported aromatic SMILES handling, and conservative unsupported-ring behavior.
- Propagates invalid imported-aromatic representation and localization-budget
  failures through `SanitizeError::Normalization` while retaining
  whole-pipeline rollback.

## Tests

- Unit tests cover parse-without-sanitize behavior, every option combination,
  installed ring-basis reuse, transient aromaticity ring computation,
  source-stereo normalization regardless of the stereo-perception option,
  idempotence, coordinate-only stereo staying outside sanitization, and exact
  rollback after normalization, valence, aromaticity, or stereo failure.
- A complete represented-state snapshot regression covers atoms, bonds,
  stable topology layout, stereo and source history, properties, and
  conformers across valence -> rings -> aromaticity, excluding only installed
  perception and its test-only mirrors.

## Benchmarks

- RDKit-generated goldens compare sanitized atom state for external PubChem fixtures.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Full RDKit sanitization parity, exact CIP assignment, coordinate-derived
  stereo assignment, cleanup transforms, and organometallic handling.

## Revision Notes

- v1: Explicit sanitization pipeline.
- v2: Validated through the corrected valence, ring, and aromaticity passes.
- v3: Add RDKit-like oxyhalogen cleanup and pass PubChem-100 through the corrected valence and aromaticity stack.
- v4: Incorporate broader pubchem-1k-driven valence and aromaticity behavior; pubchem-1k remains pending on fused aromatic bond selection and remaining valence-table coverage.
- v5: Make sanitization transactional and define fresh/stale state outcomes for every option combination.
- v6: Sanitize imported aromatic SMILES with corrected aromatic valence and atom-contribution aromaticity behavior.
- v7: Accept explicit ring-work limits and preserve transactional rollback on ring resource errors.
- v8: Move the public small-molecule sanitizer API under the `perception` facade.
- v9: Add PubChem-100k as required broad-corpus external-parity evidence.
- v10: Normalize pyrrolic aromatic nitrogen donor hydrogens to RDKit-style sanitized `nH` atom state.
- v11: Add stereo perception as an explicit sanitizer stage with options,
  reporting, freshness-state handling, and transactional rollback on stereo
  perception issues.
- v12: Accept representable Molfile double-bond either marks by assembling
  explicit unknown double-bond stereo, align conflicting multi-wedge handling
  with RDKit's non-fatal warning semantics, and gate oxyhalogen charge cleanup
  on a general terminal-oxygen valence pattern.
- v13: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v14: Inherit the exact RDKit valence-table cleanup and unified transactional
  aromaticity engine, including structured imported-aromatic matching limits.
- v15: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v16: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v17: Propagate transactional `ValenceError` failures and remove the redundant
  empty-success valence field from `SanitizeReport`.
- v18: Propagate structured transactional stereo perception errors directly,
  retain successful created-element IDs and warnings, and remove fatal-issue
  inspection of successful reports.
- v19: Define sanitization as transactional normalization plus perception,
  reuse its installed ring basis during aromaticity, and remove redundant ring
  count output from `SanitizeReport`.
- v20: Delegate representation cleanup to the first-class normalization layer
  and propagate focused normalization failures transactionally.
- v21: Delegate imported aromatic localization to normalization, remove the
  transitional perception restore, and call purely perceptual aromaticity on
  localized ordinary bond orders.
- v22: Delegate source-declared stereo assembly and diagnostics to
  normalization, expose normalization sidecar output, and leave the sanitizer's
  stereo step responsible only for true perception/validation work.
- v23: Remove post-aromaticity nitrogen/H representation feedback and enforce
  the one-way represented-molecule-to-perception-state boundary across the
  default discrete chemical passes.
