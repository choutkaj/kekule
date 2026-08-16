# RDKit-like Valence Perception

## Summary

Provide conservative valence perception for normalized small-molecule graphs.

## Behavior/API

- Exposes `perception::valence::{ValenceModel, ValenceOptions, ValenceError,
  ValenceIssue, perceive_valence, perceive_valence_with_options}`.
- `perceive_valence` and `perceive_valence_with_options` return
  `Result<(), ValenceError>`; successful execution has no separate report.
- Requires normalized represented chemistry with localized ordinary bond
  orders and computes explicit valence from those orders plus represented
  explicit hydrogens.
- Reports every remaining `BondOrder::Aromatic` as
  `ValenceIssue::UnsupportedBondOrder` and returns before calculating or
  installing any atom assignment.
- Assigns implicit hydrogens when a common allowed valence can be selected.
- Does not read installed ring or semantic aromaticity state. Valence can run
  first, and aromatic systems derive their implicit hydrogens directly from
  the normalized localized representation.
- Returns unsupported elements or valence excesses as structured
  `ValenceIssue` values in `ValenceError` instead of silently accepting them.
- Defaults to strict behavior. When strict issues are present, perception
  leaves the molecule's complete previous `PerceptionState` unchanged.
  `ValenceOptions { strict: false }` still computes and installs assignments
  while suppressing unsupported-element and valence-excess issues for
  inspection workflows; the normalized-input preflight remains mandatory.
  Default perception uses strict mode.
- Exact installed state is publicly distinguishable as absent, model-neutral,
  or `ValenceModel::RdkitLike` and can be detached and transactionally restored
  with every atom-wise implicit-H assignment.

## Implementation Notes

- The current model uses RDKit's periodic-table allowed-valence entries,
  including the unrestricted `-1` sentinel, charge-adjusted isoelectronic
  lookup, and the P/S/As/Se hypervalent-anion adjustments used by RDKit.
- Preserves RDKit's historical acceptance of two-coordinate hydride.
- Neutral alkali and alkaline-earth atoms receive the implicit hydrogens implied
  by those allowed valences; there is no corpus-specific electropositive-atom
  suppression rule.
- Radical electrons participate in target-valence selection, and explicit
  valence/count arithmetic uses `usize` so large malformed graphs return
  structured issues instead of truncating or panicking.
- Benzene, heteroaromatics, charged aromatic atoms, and aromatic radicals use
  the same localized bond-order and charge-adjusted periodic-table path as
  other represented chemistry; there is no aromaticity-specific valence
  target.
- The allowed-valence table is confined to valence perception. Molfile
  interpretation does not call it to manufacture represented hydrogen
  carriers.
- The pass stages every implicit-hydrogen assignment and all strict issues
  before mutation. It installs semantic valence state only on success; failure
  diagnostics remain sidecar error information.
- A representation preflight collects aromatic source bond IDs in stable order
  and returns before the normal atom-valence loop, preserving the complete
  prior perception state even in permissive mode.

## Tests

- Unit tests cover neutral organics, charged species, unsupported targets,
  exceeded valence, strict/permissive behavior, and exact preservation of the
  complete prior perception state on strict failure.
- Focused normalized-input regressions cover benzene, pyridine, represented
  pyrrolic `[nH]`, furan, thiophene, cationic and anionic aromatics, an aromatic
  carbon radical, fused naphthalene, and an unusual dye fixture. They run
  valence before rings/aromaticity and then verify the expected downstream
  semantic aromatic state.
- Raw interpreted benzene verifies that every aromatic source bond is rejected
  without implicit-H publication, that a nonempty previous `PerceptionState`
  is preserved exactly, and that normalized benzene still receives one
  implicit hydrogen per carbon.

## Benchmarks

- RDKit-generated goldens compare valence status, explicit valence, and implicit hydrogen assignments for external PubChem fixtures.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Imported-aromatic localization, query atoms, bond-order-dependent
  organometallic interpretation, valence tautomer handling, and default
  perception orchestration.

## Revision Notes

- v1: Conservative valence perception.
- v2: Benchmark contract narrowed to valence-specific outputs and matched the RDKit-backed `smoke` corpus.
- v3: Add corpus-driven RDKit-compatible valence cases for charged halides, boron anions, alkali counterions, hypervalent halogens, and simple mercury salts.
- v4: Expand corpus-driven RDKit-compatible valence cases for PubChem-100 salts, silicon, phosphonium, and selected metal centers.
- v5: Generalize pubchem-1k-driven valence handling for transition-metal coordination, group-14/group-15 heavy elements, oxonium centers, chalcogens, and radicals; pubchem-1k still requires further table coverage.
- v6: Add aromatic imported-SMILES valence targets so lowercase aromatic systems sanitize with RDKit-like hydrogen counts.
- v7: Move the public expert API under the `perception::valence` facade.
- v8: Add PubChem-100k as required broad-corpus external-parity evidence.
- v9: Expand RDKit-like simple-ion and main-group valence support for PubChem salts while leaving actinide and coordination-heavy cases as structured unsupported chemistry.
- v10: Allow isolated unsupported atoms as zero-valence spectators so disconnected PubChem salt fragments do not block sanitization of descriptor-bearing organic components.
- v11: Add strict/permissive options, charge-adjusted isoelectronic valence
  lookup, unrestricted-valence elements, implicit-hydrogen suppression for
  electropositive centers, and reuse the same allowed-valence table for
  Molfile tetrahedral hydrogen-carrier preservation.
- v13: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v14: Replace corpus-era charged-element exceptions with the exact RDKit
  fixed/unrestricted valence table and isoelectronic rules, restore RDKit-like
  hydrogens for neutral electropositive atoms, include radical electrons, and
  widen valence accounting to graph-sized integers.
- v15: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v16: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v17: Expose lossless installed-valence section inspection and canonical
  reconstruction, including the distinct model-neutral state.
- v18: Replace the empty-success report with `Result<(), ValenceError>` and
  install semantic valence state transactionally only after strict validation
  succeeds.
- v19: Require normalized localized represented chemistry, remove implicit-H
  dependence on installed semantic aromaticity, and verify valence-first
  aromatic workflows across carbon, heteroatom, charged, radical, fused, and
  unusual fixtures.
- v20: Enforce the normalized-input precondition by collecting every remaining
  aromatic source bond as a structured valence issue before atom assignments,
  with exact transactional rollback in strict and permissive modes.
