# Reference Validation Harness

## Summary

Provide repeatable infrastructure for comparing Rust behavior against generated
or manually curated external-reference golden data.

## Behavior/API

- Exposes `cargo xtask validate --feature FEATURE_ID|all --corpus CORPUS_ID|all [--update] [--jobs N]`.
- Defaults omitted `--corpus` to `smoke`.
- Discovers corpus manifests under `validation/corpora/<corpus-id>/features/<feature-id>.toml`.
- Verifies listed fixture paths exist.
- Requires one deterministic gzip golden under the corpus `golden/<feature-id>/` directory for each listed fixture.
- Compares normalized Rust implementation output against each golden file's `expected` payload.
- Compares fixtures in parallel by default using all available processors; `--jobs N` bounds validation to a fixed worker count.
- Prints compact progress output for the overall target set and for fixture comparison within each feature/corpus target.
- Accepts only declared implementation-vs-golden comparison manifests for required implemented validation.
- Verifies golden `feature_id`, `corpus_id`, `fixture_path`, current fixture SHA-256, and reference tool/version metadata before comparing payloads.
- Records content-addressed pass evidence over manifests, source locks, fixtures, goldens, Rust source, Cargo manifests, feature metadata, and reference generator/environment files when the reference tool is generator-backed.
- Supports manually curated semantic references when the manifest
  `reference_tool` ends in `-manual-semantic`; these are evidence-backed by
  pinned source fixtures, manifest metadata, committed goldens, and Rust
  implementation sources rather than by local reference generator files.
- Canonicalizes UTF-8 text inputs to LF before hashing so evidence is stable across Windows and Linux checkouts; binary inputs remain byte-exact.
- Normalizes representation-only graph differences such as undirected bond endpoint orientation, bond array order, and ring atom order before comparison.
- Treats non-applicable feature/corpus combinations as skips and missing required manifests as errors.
- Exposes `cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]`.
- Keeps ordinary validation read-only; `--update` clears selected stale passes before running, records evidence for successful selected targets, records fixture-level failure summaries for failed comparisons, synchronizes overall `validated`, and regenerates the dashboard.
- Supports corpus pack member checks using PubChem defaults or corpus-declared SDF member properties and SMILES title prefixes.
- Preserves all records in SDF pack fixtures for stereo implementation comparisons.
- Emits structured semantic JSON for stored axis stereo validation issues and
  includes assigned axis `M`/`P` descriptors in bond descriptor maps.

## Implementation Notes

- RDKit reference generators live under `validation/reference/rdkit/`.
- Biopython reference generators live under `validation/reference/biopython/`.
- `*-manual-semantic` manifests use externally pinned fixtures and manually
  reviewed semantic goldens; they do not require local reference generator
  files. Current semantic reference labels include `pubchem-manual-semantic`
  and `enamine-manual-semantic`.
- Golden data should be normalized JSON and include reference tool versions.
- Corpus descriptors and feature manifests use typed TOML; source selection and checksums live in `sources.lock.json`.
- Source pack records may declare `member_id_property` for SDF packs or `member_title_prefix` for SMILES packs when the corpus does not use PubChem CID metadata.
- Status evidence records fixture and comparison counts, reference versions, the manifest SHA-256, a versioned evidence input list, evidence SHA-256, and validation time.
- Evidence is considered current only when recomputing it from the current checkout produces the stored schema version and hash.
- Failed comparison status records fixture count, successful comparison count, failed count, and the first fixture-level failure without recording pass evidence.
- Evidence schema v2 includes cross-platform text line-ending normalization.
- Repeated `--update` runs preserve timestamps when the evidence hash is unchanged.
- The validation command uses the Rust implementation only; RDKit, Biopython, or manually pinned external sources are used to generate or curate goldens, not to run validation.
- Fixture comparison uses a bounded worker pool while status writes and dashboard regeneration remain single-threaded.
- Progress output uses plain ASCII bars and throttled checkpoint updates so it stays readable in terminals and captured logs.
- Reference tools are never Rust runtime dependencies.

## Validation

- Current coverage is infrastructure unit-test based plus live corpus comparisons against committed external-source goldens.
- Passing comparisons are evidence for the compared behavior; failing comparisons identify implementation gaps and should not be papered over.

## Out Of Scope

- Chemistry algorithms.
- Runtime RDKit or Biopython dependencies.
- Regenerating all goldens by default.

## Revision Notes

- v1: Manifest discovery, fixture path checks, and reference generator conventions.
- v2: Implementation-vs-golden comparisons for committed per-feature golden JSON.
- v3: Named corpora, all-feature/all-corpus selection, generated evidence status, and dashboard synchronization.
- v4: Corpus-owned layout, typed TOML, compressed goldens, source locks, and corpus integrity checks.
- v5: Content-addressed validation evidence, strict golden metadata checks, non-empty comparison enforcement, and failure-safe selected status updates.
- v6: Make evidence hashes portable across LF and CRLF working trees.
- v7: Align SMILES semantic comparison with RDKit aromatic carbonyl valence and aromatic nH no-implicit handling.
- v8: Add large-corpus pack member metadata for PL-REX, Enamine Diversity, and PubChem-100k validation wiring.
- v9: Parallelize fixture comparison by default and add `--jobs N` for bounded validation runs.
- v10: Preserve fixture-level failure summaries in corpus status so the dashboard can show compact non-passing counts.
- v11: Add clean overall and per-target fixture progress bars to `cargo xtask validate`.
- v12: Allow manually curated `pubchem-manual-semantic` validation manifests without local reference generator files.
- v13: Generalize manual semantic reference manifests to named
  `*-manual-semantic` sources and let broad stereo perception validation record
  per-record sanitize errors instead of aborting an entire fixture.
- v14: Preserve every record from SDF pack fixtures when running stereo
  implementation comparisons, enabling descriptor-bearing ligand packs such as
  PL-REX to validate CIP behavior.
- v15: Add semantic output support for stored-axis validation issues and axis
  bond CIP descriptors.
