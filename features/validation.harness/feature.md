# Reference Validation Harness

## Summary

Provide repeatable infrastructure for comparing Rust behavior against reference-generated golden data.

## Behavior/API

- Exposes `cargo xtask validate --feature FEATURE_ID|all --corpus CORPUS_ID|all [--update]`.
- Defaults omitted `--corpus` to `tiny` for compatibility.
- Discovers corpus manifests under `validation/corpora/<corpus-id>/features/<feature-id>.toml`.
- Verifies listed fixture paths exist.
- Requires one deterministic gzip golden under the corpus `golden/<feature-id>/` directory for each listed fixture.
- Compares normalized Rust implementation output against each golden file's `expected` payload.
- Accepts only declared implementation-vs-golden comparison manifests for required implemented validation.
- Verifies golden `feature_id`, `corpus_id`, `fixture_path`, current fixture SHA-256, and reference tool/version metadata before comparing payloads.
- Records content-addressed pass evidence over manifests, source locks, fixtures, goldens, Rust source, Cargo manifests, feature metadata, and reference generator/environment files.
- Normalizes representation-only graph differences such as undirected bond endpoint orientation, bond array order, and ring atom order before comparison.
- Treats non-applicable feature/corpus combinations as skips and missing required manifests as errors.
- Exposes `cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]`.
- Keeps ordinary validation read-only; `--update` clears selected stale passes before running, records evidence only for successful selected targets, synchronizes overall `validated`, and regenerates the dashboard.

## Implementation Notes

- RDKit reference generators live under `validation/reference/rdkit/`.
- Biopython reference generators live under `validation/reference/biopython/`.
- Golden data should be normalized JSON and include reference tool versions.
- Corpus descriptors and feature manifests use typed TOML; source selection and checksums live in `sources.lock.json`.
- Status evidence records fixture and comparison counts, reference versions, the manifest SHA-256, a versioned evidence input list, evidence SHA-256, and validation time.
- Evidence is considered current only when recomputing it from the current checkout produces the stored schema version and hash.
- Repeated `--update` runs preserve timestamps when the evidence hash is unchanged.
- The validation command uses the Rust implementation only; RDKit and Biopython are used to generate goldens, not to run validation.
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
