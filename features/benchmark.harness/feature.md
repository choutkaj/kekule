# Optional Reference Benchmark Harness

## Summary

Provide convenient, repeatable comparisons between the Rust implementation and
the checked-in RDKit, Biopython, DSSP, and manually reviewed reference goldens.
Benchmark observations are informational: they do not determine feature health,
release status, CI success, or whether the repository is ready for development.

## Behavior/API

- Exposes `cargo xtask benchmark --feature FEATURE_ID|all [--corpus CORPUS_ID|all] [--fixture PATH] [--accept-implementation-goldens] [--jobs N]`.
- With no `--corpus`, selects only manifest-backed PubChem-1k and PDB-100
  targets. Explicit `--corpus all` also selects registered broad corpora.
- Discovers availability exclusively from
  `benchmarks/corpora/<corpus-id>/features/<feature-id>.toml`.
- Skips absent combinations and planned features for selectors containing
  `all`. A concrete feature/corpus request without a manifest is an error.
- Requires deterministic gzip goldens for every selected fixture and compares
  normalized implementation output with each golden's `expected` payload.
- Runs fixture comparisons in a bounded worker pool; `--jobs N` overrides the
  automatic maximum of four workers.
- A difference or execution error produces a nonzero exit code for the caller
  who explicitly requested the benchmark.
- Every selected target records an observation before return, including
  matches, differences, partial fixture runs, and execution errors, then
  regenerates the dashboard.
- Stores tracked observations in corpus-owned `results.toml` files with
  outcome, scope, counts, first detail, reference versions, timestamp,
  normalized manifest digest, and identity-neutral input digest.
- Exposes `cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]`
  as an optional integrity utility. `--require-data` additionally checks ignored
  local fixture bytes and pack membership.
- Allows explicit implementation-golden acceptance only for one concrete
  `*-manual-semantic` feature/corpus target. It cannot replace RDKit- or
  Biopython-generated goldens.

## Implementation Notes

- RDKit reference tools live under `benchmarks/reference/rdkit/`; Biopython and
  DSSP reference tools live under `benchmarks/reference/biopython/`.
- Corpora, source locks, manifests, and goldens remain provenance-pinned.
  Source-lock entry order and categories are canonical membership; builders
  reuse them unless an explicit reselection mode names a new corpus version.
- The benchmark input digest covers manifests, corpus descriptors, source
  locks, fixtures, goldens, implementation/comparison source, reference
  generators/environments, and normalized external dependency lock entries.
- Core implementation sources use package-neutral
  `implementation/core/src/...` digest labels, so renaming the package
  directory does not change benchmark identity by itself.
- Digest inputs exclude absolute checkout paths, repository and package
  identity, workspace-local package names, feature status/docs, and timestamps.
- UTF-8 inputs are normalized to LF before hashing; binary inputs remain
  byte-exact.
- `results.toml` snapshots are observations, not freshness or health claims.
  Historical validation-status v2 snapshots are explicitly labeled legacy;
  the next benchmark replaces each selected target with the current schema.
- Representation-only graph differences such as undirected bond orientation,
  bond-array order, and ring-atom order are normalized before comparison.
- Sanitized SMILES semantics derive RDKit-comparable aromatic-nitrogen
  hydrogen/valence output from Kekule's represented explicit-H plus perceived
  implicit-H layers; the benchmark adapter does not require perception to
  rewrite represented atoms.
- Reference tools are never Rust runtime dependencies.

## Tests

- xtask regressions cover optional manifest discovery, small-corpus defaults,
  explicit broad selection, missing-manifest errors, planned-feature exclusion,
  result recording, partial/error observations, legacy replacement, neutral
  dashboard rendering, and removed validation schema/CLI names.
- Digest regressions prove checkout/repository/package identity is ignored and
  material implementation, fixture, golden, manifest, generator, and external
  dependency changes are detected.
- Corpus checks retain source-lock, nested-prefix, provenance, gzip,
  comparison-normalization, and manual-golden safeguards.
- SMILES semantic regressions cover represented `[nH]` and a perceived
  aromatic-nitrogen hydrogen whose reference-facing valence is derived without
  sanitizer representation feedback.

## Benchmarks

- Manifest-backed corpora provide optional external-reference coverage. Run a
  benchmark when its information is useful; routine development and release
  checks do not invoke it.

## Out Of Scope

- Chemistry algorithms.
- Runtime RDKit, Biopython, or DSSP dependencies.
- Regenerating or accepting goldens by default.
- Treating a benchmark snapshot as feature health or a release gate.

## Revision Notes

- v1-v40: Evolved the former required external-parity harness,
  corpus layout, evidence snapshots, parallel comparison, and integrity tools.
- v41: Recast the required external-parity layer as optional benchmarking, remove the
  validation CLI and required-manifest schema, default to PubChem-1k/PDB-100,
  record neutral match/difference/error observations, migrate legacy snapshots,
  and make input digests independent of repository identity.
- v42: Discover implementation sources from the Kekule core crate while
  labeling them through a package-neutral digest namespace.
- v43: Derive reference-facing aromatic-nitrogen hydrogen and valence semantics
  from separated represented and perceived state after removal of sanitizer
  representation feedback.
- v44: Migrate the former sanitizer benchmark identity and Kekule call paths to
  the explicit default perception feature while preserving pinned external
  reference schemas and goldens.
