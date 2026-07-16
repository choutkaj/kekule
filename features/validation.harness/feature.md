# Reference Validation Harness

## Summary

Provide repeatable infrastructure for comparing Rust behavior against generated
or manually curated external-reference golden data.

## Behavior/API

- Exposes `cargo xtask validate --feature FEATURE_ID|all --corpus CORPUS_ID|all [--fixture PATH] [--update] [--jobs N]`.
- Defaults omitted `--corpus` to `all`.
- Discovers corpus manifests under `validation/corpora/<corpus-id>/features/<feature-id>.toml`.
- Verifies listed fixture paths exist.
- Requires one deterministic gzip golden under the corpus `golden/<feature-id>/` directory for each listed fixture.
- Compares normalized Rust implementation output against each golden file's `expected` payload.
- Compares fixtures in parallel by default using all available processors; `--jobs N` bounds validation to a fixed worker count.
- Prints compact progress output for the overall target set and for fixture comparison within each feature/corpus target.
- Accepts only declared implementation-vs-golden comparison manifests for
  required validation of available features.
- Declares corpus availability in typed metadata; local-only corpora remain
  explicitly runnable with a concrete `--corpus` selector but are rejected from
  feature `validation_required` lists because a clean checkout cannot recompute
  their fixture evidence. `--corpus all`, including the default when the flag is
  omitted, selects every required target plus every manifest-backed target
  whose feature status is `experimental`, `supported`, or `deprecated`, across
  all registered corpora, including local-only large corpora. `planned`
  features are excluded because they have no usable implementation.
- Verifies golden `feature_id`, `corpus_id`, `fixture_path`, current fixture SHA-256, and reference tool/version metadata before comparing payloads.
- Records content-addressed pass evidence over manifests, source locks, fixtures, goldens, Rust source, Cargo manifests, feature metadata, and reference generator/environment files when the reference tool is generator-backed.
- Supports manually curated semantic references when the manifest
  `reference_tool` ends in `-manual-semantic`; these are evidence-backed by
  pinned source fixtures, manifest metadata, committed goldens, and Rust
  implementation sources rather than by local reference generator files.
- Canonicalizes UTF-8 text inputs to LF before hashing so evidence and the
  separately recorded manifest hash are stable across Windows and Linux
  checkouts; binary inputs remain byte-exact.
- Normalizes representation-only graph differences such as undirected bond endpoint orientation, bond array order, and ring atom order before comparison.
- Treats non-applicable feature/corpus combinations as skips and missing required manifests as errors.
- Exposes `cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]`.
- Always verifies the committed 20-case smoke fixture bytes, even when
  `--require-data` is omitted; the flag remains necessary for ignored larger
  corpora.
- Keeps ordinary validation read-only; `--update` clears selected stale passes before running, records evidence for successful selected targets, records fixture-level failure summaries for failed comparisons, and regenerates the dashboard.
- Allows an independently reviewed implementation-semantic change to be
  accepted explicitly with `--accept-implementation-goldens`, but only for one
  concrete feature/corpus whose reference tool is `*-manual-semantic`. This
  mode cannot replace RDKit- or Biopython-generated goldens and does not update
  validation status or the dashboard.
- Allows `--fixture PATH` to isolate one declared fixture for diagnosis or
  reviewed manual-golden acceptance without changing status evidence.
- Supports corpus pack member checks using PubChem defaults or corpus-declared SDF member properties and SMILES title prefixes.
- Preserves all records in SDF pack fixtures for stereo implementation comparisons.
- Emits structured semantic JSON for stored axis stereo validation issues and
  includes assigned axis `M`/`P` descriptors in bond descriptor maps.
- Keeps non-isomeric SMILES parse/write comparison filtering synchronized with
  the RDKit reference generator; stereo-bearing syntax is validated by the
  dedicated representation, perception, CIP, and isomeric-SMILES features.

## Implementation Notes

- RDKit reference generators live under `validation/reference/rdkit/`.
- Biopython reference generators live under `validation/reference/biopython/`.
- `*-manual-semantic` manifests use externally pinned fixtures and manually
  reviewed semantic goldens; they do not require local reference generator
  files. Current semantic reference labels include `pubchem-manual-semantic`
  and `enamine-manual-semantic`.
- Golden data should be normalized JSON and include reference tool versions.
- Corpus descriptors and feature manifests use typed TOML; source selection and checksums live in `sources.lock.json`.
- The smoke corpus data directory is intentionally checked in and exempt from the default ignore rule for larger generated corpus data.
- Every registered non-smoke corpus declares `local_only = true` because its
  data directory is ignored. These corpora are excluded from
  `validation_required`, but are selected by explicit or default `--corpus all`
  when a feature manifest exists. Their local data must therefore be present for
  an all-corpus run.
- Source pack records may declare `member_id_property` for SDF packs or `member_title_prefix` for SMILES packs when the corpus does not use PubChem CID metadata.
- Status evidence records fixture and comparison counts, reference versions,
  the line-ending-normalized manifest SHA-256, a versioned evidence input list,
  evidence SHA-256, and validation time.
- Evidence is considered current only when recomputing it from the current checkout produces the stored schema version and hash.
- Failed comparison status records fixture count, successful comparison count, failed count, and the first fixture-level failure without recording pass evidence.
- Evidence schema v2 includes cross-platform text line-ending normalization.
- Repeated `--update` runs preserve timestamps when the evidence hash is unchanged.
- The validation command uses the Rust implementation only; RDKit, Biopython, or manually pinned external sources are used to generate or curate goldens, not to run validation.
- The manual GitHub validation workflow likewise checks committed fixtures and
  goldens without installing unused RDKit or Biopython runtime environments.
- Manual workflow inputs are passed through environment variables rather than
  interpolated into shell source.
- Manual workflow log capture enables shell pipeline failure propagation, so a
  failed corpus check or validation command cannot be masked by `tee`.
- The main CI release matrix exercises all workspace features against the
  committed dependency lockfile.
- Byte-exact fixture and corpus checks stream file contents through SHA-256
  instead of loading entire broad-corpus packs into memory; SHA-256 is optimized
  in development builds so local integrity checks remain practical.
- One validation command caches exact and line-ending-normalized hashes for
  immutable inputs, so features sharing broad-corpus fixtures verify the same
  bytes without repeatedly reading them from disk.
- Development builds optimize only the third-party hashing, gzip, and JSON
  stack used by the harness, keeping project code debuggable while making large
  golden verification practical.
- Corpus integrity checks stream-decode golden JSON into a syntax sink instead
  of materializing each complete expected-value tree when only structural
  validity is required.
- Independent compressed goldens are syntax-checked in parallel with stable
  path-ordered error reporting.
- Manual-semantic golden acceptance is an explicit review operation, not a way
  to make an unexplained failure pass; the motivating chemistry or file-format
  behavior must be independently checked before accepting the snapshot.
- Manual-semantic acceptance uses the same bounded worker-count policy as
  comparison and writes deterministic gzip streams with zero timestamps.
- Fixture comparison uses a bounded worker pool while status writes and dashboard regeneration remain single-threaded.
- Progress output uses plain ASCII bars and throttled checkpoint updates so it stays readable in terminals and captured logs.
- Reference tools are never Rust runtime dependencies.

## Validation

- Current coverage is infrastructure unit-test based plus live corpus comparisons against committed external-source goldens.
- Passing comparisons are evidence for the compared behavior; failing comparisons identify implementation gaps and should not be papered over.
- Unit coverage rejects local-only corpus IDs in `validation_required` while
  verifying that explicit and default all-corpus selection includes every
  available manifest-backed local corpus and excludes planned features.

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
- v16: Restore RDKit-generator parity for the non-isomeric SMILES validation
  subset after the stereo validation expansion.
- v17: Add fixture-isolated validation, report every failing fixture, and add a
  guarded acceptance command for independently reviewed `*-manual-semantic`
  implementation goldens.
- v18: Make smoke a fixed, fully committed 20-case suite and always verify its
  source fixtures during corpus integrity checks.
- v19: Add typed local-only corpus metadata, reject every ignored corpus from
  repository-wide required evidence, and retain them for explicit runs.
- v20: Normalize the separately recorded manifest hash so evidence generated on
  Windows remains current in Linux CI checkouts.
- v21: Display every registered corpus on the dashboard while keeping
  local-only evidence optional and outside required routine parity checks.
- v22: Remove the global feature validation boolean; corpus status files and
  dashboard cells now carry the full parity result without metadata syncing.
- v23: Make explicit and omitted `--corpus all` select every applicable
  registered corpus, including large local-only corpora.
- v24: Select manifest-backed validation targets from explicit feature release
  status, including experimental and deprecated implementations while
  excluding planned features.
- v25: Remove unused reference-environment setup from the manual validation
  workflow so it matches the harness's committed-golden execution model.
- v26: Pass manual workflow inputs through environment variables to avoid
  treating dispatch values as shell source.
- v27: Enable explicit pipeline failure propagation for logged manual
  validation commands.
- v28: Run the main release checks with all workspace features and the committed
  dependency lockfile.
- v29: Stream byte-exact corpus hashes, remove a duplicate evidence-buffer copy,
  and optimize SHA-256 in development builds for bounded, practical broad-corpus
  checks.
- v30: Cache exact and evidence-normalized input hashes for the duration of one
  validation run, eliminating repeated broad-corpus reads across feature
  targets without weakening per-input evidence.
- v31: Optimize the harness's third-party hashing, gzip, and JSON dependencies
  in development builds so broad validation spends time on chemistry rather
  than debug-mode evidence decoding.
- v32: Validate compressed golden JSON as a stream during corpus integrity
  checks, avoiding unnecessary full expected-tree allocation.
- v33: Parallelize independent golden syntax checks while preserving
  deterministic path-ordered failures.
- v34: Normalize expected and implementation JSON in place during comparison,
  eliminating two full-tree clones per parallel fixture and substantially
  reducing peak memory for broad structural-output validation.
- v35: Cap the automatic fixture-worker default at four so machines with many
  logical CPUs do not multiply large JSON comparison memory unexpectedly;
  explicit `--jobs N` remains available for provisioned validation hosts.
