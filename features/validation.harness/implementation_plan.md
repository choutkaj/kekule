# Implementation Plan

## Feature ID

`validation.harness`

## Goal

Provide repeatable repository infrastructure for comparing Rust behavior against reference-generated golden data. RDKit and Biopython are reference tools only; they must not become runtime dependencies of the Rust library.

The harness should make it clear when a feature has local unit coverage only, reference golden data, or documented manual validation.

## Public API

The user-facing API is command-line oriented:

- `cargo xtask validate --feature FEATURE_ID`
- Optional future flags:
  - `--update-golden` to regenerate reference golden files.
  - `--reference rdkit|biopython|manual` to select a reference backend.
  - `--fixture FIXTURE_NAME` to narrow validation.

Initial behavior:

- Reject unknown feature IDs.
- Discover validation manifests under `validation/features/<feature-id>/`.
- Report when a feature has no reference validation configured.
- Run configured fixture comparisons and fail on mismatches.
- Print reference tool versions used to produce or compare golden data.

No Rust library API should be added for this feature.

## Internal Modules Touched

Expected scope:

- `crates/xtask/src/main.rs` or a small `xtask` module split if validation grows.
- `validation/` directory structure for manifests, fixtures, generated goldens, and helper scripts.
- Optional reference scripts under `validation/reference/rdkit/` and `validation/reference/biopython/`.
- Feature docs under `features/validation.harness/`.

Do not add RDKit, Biopython, Python bindings, or chemistry reference crates to `crates/molecules`.

## Data Model

Recommended filesystem layout:

```text
validation/
  features/
    <feature-id>/
      validation.toml
      fixtures/
      golden/
      README.md
  reference/
    rdkit/
    biopython/
```

`validation.toml` should describe:

- `feature_id`
- `reference_tool`
- `reference_version`
- `fixtures`
- `golden_files`
- `comparison_mode`
- `notes`

Golden files should be normalized JSON. Each golden file should include:

- Feature ID.
- Fixture ID.
- Reference tool and version.
- Input source and checksum when useful.
- Normalized expected output.
- Generation timestamp only if it does not make comparison nondeterministic.

## Algorithm Outline

1. Resolve and validate `--feature`.
2. Locate `features/<feature-id>/feature.toml`.
3. Locate optional `validation/features/<feature-id>/validation.toml`.
4. If no validation manifest exists, report that reference validation is not configured and exit successfully unless a strict flag is later added.
5. Parse the manifest and verify all referenced fixture and golden paths exist.
6. For comparison mode, run the Rust behavior extractor for each fixture.
7. Load normalized golden JSON.
8. Compare normalized JSON exactly unless a manifest declares allowed tolerance for numeric fields.
9. Print a per-fixture summary.
10. Exit nonzero on missing fixtures, malformed manifests, reference command failures, or mismatches.

Golden generation should be an explicit mode separate from comparison.

## Tests

Add tests for:

- Known feature lookup succeeds.
- Unknown feature lookup fails.
- Missing validation manifest produces the documented non-strict message.
- Malformed validation manifest fails.
- Missing fixture path fails.
- Missing golden path fails.
- Exact JSON match passes.
- JSON mismatch fails with useful context.
- Tool version metadata is required when reference-generated goldens are present.
- Numeric tolerance comparison works only for fields declared tolerant.
- `--update-golden` cannot run accidentally during check mode.

Prefer unit tests for manifest parsing and comparison helpers, plus one small integration-style test using temporary validation directories.

## Reference Validation

This feature defines how reference validation works, so its own validation should be infrastructure-focused:

- Unit tests for manifest parsing and comparison.
- A tiny self-test fixture that compares deterministic JSON without RDKit or Biopython.
- Manual validation notes showing that a missing manifest is reported clearly.

Do not mark `validated = true` until there is test or manual evidence for the harness behavior itself.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature validation.harness
```

## Risks

- Reference tool versions may drift and change golden output.
- Golden regeneration can hide regressions if it is not explicit and reviewed.
- JSON normalization rules can accidentally erase meaningful differences.
- Python environment setup can make validation hard to reproduce.
- Coupling the Rust library to reference tools would violate architecture rules.

## Edge Cases

- Feature exists but has no validation manifest.
- Manifest references a feature ID that does not match its directory.
- Fixture names collide.
- Golden JSON contains unordered collections.
- Reference output contains platform-dependent floating point formatting.
- RDKit or Biopython is not installed.
- A feature has manual validation only.

## Explicitly Out of Scope

- Implementing chemistry algorithms.
- Adding RDKit or Biopython as Rust runtime dependencies.
- Automatically marking features validated.
- Regenerating all goldens by default.
- CI provider-specific workflow configuration beyond commands that CI can run.
