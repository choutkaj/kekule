# Contributing

This repository is organized around feature-scoped development.

## Before implementing

1. Pick or create a feature directory under `features/`.
2. Fill in or update schema-v2 `feature.toml` and canonical `feature.md`.
3. Confirm how the feature interacts with the central architecture.

## Before submitting

```bash
cargo fmt --all
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo test --workspace --doc
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature all --corpus tiny
```

Reference validation should be run when the feature has RDKit or Biopython golden data.

Changes to `.github/`, validation generators, source locks, goldens, corpus
descriptors, or feature metadata require explicit owner review. `CODEOWNERS`
records that policy, but it is not enforced until branch protection or a
ruleset is enabled.
