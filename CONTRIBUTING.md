# Contributing

This repository is organized around feature-scoped development.

## Before implementing

1. Pick or create a feature directory under `features/`.
2. Fill in or update `feature.toml`, `spec.md`, `algorithm.md`, and `validation.md`.
3. Confirm how the feature interacts with the central architecture.

## Before submitting

```bash
cargo fmt --all
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
```

Reference validation should be run when the feature has RDKit or Biopython golden data.
