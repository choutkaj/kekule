## Feature ID

<!-- Example: core.graph -->

## Summary

## Validation

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --workspace --all-targets -- -D warnings`
- [ ] `cargo test --workspace`
- [ ] `cargo test --workspace --doc`
- [ ] `cargo xtask dashboard --check`
- [ ] `cargo xtask skills --check`
- [ ] `cargo xtask corpus check --corpus tiny --require-data`
- [ ] `cargo xtask validate --feature all --corpus tiny`
- [ ] Reference validation, if applicable

## Notes

## Commands Not Run

<!-- List every omitted applicable command and the reason. -->

## Release-Sensitive Files

- [ ] Changes under `.github/`, validation generators/locks/goldens, corpus
      descriptors, or feature metadata received owner review.
