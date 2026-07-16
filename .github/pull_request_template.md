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
- [ ] Targeted `cargo xtask corpus check --corpus <id> --require-data`, if applicable
- [ ] Targeted `cargo xtask validate --feature <id> --corpus <id>`, if applicable
- [ ] Broad or full reference validation, if release scope requires it

## Notes

## Commands Not Run

<!-- List every omitted applicable command and the reason. -->

## Release-Sensitive Files

- [ ] Changes under `.github/`, validation generators/locks/goldens, corpus
      descriptors, or feature metadata received owner review.
