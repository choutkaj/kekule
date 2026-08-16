## Feature ID

<!-- Example: core.graph -->

## Summary

## Tests

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --workspace --all-targets -- -D warnings`
- [ ] `cargo test --workspace`
- [ ] `cargo test --workspace --doc`
- [ ] `cargo xtask dashboard --check`
- [ ] Targeted `cargo xtask corpus check --corpus <id> --require-data`, if applicable

## Optional Benchmarks

- [ ] Targeted `cargo xtask benchmark --feature <id> --corpus <id>`, if deliberately run
- [ ] Broad benchmark, if deliberately run

<!-- Benchmark results are informational and are never required for PR health. -->

## Notes

## Commands Not Run

<!-- List every omitted applicable command and the reason. -->

## Release-Sensitive Files

- [ ] Changes under `.github/`, benchmark generators/locks/goldens, corpus
      descriptors, or feature metadata received owner review.
