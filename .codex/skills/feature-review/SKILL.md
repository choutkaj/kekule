---
name: feature-review
description: Independently audit molecules feature work for architecture compliance, correctness, validation claims, stale docs, metadata drift, tests, and dashboard or skill synchronization.
---

# Feature Review

Use this skill for independent audit, not builder-mode implementation.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and `features/<feature-id>/`.
2. Read `feature.toml`, `feature.md`, related validation manifests, and relevant code.
3. Keep the review scoped to the feature and direct infrastructure support.

## Review Focus

Lead with findings ordered by severity. Check:

- Architecture compliance.
- Correctness and edge cases.
- Mutation and perception invalidation, when relevant.
- Parsing versus sanitization/perception boundaries, when relevant.
- Validation requirements, per-corpus generated evidence, manifest hashes, and whether overall `validated = true` is derived correctly.
- Metadata and `feature.md` sync.
- Dashboard generation and skill workflow sync.
- Tests and missing coverage.

## Required Checks

When feasible, run or verify:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask validate --feature <feature-id> --corpus tiny
cargo xtask validate --feature all --corpus all
cargo xtask corpus check --corpus all
```

Report any checks that were not run.

## Output

Use code-review style:

1. Findings first, with file and line references.
2. Open questions or assumptions.
3. Short summary only after findings.
