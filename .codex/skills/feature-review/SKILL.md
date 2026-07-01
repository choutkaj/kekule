---
name: feature-review
description: Independently audit molecules feature work for architecture compliance, correctness, validation claims, tests, resource safety, stale docs, metadata drift, and dashboard or skill synchronization.
---

# Feature Review

Use this skill for independent audit, not builder-mode implementation.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and every affected `features/<feature-id>/` directory.
2. Read `feature.toml`, `feature.md`, related validation manifests, relevant code, and public examples.
3. Keep the review scoped to the selected feature IDs and direct infrastructure support.

## Review focus

Lead with findings ordered by severity. Check:

- Architecture compliance and public API shape.
- Correctness, edge cases, and panic safety.
- Regression coverage for the defect or contract being changed.
- Mutation invalidation and transactional failure behavior.
- Parsing versus sanitization, macromolecule validation, and preparation boundaries.
- Structured errors for untrusted input.
- Non-lossy parser/writer round trips and explicit rejection of unsupported chemistry.
- Algorithmic resource limits and stack safety where relevant.
- Validation claims, corpus evidence, manifest hashes, and derived `validated = true` status.
- Whether comparisons, asserted fields, tests, or goldens were weakened merely to obtain a pass.
- Synchronization of feature metadata, feature docs, dashboard, and skills.
- Commands not run and unsupported claims in the PR description.

## Required checks

When feasible, run or verify:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature <feature-id> --corpus tiny
```

Report every check that was not run.

## Output

Use code-review style:

1. Findings first, with file and line references.
2. Acceptance-gate status for the selected feature IDs.
3. Open questions or assumptions.
4. Short summary only after findings.
