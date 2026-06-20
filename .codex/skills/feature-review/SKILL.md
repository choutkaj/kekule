---
name: feature-review
description: Independently audit molecules feature or remediation work for architecture compliance, correctness, validation claims, stale docs, metadata drift, tests, resource safety, and dashboard or skill synchronization.
---

# Feature Review

Use this skill for independent audit, not builder-mode implementation.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and every affected `features/<feature-id>/` directory.
2. Read `feature.toml`, `feature.md`, related validation manifests, and relevant code.
3. When the pull request references `fixes.md`, read the complete numbered stage and verify every acceptance gate.
4. Keep the review scoped to the feature or one remediation stage and its direct infrastructure support.

## Review focus

Lead with findings ordered by severity. Check:

- Architecture compliance.
- Correctness and edge cases.
- Whether a regression test demonstrates the reported defect.
- Mutation, perception invalidation, and transactional failure behavior.
- Parsing versus sanitization/perception boundaries.
- Panic safety and structured errors for untrusted input.
- Non-lossy parser/writer round trips and explicit rejection of unsupported chemistry.
- Algorithmic resource limits and stack safety where relevant.
- Validation claims, per-corpus generated evidence, manifest hashes, and whether overall `validated = true` is derived correctly.
- Whether comparisons, asserted fields, tests, or goldens were weakened merely to obtain a pass.
- Metadata and `feature.md` synchronization.
- Dashboard generation and skill workflow synchronization.
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

For a remediation PR, verify the complete stage-specific gate in `fixes.md`. When available and required:

```bash
cargo xtask corpus check --corpus all --require-data
cargo xtask validate --feature all --corpus all
```

Report every check that was not run.

## Output

Use code-review style:

1. Findings first, with file and line references.
2. Acceptance-gate status for the selected feature or `fixes.md` stage.
3. Open questions or assumptions.
4. Short summary only after findings.
