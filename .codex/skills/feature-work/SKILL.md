---
name: feature-work
description: Add, research, plan, implement, maintain, and remediate molecules features through the canonical feature and audit workflow. Use for feature metadata, feature.md, implementation, validation fixtures, dashboard updates, or a numbered fixes.md stage.
---

# Feature Work

Use this skill for builder-mode feature work in `molecules`.

## Modes

### Ordinary feature work

```text
add -> optional research -> plan -> implement
```

Start from one canonical feature ID and keep the change scoped to that feature plus direct infrastructure support.

### Audit remediation

When the request names `fixes.md` or an audit stage, use the `audit-remediation` skill. A remediation stage may span multiple features, but only one numbered stage may be implemented per branch and pull request.

Review is separate: use `feature-review` for independent audit.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and the relevant feature directory.
2. For remediation, also read `fixes.md` and the complete selected stage.
3. Identify every affected feature ID before editing.
4. Treat `feature.toml` as machine-readable truth and `feature.md` as human-readable truth.
5. Add a regression test for the defect or contract change before closing the work.

## Metadata

Feature metadata uses schema v2:

- `id`
- `title`
- `area`
- `version`
- `implemented`
- `validated`
- `description`
- `depends_on`
- `validation_required`

Do not use `priority`, `status`, or `last_ai_review`.

Increment `version` only when the feature's behavior, public API, or validation contract intentionally changes. Do not increment for typo fixes or behavior-preserving file moves.

Set `implemented = true` only when implementation is complete. Declare broad-validation corpora in `validation_required`. Do not hand-set corpus results: `cargo xtask validate ... --update` records passing evidence and synchronizes overall `validated`.

Feature IDs and titles must describe canonical long-term capabilities, not maturity levels. Use `version`, `implemented`, `validated`, the Validation section, and Revision Notes to describe maturity, partial coverage, or missing goldens.

Molecular validation fixtures must be externally supplied, not invented toy systems. Keep corpus descriptors, source locks, inputs, feature manifests, goldens, and evidence under `validation/corpora/<corpus-id>/`. Record source URLs and checksums in `sources.lock.json`, and generate molecular golden data only with the declared reference software.

The `tiny` corpus is a fast wiring and regression tier, not broad validation by itself. Use declared PubChem, PDB, PL-REX, and Enamine corpora where applicable. Plain validation is read-only; use `--update` only after successful implementation-versus-golden comparison should become committed evidence.

## Feature docs

Every feature must have `features/<feature-id>/feature.md` with:

- Summary
- Behavior/API
- Implementation Notes
- Validation
- Out Of Scope
- Revision Notes

Keep the file concise and current. Do not recreate stale phase-specific planning, algorithm, specification, or validation documents.

## Implementation constraints

- One core molecular graph is shared by `SmallMolecule` and `MacroMolecule`.
- Biomolecular labels belong in `BioHierarchy`, not core `Atom`, unless chemically general.
- Parsing and chemical perception remain separate.
- Topology or chemistry-relevant mutation invalidates computed perception state.
- Failed transactional operations leave their input unchanged.
- Parsers return structured errors rather than panic.
- Writers reject representations they cannot encode faithfully.
- RDKit and Biopython are reference tools only, not Rust runtime dependencies.
- Do not weaken normalized comparisons, remove asserted fields, or regenerate goldens merely to obtain a pass.

## Checks

Run applicable checks before handoff:

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

For remediation, run the complete stage gate in `fixes.md`. When available and required:

```bash
cargo xtask corpus check --corpus all --require-data
cargo xtask validate --feature all --corpus all
```

If metadata changes, run `cargo xtask dashboard` before `cargo xtask dashboard --check`.

Report every command that was not run and why.

## Commit attribution

End every commit message with:

```text
Co-authored-by: Codex <noreply@openai.com>
```
