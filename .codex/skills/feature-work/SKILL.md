---
name: feature-work
description: Add, plan, implement, or maintain molecules features through the canonical feature workflow. Use for feature metadata, feature.md, implementation, validation fixtures, or dashboard updates.
---

# Feature Work

Use this skill for builder-mode work in `molecules`.

Canonical flow: add -> optional research -> plan -> implement.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and the relevant `features/<feature-id>/` directory.
2. Identify the canonical feature ID and any directly affected dependent feature IDs.
3. Treat `feature.toml` as machine-readable truth and `feature.md` as human-readable truth.
4. Keep the implementation scoped to the feature plus direct infrastructure support.
5. Add or update a regression test for each defect, behavior contract, or public API change.

## Feature metadata

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

Increment `version` only when behavior, public API, or validation contract intentionally changes. Set `implemented = true` only when implementation is complete. Let `cargo xtask validate ... --update` record passing evidence and synchronize `validated`.

Feature IDs and titles describe long-term capabilities, not temporary maturity levels. Use `version`, `implemented`, `validated`, Validation, and Revision Notes to describe partial coverage or missing goldens.

## Feature docs

Every feature has `features/<feature-id>/feature.md` with:

- Summary
- Behavior/API
- Implementation Notes
- Validation
- Out Of Scope
- Revision Notes

Keep feature docs concise and current. Do not recreate stale phase-specific planning, algorithm, specification, or validation documents.

## Implementation guardrails

- Follow `ARCHITECTURE.md` for API shape and module boundaries.
- Keep `Molecule` as the raw graph kernel shared by `SmallMolecule` and `MacroMolecule`.
- Keep parsing, sanitization, macromolecule validation, and preparation separate.
- Keep biomolecular labels in `SmcraHierarchy`, not core `Atom` or `Bond`, unless chemically general.
- Keep small-molecule and macromolecule sanitization options/reports/errors separate.
- Mutations that affect topology or interpreted chemistry must invalidate computed state.
- Failed transactional operations must leave inputs unchanged.
- Parsers return structured errors rather than panics.
- Writers reject chemistry they cannot encode faithfully.
- RDKit and Biopython are reference tools only, not Rust runtime dependencies.

## Validation data

Molecular validation fixtures must be externally supplied and provenance-pinned. Keep corpus descriptors, source locks, inputs, feature manifests, goldens, and evidence under `validation/corpora/<corpus-id>/`.

The `smoke` corpus is a fast wiring/regression tier, not broad validation by itself. During routine coding work, primarily validate against `smoke`, `pubchem-100`, `pubchem-1k`, `pdb-10`, `pdb-100`, and `pl-rex` where applicable. Treat `pubchem-100k` and `enamine-diversity` as large, occasional validation runs; do not run them routinely unless the user asks for them or the change clearly needs large-corpus coverage. Plain validation is read-only; use `--update` only after implementation-versus-golden comparison passes and should become committed evidence.

## Checks

Run applicable checks before handoff:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus smoke --require-data
cargo xtask validate --feature <feature-id> --corpus smoke
```

If metadata changes, run `cargo xtask dashboard` before `cargo xtask dashboard --check`.

Report every command that was not run and why.
