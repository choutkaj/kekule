---
name: feature-work
description: Add, plan, implement, or maintain molecular features through the canonical feature workflow. Use for feature metadata, feature.md, implementation, validation fixtures, or dashboard updates.
---

# Feature Work

Use this skill for builder-mode work in `molecular`.

Canonical flow: add -> optional research -> plan -> implement.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and the relevant `features/<feature-id>/` directory.
2. Identify the canonical feature ID and any directly affected dependent feature IDs.
3. Treat `feature.toml` as machine-readable truth and `feature.md` as human-readable truth.
4. Keep the implementation scoped to the feature plus direct infrastructure support.
5. Add or update a regression test for each defect, behavior contract, or public API change.

## Feature metadata

Feature metadata uses schema v5:

- `id`
- `title`
- `area`
- `domains`
- `version`
- `status`
- `description`
- `depends_on`
- `validation_required`

Do not use `priority`, `implemented`, `last_ai_review`, or the removed global
`validated` flag. `status` is one of:

- `planned`: tracked design intent with no usable implementation.
- `experimental`: usable implementation whose contract may still change.
- `supported`: release-quality contract intended for ordinary users.
- `deprecated`: usable compatibility surface scheduled for removal.

`domains` is a non-empty list containing `small-molecule`, `macromolecule`, or
`infrastructure`. Shared graph, API, unit, and modeling features declare both
chemistry domains; infrastructure cannot be combined with a chemistry domain.

`depends_on` declares semantic feature prerequisites, not source-file imports.
Dependency IDs must exist and form a directed acyclic graph. Duplicate and
self-dependencies are invalid. A `supported` feature may depend only on
`supported` features; an `experimental` feature may depend on `experimental`
or `supported` features; a `deprecated` feature may depend on any implemented
feature; a `planned` feature may depend on any registered feature.

Increment `version` only when behavior, public API, or validation contract
intentionally changes. Set `status = "supported"` only when the public contract
is release-quality. Let `cargo xtask validate ... --update` record per-corpus
parity evidence.

Feature IDs and titles describe long-term capabilities, not temporary maturity
levels. Use `version`, `status`, the Validation section, per-corpus evidence,
and Revision Notes to describe partial coverage or missing goldens.

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

Use `pubchem-1k` and `pdb-100` as the normal required baselines where their external parity contracts apply. Treat `pubchem-100k`, `enamine-diversity`, `pdb-1000`, and domain-specific `pl-rex` as deliberate broader runs. Historical smoke sets are internal regression fixtures, not validation corpora. Plain validation is read-only; use `--update` only after implementation-versus-golden comparison passes and should become committed evidence.

## Checks

Run applicable checks before handoff:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus <corpus-id> --require-data
cargo xtask validate --feature <feature-id> --corpus <corpus-id>
```

If metadata changes, run `cargo xtask dashboard` before `cargo xtask dashboard --check`.

Report every command that was not run and why.
