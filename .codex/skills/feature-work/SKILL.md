---
name: feature-work
description: Add, research, plan, implement, and maintain molecules features through the canonical feature workflow. Use when working on feature metadata, feature.md, implementation, validation fixtures, dashboard updates, or feature lifecycle changes in this repository.
---

# Feature Work

Use this skill for builder-mode feature work in `molecules`.

The workflow is:

```text
add -> optional research -> plan -> implement
```

Review is separate: use `feature-review` for independent audit.

## Start

1. Read `ARCHITECTURE.md`, `AGENTS.md`, and the feature directory.
2. Identify or create the feature ID under `features/`.
3. Keep changes scoped to the requested feature and direct infrastructure support.
4. Treat `feature.toml` as machine-readable truth and `feature.md` as human-readable truth.

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

Increment `version` only when the feature's behavior, public API, or validation contract intentionally changes. Do not increment for typo fixes.

Set `implemented = true` only when implementation is complete. Declare the corpora required for broad validation in `validation_required`. Do not hand-set corpus results: `cargo xtask validate ... --update` records passing evidence and synchronizes overall `validated`.

Feature IDs and titles must describe canonical long-term capabilities, not maturity levels. Do not encode incomplete implementation status in feature IDs or titles. Use `version`, `implemented`, `validated`, the Validation section, and Revision Notes to describe maturity, partial coverage, or missing goldens.

Molecular validation fixtures must be externally supplied, not invented toy systems. Record source URL and checksum provenance in each corpus manifest, and generate molecular golden data only with the declared reference software.

The `tiny` corpus is a fast wiring and regression tier, not broad validation by itself. Use the declared PubChem, PDB, PL-REX, and Enamine corpora where applicable. Plain validation is read-only; use `--update` only after a successful implementation-vs-golden comparison should become committed evidence.

## Feature Docs

Every feature must have `features/<feature-id>/feature.md` with these sections:

- Summary
- Behavior/API
- Implementation Notes
- Validation
- Out Of Scope
- Revision Notes

Keep the file concise and current. Do not recreate stale phase-specific docs such as separate planning, algorithm, spec, or validation markdown files.

## Implementation

Preserve project architecture:

- One core molecular graph is shared by `SmallMolecule` and `MacroMolecule`.
- Biomolecular labels belong in `BioHierarchy`, not core `Atom`, unless chemically general.
- Parsing and chemical perception remain separate.
- Topology or chemistry-relevant mutation invalidates computed perception state.
- RDKit and Biopython are reference tools only, not Rust runtime dependencies.

## Checks

Run relevant checks before handoff:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask validate --feature <feature-id> --corpus tiny
cargo xtask validate --feature all --corpus all
```

If metadata changes, run `cargo xtask dashboard` before `cargo xtask dashboard --check`.
