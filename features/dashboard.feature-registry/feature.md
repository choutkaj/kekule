# Feature Registry and Dashboard

## Summary

Keep feature metadata as the machine-readable source of truth and generate a deterministic dashboard from it.

## Behavior/API

- Exposes `cargo xtask features`.
- Exposes `cargo xtask dashboard` and `cargo xtask dashboard --check`.
- Exposes corpus-aware validation with feature and corpus `all` selectors.
- Exposes `cargo xtask skills --check` for repo-local feature workflow skill checks.

## Implementation Notes

- Feature metadata schema v2 requires `id`, `title`, `area`, `version`, `implemented`, `validated`, `description`, `depends_on`, and `validation_required`.
- Deprecated metadata keys `priority`, `status`, and `last_ai_review` are rejected.
- Each tracked feature directory must include `feature.md`.
- The dashboard renders a generated HTML table with feature metadata and one column per known corpus.
- Dashboard cells and the cached overall `validated` flag report structurally valid recorded pass evidence so generation is deterministic on clean checkouts without ignored large-corpus fixtures.
- `cargo xtask validate` and `cargo xtask features` remain the authority for whether recorded evidence is current for the files available in the checkout.
- Per-feature evidence is read from each corpus-owned `status.toml`.
- Dashboard generation rejects drift between generated evidence and the cached overall `validated` field.
- Boolean dashboard values render as check and cross marks, with compact failure counts only for recorded fixture-level validation failures.
- The dashboard omits the redundant overall `validated` column, uses centered rotated compact headers, includes corpus case counts from `corpus.toml`, and supports client-side column sorting.
- `features/DASHBOARD.html` remains the generated source artifact and is published as the GitHub Pages dashboard.

## Validation

- Current coverage is `xtask` unit-test based.
- CI should run formatting, clippy, workspace tests, dashboard check, and skill check.
- This feature does not require chemistry reference data.

## Out Of Scope

- Chemistry implementation or validation.
- Pulling feature metadata from external services.
- Automatically marking features implemented.

## Revision Notes

- v1: Registry, dashboard, validation command, schema v2, and repo-skill check behavior.
- v2: Corpus requirements, generated validation evidence, and corpus dashboard columns.
- v3: Corpus-owned evidence and validation layout.
- v4: Switch generated dashboard from Markdown to sortable HTML with compact rotated headers.
- v5: Make dashboard generation portable across clean checkouts while preserving content-addressed freshness checks in validation and feature listing.
- v6: Publish the HTML dashboard through GitHub Pages, show corpus counts in headers, center rotated labels, and surface compact fixture failure counts.
