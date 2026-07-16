# Feature Registry and Dashboard

## Summary

Keep feature metadata as the machine-readable source of truth and generate a deterministic dashboard from it.

## Behavior/API

- Exposes `cargo xtask features`.
- Exposes `cargo xtask dashboard` and `cargo xtask dashboard --check`.
- Exposes corpus-aware validation with feature and corpus `all` selectors.
- Exposes `cargo xtask skills --check` for repo-local feature workflow skill checks.

## Implementation Notes

- Feature metadata schema v3 requires `id`, `title`, `area`, `version`, `implemented`, `description`, `depends_on`, and `validation_required`.
- Deprecated metadata keys `priority`, `status`, and `last_ai_review`, plus the
  removed global `validated` key, are rejected.
- Each tracked feature directory must include `feature.md`.
- The dashboard renders a generated HTML table with feature metadata and one column per known corpus.
- Dashboard cells report structurally valid recorded per-corpus parity evidence so generation is deterministic on clean checkouts without ignored large-corpus fixtures.
- `cargo xtask validate` is the authority for checking parity against the current checkout; `cargo xtask features` lists feature metadata without deriving a global validation result.
- Per-feature evidence is read from each corpus-owned `status.toml`.
- Manifest-backed local-only corpus cells display their recorded optional
  evidence without contributing to routine required parity checks.
- Boolean dashboard values render as check and cross marks, with compact failure counts only for recorded fixture-level validation failures.
- Required validation with no current recorded status, stale or incomplete evidence, or no fixture-level failure count renders as an unknown `?` marker rather than a confirmed failure.
- The feature schema has no global validation boolean; the dashboard presents the per-corpus matrix directly, uses centered rotated compact headers, includes corpus case counts from `corpus.toml`, and supports client-side column sorting.
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
- v7: Distinguish uncertain required validation from confirmed fixture failures with a question-mark dashboard marker.
- v8: Render all registered corpora, including local-only optional validation
  tiers, as dashboard columns.
- v9: Remove the global `validated` metadata flag and treat per-corpus parity
  evidence as the complete validation state.
