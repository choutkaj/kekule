# Feature Registry and Dashboard

## Summary

Keep feature metadata as the machine-readable source of truth and generate a deterministic dashboard from it.

## Behavior/API

- Exposes `cargo xtask features`.
- Exposes `cargo xtask dashboard` and `cargo xtask dashboard --check`.
- Exposes `cargo xtask validate --feature FEATURE_ID`.
- Exposes `cargo xtask skills --check` for repo-local feature workflow skill checks.

## Implementation Notes

- Feature metadata schema v2 requires `id`, `title`, `area`, `version`, `implemented`, `validated`, `description`, and `depends_on`.
- Deprecated metadata keys `priority`, `status`, and `last_ai_review` are rejected.
- Each tracked feature directory must include `feature.md`.
- The dashboard renders Feature, Title, Area, Version, Implemented, and Validated.
- Boolean dashboard values render as check and cross marks.

## Validation

- Current coverage is `xtask` unit-test based.
- CI should run formatting, clippy, workspace tests, dashboard check, and skill check.
- This feature does not require chemistry reference data.

## Out Of Scope

- Chemistry implementation or validation.
- Web dashboards.
- Automatically marking features implemented or validated.
- Pulling feature metadata from external services.

## Revision Notes

- v1: Registry, dashboard, validation command, schema v2, and repo-skill check behavior.
