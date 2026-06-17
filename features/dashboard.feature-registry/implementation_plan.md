# Implementation Plan

## Feature ID

`dashboard.feature-registry`

## Goal

Make feature metadata under `features/*/feature.toml` the source of truth for repository feature status, and generate `features/DASHBOARD.md` deterministically from that metadata.

This feature is infrastructure only. It should not change chemistry APIs, validation semantics, parser behavior, or feature implementation status outside its own metadata when the feature is truly complete.

## Public API

The user-facing API is the `xtask` command surface:

- `cargo xtask features`
- `cargo xtask dashboard`
- `cargo xtask dashboard --check`
- `cargo xtask validate --feature FEATURE_ID`

Expected behavior:

- `cargo xtask features` lists discovered features in deterministic order.
- `cargo xtask dashboard` rewrites `features/DASHBOARD.md` from feature metadata.
- `cargo xtask dashboard --check` exits successfully when the committed dashboard matches generated output.
- `cargo xtask dashboard --check` exits nonzero with a clear message when the dashboard is stale.
- `cargo xtask validate --feature FEATURE_ID` rejects unknown features and accepts known feature IDs.

No library API should be added for this feature.

## Internal Modules Touched

Expected scope:

- `crates/xtask/src/main.rs`.
- `features/registry.toml`, only if a registry-level schema version or settings field is needed.
- `features/DASHBOARD.md`, generated only by `cargo xtask dashboard`.
- Unit tests for pure parsing/rendering helpers if helpers are split into testable functions.
- Optional integration-style tests using temporary feature directories if the current single-file `xtask` layout remains simple enough.

Do not hand-edit `features/DASHBOARD.md` content. Always generate it.

## Data Model

Represent each feature with required metadata:

- `id: String`
- `title: String`
- `area: String`
- `priority: String`
- `status: String`
- `implemented: bool`
- `validated: bool`
- `last_ai_review: String`
- `description: String`
- `depends_on: Vec<String>`

The first implementation may continue to parse the simple current TOML shape manually, but it should validate types enough to avoid silently accepting malformed values. If parsing grows beyond the current limited shape, switch to a TOML parser in `xtask` rather than extending ad hoc string handling too far.

Recommended allowed values:

- `priority`: `P0`, `P1`, `P2`, `P3`.
- `status`: `planned`, `implemented`, `validated`, `deferred`, `blocked`.
- `implemented`: boolean.
- `validated`: boolean.

The generated dashboard should display:

- Feature ID.
- Title.
- Area.
- Priority.
- Status.
- Implemented yes/no.
- Validated yes/no.
- Last AI review date.

## Algorithm Outline

1. Read `features/` directory entries.
2. Skip non-directories and template or hidden feature directories whose names start with `_`.
3. For each feature directory, require `feature.toml`.
4. Parse feature metadata into a typed `Feature` struct.
5. Validate required keys and basic field values.
6. Validate that the `id` matches the feature directory name.
7. Validate that each `depends_on` entry points to an existing feature ID.
8. Sort features by ID for stable output.
9. Render `features/DASHBOARD.md` with deterministic markdown.
10. In write mode, write the rendered dashboard.
11. In check mode, compare existing dashboard text byte-for-byte against rendered output and fail if stale.

Keep check mode read-only.

## Tests

Add tests for:

- Parsing a valid `feature.toml` with all required fields.
- Rejecting missing required fields.
- Rejecting malformed booleans for `implemented` and `validated`.
- Rejecting unknown `priority` values.
- Rejecting unknown `status` values.
- Rejecting an `id` that does not match its directory name.
- Rejecting dependency IDs that are not present.
- Skipping `_template` and other underscore-prefixed directories.
- Sorting features by ID.
- Rendering yes/no from booleans.
- Rendering stable markdown for a small two-feature fixture.
- `dashboard --check` succeeds when generated output matches.
- `dashboard --check` fails when generated output differs.
- `validate --feature` succeeds for known features and fails for unknown features.

If tests use temporary directories, keep fixture content small and avoid touching the real repository dashboard.

## Validation

This feature does not require RDKit, Biopython, or chemistry golden data.

Validation is local infrastructure validation:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature dashboard.feature-registry
```

Manual validation should include changing a temporary or fixture feature metadata file and confirming dashboard check mode reports staleness without writing the dashboard.

Leave `validated = false` unless there is documented manual validation or test evidence accepted as validation for this infrastructure feature.

## Risks

- Ad hoc TOML parsing may silently accept malformed metadata if the schema expands.
- Byte-for-byte dashboard checks can fail from newline normalization if file writing is inconsistent across platforms.
- A generated dashboard can drift if metadata updates skip `cargo xtask dashboard`.
- Dependency validation must avoid cycles only if future workflow needs topological ordering; current dashboard rendering does not require cycle detection.
- Treating feature metadata as implementation evidence would be misleading; metadata should only reflect real implementation and validation status.

## Edge Cases

- Feature directories without `feature.toml` should be ignored or reported consistently; the initial recommendation is to ignore non-feature directories and require metadata only for feature directories intended to be tracked.
- Empty `features/` should generate a dashboard with the header and no rows.
- Unknown boolean text should fail rather than render `unknown`.
- Descriptions containing markdown table separators should be excluded from the dashboard or escaped if later displayed.
- Duplicate feature IDs should fail even if they come from different directories.
- Missing `features/DASHBOARD.md` in check mode should fail with a clear stale or missing-file message.
- `depends_on = []` should parse as an empty dependency list.

## Explicitly Out of Scope

- Chemistry implementation or validation.
- Editing feature implementation statuses for unrelated features.
- Hand-editing generated dashboard rows.
- Building a web dashboard.
- Pulling feature data from GitHub, Slack, or external services.
- Automatically marking features implemented or validated.
