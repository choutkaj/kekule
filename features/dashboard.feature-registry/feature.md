# Feature Registry and Dashboard

## Summary

Keep feature metadata as the machine-readable source of truth and generate a deterministic dashboard from it.

## Behavior/API

- Exposes `cargo xtask features`.
- Exposes `cargo xtask dashboard` and `cargo xtask dashboard --check`.
- Exposes corpus-aware validation with feature and corpus `all` selectors.
- Exposes `cargo xtask skills --check` for repo-local feature workflow skill checks.

## Implementation Notes

- Feature metadata schema v5 requires `id`, `title`, `area`, `domains`,
  `version`, `status`, `description`, `depends_on`, and
  `validation_required`.
- `status` uses the release vocabulary `planned`, `experimental`, `supported`,
  and `deprecated`. The removed `implemented` and global `validated` booleans,
  plus deprecated `priority` and `last_ai_review`, are rejected.
- `depends_on` declares semantic feature prerequisites. Dependency IDs must
  exist and form a directed acyclic graph; duplicate dependencies,
  self-dependencies, and cycles are rejected.
- Status compatibility is enforced across graph edges: `supported` features
  require only `supported` prerequisites; `experimental` features require
  `experimental` or `supported` prerequisites; `deprecated` features may use
  any implemented prerequisite; `planned` features may name any registered
  prerequisite.
- Each tracked feature directory must include `feature.md`.
- The dashboard renders separate generated HTML tables for small molecules,
  macromolecules, and infrastructure. Shared chemistry foundations appear in
  both chemistry tables.
- The dashboard procedurally renders the complete dependency DAG as a
  deterministic SVG. Arrows run from prerequisites to dependents, and columns
  are assigned from dependency depth. The feature tables precede the graph so
  the primary release and validation overview appears first.
- Small-molecule and PDB-derived corpora are selected from typed corpus
  `kind` metadata. Unregistered internal smoke sets never become dashboard columns.
- Each chemistry section displays the exact reference codebase version found
  in its manifests, plus supplemental semantic reference labels where
  applicable. Individual parity-cell tooltips identify their exact reference.
- Dashboard cells report recorded per-corpus parity evidence only when a current feature manifest exists for that corpus, so removed manifests cannot leave ghost status cells.
- `cargo xtask validate` is the authority for checking parity against the current checkout; `cargo xtask features` lists feature metadata without deriving a global validation result.
- Per-feature evidence is read from each registered corpus-owned `status.toml`; status entries without a current manifest are ignored and pruned on the next update.
- Manifest-backed local-only corpus cells may represent required baseline evidence or optional broad evidence according to each feature's `validation_required` metadata.
- Release statuses render as labeled, color-coded pills. Compact failure counts
  are reserved for recorded fixture-level validation failures.
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
- Automatically assigning or promoting feature release statuses.

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
- v10: Split the dashboard into small-molecule, macromolecule, and
  infrastructure tables; add explicit feature-domain metadata and display
  versioned external reference information.
- v11: Replace the implementation boolean with explicit release statuses,
  enforce the feature dependency DAG and status compatibility, and render the
  generated graph in the HTML dashboard.
- v12: Render only registered public corpora, allow required local-only baselines, and ignore or prune status evidence that has no current feature manifest.
- v13: Place the generated feature dependency graph after all feature tables.
