# Feature Registry and Dashboard

## Summary

Keep feature metadata as the machine-readable source of release truth and
generate a deterministic dashboard that displays optional benchmark
availability and observations without folding them into feature health.

## Behavior/API

- Exposes `cargo xtask features`.
- Exposes `cargo xtask dashboard` and `cargo xtask dashboard --check`.
- Exposes `cargo xtask skills --check` for repository-local workflow checks.
- Discovers benchmark availability only from manifests under
  `benchmarks/corpora/<corpus>/features/`.
- Allows features with no benchmark manifest; missing manifests cannot prevent
  feature listing or dashboard generation.

## Implementation Notes

- Feature metadata requires `id`, `title`, `area`, `domains`, `version`,
  `status`, `description`, and `depends_on`. The removed
  `validation_required` key has no replacement.
- Release status uses `planned`, `experimental`, `supported`, and `deprecated`.
  Removed global `implemented`/`validated` booleans and deprecated metadata
  keys are rejected.
- Dependencies must exist, be acyclic, and satisfy release-status
  compatibility.
- Every tracked feature directory includes `feature.md` with a required
  `Tests` section; `Benchmarks` is optional.
- The dashboard renders separate small-molecule, macromolecule, and
  infrastructure tables followed by the complete deterministic dependency DAG.
- Benchmark cells are neutral observations: `available`, `last match`, `last
  differences`, or `last error`, including timestamp and scope in tooltips.
  They never alter or summarize release status.
- Results are read from corpus-owned `results.toml` files only when a current
  manifest exists, preventing removed manifests from leaving ghost cells.
- `features/DASHBOARD.html` is generated and must not be edited by hand.

## Tests

- xtask tests cover feature schema and dependency checks, optional manifest
  discovery, neutral benchmark cells, legacy/current result rendering, stable
  generation, and platform-independent line endings.
- Routine checks run formatting, clippy, workspace tests, documentation,
  dashboard consistency, and skill consistency without external benchmark
  data.

## Out Of Scope

- Chemistry implementation or external-reference comparison.
- Pulling feature metadata from external services.
- Automatically promoting feature release statuses.
- Treating benchmark results as pass/fail repository health.

## Revision Notes

- v1-v13: Established the typed registry, release statuses, dependency graph,
  corpus matrix, and generated sortable dashboard.
- v14: Remove required benchmark metadata and health semantics; discover
  optional manifests directly and render neutral availability and historical
  match/difference/error observations.
