# Release Readiness

Status inspected on June 21, 2026.

## Verified Repository State

- The repository is private and its default branch is `main`.
- GitHub Actions workflows `CI` and `Reference validation` are active.
- Default workflow token permissions are read-only and workflows cannot
  approve pull requests.
- The ordinary CI workflow runs formatting, all-target checks, clippy, unit and
  doc tests, rustdoc, dashboard/skill checks, tiny corpus integrity, tiny
  reference comparison, and bounded parser fuzz smoke.
- A weekly workflow runs longer bounded campaigns for all four parser fuzz
  targets.
- The only repository collaborator returned by GitHub is the owner.

## Owner-Managed Policy Blockers

- `main` is not protected. Required status checks, pull-request review, and
  path-sensitive controls for generators, source locks, goldens, and workflows
  are not enforced. Tracked in
  [#13](https://github.com/choutkaj/molecules/issues/13).
- The private repository plan rejected branch-protection and ruleset API reads
  with HTTP 403. Protection cannot be claimed until the owner enables a plan or
  repository visibility that supports it and the settings are inspected again.
- No license is selected. The owner must choose one before any public release,
  after which canonical license text, Cargo/README metadata, dependency
  compatibility, and corpus/golden redistribution terms must be reviewed.
  Tracked in [#15](https://github.com/choutkaj/molecules/issues/15).

## Data And Scheduled Validation Blockers

Tiny corpus fixtures are committed and run in ordinary CI. PubChem 100/1000
and PDB 10/100 passed locally against pinned locks and committed goldens, but
their fixture data is intentionally ignored.

The PDB locks contain directly retrievable per-file URLs. PubChem locks record
hashes for extracted records whose source URLs point to very large bulk
archives, so a clean runner currently has no bounded pinned-artifact retrieval
path. PL-REX and Enamine descriptors remain `ready = false` and have no feature
manifests. A scheduled large-corpus workflow would therefore be misleading
until an owner-approved artifact store or deterministic bounded fetch format is
available. Tracked in [#14](https://github.com/choutkaj/molecules/issues/14).

Release requires one of:

1. A provenance-pinned artifact store that reconstructs every ignored fixture
   and passes `cargo xtask corpus check --corpus all --require-data`.
2. Committed redistributable fixtures after license and size review.
3. Explicitly narrowing required release corpora in feature metadata, with an
   architecture and validation-contract review.

## Audit Closure

The original high-severity implementation findings are covered by Stages 1
through 8 and their regression/reference checks. Remaining release blockers
are governance, license, and reproducible large-corpus distribution rather
than silent chemistry or parser deferrals.
