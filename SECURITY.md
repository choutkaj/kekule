# Dependency security and maintenance

The dependency audit checks both the workspace and fuzz lockfiles on dependency
pull requests, dependency changes merged to `main`, weekly, and on manual
dispatch. Vulnerability advisories fail the job. Informational advisories remain
visible in the logs and require maintenance review; they are not suppressed by
an ignore list. Scheduled audits only report CI status and do not create issues
or request repository write permissions.

Before releasing, run `cargo audit --file Cargo.lock` and
`cargo audit --file fuzz/Cargo.lock --no-fetch` against a freshly fetched RustSec
database, review informational warnings, and run the repository's documented
format, check, Clippy, test, documentation, and package gates. CI also checks the
declared Rust 1.89 minimum against all workspace targets and features.

## Tracked maintenance concern

As reviewed on 9 September 2026, the pinned DREIDING dependency tree brings in
`paste` 1.0.15 through its numerical libraries. The
[RustSec advisory RUSTSEC-2024-0436](https://rustsec.org/advisories/RUSTSEC-2024-0436.html)
reports that this compile-time procedural macro is unmaintained; it does not
report a vulnerability. It remains an accepted maintenance warning while those
upstream dependencies require it, and must be reconsidered on each DREIDING
dependency update or release review.

Prefer an upstream maintained replacement when the dependency chain supports
one. Changes to the pinned scientific dependencies must preserve their numerical
contracts and pass the relevant regression and scientific-reference checks.
Do not hide the warning, remove asserted scientific fields, or regenerate
expected outputs just to accommodate a dependency update. A new vulnerability
advisory requires its own assessment and is not covered by this maintenance
decision.
