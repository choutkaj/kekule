# Security policy

## Reporting a vulnerability

Report vulnerabilities privately through GitHub's
[private vulnerability reporting](https://github.com/choutkaj/kekule/security/advisories/new)
(the **Report a vulnerability** button on the repository's Security tab). Do not
open a public issue, pull request or discussion for a suspected vulnerability.

Include the affected crate and version, the input or call sequence that
triggers the problem, and what happens (panic, hang, excessive memory, wrong
result). Attach reproducing inputs as files, but only structures you are
allowed to share.

You can expect an acknowledgement within a week. Fixes are released as a patch
of the latest minor version, and confirmed vulnerabilities are published as a
GitHub security advisory and submitted to the
[RustSec advisory database](https://rustsec.org/).

## Supported versions

Kekule is pre-1.0. Only the latest 0.x minor release of each crate receives
security fixes; older minor versions are not patched.

| Crate | Supported |
| --- | --- |
| `kekule`, `kekule-traj`, `kekule-openff`, `kekule-potentials` | Latest 0.x minor (currently 0.4) |
| `kekule-openff-ash` (model data only) | Latest release (currently 0.3) |

## Scope

The parsers are written to accept untrusted input: SMILES and CXSMILES, SMARTS,
Molfile V2000 and V3000, SDF, mmCIF, PDB, OFFXML and NAGL model bundles, and
XYZ, DCD, TRR and XTC trajectories. A panic, unbounded memory or time, or
memory unsafety on any such input is in scope, as is input within the
documented resource limits that still exhausts resources. All crates forbid
`unsafe` code.

Out of scope:

- Resource use that a caller asked for by raising or disabling a documented
  limit.
- Scientifically wrong results that need no adversarial input; report these as
  ordinary issues.
- Vulnerabilities in dependencies that do not affect Kekule's use of them;
  report those upstream. The dependency audit below tracks them.
- The benchmark tooling, fuzz harness and development scripts, which are not
  published.

## Dependency security and maintenance

The dependency audit checks both the workspace and fuzz lockfiles on dependency
pull requests, dependency changes merged to `main`, weekly, and on manual
dispatch. Vulnerability advisories fail the job. Informational advisories remain
visible in the logs and require maintenance review; they are not suppressed by
an ignore list. Scheduled audits only report CI status and do not create issues
or request repository write permissions.

Before releasing, run `cargo audit --file Cargo.lock` and
`cargo audit --file fuzz/Cargo.lock --no-fetch` against a freshly fetched RustSec
database, review informational warnings, and run the repository's documented
format, check, Clippy, test, documentation, and package gates (see
`RELEASING.md`). CI also checks the declared Rust 1.89 minimum against all
workspace targets and features, and runs `cargo deny check` against the
licence, advisory, ban and source policy in `deny.toml`.

### Tracked maintenance concerns

None are currently tracked. The unmaintained `paste` warning
([RUSTSEC-2024-0436](https://rustsec.org/advisories/RUSTSEC-2024-0436.html))
reviewed on 9 September 2026 came only from the DREIDING dependency tree, which
was removed together with the DREIDING potential.

When a new informational advisory is accepted, prefer an upstream maintained replacement when the dependency chain supports
one. Changes to the pinned scientific dependencies must preserve their numerical
contracts and pass the relevant regression and scientific-reference checks.
Do not hide the warning, remove asserted scientific fields, or regenerate
expected outputs just to accommodate a dependency update. A new vulnerability
advisory requires its own assessment and is not covered by this maintenance
decision.
