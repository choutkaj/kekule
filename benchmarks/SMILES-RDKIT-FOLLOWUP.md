# Four canonical SMILES alignment changes, 2026-09-28

This follows [the initial SMILES alignment](SMILES-VALIDATION.md). The target
remains RDKit 2026.03.3 with the existing exact CXSMILES contract. Neither the
independent reference observations nor the identity comparison was changed.

## Implemented rules

1. Formal charge uses RDKit's unsigned 32-bit comparison ordering. For example,
   acetate now writes `CC(=O)[O-]` rather than `CC([O-])=O`.
2. Canonical ring closures follow DFS back-edge discovery order. Closure bond
   symbols are emitted at the closing endpoint, and a number closed on an atom
   cannot be reused to open another ring on that same atom. Tetrahedral spelling
   uses the actual closure order. Ordinary-mode ring-label capacity is retained.
3. Refinement includes tetrahedral configuration when carrier classes are
   distinct, enhanced-group kind and member classes. Complete bounded graph and
   stereo certificates still resolve remaining ties. AND/OR groups normalize
   the lowest-ranked member in its emitted carrier frame to `@`, rather than
   normalizing whichever member traversal reaches first. CX emits OR groups
   before AND groups, with independently canonical group numbers. Legacy
   relative groups keep their existing first-emitted-member policy and `r`
   semantics; they are not silently converted into AND or absolute groups.
4. All SMILES modes emit redundant `^1`, `^2` and `^5` CX radical annotations
   for one, two and three electrons. Atom indices are assigned after traversal
   and offset after component ordering. Unspecified spin remains unspecified;
   explicit spin still fails instead of being discarded. Higher occupancies
   retain their bracket encoding because these CX codes do not represent them.

Refinement constructs each enhanced-group signature once per round, rather
than copying a whole group into every member's key. Its additional stereo and
group visits count against the existing refinement budget. No RDKit runtime
dependency or changes to the public molecular symmetry algorithm were added.

## Exact strings on the same full corpus

| Corpus | Cases | Initial alignment | Four-rule follow-up | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 83,614 | 96,480 | 9 |
| Enamine | 50,240 | 44,663 | 48,190 | 0 |
| Smoke | 8 | 7 | 7 | 0 |
| Total | 150,248 | 128,284 (85.38%) | 144,677 (96.29%) | 9 |

The same 2,008-case deterministic sample improves from 1,735/2,008 (86.40%)
to 1,932/2,008 (96.22%): PubChem 965/1,000, Enamine 960/1,000, smoke 7/8.
The nine full-corpus errors remain shared valence rejections. There are no
native-only failures, including the harness's numbering and fixed-point checks.

This improvement is not monotonic case by case: 16,990 strings change from
different to equal, 597 from equal to different, 127,687 stay equal, and 4,965
stay different. The independent expected observations were checked unchanged
for every paired case. All remaining mismatches remain failures in the report.

## Identity audit

The ordinary identity harness still agrees on 2,007/2,008 sampled canonical
outputs, retaining the same `Z7371431020` legacy-`r` mismatch.

An additional paired audit reread all 150,239 successfully emitted full-corpus
records through the existing RDKit `writer_value('io.smiles.canonical', ...)`
adapter and compared them with the independent canonical CX expectations,
including titles. It checked both the preceding and new implementations; an
identical emitted observation reused the same deterministic reader result.
This audit does not substitute normalized strings into the exact-text metric.

Both implementations have 150,175 identity agreements, 64 identity mismatches,
and nine upstream valence errors. Of the successful records, 150,172 retain
agreement, 61 retain disagreement, three gain agreement and three lose it.
The three new mismatches are Enamine `Z9903282369`, `Z9903282301` and
`Z8158346544`; all involve legacy `r` plus explicit enhanced groups. Native
relative relationships are retained, but RDKit's different interpretation
makes the chosen absolute representative observable to its identity reader.

Of the 64 identity mismatches, 62 involve legacy `r`. The other two, Enamine
`Z9082652630` and `Z9237738672`, are pre-existing enhanced spiro-stereo cases.
For `Z9082652630`, the new raw string exactly matches RDKit's stored string,
but rereading that string through RDKit changes its canonical stereo spelling.
Thus exact-text success and this identity check can disagree even when the
emitted text equals the independent reference byte for byte. These cases are
retained, not reclassified or hidden.

## Remaining text differences

A lexical breakdown of the 5,562 mismatches finds 2,456 other body differences,
2,334 ring-label differences, 675 stereo-mark differences, 40 combined ring and
stereo differences, and 57 CX-only differences. These are diagnostics, not
proof of a common cause; body categories may additionally differ in CX fields.
Remaining ranking and tie-breaking choices can also change which ring labels
appear even where the atom and branch token sequence looks the same.

## Validation and reproduction

Reports and detailed observations are local under `target/smiles-four-rules/`:
`final-full.json`, `final-sample.json`, `final-identity-sample.json`,
`final-full.identity-audit.json`, and `paired-text-analysis.json`. The identity
audit's script is retained there as `audit_identity.py`. Bulk data and reports
remain untracked, as required by the benchmark policy.

```text
cargo benchmark --feature io.smiles.text.canonical --dataset all --jobs 4
cargo benchmark --feature io.smiles.text.canonical --dataset all --limit 1000 --jobs 4
cargo benchmark --feature io.smiles.canonical --dataset all --limit 1000 --jobs 4 --goldens IDENTITY_REFERENCES --writer-python PATH
```

Focused regressions assert RDKit charge, fused-ring, symmetric-stereo,
enhanced-group and radical strings; they also assert fixed points. Randomized
atom/bond/endpoint permutations now retain and exercise enhanced groups.
Radical regressions check electron-count/spin preservation and rejection of
explicit spin. Workspace tests pass: 1,477 passed, three pre-existing ignored.

Formatting, workspace/all-target/all-feature checking, warnings-denied Clippy,
warnings-denied documentation, Rust 1.89 checking, fuzz-target compilation,
potentials tests/documentation without default features, foundational package
verification and companion package-file listings were checked. The 84 Python
reference/dashboard tests and JavaScript dashboard tests passed. Final command
logs are retained beside the reports. An intermediate concurrent build caused
stale-artifact doctest failures; the final serialized validation run supersedes it.

Linux CI and fuzz execution were not run on this Windows host. Full companion
package verification requires the foundational crate published on crates.io;
the package file lists are checked instead, as in repository CI. No README,
benchmark assertions or independent reference payloads were changed by this
follow-up.

## Cleanup audit

The subsequent cleanup consolidates the implementation without changing its
canonical ordering policies. Runtime SMILES code contains no corpus IDs,
molecule lookup tables, expected-output substitutions or benchmark-result
exceptions. Exact comparison still compares the emitted text directly, and
identity remains a separate measurement.

Tree planning and emission share one neighbor-ordering function and one ring
membership calculation. Stereo representative selection and emission share
the same child/carrier ordering. The group-normalization prepass runs only
when AND/OR reference centers exist. Canonical export relies on the published
nonempty, connected `Molecule` invariant; component ordering remains owned by
topology export. The unused component-search helper and an unnecessary atom
clone were removed.

Initial writer classes are named separately from chemical symmetry classes.
Stereo refinement uses named relationship variants instead of numeric category
tests. These remain serialization preferences; the complete graph/stereo
certificate resolves ties. Legacy relative relationships retain their distinct
semantics. Internal live-ID lookups now expose invariant violations instead of
silently omitting data.

The unchanged quadratic admission guard is documented as an input complexity
score, `2*n*(n+2*m)`, independently of metered search/refinement work. Its error
diagnostic now says "input complexity" instead of the obsolete "candidate
traversal" wording. A boundary regression checks admission/rejection, and
the stereo permutation regression also exercises deleted atom and bond slots.

The post-cleanup full comparison has zero changed actual or expected
observations, statuses or differences across 151,992 paired harness rows:
150,248 applicable cases plus 1,744 not-applicable rows. Agreement remains
144,677/150,248 (96.29%), with the same 5,562 differences and nine shared
errors. The comparison exits nonzero for these retained mismatches; this is
not a new cleanup failure. Because every emitted observation is unchanged,
the prior full identity audit also applies without rerunning the external
reader. Reports, the strict pairing script and validation logs are under
`target/smiles-cleanup/` (`full.json`, `paired-output-check.json`,
`compare_outputs.py` and `checks.json`).

Post-cleanup validation passed: formatting; workspace/all-target/all-feature
check and warnings-denied Clippy; 1,478 workspace tests including doctests
(three pre-existing ignored); warnings-denied documentation; Rust 1.89 checks;
fuzz-target compilation; no-default-feature potentials tests and documentation;
core package build verification; companion and benchmark package listings;
84 Python reference/dashboard tests; JavaScript dashboard tests; and
`git diff --check`. Exact command arguments are in the local `validate.ps1`;
the separate workspace check is logged in `check.log`.

Linux CI and `cargo +nightly fuzz run` were not run on this Windows host.
Full `cargo package -p kekule-potentials`, `-p kekule-traj` and
`-p kekule-openff` build verification remains deferred until the foundational
crate is published; their `--list` checks follow repository CI. No additional
identity-reader run was needed because the complete paired emitted
observations are identical. README and independent reference payloads were
not modified by this cleanup.
