# Benchmark validation — 2026-09-16

The benchmark integrity checks pass. Scientific comparisons continue to expose
Kekule/reference differences; the goldens were not adjusted to make Kekule pass.

## Changes verified

- Streamed golden storage with schema, checksum, input-lock, reference-version,
  ordering and completeness checks; immutable generation and partial-run reports.
- Explicit applicability/error accounting, raw structural/numerical differences,
  declared reference precision, panic isolation and reference process deadlines.
- Correct CXSMILES handling, requested writer formats, the MOL model title
  contract, DSSP beta-partner identities and omega observations.
- Normal production imports, shared query inputs, fewer record wrappers and
  clones, and error propagation that retains causes.
- Removal of 4,751 retired golden files and unused indexes/helpers. The removed
  assets occupied 822,413,805 bytes including a few cached Python files; untracked
  bytecode was backed up before removal. Tracked assets remain in Git history.
- Correction of 2,881 Enamine CXSMILES cases across 17 features, preserving every
  unaffected golden line. Expanded DSSP references preserve every previously
  asserted successful value. Recompression saved 15,285,637 bytes overall.
- Provenance validation for optional trajectory comparisons, including input and
  artifact hashes, units, tolerance, versions and dimensions.

## Checks

| Command or check | Result |
| --- | --- |
| `cargo fmt --all -- --check` | Passed |
| `cargo check --workspace --all-targets --all-features --locked --offline` | Passed |
| `cargo +1.89.0 check --workspace --all-targets --all-features --locked --offline` with a separate target directory | Passed |
| `cargo clippy --workspace --all-targets --all-features --locked --offline -- -D warnings` | Passed |
| `cargo test --workspace --all-features --locked` | Passed, including documentation tests |
| `cargo test -p kekule-bench --locked` after final benchmark changes | 41 unit tests and 2 integration tests passed |
| `cargo test -p kekule-traj --example trajectory_periodic_reference --locked` | Provenance regression passed |
| `cargo doc --workspace --all-features --no-deps --locked --offline` with `RUSTDOCFLAGS=-D warnings` | Passed |
| `cargo package -p kekule --locked --allow-dirty --offline` | Packaged and verified |
| `cargo package -p PACKAGE --locked --allow-dirty --offline --list`, for kekule-potentials, kekule-traj and kekule-bench | Passed file-list checks |
| RDKit / Biopython / shared reference-runner Python tests | 16 / 6 / 8 passed |
| External input verification through `data.py` | All 1,494 files passed |
| All dataset/feature pairs, one source ID per dataset | All 125 golden files and 5,741,750 records passed storage/schema validation |
| `git diff --check` | Passed |

Python checks used the pinned reference environments. The Biopython check needed
its environment's library path and one BLAS thread on this Windows installation;
an initial unactivated invocation stalled in NumPy and was stopped, then rerun
successfully with that configuration. No expectations were changed for this.

Full companion package builds were not run: the repository's release workflow
uses their file-list checks until foundational dependency publication is ready.
The unpublished benchmark package was checked by file list. Linux fuzzing
and unrelated optional-feature matrices were not rerun; no runtime implementation
or feature configuration was changed. Full scientific evaluation of millions of
Kekule cases and the external PDB/XTC trajectory profile were not run; the former
was sampled, and no new external trajectory/reference export was prepared.

## Local-only data distribution

The 100 full-dataset golden archives stay on the local filesystem and are ignored
by Git. The 25 smoke archives total 141,733 bytes; they and all 125 provenance
manifests remain tracked. Full input datasets were already local and ignored.
Cargo package rules likewise exclude full input and golden payloads. A regression
stages a temporary checkout and verifies that current and future bulk datasets
remain local while smoke inputs, goldens, source locks and manifests are included.
The unpublished commit was rewritten to omit the full golden archives before
pushing, with the original commit preserved in a local backup ref. No shared
history was rewritten, and no external data upload was configured.

All 125 local archives were rehashed after the history/index changes and remain
byte-for-byte identical to their preceding inventory and tracked manifests. The
staged checkout was exported without bulk files: all 43 benchmark tests passed,
and all 25 smoke features reproduced the results below. This machine's inherited
Cargo Git-source patch was redirected to the exported source for that build.
The package file list contains 25 smoke archives and 125 manifests, with no bulk
inputs or golden archives. Formatting, workspace check, clippy and Rustdoc also
passed after the distribution changes.

## Dashboard validation

`dashboard.py` builds a self-contained HTML snapshot from explicit comparison
summaries and the tracked reference catalogue. It reads no golden payloads or
case records and exports no local paths or raw error messages. Runs stay separate;
selection scope, incomplete runs, stale provenance and unrun pairs remain visible.
Fixed snapshots stay under the ignored `target/` directory. Automatic run history
and its live page are kept under the ignored `benchmarks/runs/` directory.

Validation for this dashboard addition:

- All 21 standard-library Python regressions passed, including inconsistent
  counts, incomplete runs, overlapping error origins, changed provenance, duplicate
  imports/rows, safe embedded JSON, preservation of input reports, corpus sizes
  and formats, recorded-time ordering, legacy filenames, invalid-report isolation,
  empty/incomplete history, the combined agreement display and simplified tables. These tests
  are now included in CI without invoking scientific reference software.
- JavaScript syntax and in-memory DOM checks passed for 125 cells, filtering,
  empty search results, run switching, incomplete-run provenance, corpus rows,
  combined agreement bar widths
  and downloading the sanitized plotted data. The tracked Node regression also
  verifies live updates, unchanged-history selection, newest-run selection, empty
  history and retry after a missing local data file. This is not a visual browser test.
- Rust formatting, workspace check, clippy, Rustdoc, 45 benchmark tests and the
  benchmark package file-list check passed. The package still omits bulk data.
- Fresh full-smoke and one-ID-per-dataset comparisons completed and reproduced
  the totals below. Both correctly returned exit status 1 for scientific
  disagreements. The initial debug sample was stopped after preserving its
  incomplete checkpoint, then completed using the optimized executable.
- Automated visual preview could not run: the browser tool's URL security policy
  blocked local `file:` navigation. No alternate browser route was attempted.

The automatic page is `benchmarks/runs/index.html`; comparison commands refresh
its sanitized `dashboard-data.js` feed after success, disagreement or a handled
execution error. Each report records `started_at_unix_ms`. Updates are serialized
with an OS file lock and published by atomic replacement. Explicit `--output`
summaries join the history without copying their potentially large case files.
Missing Python leaves the report and scientific exit status intact and emits a
warning. Integration tests cover these cases against the checked-in smoke data.
The Git staging and package checks exclude run history and `.dashboard-python`.
An initially stale Cargo test executable created temporary missing-golden reports
in the default history. The benchmark build cache was cleared, all benchmark tests
were rebuilt and passed with isolated history, and only those test artifacts were
removed. Workspace check, clippy and Rustdoc were repeated after the clean build.

The existing Enamine comparison `run-32348-1789591874492654300.json` and the two
earlier dashboard summaries were imported without modifying the source reports;
the large case records remain at their original local paths. The earlier fixed
snapshot remains at `target/benchmark-dashboard/index.html`.
Full runtime/MSRV test matrices and runtime package builds were not repeated for
these benchmark reporting changes; their preceding audit results are above.
No runtime implementation, comparison contract, reference adapter or golden was
changed. Full-corpus Kekule evaluation and external trajectory profiling remain
outside this dashboard validation. No website publication was performed.

## Scientific results

The full smoke corpus measures 452 applicable cases across all 25 features:
**356 agree, 92 disagree, and 4 report Kekule errors**. Another 126 selected
ID/feature combinations lack the supplied format and are counted separately.
The four errors are two ordinary SMILES writer rejections of stereo input and
two unsupported SMARTS cases. No reference or schema errors occurred.

The one-ID-per-dataset pass completes all 125 feature/dataset pairs:
**92 agree, 35 disagree, and 1 error**, across 128 applicable cases, with 41
missing-format combinations. These small samples are validation of the harness,
not estimates of whole-corpus correctness. Both scientific commands correctly
return a failing exit status while recording `complete: true`.

Reports in the working directory:

- `target/benchmark-audit-final-smoke.json` and its `.cases.jsonl` companion.
- `target/benchmark-audit-fixed-all-sample.json` and its `.cases.jsonl` companion.
- `target/benchmark-audit-20260916/` contains verification logs and the exact
  reference-update/removal inventories.

Remaining disagreements require scientific interpretation. For example, the
smoke DSSP case matches all structural fields but retains numerical deviations
in independently calculated omega angles, including Biopython coordinate
precision effects. Descriptor masses retain differences between element/isotope
constants. These are not silently turned into agreement. Historical PubChem
selection bias, incomplete API/query coverage, and missing original generator
fingerprints for imported goldens remain explicit limitations in
[GUIDE.md](GUIDE.md) and [GOLDENS.md](GOLDENS.md).
