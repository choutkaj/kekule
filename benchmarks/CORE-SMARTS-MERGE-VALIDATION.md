# Core SMARTS integration validation

Validated on 2026-09-20 after integrating main at `e4f8054b` (PR #183).
This record supersedes the earlier benchmark-only implementation checks.
`README.md` is unchanged.

## Repository checks

The following passed on the combined tree. Rust commands used
`--target-dir target/smarts` unless noted; MSRV used `target/mdl-msrv`.

- `cargo fmt --all -- --check`
- `cargo check --workspace --all-targets --all-features --locked`
- `cargo +1.89.0 check --workspace --all-targets --all-features --locked -j2`
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`
- `cargo test --workspace --all-features --locked -j2`
- `cargo test --workspace --all-features --doc --locked`
- `cargo doc --workspace --all-features --no-deps --locked`, with `RUSTDOCFLAGS=-D warnings`
- `cargo test -p kekule-potentials --no-default-features --locked`
- `cargo doc -p kekule-potentials --no-default-features --no-deps --locked`, with warnings denied
- `cargo package -p kekule --allow-dirty --locked`, including package compilation
- Package file-list checks for `kekule-potentials`, `kekule-traj` and `kekule-bench`, and packaged license-copy comparisons
- `cargo check --manifest-path fuzz/Cargo.toml --bins --locked -j2` under Ubuntu WSL, with `CARGO_TARGET_DIR=target/linux-fuzz` (resolved to the workspace's absolute path)
- `python -m unittest discover -s benchmarks -p 'test_*.py'` (44 tests)
- `python -m unittest discover -s benchmarks/reference/rdkit -p 'test_*.py'` (47 tests, pinned RDKit environment)
- `python -m unittest discover -s benchmarks/reference -p test_runner.py` (8 tests)
- `node benchmarks/test_dashboard.cjs`

The first Windows test build exceeded the MSVC PDB limit. Subsequent builds
disabled development/test debug symbols and incremental compilation and used
two build jobs. These were environment overrides, not repository settings.
The full suite exposed an obsolete merge restriction on negated atom stereo;
removing that restriction restored the existing Boolean-stereo regression.
The final full suite, clippy and core package verification passed after that fix.

Full companion-crate package verification was not run: repository CI defers it
until the foundational core release is available on crates.io. Package file
lists and licenses were checked instead. Long fuzz campaigns were not rerun;
all fuzz targets compiled, and earlier bounded SMARTS fuzz evidence remains
documented. The entire historical MDL corpus was not rerun for this merge.

## External comparisons

The reference environment pins RDKit **2026.03.3** and Python **3.11**; Rust has
no runtime dependency on RDKit or OpenFF. These development comparisons are
not routine CI gates.

- `query.smarts`, all seven registered corpora: **150,765 exact agreements,
  zero disagreements, one resource error, 1,226 not-applicable cases**.
  CID 447702 retains its bounded-search error in both engines. The 518 RDKit
  query patterns contribute 8,288 query/target pairs, 136 positive pairs and
  265 mappings. This limited target panel does not establish universal
  conformance. See [results](reports/query-smarts-merged-results.json) and
  [SMARTS validation](QUERY-SMARTS-VALIDATION.md).
- `algo.aromaticity.mdl`, the newly merged RDKit structure corpus: **46/50
  exact agreements, zero flag disagreements, four import errors**. All inputs
  and errors remain recorded: one unsupported V3000 CFG and three generic
  R/R1 labels. See [results](reports/mdl-structures-results.json),
  [errors](reports/mdl-structure-errors.json) and [MDL validation](MDL-VALIDATION.md).
  The earlier full MDL report still contains an unresolved discrepancy and is
  retained as historical evidence.

SMARTS/MDL compatibility metadata was explicitly migrated to main's contract
3 because these observations contain none of its changed radical/spin/DSSP
fields. All twelve preexisting payload digests, generator identities, source
locks, reference identities and case counts were audited unchanged. Manifests
record their old contracts in `origin`; historical reports retain their old
identities. Only the newly added MDL structure baseline was freshly generated.
