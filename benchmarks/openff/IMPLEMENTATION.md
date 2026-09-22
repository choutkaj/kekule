# Native OpenFF implementation checks

The [crate contract](../../crates/kekule-openff/CONTRACT.md) defines the supported
runtime behavior. This directory contains optional reference tools, never Rust
runtime dependencies. The original audit/report bytes remain unchanged; the
user-approved charge comparison cutoff is now **0.00005 e**.

## Reproduce

Use the versions in `environment.yml`. Generate the data-only model and build
the native observer:

```text
micromamba run -p target/openff-reference python benchmarks/openff/export_model.py target/openff-ash
cargo build -p kekule-bench --bin openff_parameterize --release --locked
python benchmarks/openff/compare_native.py --binary target/release/openff_parameterize.exe --model target/openff-ash --output target/openff-native.json
micromamba run -p target/openff-reference python benchmarks/openff/reference_parameters.py --binary target/release/openff_parameterize.exe --model target/openff-ash --output target/openff-parameters.json
python benchmarks/openff/validate_lookup.py --binary target/release/openff_parameterize.exe --model target/openff-ash --output target/openff-lookup.json
micromamba run -p target/openff-reference python benchmarks/openff/classify_lookup.py target/openff-lookup.json target/openff-lookup-classified.json
```

Omit `.exe` outside Windows. `compare_native.py` and `validate_lookup.py` use only
standard-library Python. The latter deliberately exits nonzero when any stored
key disagrees, even if the current official toolkit also rejects that entry.
`classify_lookup.py` records official outcomes without editing or excluding
native failures. No lookup table, original charge array or original golden was
changed to obtain agreement.

The extended cases come from the original 23 externally sourced audit inputs
and ten additional PubChem responses. `supplementary-sources.lock.json` pins
the response URLs and bytes. `fetch_cases.py` is an explicit acquisition tool;
ordinary tests never fetch or refresh fixtures. The supplement adds 3-/4-membered
rings, cubane, norbornane, fused rings, cis/trans alkenes, glucose, cysteine and
dimethyl sulfoxide. Each is also tested in reversed input atom order.

`native-reference.json.gz` archives all 66 live observations, including complete
native output, official feature/inference arrays, official parameter values,
parameter multiplicities and improper energy comparisons at three deterministic
random coordinate sets. Improper tuples are compared with the center first,
retaining all three trefoil terms; opposite outer-atom orientation is checked
by evaluating the emitted terms, not by dropping duplicates.

`lookup-identity.json.gz` archives the exhaustive stored-key scan and official
classification. `implementation-reports.lock.json` pins both reports and the
source files used for the implementation run; these hashes record that run,
not a requirement that future source versions stay unchanged.

## Recorded outcomes (2026-09-22, Windows CPU)

* Original native comparison: 23/23 pass, exact labels and InChI identifiers;
  maximum charge difference `1.22e-7 e`.
* Extended reference: 66/66 pass. Maximum forced neural inference difference
  `2.54e-7 e`; maximum complete system charge difference `3.00e-5 e`. All feature
  columns agree within `1e-6`; numeric parameters and improper energies agree
  within `atol=1e-10, rtol=1e-12`.
* Stored lookup identities: 13,235/13,944 exact matches. The 709 recorded
  disagreements comprise 689 inputs rejected by the current official toolkit,
  19 keys on which current native and official identifiers agree but differ
  from the stored key, and one native InChI rejection accepted by the toolkit:
  entry 4548, `[H:5][O:3][N+:2]#[S:1][H:4]`. The last is outside Ash's inference
  domain due to its sulfur triple bond. This is a retained compatibility limit.
* Exporting the checkpoint twice produced byte-identical JSON and weight files.

These checks establish a working port on this corpus, not exhaustive SMIRNOFF,
normalization, resonance, atom-remapping or large-system performance coverage.
Reference mismatches uncovered and fixed tetrahedral parity, improper tuple
ordering and neutral valence-five N/P input normalization. Focused regressions
and the full original feature/label corpus now run in native tests.

The later [110-molecule robustness pass](ROBUSTNESS.md) adds protein-sized cases,
frozen independent OpenMM energies, a `potentials` benchmark consumer, and
regressions for isotope charge lookup, bounded model loading, and SMARTS scaling.
Its raw baseline failures and remaining identifier limits are preserved separately.

## Maintenance validation

Passed during implementation:

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo +1.89.0 check --workspace --all-targets --all-features --locked
cargo doc --workspace --all-features --no-deps --locked
cargo package -p kekule-openff --allow-dirty --locked --list
cargo package -p kekule-bench --allow-dirty --locked --list
python -m unittest discover -s benchmarks/openff -p test_audit.py
git diff --check
```

Workspace tests included 1,397 passing tests/doctests and the explicitly ignored
model-dependent test; that model test was separately run successfully with
`KEKULE_OPENFF_MODEL` pointing to the absolute bundle path. Two additional
scaling/unit regressions were subsequently checked with the complete new-crate
test suite and Clippy. Offline audit/provenance checks contain 14 passing tests.
Rustdoc ran with `RUSTDOCFLAGS=-D warnings`; workspace tests used
`RUST_TEST_THREADS=1` on Windows. The runnable ethanol example was exercised.

Not run: a separate workspace `--doc` test command (already included in workspace
tests); full companion package verification/publication (the workspace's core
0.2.1 dependency is not yet published, so CI also checks companion file lists);
unchanged core/other companion package builds; unrelated no-default-feature
potential checks, fuzz targets and external benchmark suites. Linux execution
and performance benchmarks were not performed locally. CI includes the new
crate in Linux/MSRV/Windows workspace checks and its package/license file checks.
