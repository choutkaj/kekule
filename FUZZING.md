# Parser Fuzzing

Parser fuzz targets live in the standalone `fuzz/` package so fuzz-only dependencies do not enter the runtime workspace graph.

Run a bounded smoke campaign from the repository root:

```bash
cargo install cargo-fuzz --locked
cargo +nightly fuzz run mol_v2000 -- -runs=256 -max_len=4096 -seed=1
cargo +nightly fuzz run sdf_v2000 -- -runs=256 -max_len=4096 -seed=2
cargo +nightly fuzz run smiles -- -runs=256 -max_len=4096 -seed=3
cargo +nightly fuzz run mmcif -- -runs=256 -max_len=4096 -seed=4
```

Longer manual campaigns can omit `-runs` and raise `-max_len`. Seed inputs are committed under `fuzz/corpus/<target>/`. Crashing inputs are written under ignored `fuzz/artifacts/`; preserve and add any reproducer as a focused regression test before fixing it.
