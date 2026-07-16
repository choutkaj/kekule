# Parser Fuzzing

Parser fuzz targets live in the standalone `fuzz/` package so fuzz-only dependencies do not enter the runtime workspace graph.

Run a bounded smoke campaign from the repository root:

```bash
cargo install cargo-fuzz --locked
cargo +nightly fuzz run mol_v2000 -- -runs=256 -max_len=4096 -seed=1
cargo +nightly fuzz run mol_v3000 -- -runs=256 -max_len=4096 -seed=2
cargo +nightly fuzz run sdf_v2000 -- -runs=256 -max_len=4096 -seed=3
cargo +nightly fuzz run smiles -- -runs=256 -max_len=4096 -seed=4
cargo +nightly fuzz run mmcif -- -runs=256 -max_len=4096 -seed=5
```

Longer manual campaigns can omit `-runs` and raise `-max_len`. Seed inputs are committed under `fuzz/corpus/<target>/`. Crashing inputs are written under ignored `fuzz/artifacts/`; preserve and add any reproducer as a focused regression test before fixing it.

The `Scheduled parser fuzzing` workflow runs each target for up to 30 minutes
every Monday and uploads crash artifacts on failure. Download an artifact
before the workflow retention period expires, reproduce it locally with
`cargo +nightly fuzz run <target> <artifact>`, and commit a minimized input only
when its redistribution terms permit it. Never commit inputs containing
secrets, private structures, or unreviewed third-party data.

Fuzzing demonstrates the explored executions, not parser correctness or
unbounded-input safety. The targets cap generated input length; the public
Molfile, SDF, SMILES, and mmCIF parsers also enforce their documented runtime
limits.
