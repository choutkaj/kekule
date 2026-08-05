# Fuzzing

Fuzz targets live in the standalone `fuzz/` package so fuzz-only dependencies
do not enter the runtime workspace graph.

List the registered targets and run a bounded campaign from the repository root:

```bash
cargo install cargo-fuzz --locked
cargo +nightly fuzz list
cargo +nightly fuzz run <target> -- -runs=256 -max_len=4096 -seed=1
```

The input-parser targets are `mol_v2000`, `mol_v3000`, `sdf_v2000`, `smiles`,
and `mmcif`. The trajectory targets are `trajectory_detection`,
`trajectory_xyz`, `trajectory_dcd`, `trajectory_trr`, and `trajectory_xtc`.

Longer manual campaigns can omit `-runs` and raise `-max_len`. Seed inputs are
committed under `fuzz/corpus/<target>/`. Crashing inputs are written under
ignored `fuzz/artifacts/`; preserve and add any reproducer as a focused
regression test before fixing it.

The manual bounded-fuzzing workflow runs one selected target for five minutes
and uploads crash artifacts on failure. The manually dispatched CI workflow
runs a 256-input smoke campaign for every target reported by
`cargo +nightly fuzz list`. Download an artifact before the workflow retention
period expires, reproduce it locally with `cargo +nightly fuzz run <target>
<artifact>`, and commit a minimized input only when its redistribution terms
permit it. Never commit inputs containing secrets, private structures, or
unreviewed third-party data.

Fuzzing demonstrates the explored executions, not parser correctness or
unbounded-input safety. The targets cap generated input length; the public
Molfile, SDF, SMILES, mmCIF, and trajectory readers also enforce their
documented runtime limits.
