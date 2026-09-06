# Trajectory foundation validation

This optional scientific check uses externally supplied trajectories. MDTraj and
MDAnalysis are development references only, never Rust runtime dependencies. The
comparison is not a routine CI or release gate.

## Reproduction

Install `MDAnalysis==2.9.0` and `mdtraj==1.11.1.post1` in a separate Python
environment. Supply a matching PDB/XTC pair containing one connected molecule.
The PDB provides atom order and bond connectivity for numerical reconstruction;
this does not establish chemical bond orders. The script records input SHA-256
hashes and reference-library versions alongside the exported coordinates.

```text
python benchmarks/reference/trajectory/export_periodic.py --topology SYSTEM.pdb --trajectory INPUT.xtc --output target/periodic-reference
cargo run -p kekule-traj --release --example trajectory_periodic_reference -- target/periodic-reference/topology.txt INPUT.xtc target/periodic-reference
```

All coordinate components must agree within the predeclared tolerance of
`0.0001 nm`. Field presence and non-position metadata are checked exactly, as is
agreement between loaded and streaming transformations. A comparison failure must
be investigated; do not regenerate references or increase tolerances to conceal it.

## Recorded external check, 2026-09-06

The MDTraj lambda phage lysozyme example contains 51 frames, 2,504 atoms, and
2,529 bonds, with the protein crossing periodic boundaries. Inputs came from
MDTraj revision `c2e0f5a4dcc207dd8d2b63857584fd3b51928518`:

- [1am7_protein.pdb](https://raw.githubusercontent.com/mdtraj/mdtraj/c2e0f5a4dcc207dd8d2b63857584fd3b51928518/tests/data/1am7_protein.pdb), SHA-256 `cf3e4b921d5f42fbcbf9a0669f37c3fb8e012e05a73985d2e855c11dc1c01d88`.
- [1am7_uncorrected.xtc](https://raw.githubusercontent.com/mdtraj/mdtraj/c2e0f5a4dcc207dd8d2b63857584fd3b51928518/tests/data/1am7_uncorrected.xtc), SHA-256 `cdc31fb1909c99a2b5d9ba1f2c21c322e1d7cdc536955b52ded90614c719831f`.

The files and generated references are external scratch inputs, not bundled test
fixtures. This check used Python 3.12, NumPy 2.5.3, and the pinned versions above.

| Operation | Independent reference | Maximum component deviation (nm) |
|---|---|---:|
| Raw XTC coordinates | MDTraj reader | 0 |
| Make molecules whole | MDTraj `make_molecules_whole` | 2.655934565609641e-7 |
| Image molecules | MDTraj `image_molecules`, whole protein as explicit anchor | 8.282691874583747e-6 |
| Temporal unwrap | MDAnalysis sequential `NoJump` | 9.179115298962873e-7 |

Loaded/streaming equality and all non-position metadata checks passed. This profile
validates one real orthorhombic system; focused Rust regressions separately cover
skew, rotated, and partially periodic cells, rings, and ambiguous crossings.

The first run exposed an XTC decoder defect: a zero run-change flag incorrectly
reset the preceding run length. The independent xdrfile decoder retains it. A
hand-packed bitstream regression reproduces this case without relying on Kekule's
encoder. Fixing the state lifetime made the real trajectory decode exactly, without
relaxing compressed-data validation.

## Informational performance check

```text
cargo bench -p kekule-traj --bench frame_workflow -- SYSTEM.cif
```

The benchmark repeats the coordinates and annotations of an external mmCIF model
64 times to isolate container access and fitting costs. It does not synthesize a
molecular fixture or model physical trajectory sampling. It is informational,
without thresholds or CI gates.

On Windows with the release profile, using the supplied RCSB `1CRN.cif` (327
atoms), one million random stored-frame accesses averaged approximately 103 ns
before the refactor and 6.4 ns after it. Iteration averaged 101 ns and 0.39 ns per
frame respectively; the latter loop can be heavily optimized after validation is
removed from reads. Superposing all 64 frames averaged 3.04 ms before and 3.33 ms
after, providing no evidence of a fitting speedup in this small single-run check.
The baseline used a snapshot of the pre-refactor sources with the same benchmark.
Ordinary superposition also no longer allocates a discarded report vector.
