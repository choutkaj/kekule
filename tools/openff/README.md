# OpenFF data provenance

These scripts produced the OpenFF data shipped in `kekule-openff` and
`kekule-openff-ash`. Python, PyTorch and the OpenFF packages are exporter
dependencies only, never Rust runtime dependencies.

| File | Purpose |
| --- | --- |
| `environment.yml` | Pinned conda-forge environment the exports were made in |
| `export_model.py` | Converts a checksum-verified NAGL checkpoint to JSON metadata and float32 weights |
| `package_ash.py` | Stores the exported Ash weights as zlib-compressed byte planes in `crates/kekule-openff-ash/data` |
| `export_normalizations.py` | Extracts `crates/kekule-openff/data/normalizations.json` from the pinned NAGL source in `upstream/` |
| `export_test_bundles.py` | Exports the two bundles that the ignored `kekule-openff` model tests read |
| `sources.lock.json` | URL, revision and SHA-256 of each upstream source and shipped licence |
| `models.lock.json` | Checkpoint and bundle fingerprints of both NAGL models |
| `test_provenance.py` | Offline check that the shipped data still matches these records |

`upstream/openff.py` is an unmodified file from
[openff-nagl](https://github.com/openforcefield/openff-nagl) (MIT,
`crates/kekule-openff/data/LICENSE-nagl`).

## Reproduce

Run from the repository root:

```sh
micromamba create -y -p target/openff-reference -f tools/openff/environment.yml
micromamba run -p target/openff-reference python tools/openff/export_model.py target/ash-export
python tools/openff/package_ash.py target/ash-export crates/kekule-openff-ash/data
micromamba run -p target/openff-reference python tools/openff/export_normalizations.py
python -m unittest discover -s tools/openff -p "test_*.py"
```

The ignored two-model tests need both exported bundles:

```sh
micromamba run -p target/openff-reference python tools/openff/export_test_bundles.py target/openff-models
KEKULE_OPENFF_MODELS=target/openff-models cargo test -p kekule-openff --test models -- --ignored
```

The numerical validation against OpenFF and OpenMM that used to live in
`benchmarks/openff` was retired; its report and scripts remain in git history.
