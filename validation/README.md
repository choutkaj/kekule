# Validation

Validation has two top-level areas:

- `corpora/<corpus-id>/` owns the corpus descriptor, pinned source lock, local inputs, feature manifests, compressed reference goldens, and generated evidence.
- `reference/` owns RDKit and Biopython acquisition, corpus construction, and golden-generation tooling.

All registered corpus `data/` directories are generated locally and ignored. `sources.lock.json` pins every selected external record, URL, checksum, category, and generated pack. Locks, manifests, deterministic goldens, and generated status evidence are repository artifacts.

```bash
cargo xtask corpus check --corpus all
cargo xtask corpus check --corpus pubchem-1k --require-data
cargo xtask corpus check --corpus pdb-100 --require-data
cargo xtask validate --feature <feature-id> --corpus <corpus-id>
```

The normal required baselines are `pubchem-1k` for applicable small-molecule features and `pdb-100` for applicable macromolecular features. Broad deliberate runs add `pubchem-100k`, `enamine-diversity`, `pl-rex`, and `pdb-1000` where manifests exist.

To diagnose one declared fixture without changing generated status, use `--fixture`. If an implementation-semantic contract has changed and the new behavior has been independently reviewed, its snapshot can be accepted with `--accept-implementation-goldens`. Acceptance is restricted to one concrete feature/corpus using a `*-manual-semantic` reference and is deliberately separate from `--update`; generator-backed RDKit and Biopython goldens cannot be replaced this way.

The following command runs every required or implemented manifest-backed feature/corpus comparison across all registered corpora and updates recorded evidence plus the dashboard. Omitting `--corpus` has the same behavior. Every referenced local corpus dataset must already be present.

```bash
cargo xtask validate --feature all --corpus all --update
cargo xtask dashboard --check
```

The automatic worker count is capped at four to bound memory. Provisioned hosts may override it explicitly, for example with `--jobs 8`.

Molecular fixtures must be externally supplied. RDKit and Biopython are reference-only tools and must not become Rust runtime dependencies.