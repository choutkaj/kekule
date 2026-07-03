# Validation

Validation has two top-level areas:

- `corpora/<corpus-id>/` owns the corpus descriptor, pinned source lock, local inputs,
  feature manifests, compressed reference goldens, and generated evidence.
- `reference/` owns RDKit and Biopython acquisition and golden-generation tooling.

Large corpus `data/` directories are ignored. `sources.lock.json` pins every selected external
record, URL, checksum, category, and generated pack. The `smoke` corpus is committed in full.

```bash
cargo xtask corpus check --corpus all
cargo xtask corpus check --corpus all --require-data
cargo xtask validate --feature all --corpus smoke
```

The following command will run a full validation of all features on all corpora, and update the feature dashboard.
```bash
cargo xtask validate --feature all --corpus all --update
```

After it finishes, a good sanity check is:
```bash
cargo xtask dashboard --check
```

By default, the xtask validate command will use all available processors. You can override it with for example:
```bash
cargo xtask validate --feature all --corpus all --update --jobs 8
```

Molecular fixtures must be externally supplied. RDKit and Biopython are reference-only tools and
must not become Rust runtime dependencies.
