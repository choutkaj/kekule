# Validation

Validation has two top-level areas:

- `corpora/<corpus-id>/` owns the corpus descriptor, pinned source lock, local inputs,
  feature manifests, compressed reference goldens, and generated evidence.
- `reference/` owns RDKit and Biopython acquisition and golden-generation tooling.

Large corpus `data/` directories are ignored. `sources.lock.json` pins every selected external
record, URL, checksum, category, and generated pack. The 20-case `smoke` corpus is committed in
full, and even the default corpus check verifies all of its fixture bytes.

```bash
cargo xtask corpus check --corpus all
cargo xtask corpus check --corpus all --require-data
cargo xtask validate --feature all --corpus smoke
cargo xtask validate --feature all
```

To diagnose one declared fixture without changing generated status, use
`--fixture`. If an implementation-semantic contract has changed and the new
behavior has been independently reviewed, its snapshot can be accepted with
`--accept-implementation-goldens`. Acceptance is restricted to one concrete
feature/corpus using a `*-manual-semantic` reference and is deliberately
separate from `--update`; generator-backed RDKit and Biopython goldens cannot be
replaced this way.

The following command runs every required or implemented manifest-backed
feature/corpus parity check across all registered corpora, including the ignored
local-only large corpora, and updates the recorded evidence and feature
dashboard. Omitting `--corpus` has the same all-corpus behavior. All local
corpus data must already be present.
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
