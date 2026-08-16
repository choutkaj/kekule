# Benchmarks

This tree provides optional comparisons between Kekule and pinned external references:

- `corpora/<corpus-id>/corpus.toml` describes a corpus;
- `sources.lock.json` pins record identity, URLs, checksums, categories, and packs;
- `features/*.toml` contains benchmark manifests;
- `golden/<benchmark-id>/` contains deterministic compressed reference outputs;
- `reference/` contains RDKit, Biopython, and DSSP reproduction tools.

The `features` directory name and manifest `feature_id` field are legacy benchmark-schema vocabulary. They identify benchmark targets and do not refer to a repository feature registry.

All tracked manifests are discovered directly from the filesystem. Run one comparison or a deterministic selection of all available manifests with:

```bash
cargo xtask benchmark --benchmark io.smiles.parse --corpus pubchem-1k
cargo xtask benchmark --benchmark all --corpus pubchem-1k
cargo xtask benchmark --benchmark all --corpus all
```

Use `--fixture PATH` for a partial diagnostic run. Independently reviewed `*-manual-semantic` goldens can be accepted only with one concrete benchmark/corpus pair plus `--accept-implementation-goldens`; generated RDKit, Biopython, and DSSP goldens cannot be replaced by that option.

`cargo xtask corpus check --corpus all` validates tracked provenance, nesting, and compressed goldens. Add `--require-data` to require and byte-check ignored local fixtures and packs. `--jobs N` overrides the automatic comparison worker count, which is capped at four.

Benchmark invocations report current matches and differences without updating unrelated repository metadata. They are deliberate scientific tools, not CI or release gates. External fixtures remain provenance-pinned, and reference tools never become Rust runtime dependencies.
