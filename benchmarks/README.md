# Benchmarks

This tree provides optional external-reference comparisons:

- `corpora/<corpus-id>/` owns corpus metadata, pinned membership, local inputs,
  feature manifests, compressed goldens, and informational `results.toml`
  observations.
- `reference/` owns RDKit, Biopython, and DSSP corpus/reproduction tooling.

All registered `data/` directories are generated locally and ignored.
`sources.lock.json` pins record identity, URL, checksum, category, and pack
membership. Locks, manifests, goldens, and result observations are tracked.

```bash
cargo xtask benchmark --feature <feature-id> --corpus <corpus-id>
cargo xtask benchmark --feature all
cargo xtask benchmark --feature all --corpus all
cargo xtask corpus check --corpus all
```

Omitting `--corpus` considers only manifest-backed PubChem-1k and PDB-100
targets. Explicit `--corpus all` also considers registered broad corpora.
Absent combinations and planned features are skipped for `all` selectors; a
concrete feature/corpus request without a manifest is an error.

Every requested target records `match`, `differences`, or `error` with its
scope and timestamp, then refreshes the dashboard. A difference or error is
nonzero for that explicit invocation, but no repository workflow uses
benchmarks as a gate. Result snapshots never affect feature or release status.

Use `--fixture PATH` for a partial diagnostic run. Independently reviewed
`*-manual-semantic` goldens can be accepted only with one concrete
feature/corpus plus `--accept-implementation-goldens`; generated RDKit,
Biopython, and DSSP goldens cannot be replaced by that option.

`cargo xtask corpus check` is an optional integrity utility. Add
`--require-data` to require and byte-check ignored local fixtures and packs.
The automatic comparison worker count is capped at four; `--jobs N` overrides
it.

External fixtures remain provenance-pinned, and reference tools must never
become Rust runtime dependencies.
