# Validation

Reference validation is separated from the pure Rust library.

- `validation/features/<feature-id>/` contains feature-scoped manifests and compact fixtures.
- `validation/corpora/` describes optional large local datasets that are not committed to git.
- `validation/data/` is for local unpacked corpora and is ignored by git.
- `validation/reference/rdkit/` contains small-molecule reference generators.
- `validation/reference/biopython/` contains macromolecular reference generators.
- `validation/golden/` contains normalized JSON expectations checked by Rust tests or `xtask` commands.

Golden files should record the input, expected normalized behavior, reference implementation, and reference version.

The checked-in fixtures are intentionally small and durable. Some fixtures may exercise behavior that
the current prototype does not fully support yet; that is acceptable. The fixture set is a long-term
proving ground, not a snapshot tailored to the current implementation.

Large corpora should be treated as optional validation tiers:

- normal CI: checked-in tiny fixtures and any committed golden files
- local corpus smoke: sampled files from `validation/data/`
- full corpus regression: complete local datasets, usually outside normal CI

Do not mark a feature as validated merely because a manifest or fixture exists. `validated = true`
requires reference-generated golden data or documented manual validation evidence.
