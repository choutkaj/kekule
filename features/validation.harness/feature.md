# Reference Validation Harness

## Summary

Provide repeatable infrastructure for comparing Rust behavior against reference-generated golden data.

## Behavior/API

- Exposes `cargo xtask validate --feature FEATURE_ID`.
- Discovers optional validation manifests under `validation/features/<feature-id>/validation.toml`.
- Verifies listed fixture paths exist.
- Reports clearly when a feature has no reference validation manifest.

## Implementation Notes

- RDKit reference generators live under `validation/reference/rdkit/`.
- Biopython reference generators live under `validation/reference/biopython/`.
- Golden data should be normalized JSON and include reference tool versions.
- Reference tools are never Rust runtime dependencies.

## Validation

- Current coverage is infrastructure unit-test based.
- Existing manifests are fixture-discovery only until reference goldens are committed and compared.

## Out Of Scope

- Chemistry algorithms.
- Runtime RDKit or Biopython dependencies.
- Automatically marking features validated.
- Regenerating all goldens by default.

## Revision Notes

- v1: Manifest discovery, fixture path checks, and reference generator conventions.
