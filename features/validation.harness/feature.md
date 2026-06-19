# Reference Validation Harness

## Summary

Provide repeatable infrastructure for comparing Rust behavior against reference-generated golden data.

## Behavior/API

- Exposes `cargo xtask validate --feature FEATURE_ID`.
- Discovers optional validation manifests under `validation/features/<feature-id>/validation.toml`.
- Verifies listed fixture paths exist.
- Requires one golden JSON file under `validation/features/<feature-id>/golden/` for each listed fixture.
- Compares normalized Rust implementation output against each golden file's `expected` payload.
- Normalizes representation-only graph differences such as undirected bond endpoint orientation, bond array order, and ring atom order before comparison.
- Reports clearly when a feature has no reference validation manifest.

## Implementation Notes

- RDKit reference generators live under `validation/reference/rdkit/`.
- Biopython reference generators live under `validation/reference/biopython/`.
- Golden data should be normalized JSON and include reference tool versions.
- The validation command uses the Rust implementation only; RDKit and Biopython are used to generate goldens, not to run validation.
- Reference tools are never Rust runtime dependencies.

## Validation

- Current coverage is infrastructure unit-test based plus live `cargo xtask validate --feature ...` comparisons against committed external-source goldens.
- Passing comparisons are evidence for the compared behavior; failing comparisons identify implementation gaps and should not be papered over.

## Out Of Scope

- Chemistry algorithms.
- Runtime RDKit or Biopython dependencies.
- Automatically marking features validated.
- Regenerating all goldens by default.

## Revision Notes

- v1: Manifest discovery, fixture path checks, and reference generator conventions.
- v2: Implementation-vs-golden comparisons for committed per-feature golden JSON.
