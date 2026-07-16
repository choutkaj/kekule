# Changelog

All notable changes to this project will be documented in this file.

The project follows Cargo semantic-versioning conventions. During the `0.x`
series, breaking public API changes increment the minor version.

## 0.1.0 - 2026-07-16

Initial release.

- Stable-ID molecular graph kernel with explicit perception state.
- Separate small-molecule and macromolecule domain boundaries.
- Staged SMILES, Molfile, SDF, and mmCIF parsing and interpretation.
- Configurable parser resource limits and structured rejection of malformed or
  unsafe record boundaries.
- Qualified biomolecular hierarchy and mmCIF provenance.
- Fixed-topology molecular models and the DREIDING adapter.
- Explicit sanitization, hydrogen normalization, query, substructure,
  canonicalization, and modelling workflows.
- Bounded parser fuzz smoke tests in CI and longer scheduled campaigns.
