# Changelog

All notable changes to this project will be documented in this file.

The project follows Cargo semantic-versioning conventions. During the `0.x`
series, breaking public API changes increment the minor version.

## Unreleased

- Begin the deliberate 0.2.0 topology-centered API transition: promote an
  immutable coordinate-free `Topology`, separate configurations from topology,
  add finite ensembles and streaming-first trajectories, bind selections and
  prepared potentials to exact topology identity, and remove the obsolete
  model-topology vocabulary once consumers are migrated.
- Reclassify RDKit, Biopython, and DSSP external parity from a required
  repository gate to optional benchmarking; preserve corpora and goldens, add
  neutral result observations, and remove repository identity from benchmark
  input digests and corpus selection identifiers.
- Hard-rename the repository and publishable crates from `molecules` to
  `molecular`, including the DREIDING adapter, Rust import paths, benchmark
  tooling, generated writer provenance, and project branding.

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
- Generated feature/corpus benchmark matrix with PubChem-1k and PDB-100
  defaults, broader PubChem/Enamine/PL-REX/PDB-1000 tiers, and reproducible
  RDKit/Biopython/mkdssp reference artifacts.
