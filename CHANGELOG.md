# Changelog

All notable changes to Kekule are documented in this file.

## [Unreleased]

### Added

- Add `kekule_traj::io::read_trajectory` and `read_trajectory_with_options` to load
  complete trajectories through the existing streaming codecs, preserving decoded
  frame state, topology sharing, and validation.
- Add a trajectory workflow example that loads an mmCIF topology, prints frame
  information, and aligns to the first frame. Streaming readers remain public
  for processing files without loading every frame into memory.

## [0.2.1] - 2026-09-01

This compatible workspace release adds canonical molecule and residue
classification and uses it for ordinary mmCIF entity planning.

### Added

- Add topology-owned `MoleculeClass` and `ResidueClass` values with automatic,
  conservative inference during topology publication.
- Add definition-, instance-, residue-, builder-override-, and typed-selection
  APIs for canonical classification.
- Preserve classification through definition reuse and append-only topology
  transformations while re-inferring it for structural subsets.

### Changed

- Derive ordinary mmCIF polymer, water, and non-polymer entity kinds from
  canonical topology classification; explicit expert classifications and
  source interpretation reports remain authoritative overrides.
- Keep carbohydrate projection conservative: a single carbohydrate residue is
  written as a non-polymer, while multi-residue carbohydrates require explicit
  or source-preserved mmCIF entity semantics.
- Simplify the combined SDF and mmCIF README workflow so generic models no
  longer need a manually assembled entity-classification sidecar.

## [0.2.0] - 2026-08-30

This release establishes the canonical object model described in
`ARCHITECTURE.md`. It contains breaking API changes throughout the workspace.

### Breaking changes

- Make every published `Molecule` one non-empty connected chemical graph.
  Disconnected salts, complexes, and systems are represented by multiple
  molecule instances in a `Topology`.
- Move chain, residue, and atom-site hierarchy ownership to `Topology`, with
  topology-qualified atom and bond identities and deterministic dense ordering.
- Replace parallel annotation containers with the unified `Properties`,
  `PropertyTable`, `PropertyColumn`, `PropertyKey`, and `PropertyValue` APIs.
- Separate format parsing from canonical chemical interpretation. SMILES,
  Molfile, SDF, and mmCIF now expose format-specific documents and
  interpretation/projection APIs.
- Replace subsystem-specific model units with one library-wide canonical unit
  system based on nanometers, daltons, picoseconds, kilojoules per mole,
  elementary charge, kelvin, and radians.
- Refactor geometry-bearing objects around shared immutable `Topology` values,
  topology-free dense storage, and explicit `Model`, `Ensemble`, and
  `Trajectory` realization views.

### Added

- Canonical SMILES, Molfile, SDF, and mmCIF writing APIs, including streaming
  writers and explicit format-loss/error reporting.
- Public molecule build/edit transactions that enforce graph publication
  invariants.
- Topology construction, hierarchy navigation, selection, slicing, and
  operation-specific source-to-target correspondence.
- Configurable rotatable-bond detection.
- Richer trajectory streaming, indexing, slicing, RMSD, superposition, and
  DCD, TRR, XTC, and XYZ codec workflows in `kekule-traj`.
- Shared structural-view integration for `kekule-potentials`.

### Changed

- Aromatic source notation is localized before molecule publication;
  aromaticity remains derived perception state.
- Source stereochemistry is normalized into canonical graph stereo, and CIP,
  valence, ring, and aromaticity perception updates are transactional.
- Hydrogen declarations, mmCIF connectivity and hierarchy reconstruction, and
  Molfile/SDF component handling are stricter and more explicit.
- Crate documentation, examples, package metadata, and public API regression
  coverage now describe the canonical workflows.

### Migration notes

- Use format namespaces such as `kekule::smiles`, `kekule::molfile`,
  `kekule::sdf`, and `kekule::mmcif` for parsing, interpretation, and writing.
- Use `to_molecules()` for source scopes that can contain disconnected
  components; no conversion silently chooses a main component.
- Build systems with `Topology::from_molecule`,
  `Topology::from_molecules`, or `TopologyBuilder`, then combine a shared
  topology with `Positions` through `Model::new`.
- Access hierarchy through `Topology` or topology-bound molecule/model views.
- Store extensible annotations at the narrowest valid owner scope through the
  unified property APIs.

## [0.1.0] - 2026-08-05

- Initial release of `kekule`, `kekule-traj`, and `kekule-potentials`.

[0.2.1]: https://github.com/choutkaj/kekule/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/choutkaj/kekule/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/choutkaj/kekule/releases/tag/v0.1.0
