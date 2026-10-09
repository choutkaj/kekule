# Changelog

All notable changes to Kekule are documented in this file.

## [Unreleased]

### Changed

- **Breaking:** `Topology` stores its dense atom order instead of deriving it
  from instance order, and instances need not be contiguous. Format
  interpretation keeps source atom-row order (mmCIF `_atom_site` rows,
  Molfile/SDF atom blocks, SMILES atom order). `TopologyBuilder::set_atom_order`
  and `ModelBuilder::set_atom_order` accept explicit permutations. Editors,
  subsets, and instance filters keep surviving atoms in source order and append
  new atoms. Qualified-ID and dense-index lookups are constant time, and
  `TopologyBuilder`/`ModelBuilder` `atom_ids` and `bond_ids` return slices.
- **Breaking:** mmCIF writers emit `_atom_site` rows in dense atom order instead
  of hierarchy order, so written structures stay aligned with trajectories
  written in the same order.
- **Breaking:** published molecules have dense IDs. `MoleculeEditor::finish`
  renumbers atoms, bonds, stereo elements, and stereo groups after deletions;
  `finish_with_correspondence` returns the draft-to-published
  `MoleculeCorrespondence` (renamed from `MoleculeAppendMapping`), and
  `RemoveHydrogensReport::correspondence` maps input IDs to the result. Remove
  `Molecule::stereo_group_slots`, `stereo_group_slot_count`, and
  `MoleculeEditor::append_stereo_group_tombstone`; rename `RingMembership`
  slot-flag accessors to `from_flags`, `atom_flags`, and `bond_flags`.
- **Breaking:** rebuild `kekule-potentials` around a backend-independent
  `Potential` trait over `ModelView`, with validated `Energy` decompositions and
  `Evaluation` gradients. Potentials evaluate through `&self` and are `Send + Sync`.
- **Breaking:** move potential evaluation and minimization out of `kekule`; the
  `kekule::modeling` module, `HarmonicBondPotential`, and steepest-descent
  `minimize` are removed.
- **Breaking:** remove the DREIDING potential and its `dreid-forge`/`dreid-kernel`
  dependencies, which also removes the unmaintained `paste` advisory.
- **Breaking:** `kekule-openff` interactions and per-atom vdW entries share one
  `Arc` parameter allocation per rule, and nonbonded methods are typed
  `VdwMethod`/`ElectrostaticsMethod` values instead of strings.

### Added

- Add `kekule_potentials::openff::OpenFfPotential`, which evaluates OpenFF
  energies, gradients, and per-component gradients from a `ParameterizedTopology`
  in vacuum without cutoffs, with optional replacement charges.
- Add L-BFGS `minimize` and `minimize_with_observer` with a strong-Wolfe line
  search, per-step displacement bound, singular-trial backtracking, and
  `Minimization::to_model`.
- The optional OpenFF energy and gradient benchmark now exercises the public
  `OpenFfPotential` instead of a benchmark-local evaluator.

- Add `io::write_trajectory` and `write_trajectory_with_options` for atomic saving
  through the strict codecs, with extension inference, explicit format/precision
  options, metadata preservation checks, and overwrite protection.
- Add ordered `Trajectory::select_frames` for ranges, strides, reordering, and
  repeated indices; preserve original time, step, properties, and shared topology.
- Add complete `TrajectoryFrameView::to_frame`, validated `replace_frame`, and a
  restricted `TrajectoryFrameMut` editor whose setters preserve dense dimensions.
- Add reusable `analysis::FrameSuperposer`, `periodic::MoleculeImager`, and stateful
  `periodic::TrajectoryUnwrapper` for streaming the same transformations used by
  loaded trajectories. Unwrapping retains state across chunks, checks consecutive
  source indices, supports explicit reset, and rolls back on failure.
- Add an optional external-trajectory comparison against pinned MDTraj and
  MDAnalysis versions, and an informational frame-access/superposition benchmark.
- Add `kekule_traj::io::read_trajectory` and `read_trajectory_with_options` to load
  complete trajectories through the existing streaming codecs, preserving decoded
  frame state, topology sharing, and validation.
- Add a trajectory workflow example that loads an mmCIF topology, prints frame
  information, and aligns to the first frame. Streaming readers remain public
  for processing files without loading every frame into memory.
- Add infallible `AtomSelection::all(&topology)` with authoritative dense order
  and exact topology sharing.
- Add `Trajectory::make_molecules_whole`, `image_molecules`, and temporal `unwrap`,
  with copy-returning and transactional `_in_place` variants. Support triclinic,
  rotated, and partially periodic cells; validate bonded ring closure and reject
  ambiguous temporal crossings. Imaging uses explicitly selected anchor molecules.
- Add explicit in-place superposition and opt-in superposition reports, and update
  the runnable workflow example to use ordinary alignment without policy options.

### Fixed

- Trajectories read with an interpreted structure topology assign coordinates
  to the right atoms when the source interleaves molecules, for example a
  covalently linked ligand listed after water. Dense order used to follow
  molecule instances, so such files silently swapped coordinates.
- `kekule-openff` places per-atom vdW parameters and charges in dense atom
  order; instance order misassigned them when instances were not contiguous.
- Fix compressed XTC decoding when a frame reuses a preceding nonzero coordinate
  run length. Valid trajectories from external GROMACS tooling no longer fail
  mixed-radix bounds checks; existing corruption checks remain intact.
- Retain underlying frame, property, position, unit, and model errors through
  `std::error::Error::source`; model failures are no longer flattened into topology
  mismatch strings.

### Changed

- Borrow validated stored frames without rescanning dense properties on every
  access, and avoid allocating discarded per-frame superposition report vectors.
- Temporal unwrapping rejects decreasing available times, including across frames
  without timestamps. Equal times remain allowed.

### Breaking changes

- `superpose_to_frame` and `superpose_to_frame_with_options` now return a new
  trajectory, leaving their source unchanged. Use the corresponding `_in_place`
  methods for mutation or `_with_report` methods for `(trajectory, report)` results.
  `Trajectory` is marked `must_use` to diagnose accidentally discarded copies.
- Kabsch fitting, trajectory superposition, and RMSD now use stored Cartesian
  coordinates by default, including periodic frames. Strict periodic rejection is
  still available through explicit options. Molecular reconstruction, imaging, and
  temporal unwrapping are separate preprocessing operations.

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
