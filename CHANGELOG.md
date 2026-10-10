# Changelog

All notable changes to Kekule are documented in this file.

## [0.3.0] - 2026-10-10

This release adds the `kekule-openff` and `kekule-openff-ash` crates for OpenFF
parameterization, rebuilds `kekule-potentials` around OpenFF energies and L-BFGS
minimization, and moves trajectories and alignment into `kekule`. It contains
breaking API changes throughout the workspace.

### Changed

- **Breaking:** An atom's implicit hydrogen count is either fixed or wholly
  inferred. `HydrogenDeclaration { Infer { explicit }, Fixed }` becomes
  `ImplicitHydrogens { Inferred, Fixed }`, and the mixed "stored count plus
  inference" state is removed. `Molecule`, `MoleculeEditor`, and topology atom
  views report `explicit_hydrogens`, `implicit_hydrogens`, and
  `total_hydrogens`; `implicit_hydrogens` is now the complete non-graph count
  instead of the perceived part, which `Perception::inferred_hydrogens` reports.
  `AddHydrogensOptions::explicit_only` becomes `fixed_only`,
  `AddedHydrogenOrigin::ExplicitCount` and `Implicit` become `Fixed` and
  `Inferred`, `HydrogenCountAdjustment` reports the resulting `hydrogens`
  declaration, and `HydrogenTransformError::HydrogenCountNotPreserved` is
  removed. Installing perception that infers hydrogens on a fixed count fails
  with `PerceptionInstallError::InferredHydrogensOnFixedAtom`. V2000 and V3000
  writers no longer have an unencodable hydrogen state. See
  [hydrogen semantics](docs/hydrogen-semantics.md) for migration.
- **Breaking:** `Trajectory`, `TrajectoryFrame`, `Velocities`, and `Forces` move
  from `kekule-traj` into `kekule::structure`, next to `Ensemble`. The two
  collections stay distinct types over one private store and share `new`,
  `from_items`, `get`/`get_mut`/`iter`, `push`, `replace`, `remove`,
  `replace_positions`, `select`, `subset`, `perceive`, `into_parts`, and
  `into_items`; item views are `EnsembleMemberView`/`EnsembleMemberMut` and
  `TrajectoryFrameView`/`TrajectoryFrameMut`. `select` replaces
  `Trajectory::slice` and copies items in any order, including ranges, strides,
  and repeated indices. `Ensemble::from_models` consumes models.
  `Trajectory::into_ensemble` projects frames onto equally weighted members
  without copying conformations. `TrajectoryError` wraps core
  `ConformationError` and `RealizationError`.
- **Breaking:** every `EnsembleMember` carries a finite, strictly positive
  relative weight: `EnsembleMember::new(conformation, weight)` is fallible,
  `weight()` returns `f64`, and `set_weight` takes `f64`. Models, mmCIF
  coordinate models, and trajectory frames become equally weighted members
  (weight `1.0`). `Ensemble::normalize_weights` keeps weight ratios and cannot
  overflow; `RealizationError::MissingWeight` and `ZeroTotalWeight` are removed.
- **Breaking:** realization payloads hold a `Conformation`: dense positions,
  optional cell, typed occupancies and B factors, and `RealizationProperties`,
  addressed by `TopologyAtomIndex`. `Model` is a topology plus one conformation;
  `EnsembleMember` adds a weight and `TrajectoryFrame` adds velocities, forces,
  time, and step. Conformation state is edited through the `ConformationMut`
  guard from `conformation_mut()` on models, payloads, and collection items.
  Kernels accept any `AsModelView` source (models, collection items, frame
  buffers).
- **Breaking:** typed, scoped properties. `Properties` is replaced by
  `MoleculeProperties`, `TopologyProperties`, `RealizationProperties`, and
  collection `OwnerProperties`, built from `PropertyTable<R>` columns whose row
  type names the domain (`AtomId`, `TopologyAtomIndex`, `ChainId`, ...). Writes go
  through `properties_mut()` guards with whole-column, single-cell, and
  transactional batch updates. Keys are no longer reserved: occupancy and B
  factors are typed conformation fields, and a generic `occupancy` property is
  independent of them. A draft editor rejects values on deleted slots with
  `PropertyError::RemovedRow`. Detached payloads have atom rows only and gain
  bond rows when bound to a collection.
- **Breaking:** superposition and RMSD move from `kekule-traj` to every
  collection in `kekule::alignment`. `Trajectory::superpose_to_frame`,
  `rmsd_to_frame`, and `aligned_rmsd_to_frame` become `superpose`, `rmsd`, and
  `aligned_rmsd` (each still with `_with_options`), which take a `Reference`
  (an item index or any borrowed view) and `FitAtoms` (a selection or an
  `AtomCorrespondence`). Superposition stays in place and transactional and
  returns a `SuperpositionReport`. One `AlignmentOptions` (`Weighting`,
  `PeriodicPolicy`) serves fitting and RMSD. Collection-wide failures are
  reported once instead of per item.
- **Breaking:** fitting and RMSD use stored Cartesian coordinates by default,
  including for periodic frames, which were previously rejected.
  `PeriodicPolicy::RejectPeriodic` restores rejection. Molecular
  reconstruction, imaging, and temporal unwrapping are separate preprocessing
  operations in `kekule_traj::periodic`.
- **Breaking:** topology reads return views. `Topology::atom`, `bond`,
  `molecule`, `chain`, `residue`, and `atom_site` return `Option` views
  (`AtomView`, `BondView`, ...) with neighbors, hydrogens, aromaticity,
  hierarchy, and static properties; `Model::atom` adds the atom's position,
  occupancy, B factor, and realization properties. Selection predicates take one
  view. Topology-level hydrogen counters and lookup forwarders on `Model` are
  removed. `Molecule` and `MoleculeEditor` share their read API through `Deref`
  to `Graph`.
- **Breaking:** one writer pair per format: `write(source, options)` and
  `write_to(writer, source, options)` for SMILES (`SmilesWriteOptions::ordinary`,
  `isomeric`, `canonical`), Molfile (`MolfileSource`), SDF (records from models,
  collection items, or interpreted records), and mmCIF (`MmcifBlockSource` with
  optional classifications or reports; multiple blocks are suffixed `_1`, `_2`,
  ...). Version- and source-specific writer functions are removed.
- **Breaking:** constructors that take chemistry take ownership:
  `add_molecule_definition`, `add_molecule`, `Topology::from_molecule`,
  `Topology::from_molecules`, and `Model::from_molecule` take `Molecule` values
  instead of references. Substructure search uses `find_match`, `find_matches`,
  `visit_matches`, and `find_topology_matches` (with `_with_options`) in place
  of `find_substructure_match` and `find_substructure_matches`; a match limit is
  a resource error, never a truncation.
- **Breaking:** editor publication no longer rewrites oxohalogens silently.
  Call `MoleculeEditor::normalize_oxohalogens` first; `finish` rejects
  unnormalized ones with `UnnormalizedOxohalogen`. Format interpretation still
  produces normalized molecules.
- **Breaking:** perception keeps a topology's layout identity. `Topology::perceived`
  and model, ensemble, and trajectory `perceive` share the published layout
  instead of copying it. Selections, frame buffers, readers and writers,
  alignment, measurements, geometry edits, reductions, periodic tools, subsets,
  and potentials accept any snapshot sharing their layout, checked with the new
  `Topology::shares_layout`, so perception no longer detaches them. Independently
  published equal topologies stay incompatible, and prepared substructure
  targets still bind their exact snapshot. Selection equality uses layout
  identity. `Topology` hierarchy and property accessors are no longer `const fn`.
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
  `MoleculeCorrespondence`, and `RemoveHydrogensReport::correspondence` maps
  input IDs to the result. Remove `Molecule::stereo_group_slots`,
  `stereo_group_slot_count`, and `MoleculeEditor::append_stereo_group_tombstone`;
  rename `RingMembership` slot-flag accessors to `from_flags`, `atom_flags`, and
  `bond_flags`.
- **Breaking:** rebuild `kekule-potentials` around a backend-independent
  `Potential` trait over `ModelView`, with validated `Energy` decompositions and
  `Evaluation` gradients. Potentials evaluate through `&self` and are `Send + Sync`.
- **Breaking:** move potential evaluation and minimization out of `kekule`; the
  `kekule::modeling` module, `HarmonicBondPotential`, and steepest-descent
  `minimize` are removed.
- **Breaking:** remove the DREIDING potential and its `dreid-forge`/`dreid-kernel`
  dependencies, which also removes the unmaintained `paste` advisory.
- Default perception has no molecule size limit. `Molecule::perceive`,
  `perceive_ring_set`, and `perceive_aromaticity` scale their ring and
  aromaticity work bounds with the molecule (`RingPerceptionOptions::for_graph`,
  `AromaticityOptions::for_graph`); explicitly supplied options stay absolute.
  Smallest-ring searches traverse only ring bonds with workspace proportional to
  what they visit, so ring perception no longer costs rings times molecule size:
  a 70,000-atom polymer with 10,000 rings perceives in well under a second. Ring
  sets are unchanged.
- Borrow validated stored frames without rescanning dense properties on every
  access.

### Added

- The `kekule-openff` crate: pure-Rust SMIRNOFF parameter assignment with
  native NAGL partial charges. It bundles the OpenFF Rosemary force field
  (`ForceField::rosemary`) and, through the default `ash` feature, the Ash
  charge model from the new `kekule-openff-ash` data crate (`NaglModel::ash`),
  so a complete parameterization needs no Python, C toolchain, external files,
  or network access. `ForceField::parameterize` assigns a shared topology once
  per molecule definition and `parameterize_molecule` takes one molecule; both
  take a `ChargeMethod` (a NAGL model, or `ChargeMethod::LibraryOnly`). Failures
  are `Error` values with a stable `ErrorKind`. Custom OFFXML within the
  supported SMIRNOFF subset and other schema-2 NAGL bundles are accepted, and
  `kekule_openff::diagnostics` exposes NAGL features, raw inference, and lookup
  keys. There is no molecule size limit. The crate's `CONTRACT.md` records the
  supported subset, assignment semantics, error kinds, and validation. Both
  crates are licensed `(MIT OR Apache-2.0) AND CC-BY-4.0` because they embed
  CC BY 4.0 OpenFF data.
- NAGL lookup selects a stored entry only when the input is exactly that
  entry's molecule (isotopes ignored). Other bond-order or charge-placement
  forms that upstream's fixed-H InChI lookup also merges use inference instead;
  for Ash these are exotic, mostly with formal charges of magnitude 2-6.
- Add `kekule_potentials::openff::OpenFfPotential`, which evaluates OpenFF
  energies, gradients, and per-component gradients from a `ParameterizedTopology`
  in vacuum without cutoffs, with optional replacement charges.
- Add L-BFGS `minimize` and `minimize_with_observer` with a strong-Wolfe line
  search, per-step displacement bound, singular-trial backtracking, and
  `Minimization::to_model`.
- Add `kekule_traj::periodic::make_molecules_whole`, `image_molecules`, and
  temporal `unwrap`, which act in place on a `Trajectory`. They support
  triclinic, rotated, and partially periodic cells, validate bonded ring
  closure, and reject ambiguous temporal crossings. Imaging uses explicitly
  selected anchor molecules. Unwrapping rejects decreasing available times,
  including across frames without timestamps; equal times are allowed.
- Add streaming `periodic::MoleculeImager`, `periodic::TrajectoryUnwrapper`, and
  `analysis::FrameSuperposer`, which apply the same transformations to a
  `FrameBuffer` and take the caller's frame index for diagnostics. Unwrapping
  retains state across chunks, checks consecutive source indices, supports
  explicit reset, and rolls back on failure.
- Add `rmsf` and `contact_occupancy` over a `Trajectory`, where every frame
  counts once; they reject ensembles, whose statistics must be weighted.
- Add `kekule_traj::io::read_trajectory` and `read_trajectory_with_options` to
  load complete trajectories through the existing streaming codecs, preserving
  decoded frame state, topology sharing, and validation.
- Add `io::write_trajectory` and `write_trajectory_with_options` for atomic saving
  through the strict codecs, with extension inference, explicit format/precision
  options, metadata preservation checks, and overwrite protection.
- Add a trajectory workflow example that loads an mmCIF topology, prints frame
  information, and aligns to the first frame. Streaming readers remain public
  for processing files without loading every frame into memory.
- Add infallible `AtomSelection::all(&topology)` with authoritative dense order
  and exact topology sharing.
- Add an optional comparison of the periodic transformations against references
  generated with pinned MDTraj and MDAnalysis versions
  (`trajectory_periodic_reference` example), and an informational
  frame-access/superposition benchmark.

### Fixed

- `remove_hydrogens` preserves counts that valence inference cannot reproduce,
  such as `[H]S([H])([H])[H]` and `PH5`, by fixing them instead of failing with
  a count-preservation error. Aromatic `[nH]` and stereo parents collapse to
  fixed counts.
- Canonical SMILES of a perceived molecule with explicit hydrogen atoms matches
  its hydrogen-suppressed form. Collapsing those atoms used to discard the
  perceived hydrogen counts, so an aromatic NH such as pteridine-2,4-dione's
  failed with a request to perceive the molecule, and other parents were ranked
  as if they carried no hydrogens.
- Trajectories read with an interpreted structure topology assign coordinates
  to the right atoms when the source interleaves molecules, for example a
  covalently linked ligand listed after water. Dense order used to follow
  molecule instances, so such files silently swapped coordinates.
- Fix compressed XTC decoding when a frame reuses a preceding nonzero coordinate
  run length. Valid trajectories from external GROMACS tooling no longer fail
  mixed-radix bounds checks; existing corruption checks remain intact.
- Retain underlying frame, property, position, unit, and model errors through
  `std::error::Error::source`; model failures are no longer flattened into topology
  mismatch strings.
- `Molecule::remove_hydrogens` keeps a double-bond reference hydrogen as a graph
  atom, reported as `UnsupportedStereoRole`, when its endpoint has another
  hydrogen, as on a terminal `=CH2`. Collapsing it left an implicit reference
  that named neither hydrogen: CIP assignment then failed with
  `UnresolvedPriority` instead of skipping the bond as nonstereogenic, and
  `add_hydrogens` rejected the molecule. This affected 3D structures such as
  acrylamide atropisomers read with explicit hydrogens.

### Removed

- The unpublished `kekule-bench` benchmark layer under `benchmarks/`: stored
  goldens, contract hashes, dashboard, run history, reference adapters and the
  OpenFF/OpenMM numerical validation. Inputs used by crate tests moved to
  `crates/kekule/tests/fixtures/corpus`, and the exporters and provenance of
  shipped OpenFF data moved to `tools/openff`. It is replaced by a scientific
  benchmark on two curated, hash-locked datasets (small molecules and wwPDB
  entries) that compares chemistry-level observations with RDKit, gemmi,
  Biotite and mkdssp nightly and triages every difference against a
  checked-in list of known differences.

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

[0.3.0]: https://github.com/choutkaj/kekule/compare/v0.2.1...v0.3.0
[0.2.1]: https://github.com/choutkaj/kekule/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/choutkaj/kekule/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/choutkaj/kekule/releases/tag/v0.1.0
