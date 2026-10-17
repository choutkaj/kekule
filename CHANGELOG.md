# Changelog

All notable changes to Kekule are documented in this file.

## [Unreleased]

### Added

- The `kekule-traj` readers, writers, and stream adapters (`XyzReader`,
  `IndexedDcdReader`, `FileTrajectoryWriter`, `MoleculeImager`,
  `TrajectoryUnwrapper`, and the others) and `kekule::core::AtomMut` and
  `BondMut` implement `Debug`. Readers and writers show their options,
  progress, and atom count, never their stream or decode buffers.
- `XtcWriteOptions::with_cell_policy` mirrors the reader's `XtcCellPolicy`:
  `ZeroMatrixAsAbsent` writes the all-zero box GROMACS uses for
  non-periodic systems, so XTC files without a cell round-trip. The default
  still requires a periodic cell.
- `smiles::write_to` writes one SMILES record to a sink as a line ending in
  a newline, completing the `write` and `write_to` pair every writer offers.
- `ComponentLibrary::amber` and `ComponentLibrary::charmm` define the residue
  and atom names of structures that AMBER and CHARMM tools write, as
  interpretation context for `ConnectivityOptions::components`. Their
  templates derive from the bundled dictionary entries as the AMBER ff14SB,
  OL15, and OL3 and the CHARMM36 libraries define each residue: protonation
  states such as `HIE`, `HID`, `ASH`, `CYX`, `HSD`, `HSE`, and `LSN`,
  charged `ASP` and `GLU`, AMBER terminal residues such as `NALA` and `CALA`
  and nucleotides such as `DA5`, the `ACE` cap, water (`WAT`, `TIP3`), and
  ions (`Na+`, `SOD`, `CLA`). CHARMM templates carry CHARMM atom names such as
  `HN`, `HB1`, `CD`, `OT1`, and `OT2`. The libraries are opt-in, since
  force-field names collide with dictionary identifiers.
- `kekule::pdb` reads PDB-format entries. `pdb::parse_str` and `parse_reader`
  type the coordinate (`ATOM`, `HETATM`), `TER`, `MODEL`, `CONECT`, `LINK`,
  and `SSBOND` records, decoding hybrid-36 numbers, and keep every other
  record verbatim. `PdbDocument` interprets one coordinate model, an ensemble
  of models, or explicit conformations through the mmCIF pipeline and reports.
  A chain's polymer ends at its `TER` record. `LINK` and `SSBOND` records are
  covalent bonds, and a `LINK` to a metal is dative. `CONECT` records
  supplement residue templates, stand in for distance inference, and state a
  double or triple bond by listing a partner repeatedly. Atoms without an
  element symbol take the element their name's alignment implies, reported as
  `ElementsInferredFromAtomNames`, and a residue number that recurs after
  other residues names a new residue, `SourceResiduePosition::RepeatedAuthor`.
- mmCIF interpretation resolves bonds the file does not state from component
  templates, as described in ARCHITECTURE.md. `ConnectivityOptions`
  (field `connectivity` of `StructureInterpretOptions` and
  `StructureEnsembleInterpretOptions`) adds caller templates ahead of the bundled
  ones and switches distance inference for residues without templates, which
  is on by default. Every published bond records its evidence in
  `InstanceProvenance::bonds` (`BondEvidence`), and the report adds
  `ResidueBondsInferred`, `HydrogenAttachedByDistance`, `ChainBreak`,
  `ProtonationFromHydrogens`, and `ResidueFragmented` issues.
- `kekule::components` holds chemical component templates in the model of the
  wwPDB Chemical Component Dictionary: named atoms with elements, formal
  charges, aromatic and leaving-atom flags, and bonds with orders and
  metal-coordination flags. `ComponentLibrary::bundled` embeds 82 unmodified
  dictionary entries (standard amino acids and nucleotides, water, common
  ions, cofactors, and crystallization additives) plus the two neutral
  histidine tautomers, and `ComponentLibrary::from_mmcif` reads further
  definitions, up to the complete dictionary. Structural interpretation will
  use these templates to resolve connectivity.
- `kekule::diagnostics` provides the shared parse and interpretation
  vocabulary: `SourceLocation` (line and optional column in the original
  input), `DiagnosticKind` (I/O, limit, syntax, invalid value, invalid
  reference, unsupported, inconsistent, chemistry), and `UnsupportedContent`.
- `mmcif::parse_reader` and `parse_reader_with_options` read from any
  `std::io::Read`, enforcing the byte limit while reading.
- Molfile documents are typed. `MolfileDocument` exposes the decoded
  connection table: `atoms()`, `bonds()`, `stereo_groups()`, `chiral_flag()`,
  and `unsupported()`, as `MolfileAtom`, `MolfileBond`, `MolfileStereoGroup`,
  and `MolfileUnsupported` values with source locations.
- `sdf::records` and `records_with_options` stream SDF records one at a time.
  Records are split at `$$$$` lines, so a malformed, oversized, or non-UTF-8
  record yields one error and reading resumes at the next record; memory is
  bounded by one record.
- `molfile::parse_reader` and `sdf::parse_reader` (with `_with_options`
  variants) read from any `std::io::Read`, enforcing the byte limit while
  reading.
- `MmcifDocument` gains `to_molecules`, `to_topology`, and `to_model`, with
  `_with_options` variants, for documents with exactly one atom-site block,
  like the existing `MmcifBlock` conveniences.
- `smiles::records` and `records_with_options` stream multi-record SMILES
  input such as `.smi` files, one record per line. Blank lines are skipped,
  a malformed, oversized, or non-UTF-8 line yields one error located at its
  input line, and reading continues with the next line.
- SMILES documents type their CXSMILES extension: `SmilesDocument::cx_fields()`
  lists `SmilesCxField` values (radicals, enhanced stereo groups, the
  relative flag, and unsupported fields as written) with their spans.
  `SmilesInterpretOptions` with `UnsupportedContent::Report` omits
  unsupported fields and lists them in
  `SmilesInterpretation::omitted_cx_fields()`; `smiles::interpret_with_options`
  and `SmilesDocument::interpret_with_options` accept it.
- `MolfileInterpretOptions` selects the unsupported-content policy for
  `molfile::interpret_with_options`, `MolfileDocument::interpret_with_options`,
  and `SdfRecord::interpret_with_options`. With
  `UnsupportedContent::Report`, omitted content is listed in
  `MolfileInterpretation::omitted()` and
  `SdfRecordInterpretationReport::omitted()`.
- Published topology annotations can be edited without changing the layout.
  `Topology::properties_mut` edits in place, and
  `Model`, `Ensemble`, and `Trajectory` gain `topology_properties_mut`, which
  replaces a topology shared with other owners by an edited snapshot. As with
  perception, the result shares the source's layout, so selections, frame
  buffers, and prepared potentials stay usable. `Topology` implements `Clone`;
  a clone shares its source's layout.

### Changed

- **Breaking:** `Model::apply_transform` is renamed `Model::transform_atoms`:
  it moves a selection's atoms and keeps the cell. The new
  `Model::apply_rigid_transform(&RigidTransform)` moves the whole model and
  rotates its cell, as `EnsembleMemberMut` and `TrajectoryFrameMut` already
  did, so a fitted transform applied to a periodic model keeps its lattice in
  the same frame as its coordinates.
- **Breaking:** `kekule-traj` reductions accept only trajectory frames
  (ARCHITECTURE.md, "Trajectories and streaming").
  `RmsfAccumulator::observe` and `ContactOccupancyAccumulator::observe` take
  the new sealed `analysis::AsFrameView`, implemented by
  `TrajectoryFrameView` and `FrameBuffer`, instead of any `AsModelView`, so
  passing an ensemble member or a `ModelView` no longer compiles. Pass
  `&buffer` where `&buffer.as_model_view()` was passed.
- **Breaking:** `Atom` has private fields. Build one with `Atom::new` and
  `with_element`, `with_isotope`, `with_formal_charge`, `with_radical`,
  `with_hydrogens` and `with_atom_map`; read it with the same-named getters;
  change an atom in a molecule through the `AtomMut` setters (`set_element`,
  `set_isotope`, …).
  - The isotope is `Option<NonZeroU16>` and the atom map
    `Option<NonZeroU32>`, so mass number 0 and atom map 0 are not
    representable. SMILES `[0C]` and `[CH4:0]` and Molfile atom map 0 read
    as no isotope and no atom map, as isotope 0 already did, so `[CH4:0]`
    now writes `[CH4]`; Molfile writers no longer have an atom-map-0 error.
- **Breaking:** Stereo carriers have a fixed count:
  `TetrahedralStereo::carriers` is `[StereoCarrier; 4]` and
  `AxisStereo::carriers` is `[StereoCarrier; 2]`, so
  `StereoValidationIssue::InvalidAxisCarrierCount` is gone.
  `InvalidTetrahedralCarrierCount` remains for source readers that resolve
  the wrong number of neighbours. `StereoElement::group` is private behind
  `group()`; set it through the editor's stereo-group operations.
- **Breaking:** `AtomRadical::new(1, None)` records a doublet: one radical
  electron always has spin multiplicity 2, so it equals
  `AtomRadical::new(1, Some(2))`, and Molfile writers encode it as radical
  code 2 instead of rejecting it.
- **Breaking:** `Element` has a constant per element (`Element::H` …
  `Element::OG`), implements `FromStr` and `TryFrom<u8>`, and the new
  `ElementError` reports an unknown symbol or atomic number.
- **Breaking:** Topology handles follow one convention (ARCHITECTURE.md,
  "Units, topology handles, and traits"): `topology()` returns `&Topology`
  and `shared_topology()` returns `Arc<Topology>`.
  - `TopologyQueryMatch`, `ParameterizedTopology` and potentials
    `Evaluation` return `&Topology` from `topology()` and gain
    `shared_topology()`; `Potential::topology` returns `&Topology`, and
    `OpenFfPotential::shared_topology` returns the handle.
    `ModelEditor::topology()` is `topology_editor()`.
  - Constructors that store the topology take `impl Into<Arc<Topology>>`
    instead of `&Arc<Topology>`: every `AtomSelection` and `BondSelection`
    constructor, `AtomCorrespondence::{from_pairs, from_same_layout}`,
    `BondDihedral::new`, `{Distance,Angle,Dihedral}Edit::{new,
    with_moving_atoms}`, `SolvationRadii::new`, `FrameBuffer::new`,
    `MemoryTrajectoryWriter::new` and the `kekule-traj` writers; pass
    `Arc::clone(&topology)` or an owned topology. `ensure_compatible` and
    `semantic_ids` take `&Topology`.
  - `Topology::subset(self: &Arc<Self>, ..)` is the associated function
    `Topology::subset(source, selection)`, and `Topology::to_editor` takes
    `&self` and copies the topology; edit a shared snapshot without copying
    through `TopologyEditor::from_topology(Arc::clone(&topology))`.
  - `Model` and `Evaluation` equality compare topologies with
    `shares_layout` instead of `Arc` pointer identity.
- **Breaking:** `AsModelView` and `ScaleValue` are sealed. `FrameBuffer` no
  longer implements `AsModelView`; pass `buffer.as_model_view()`. The
  extension traits `TrajectoryReader`, `SeekableTrajectoryReader`,
  `TrajectoryWriter` and `Potential` document the contract implementors
  keep.
- **Breaking:** Frame indices and counts are `usize`:
  `SeekableTrajectoryReader::{frame_count, read_frame}`,
  `TrajectoryIoError::{frame, with_frame, frame_index_out_of_range}`,
  `FileTrajectoryMetadata::{declared_frame_count, indexed_frame_count}` and
  `TrajectoryIoLimits::max_frames`; `u64` remains for byte offsets.
- **Breaking:** Identifiers display in one style, `<kind> <number>`:
  `atom 5`, `bond 5`, `stereo element 5`, `stereo group 5`, `query atom 5`,
  `query bond 5`, `chain 5`, `residue 5`, `atom site 5`, `definition 5`,
  `instance 5`, `atom index 5`, `bond index 5`, `instance 2:atom 5`, and
  `edit atom 5` and so on for editor handles (were `a5`, `molecule2:a5`,
  `atom-index5`, `edit-atom5`, …).
- **Breaking:** Names follow the instance and definition vocabulary, and
  counts are `*_count()`: `Topology::{molecule, molecules}` are
  `Topology::{instance, instances}`, `AtomView::molecule` and
  `BondView::molecule` are `instance` (they return an `InstanceView`), and
  `InstanceAtomId::molecule`, `InstanceBondId::molecule` and
  `InstanceProvenance::molecule` are `instance`. `FitAtoms::len` is
  `pair_count` and `DistanceMatrix::len` is `member_count`.
- **Breaking:** Every public value with a physical dimension is a `Quantity`,
  and no public name carries a unit (ARCHITECTURE.md, "Units, topology
  handles, and traits").
  - `Element::{covalent_radius_angstrom, van_der_waals_radius_bohr,
    van_der_waals_radius_angstrom, dipole_polarizability_bohr3}` become
    `covalent_radius()`, `van_der_waals_radius(source)` and
    `dipole_polarizability()`, each `Option<Quantity<f64>>` in the tabulated
    unit (ångström, bohr, cubic bohr). `ElementReferenceData` and
    `Element::reference_data` are removed; the citations moved to the
    accessors.
  - `DsspResidue::{phi, psi, omega, alpha, kappa}` (were `*_degrees`) return
    `Option<Quantity<f64>>` in degrees, and `DsspHydrogenBond::energy()`
    (was `energy_kcal_per_mol`) returns kJ/mol.
  - `RigidTransform::new` takes the translation as a `Quantity<Vector3>` and
    `translation()` returns one; `RigidTransformError` gains `Unit` and is
    no longer `Copy` or `Eq`.
  - `PeriodicGeometry::{fractional, minimum_image, nearest_image}` take and
    `cartesian`, `minimum_image`, `nearest_image` and `shortest_translation`
    return length quantities; `PeriodicGeometryError` gains `Unit` and is no
    longer `Copy` or `Eq`. `PeriodicCell::signed_volume` returns a volume quantity.
  - `Quantity::is_close` takes its absolute tolerance as a quantity.
  - `kekule-traj`: `XtcWriteOptions::with_precision` takes and `precision()`
    returns an inverse length (`units::PER_NANOMETER`);
    `DcdWriteOptions::with_header_delta` takes `DELTA` as a time quantity,
    stored in the time policy's unit (the canonical time unit for
    `DcdTimePolicy::Absent`), and returns `Result`; `delta()` returns a
    quantity; `CoordinateEncoding::Lossy` carries `resolution` as a length
    quantity and loses `unit`.
  - `Quantity` no longer implements `Deref` or `DerefMut` and has no
    `value_mut`; read the stored value with `value()` or `into_value()`.
  - The panicking `Mul` and `Div` between units, between scalar quantities,
    and `Unit::powi` are removed; `try_mul`, `try_div` and `try_powi`
    remain, so unit composition never panics.
  - `Dimension::new([i32; 7])` is replaced by `Dimension::with_exponent`,
    which does not depend on the number of base dimensions.
  - `CANONICAL_MASS_UNIT` is `CANONICAL_MOLAR_MASS_UNIT`. New:
    `CANONICAL_AREA_UNIT`, `CUBIC_BOHR`, `CUBIC_NANOMETER`, `PER_NANOMETER`
    and `units::constants::{AVOGADRO, MOLAR_GAS_CONSTANT}`.
- **Breaking:** Every public enum is `#[non_exhaustive]` (ARCHITECTURE.md,
  "Options, results, and enums"), so matches outside the crate need a
  wildcard arm. This now includes `BondOrder`, `StereoElementKind`,
  `StereoCarrier`, `StereoDescriptor`, `StereoGroupKind`,
  `ImplicitHydrogens`, the perception model enums, `MoleculeClass`,
  `ResidueClass`, `EntityKind`, `ModelSelection`, `AltLocPolicy`,
  `FitAtoms`, `Weighting`, `Reference`, `ConnectivityCheck`,
  `MatchCompletion`, `QueryPerception`, the property value enums, the
  solvation enums, the Molfile, PDB and SMILES document enums,
  `ChargeSource`, `VdwMethod`, `ElectrostaticsMethod` and
  `MinimizationStatus`. Six enums are closed, and their documentation says
  why: `TetrahedralOrientation`, `DoubleBondOrientation`, `AxisOrientation`,
  `CorrespondenceSide`, `BaseDimension` and `DcdEndian` (the last two lose
  `#[non_exhaustive]`).
- **Breaking:** `bool` parameters and fields that select a mode are enums:
  `AddHydrogensOptions::with_fixed_only(bool)` is
  `with_scope(AddHydrogensScope::{All, FixedOnly})`,
  `ValenceOptions::with_strict(bool)` is
  `with_validation(ValenceValidation::{Strict, Permissive})`,
  `StructureInterpretOptions` and `StructureEnsembleInterpretOptions`
  `with_strict_entity_metadata(bool)` is
  `with_entity_metadata(EntityMetadataPolicy::{InferMissing, RequireDeclared})`,
  `Trajectory::validate_monotonic_time(bool)` takes
  `MissingTimePolicy::{Skip, Reject}`, and `DsspBetaPartner::parallel()` is
  `direction()` (`DsspBridgeDirection::{Parallel, Antiparallel}`). Feature
  switches such as `RotatableBondOptions::with_include_ring_bonds` stay
  `bool`.
- **Breaking:** Results, reports and their records have private fields with
  accessors and no `Default` (ARCHITECTURE.md, "Options, results, and
  enums"); a field `report.added` is now `report.added()`, and collections
  are returned as slices. This covers `AddHydrogensReport`,
  `AddedHydrogen`, `RemoveHydrogensReport`, `RemovedHydrogen`,
  `RetainedHydrogen`, `HydrogenCountAdjustment`, `CipAssignmentReport`,
  `CipAssignment`, `CipSkipped`, `CoordinateStereoResult`,
  `CoordinateStereoMaterializationReport`, `StereoCleanupReport`,
  `ModelAppendReport`, `SolvationReport`, `DsspResidueSource`,
  `DsspHydrogenBond`, `DsspBetaPartner`, `DsspSkippedResidue`,
  `SourceResidueId`, `AltLocResidue`, `AltLocDecision`, the `kekule-openff`
  parameter types (`ParameterIdentity`, `BondParameter`, `AngleParameter`,
  `TorsionTerm`, `TorsionParameter`, `VdwParameter`, `ConstraintParameter`,
  `Interaction`, `PairException`, `NonbondedSettings`, `MoleculeLabels`,
  `ChargeAssignment`), `kekule-potentials` `EnergyComponent`,
  `MinimizationStep` and `ComponentEvaluation`, and `kekule-traj`
  `TrajectoryFieldAvailability`.
  - `SourceResidueId::from_parts` builds a residue identity, for example to
    name an override.
  - `AltLocSelection` and `AtomSiteMetadata` are inputs and follow the
    options style: getters and `with_*` setters, keeping `Default`.
    `AltLocSelection::default` is `default_preference()` /
    `with_default_preference`, and `with_override(residue, label)` adds one
    override.
- **Breaking:** Options structs have private fields, `Default`, getters and
  `#[must_use]` chainable `with_*` setters (ARCHITECTURE.md, "Options,
  results, and enums"); struct literals such as
  `SubstructureMatchOptions { max_matches: 5, ..Default::default() }` become
  `SubstructureMatchOptions::default().with_max_matches(5)?`. Presets stay:
  `RotatableBondOptions::{STRICT, GENERAL}`, `SmilesWriteOptions::{ordinary,
  isomeric, canonical}`, `RingPerceptionOptions::{membership_only, for_graph,
  for_molecule}` and `AromaticityOptions::{for_graph, for_molecule}`. This
  covers every `kekule` options and limits struct, `MinimizeOptions` in
  `kekule-potentials` and `TrajectoryIoLimits` in `kekule-traj`.
  - Setters that can receive an invalid value return `Result` with the error
    the operation used to report, so an invalid options value cannot be
    built: `StereoSearchOptions::{with_max_atoms, with_max_search_states}`
    (`StereoSearchError::InvalidLimit`), `SubstructureMatchOptions::{with_max_matches,
    with_max_search_states, with_max_query_atoms, with_max_candidate_pairs}`
    (`SubstructureMatchError::InvalidLimit`),
    `DsspOptions::with_min_polyproline_stretch`
    (`DsspError::InvalidPolyprolineStretch`),
    `MmcifWriteOptions::{with_block_name, with_coordinate_precision}`,
    `SmilesParseOptions::{with_max_atoms, with_max_bonds}`, every
    `SmartsParseOptions` bound, and
    `StructureEnsembleInterpretOptions::with_model_ids` (all `FormatError`
    with the former kind and detail), and
    `MinimizeOptions::{with_gradient_tolerance, with_max_step, with_memory,
    with_max_line_search_evaluations}` (`MinimizationError`).
  - The `kekule-traj` options gain getters for every setting, and
    `DcdWriteOptions::with_cells(bool)` becomes
    `with_cell_policy(DcdCellPolicy::{Omit, Write})`.
- **Breaking:** Identity lookups return `Option` (ARCHITECTURE.md, "Lookups,
  errors, and panics"); `None` means the ID names nothing in the owner.
  Searches by label or name keep `Result`.
  - `Graph::{atom, bond, neighbors, incident_bonds, bond_between,
    stereo_element, stereo_group}` and the same `QueryGraph` lookups.
    `bond_between` returns `Option<BondId>`, `None` also when either atom is
    missing.
  - `Hierarchy::{chain, residue, atom_site}`, `TopologyBuilder::definition`,
    `TopologyEditor::{atom, bond, atom_handle, bond_handle, chain, residue,
    atom_site, chain_handle, residue_handle}`, `ModelAppend::{atom, bond,
    chain, residue, atom_site}` and `MoleculeEditor::{atom_mut, bond_mut}`.
  - `Model::position`, `ModelView::position`, `ModelBuilder::position` and
    `ModelEditor::position` return `Option<Quantity<Point3>>`.
    `Positions::position_at` is now `Positions::get`, and
    `Positions::set_position_at` is `Positions::set`.
  - Per-atom reads no longer return `Result<Option<_>>`:
    `Molecule`/`MoleculeEditor::{explicit_hydrogens, implicit_hydrogens,
    total_hydrogens, atom_is_aromatic, bond_is_aromatic, atom_is_conjugated,
    bond_is_conjugated, cip_descriptor}` return `Option`, as do
    `occupancy` and `b_factor` on models, views, members and
    `ModelEditor`. `None` covers both an unset value and an ID or index that
    names nothing.
- **Breaking:** Errors follow one policy (ARCHITECTURE.md, "Lookups, errors,
  and panics"). Which inputs fail is unchanged; callers update the error
  types they name and match.
  - Every format failure, SMARTS included, is one
    `diagnostics::FormatError` with `format()` (`diagnostics::Format`),
    `stage()` (`FormatStage::{Parse, Interpret, Write}`), `kind()`
    (`DiagnosticKind`), `location()`, `record()`, `message()`, typed
    `detail()` (`FormatErrorDetail`) and `source()`. It replaces
    `SmilesParseError`, `SmilesInterpretError`, `smiles::SmilesReadError`,
    `SmilesComponentCountError` (now `FormatErrorDetail::ComponentCount`),
    `MolfileParseError`, `MolfileInterpretError`, `MoleculeWriteError` and
    `MoleculeWriteErrorKind`, `SdfParseError`, `SdfInterpretError`,
    `SdfWriteError`, `SdfEnsembleError`, `MmcifParseError`,
    `MmcifWriteError`, `PdbParseError`, `StructureInterpretError`,
    `StructureEnsembleInterpretError` (its variants are `FormatErrorDetail`
    variants), `components::ComponentError`, and `SmartsParseError` and
    `SmartsParseErrorKind`.
  - `location()` is an `Option` for every stage. SMILES and SMARTS errors
    locate the byte column within their line; `SmilesParseError::offset()`
    and `SmartsParseError::span()` are gone (the column is `offset + 1`).
    `record()` is an `Option` and counts SDF records from 1, including the
    ensemble assembly error, which used to count from 0.
  - Writer kinds follow the writer failure taxonomy:
    `MoleculeWriteErrorKind::UnsupportedRepresentation` is
    `DiagnosticKind::Unsupported`, `MissingPerception` is the new
    `DiagnosticKind::MissingPerception`, `ResourceLimit` is
    `DiagnosticKind::Limit`, `InvalidModel` is `DiagnosticKind::Inconsistent`,
    and `Io(kind)` is `DiagnosticKind::Io` with the `std::io::Error` as
    `source()`. SMARTS `Empty` and `InvalidSyntax` are
    `DiagnosticKind::Syntax`, `ResourceLimit` is `DiagnosticKind::Limit`.
  - The mmCIF writer's typed failures are `mmcif::MmcifWriteIssue`, carried
    as `FormatErrorDetail::MmcifWrite`. Its `&'static str` fields are typed:
    `MmcifAtomField`, `MmcifMemberField` and `MmcifCellReason`;
    `UnsupportedTextValue { item }` names the mmCIF item. `InvalidModel`,
    `InvalidHierarchy` and `Io` are plain `Inconsistent` and `Io` errors with
    their cause as `source()`. The never-constructed
    `ClassificationCountMismatch` and `InvalidGroupPdb` are removed.
  - An ensemble failure inside one coordinate model is that model's error
    with `FormatErrorDetail::Model { model_id }` or
    `Conformation { selection }`, not a wrapper; the `"<model inventory>"`
    sentinel model ID is gone.
  - `Display` is `"<format> <stage> error[ in record N][ at LOCATION]:
    message"` and never repeats the source; causes that used to be formatted
    into messages, such as I/O and graph-construction errors, are `source()`.
  - Resource exhaustion is one shape: `diagnostics::LimitExceeded`, with
    `resource()` (`diagnostics::Resource`), `limit()` and `observed()`.
    It replaces `AromaticityError::ResourceLimit { limit }`,
    `RingPerceptionError::ResourceLimit { resource, observed, limit }`,
    `StereoSearchError::ResourceLimit`,
    `SubstructureMatchError::ResourceLimit` (now `LimitExceeded { limit,
    work }`), `QueryExpressionError::ResourceLimit`,
    `QueryGraphError::ResourceLimit`, `HydrogenTransformError::ResourceLimit`,
    `ResonanceError::ResourceLimit`, `CipRankingError::DepthLimitExceeded`
    and `ResourceLimitExceeded`, `CipAssignmentIssue::DepthLimitExceeded` and
    `ResourceLimitExceeded` (now `LimitExceeded { element, limit }`),
    `SymmetryError::PermutationLimit` and `SearchLimit`,
    `NormalizationError::AromaticLocalizationLimit { examined_states, .. }`,
    `DsspError::ResourceLimitExceeded` with `dssp::DsspResource`,
    `PeriodicGeometryError::ImageSearchLimit`, `SolvationError::CandidateLimit`
    and `kekule_traj::PeriodicError::ImageSearchLimit`, each now a
    `LimitExceeded` variant carrying it. A format failure of kind
    `DiagnosticKind::Limit` carries it as `FormatErrorDetail::LimitExceeded`
    and `FormatError::limit_exceeded()`.
  - Two payloads were wrong: a query graph over its 4,096-atom bound
    reported the 65,536 expression-node bound, and a SMARTS bond expression
    over its depth bound reported its node count against the node bound.
  - Default limits are constants: `DEFAULT_*` associated constants on
    `AromaticityOptions`, `RingPerceptionOptions`, `StereoSearchOptions`,
    `SubstructureMatchOptions`, `CipAssignmentOptions`, `AddHydrogensOptions`,
    `ResonanceOptions` (with `MAX_STRUCTURES`), `SymmetryLimits`,
    `DsspLimits`, `SolventOptions`, and the SMILES, SMARTS, Molfile, SDF,
    mmCIF and PDB parse options; built-in bounds are
    `query::{MAX_QUERY_RECURSION_DEPTH, MAX_QUERY_TOTAL_ATOMS,
    MAX_QUERY_TOTAL_EXPRESSION_NODES}`, `geometry::MAX_LATTICE_IMAGES`, and
    the canonical and isomeric SMILES writer bounds
    `smiles::MAX_CANONICAL_*` and `smiles::MAX_DIRECTIONAL_*`.
  - A failed `try_build` or `try_finish` returns `kekule::DraftError<D, E>`
    with `error()`, `draft()`, `into_draft()` and `into_parts()`. It replaces
    `MoleculeFinishError` (`DraftError<MoleculeEditor,
    MoleculePublicationError>`), `TopologyBuilderError`
    (`DraftError<TopologyBuilder, TopologyBuildError>`),
    `TopologyFinishError` (`DraftError<TopologyEditor, TopologyEditError>`),
    `ModelBuilderError` (`DraftError<ModelBuilder, ModelBuildError>`) and
    `ModelFinishError` (`DraftError<ModelEditor, ModelEditError>`);
    `editor()`/`builder()` are `draft()` and `into_editor()`/`into_builder()`
    are `into_draft()`. Its `Display` no longer repeats the error, which is
    its `source()`.
  - `MoleculeError::InvalidStereoReference` carries a
    `stereo::StereoValidationIssue` instead of one of 24 `&'static str`
    reasons; elements and groups being added are named by the ID they would
    receive. `StereoValidationIssue` gains the group and membership issues
    (`GroupMembershipOutsideGroupEdit`, `MissingStereoGroup`,
    `EmptyStereoGroup`, `DuplicateStereoGroupMember`,
    `MissingStereoGroupMember`, `InconsistentStereoGroupMembership`,
    `StereoElementAlreadyGrouped`, `SplitCorrelatedStereoGroup`), `Display`,
    and `#[non_exhaustive]`. `MoleculeError::UnsupportedFeature`, which was
    never constructed, and `impl Eq for MoleculeError`, unsound because the
    error can hold an `f64` property value, are removed; types that wrap it
    (`HydrogenTransformError`, `CoordinateStereoError`, `NormalizationError`
    and `SourceStereoNormalizationError`) are no longer `Eq` either.
  - `MoleculePublicationError` is the one publication error.
    `MoleculeConnectivityError` becomes `Disconnected { components }`, and
    `StereoPublicationError` becomes `InvalidStereo(StereoValidationIssue)`,
    which now names the specific issue instead of
    `InvalidElementReference`. `GraphValidationError` and the
    `InvalidGraph` variant are removed: the public API cannot corrupt graph
    storage, so publication checks it only as a debug assertion.
  - `PerceptionInstallError::InvalidDelocalization` carries a
    `perception::DelocalizationIssue` instead of a `&'static str`.
  - Other `&'static str` and stringified payloads are typed:
    `StereoSearchError::InvalidOptions` and
    `SubstructureMatchError::InvalidOptions` are `InvalidLimit { resource,
    limit }`; `QueryGraphError::InvalidStereo` carries a
    `query::QueryStereoIssue`; `ResonanceError::Materialization(String)` is
    `ContributorOutOfRange { index, count }`, `Molecule(MoleculeError)` or
    `Publication(MoleculePublicationError)`; the `AppendProperty { domain }`
    of `TopologyEditError` and `ModelEditError` is a
    `properties::PropertyDomain`; DSSP backbone atoms are
    `dssp::DsspBackboneAtom` (`AmbiguousBackboneAtom { atom }`,
    `DsspSkipReason::MissingBackboneAtoms { missing }`), and the
    single-valued `quantity` fields of `DegenerateBackboneGeometry`,
    `NonFiniteGeometry` and `CoordinateOutOfRange` are removed;
    `ConnectionResolutionReason::ConflictingSelectorValues` carries a
    `macromolecular::ConflictingSelector`.
  - Removed variants that could never be constructed:
    `AromaticityError::UnsupportedElement`,
    `ValenceIssue::UnsupportedElement` (the valence rules cover every
    element) and `DsspError::InvalidHierarchy`.
  - Wrapping errors report their cause only through `source()`: `Display`
    gives the context alone. This covers `PerceptionError`,
    `AromaticityError`, `StereoSearchError`, `CoordinateStereoError`,
    `HydrogenTransformError`, `TaggedMatchError` (now `#[non_exhaustive]`),
    `ResonanceError`, `SymmetryError`, `AlignmentError`, `ClusteringError`,
    `SolvationError`, `PeriodicBoxError`, `PeriodicCellError`,
    `PropertyError`, `PositionError`, `StructureError`, `ModelError`,
    `ModelBuildError`, `ModelEditError`, `EnsembleError`, `TrajectoryError`,
    `MeasurementError`, `BondDihedralError`, `GeometryEditError`,
    `TopologyBuildError`, `TopologyEditError` and `TopologySubsetError`.
    `HydrogenTransformError`, `TaggedMatchError` and `ResonanceError` gain
    a `source()`.
  - One error per functional area. `HierarchyLookupError` folds into
    `topology::HierarchyError` (`NotFound`, `Ambiguous`).
    `TopologyHierarchyError` and `TopologyBuildError::InvalidHierarchy` are
    removed: the one case the public API can reach, an atom site on an atom
    the topology lacks, is `TopologyBuildError::MissingHierarchyAtom { site,
    atom }`, and the other ten storage inconsistencies are debug assertions.
    `PositionError` becomes `structure::DenseArrayError`
    (`IndexOutOfRange { index, len }`, `CountMismatch`, `NonFinite { index }`,
    `Unit`), which positions, velocities, forces and per-atom scalars share;
    every `index` is the dense row. `Velocities` and `Forces` constructors
    return it instead of `StructureError`, and `StructureError`'s
    `AtomCountMismatch`, `InvalidAtomIndex`, `NonFiniteValue`,
    `NonFiniteVector` and `Position` become `StructureError::Array`. A
    non-finite position given to the model editor is
    `ModelEditError::NonFinitePosition { atom }` instead of a fabricated
    row 0, and `Model::set_atom_positions` reports the dense row instead of
    the input position. `ModelError::InvalidBondId`, never constructed, is
    removed. `CipAssignmentIssue` and `ValenceIssue` are
    `#[non_exhaustive]`.
  - `kekule_traj::TrajectoryIoError` is one struct with `kind()`
    (`TrajectoryIoErrorKind`), `operation()`, `format()`, `source_label()`,
    `frame()`, `byte_offset()`, `expected()`, `actual()`, `limit_exceeded()`,
    `field()`, `detail()`, `io_error_kind()` and `source()` (the
    `std::io::Error` or Kekule error), and public constructors for
    `TrajectoryReader` and `TrajectoryWriter` implementors. It replaces the
    `TrajectoryIoError` enum, `TrajectoryIoErrorContext`,
    `TrajectoryCodecErrorContext` and `TrajectoryCodecErrorKind`
    (`ResourceLimitExceeded` is `LimitExceeded`); `Io(_)`, `Codec(_)`,
    `Structure(_)`, `Trajectory(_)`, `TopologyMismatch` and
    `FrameIndexOutOfRange(n)` are kinds `Io`, the codec kinds,
    `InvalidStructure`, `InvalidTrajectory`, `TopologyMismatch` and
    `FrameIndexOutOfRange` with `frame()`. `UnsupportedField(&'static str)`
    is kind `UnsupportedField` with `field()` (`kekule_traj::TrajectoryField`),
    which writers now also set; the XTC and TRR property messages say
    "frame properties". The XYZ frame-byte limit carries its
    `LimitExceeded`. The never-constructed `AtomOrderMismatch`,
    `UnsupportedRandomAccess`, `MissingRequiredTopology` and
    `TrajectoryIoOperation::WriteHeader` are removed.
  - `kekule_openff::Error` returns the wrapped Kekule, I/O, JSON or OFFXML
    error from `source()` (it used to skip to that error's own source) and
    its `Display` gives `ErrorKind`'s description instead of repeating the
    wrapped message; `ErrorKind` gains `Display`. `Error::limit_exceeded()`
    gives the exceeded bound of `ErrorKind::ResourceLimit` failures.
  - `kekule_potentials`: `EvaluationError::Backend { backend, message }` is
    `Backend(Arc<dyn Error + Send + Sync>)`, built with
    `EvaluationError::backend`, with the cause as `source()`;
    `InvalidGeometry::interaction` is an `InteractionKind`;
    `MinimizationError::InvalidOptions` carries a `MinimizeOption` and
    `NumericalFailure` a `NumericalFailure`;
    `OpenFfPotentialError::UnsupportedNonbondedMethod { handler, method }` is
    `UnsupportedVdwMethod(VdwMethod)` or
    `UnsupportedElectrostaticsMethod(ElectrostaticsMethod)`;
    `InvalidParameter { parameter, message }` is `InvalidParameter {
    parameter: openff::OpenFfParameter, value }`, `ZeroTorsionDivisor`,
    `SelfPairException` or `DuplicatePairException`; and
    `AtomCountMismatch::parameter` is `array: openff::PerAtomArray`.
    `MinimizationError`, `EvaluationError` and `OpenFfPotentialError` no
    longer repeat their source in `Display`.
- **Breaking:** Public items live in the module of the domain that owns them,
  and every operation has one home. Behaviour is unchanged; callers update
  paths, names, and signatures as listed.
  - `kekule::core` is renamed `kekule::molecule`.
  - Represented stereo moves from `kekule::core` to `kekule::stereo`:
    `StereoElement`, `StereoElementId`, `StereoElementKind`,
    `TetrahedralStereo`, `DoubleBondStereo`, `AxisStereo`, `StereoCarrier`,
    `StereoGroup`, `StereoGroupId`, `StereoGroupKind`, `StereoDescriptor`,
    `TetrahedralOrientation`, `DoubleBondOrientation`, and
    `AxisOrientation`. `validate_stereo`, `StereoValidationError`, and
    `StereoValidationIssue` stay in `kekule::stereo`.
  - `kekule::stereo` splits its algorithms into submodules.
    `stereo::cip` holds `assign_cip_descriptors`, `CipAssignment`,
    `CipAssignmentError`, `CipAssignmentIssue`, `CipAssignmentOptions`,
    `CipAssignmentReport`, `CipRankingError`, `CipSkipped`, and
    `CipSkippedReason`. `stereo::coordinates` holds
    `infer_coordinate_stereo`, `materialize_coordinate_stereo`,
    `CoordinateStereoError`, `CoordinateStereoMaterializationReport`,
    `CoordinateStereoOptions`, and `CoordinateStereoResult`.
    `stereo::cleanup` holds `detect_stereo_candidates`, `cleanup_stereo`,
    `StereoCandidate`, and `StereoCleanupReport`.
  - `kekule::core::StereoPerception`, the installed CIP state, is
    `kekule::stereo::cip::CipPerception`. `stereo::StereoPerceptionOptions`
    and `stereo::StereoPerceptionError` are `stereo::StereoSearchOptions` and
    `stereo::StereoSearchError`.
  - `Perception`, `PerceptionBuilder`, `PerceptionBuildError`,
    `PerceptionComponent`, and `PerceptionInstallError` move from
    `kekule::core` to `kekule::perception`, and `MalformedRingReason` to
    `kekule::perception::rings`.
  - SMARTS is a format module: `smarts::parse_str(input, options)` replaces
    `query::parse_smarts` and `query::parse_smarts_with_options` and returns
    a `query::QueryGraph`. `SmartsParseOptions`, `SmartsParseError`, and
    `SmartsParseErrorKind` move from `kekule::query` to `kekule::smarts`.
  - Solvation is a module of its own. `kekule::solvation` holds
    `SolventOptions`, `SolvationReport`, `SolvationRadii`,
    `SolvationChargeBasis`, `SolvationError`, `WaterModel`, `PositiveIon`,
    `NegativeIon`, `BoxShape`, `PeriodicBoxOptions`, and `PeriodicBoxError`,
    formerly in `kekule::structure`. `Model::add_solvent(&options)` is
    `solvation::solvate(&mut model, &options)` and
    `Model::add_periodic_box(&options)` is
    `solvation::add_periodic_box(&mut model, &options)`; both still edit the
    model in place and are atomic.
  - The topology editor moves to `kekule::topology::edit`: `TopologyEditor`,
    `EditAtomId`, `EditBondId`, `EditChainId`, `EditResidueId`,
    `EditAtomSiteId`, `EditBond`, `EditMolecule`, `EditChain`,
    `EditResidue`, `EditAtomSite`, `EditCorrespondence`, `TopologyEditError`,
    and `TopologyFinishError`.
  - `kekule::topology::transform` is no longer a module. `TopologySubset`,
    `TopologySubsetCorrespondence`, and `TopologySubsetError` are in
    `kekule::topology`. `transform::retain_instances`,
    `transform::remove_instances`, and `TopologyTransformError` are removed;
    `Topology::subset` with `AtomSelection::for_instances` (and its
    `complement()`) selects whole instances.
  - `topology::MoleculeInstanceView` is `topology::InstanceView`.
    `topology::MoleculeInstance` is no longer public:
    `InstanceView::instance()` is removed (use `id()` and
    `definition_id()`), and `TopologyBuilder::instances()` yields
    `(MoleculeInstanceId, MoleculeDefinitionId)` pairs.
  - `kekule::canon` and `kekule::rotatable_bonds` are removed.
    `canon::atom_ranking` is `descriptors::canonical_atom_ranking`, and
    `rotatable_bonds::detect` is `descriptors::rotatable_bonds`.
    `CanonicalAtomRanking`, `RotatableBondOptions`, and `RotatableBondSet`
    move to `kekule::descriptors`.
  - `Molecule::edit`, `Model::edit`, and `Topology::edit`, which copy the
    owner into a draft, are `to_editor`; `into_editor` is unchanged.
  - Every operation takes its options as a parameter, and the
    `*_with_options` twins are removed; pass `Default::default()` for the
    former plain behaviour. This applies to `perception::valence::perceive_valence`,
    `perception::aromaticity::perceive_aromaticity`,
    `stereo::cip::assign_cip_descriptors`, `stereo::cleanup::detect_stereo_candidates`,
    `stereo::coordinates::infer_coordinate_stereo` and
    `materialize_coordinate_stereo`, `Molecule::add_hydrogens`,
    `alignment::rmsd`, `Ensemble` and `Trajectory` `superpose`, `rmsd`, and
    `aligned_rmsd`, and `Ensemble::cluster_by_radius`.
    `structure::measure::BondDihedral::new(topology, bond, options)` replaces
    `BondDihedral::new` and `with_options`.
  - `perception::rings::perceive_rings(mol, options)` replaces
    `perceive_ring_membership`, `perceive_ring_set`, and
    `perceive_ring_set_with_options`. It installs ring perception and returns
    `()`, like the other perception steps; read the result through
    `Molecule::ring_membership` and `Molecule::ring_set`. The new
    `RingPerceptionOptions::scope` field (`RingPerceptionScope::Basis` by
    default) selects membership-only perception with
    `RingPerceptionScope::Membership`, or `RingPerceptionOptions::membership_only()`,
    which ignores the bounds and cannot fail.
    `perceive_aromaticity_with_ring_options(mol, model, ring_options)` is
    `perceive_aromaticity(mol, model, AromaticityOptions { ring_options, .. })`.
  - `RingPerceptionOptions::default()` and `AromaticityOptions::default()`
    keep their fixed bounds. The bounds the former plain `perceive_ring_set`
    and `perceive_aromaticity` used, which grow with the molecule, come from
    the now public `RingPerceptionOptions::for_molecule` and
    `AromaticityOptions::for_molecule`; the two differ only for molecules
    above about 5,000 atom and bond slots.
  - Substructure search runs only on prepared targets. The free functions
    `substructure::find_match`, `find_matches`, `find_matches_with_options`,
    `visit_matches`, `visit_matches_with_options`, `find_topology_matches`,
    and `find_topology_matches_with_options` are removed.
    `PreparedTarget` and `PreparedTopologyTarget` provide
    `first(query, options)`, `all(query, options)`, and
    `visit(query, options, visitor)`, replacing `find_match`, `find_matches`,
    `find_matches_with_options`, `visit_matches`, and
    `visit_matches_with_options`; `PreparedTopologyTarget::first` is new.
    `TaggedQuery::find_matches` takes `substructure::TupleMode::All` or
    `TupleMode::Distinct` instead of a `bool`.
  - `alignment::kabsch` and `kabsch_with_options` are `alignment::fit`. The
    new free `alignment::aligned_rmsd` fits one view on one set of pairs and
    measures another, as the collection method does per item.
    `Ensemble::pairwise_rmsd` and `pairwise_rmsd_with_options` are
    `Ensemble::pairwise_fitted_rmsd(atoms, options)`.
  - `dssp::assign(model, options)` takes any `&impl AsModelView`, replacing
    `dssp::assign(view, options)` and `dssp::assign_model(&model, options)`.
  - `structure::measure::distance`, `angle`, and `dihedral` take a
    `ConnectivityCheck`, replacing their `*_with_connectivity` variants;
    `ConnectivityCheck::Unrestricted` is the former plain behaviour.
  - Format parsing takes options: `parse_str(input, options)` and
    `parse_reader(reader, options)` in `smiles`, `molfile`, `sdf`, `mmcif`,
    and `pdb`, and `records(reader, options)` in `smiles` and `sdf`. The new
    `smiles::parse_reader` reads one record within the byte limit.
  - Interpretation is a document method taking options:
    `doc.interpret(options)` and `doc.interpret_ensemble(options)` on
    `SmilesDocument`, `MolfileDocument`, `SdfRecord`, `MmcifDocument`,
    `MmcifBlock`, and `PdbDocument`, replacing the plain and
    `*_with_options` methods. The free functions `smiles::interpret`,
    `interpret_with_options`, and `interpret_base`, `molfile::interpret` and
    `interpret_with_options`, and `mmcif::` and `pdb::interpret`,
    `interpret_with_options`, `interpret_ensemble`,
    `interpret_ensemble_with_options`, and `interpret_conformations` are
    removed; call the document method instead.
  - The document shortcuts `to_molecules`, `to_topology`, and `to_model`, and
    their `*_with_options` variants, are removed from every format document;
    use `doc.interpret(options)?.into_molecules()` (or `into_topology`,
    `into_model`). `smiles::to_molecules` and `smiles::to_topology` remain
    the one-call SMILES shortcuts.
  - `kekule-traj`: `io::read_trajectory`, `open_trajectory`,
    `open_indexed_trajectory`, and `write_trajectory` take their options, and
    their `*_with_options` twins are removed. `TrajectoryWriteOptions`
    implements `Default`, which takes the format from the destination
    extension as plain `write_trajectory` did; `create_trajectory_writer`
    accepts it too. `analysis::FrameSuperposer::new(reference, atoms, options)`
    replaces `new` and `with_options`. `periodic::unwrap` is
    `periodic::unwrap_trajectory`, and `FrameBufferData` is `FrameData`.
  - `kekule-traj`: the test helpers `CoordinateFrameReader` and
    `validate_atom_order` are removed, and `MemoryTrajectoryReader` is no
    longer public; fill a `FrameBuffer` from an in-memory frame with
    `FrameBuffer::copy_from`.
  - `kekule-openff`: the `diagnostics` module is removed from the public API.
    The tests that compared lookup keys, input features, and network charges
    with reference data run inside the crate.
- **Breaking:** Every public module is an explicit list of re-exports, and
  every public item has exactly one public path; only `prelude` repeats items.
  No item is renamed and no behaviour changes, so callers only update imports.
  Items that no public path reached are now crate-private, which removes no
  public path.
  - The perception state types are exported only from their
    `kekule::perception` module, no longer from `kekule::core`:
    `valence::{ValenceModel, ValencePerception}`,
    `rings::{Ring, RingSet, RingMembership, RingBasisModel, RingBasisState,
    RingPerception}`, `aromaticity::{AromaticityModel,
    AromaticityPerception}`, `conjugation::{ConjugationModel,
    ConjugationPerception}`, and `resonance::{ResonanceGroup,
    ResonancePerception}`. `ValencePerception`, `RingPerception`, and
    `AromaticityPerception` were reachable only through `kekule::core`.
  - `kekule::rotatable_bonds::PerceptionError` is removed; use
    `kekule::perception::PerceptionError`.
  - `kekule::core::Result`, an alias of `Result<T, MoleculeError>`, is no
    longer public.
  - `kekule::prelude` holds `Molecule`, `MoleculeEditor`, `Atom`, `Bond`,
    `BondOrder`, `Element`, `ImplicitHydrogens`, `Topology`, `Model`,
    `Positions`, and `Quantity`. `AtomId`, `BondId`, and `Hierarchy` left it;
    import them from `kekule::core` and `kekule::topology`.
  - The 25 interpretation types shared by mmCIF and PDB are exported only
    from `kekule::macromolecular`, no longer from `kekule::mmcif` and
    `kekule::pdb`: `AltLocDecision`, `AltLocPolicy`, `AltLocPreference`,
    `AltLocResidue`, `AltLocSelection`, `AltLocSelectionReason`,
    `AtomProvenance`, `BondEvidence`, `BondProvenance`, `CellOmissionReason`,
    `ConnectionResolutionReason`, `ConnectivityOptions`, `EntityKind`,
    `InstanceProvenance`, `ModelSelection`, `SourceResidueId`,
    `SourceResiduePosition`, `StructureInterpretOptions`,
    `StructureInterpretation`, `StructureInterpretationReport`,
    `StructureInterpretIssue`, `StructureInterpretError`,
    `StructureEnsembleInterpretOptions`, `StructureEnsembleInterpretation`,
    and `StructureEnsembleInterpretError`.
  - `UnrepresentablePolicy` is exported only from `kekule::diagnostics`, no
    longer from `kekule::mmcif`, `kekule::molfile`, and `kekule::sdf`.
  - `MoleculeWriteError` and `MoleculeWriteErrorKind` are exported only from
    `kekule::molfile`, no longer from `kekule::smiles`, and
    `MoleculeWriteError`, `MolfileInterpretOptions`, and `MolfileWriteVersion`
    no longer from `kekule::sdf`.
  - `kekule-traj` and `kekule-openff` list their exports explicitly; none of
    their public paths changes.
- **Breaking:** `kekule-traj` indexed readers are opened directly with
  `IndexedXyzReader::new`, `IndexedDcdReader::new`, `IndexedTrrReader::new`,
  and `IndexedXtcReader::new`, which take the same arguments as the
  sequential constructors. The sequential readers' `into_indexed` is removed:
  after a sequential read it indexed the wrong frames. TRR and XTC started
  the index at the current frame, so three frames read once indexed as two
  and `read_frame(0)` returned frame 1. DCD failed on its declared frame
  count, and XYZ indexed from byte 0 of the stream instead of the position
  the reader was opened at. Each indexed reader now takes its opening
  position as frame 0.
- `kekule-traj` no longer depends on `molly`. The XTC `variant` in
  trajectory metadata no longer names the writer and reads, for example,
  `Xtc1995; first payload 92 bytes`.
- **Breaking:** `CoordinateStereoResult` and
  `CoordinateStereoMaterializationReport` gain `disagreements`, the asserted
  stereo elements whose configuration the coordinates contradict. Inference
  still fills only foci without an element and leaves asserted ones, including
  explicitly unknown configurations, unchanged; that contract is now documented.
  Unknown configurations and flat or degenerate geometry never disagree.
- **Breaking:** `MolfileInterpretationWarning` is `#[non_exhaustive]` and gains
  `AtomParityContradictsCoordinates`.
- **Breaking:** `MoleculeWriteErrorKind` says why a SMILES or Molfile writer
  failed. The new `MissingPerception` kind marks a source that must be
  perceived before writing, which used to be reported as
  `UnsupportedRepresentation`. Canonical SMILES classifies a failed hydrogen
  normalization by its cause instead of always blaming missing perception.
  `UnsupportedRepresentation` now means only that the format cannot encode the
  content. `InvalidMetadata` is removed: SDF titles and data fields have always
  been rejected as `SdfWriteError::InvalidTitle` and `InvalidDataField`.
- **Breaking:** The interpretation types shared by macromolecular formats have
  format-neutral names and a namespace of their own, `kekule::macromolecular`,
  which is their only public path. `MmcifInterpretOptions`,
  `MmcifInterpretation`, `MmcifInterpretationReport`, `MmcifInterpretIssue`,
  and `MmcifInterpretError` are now `StructureInterpretOptions`,
  `StructureInterpretation`, `StructureInterpretationReport`,
  `StructureInterpretIssue`, and `StructureInterpretError`, and the ensemble
  types `MmcifEnsembleInterpret{Options,ation,Error}` are now
  `StructureEnsembleInterpret{Options,ation,Error}`. The `Mmcif` prefix is
  dropped from `AltLocPolicy`, `AltLocPreference`, `AltLocSelection`,
  `AltLocSelectionReason`, `AltLocDecision`, `AltLocResidue`, `ModelSelection`,
  `ConnectivityOptions`, `AtomProvenance`, `InstanceProvenance`,
  `BondEvidence`, `BondProvenance`, `ConnectionResolutionReason`, and
  `EntityKind`, and `MmcifResidueId` and `MmcifResiduePosition` are now
  `SourceResidueId` and `SourceResiduePosition`. Types specific to mmCIF
  documents and writing keep their prefix.
- **Breaking:** mmCIF files without `_chem_comp_bond` rows, such as wwPDB
  archive entries, no longer fall apart into one molecule per unbonded atom:
  residue bonds come from the bundled or caller templates, and residues no
  template covers have their bonds inferred from covalent radii. A heavy atom
  whose name its residue template lacks is an error. Formal charges not stated
  by `_atom_site.pdbx_formal_charge` come from the template, so lysine,
  arginine, and histidine are cations in hydrogen-free models; model hydrogens
  adjust heteroatom charges and select histidine tautomers.
- **Breaking:** Metal coordination (`_struct_conn` type `metalc`, and template
  bonds flagged `pdbx_metal_coordination_flag`) becomes `Dative` bonds, so a
  metal site and its ligands form one molecule. The mmCIF writer emits dative
  bonds as `metalc` connections instead of rejecting them, and names water
  atoms `O`, `H1`, `H2`.
- Polymer links between residues adjacent in sequence require their linking
  atoms to lie within bonding distance; longer gaps are reported as chain
  breaks. Without `_entity_poly`, the residues' templates decide the polymer
  type: a residue without a polymer template takes its neighbour's linkage,
  and two such residues link by atom names.
- **Breaking:** `SourceResiduePosition` gains `RepeatedAuthor` for residue
  numbers that recur within a chain, and interpretation errors no longer name
  mmCIF in their messages.
- **Breaking:** mmCIF documents keep one line-normalized copy of the source
  and store every tag and value as a compact span into it, instead of one
  owned string per value. `MmcifEntry`, `MmcifItem`, `MmcifLoopTable`, and
  `MmcifValue` are now borrowed views with a lifetime: `MmcifBlock::entries()`
  and `MmcifLoopTable::row()` return iterators, values are returned by value,
  and `MmcifLoopTable::value_at` reads a cell by a column resolved once with
  `column_index`. Parsing a 400,000-row atom-site table is about four times
  faster.
- **Breaking:** `MmcifParseError` and `StructureInterpretError` report `kind()` and
  `location()` instead of `line()`. `StructureEnsembleInterpretError` gains
  `kind()`, `location()`, and an error `source()`.
- **Breaking:** Molfile and SDF record interpretation rejects unsupported
  content by default instead of ignoring it: V2000 property records other
  than `M  CHG`, `M  ISO`, and `M  RAD` (such as S-groups), atom aliases,
  values, and group abbreviations; V3000 atom options other than `CFG`,
  `CHG`, `MASS`, `RAD`, `HCOUNT`, and `VAL`, bond options other than `CFG`,
  collections other than enhanced stereo groups, blocks such as `SGROUP`, and
  other CTAB records; and text after `M  END`. Previously only V3000 atom and
  bond options failed, and the rest was listed by line in
  `MolfileInterpretationReport::ignored_record_lines`, which is removed.
  V3000 option values in parentheses or quotes are one field, so list-valued
  options such as `RGROUPS=(1 1)` are unsupported content rather than syntax
  errors.
- **Breaking:** Molfile codes are decoded while parsing. Unknown bond types,
  stereo codes that do not fit the bond order, and out-of-range mass numbers,
  `HCOUNT`, or `VAL` now fail parsing instead of interpretation, and negative
  mass numbers are rejected instead of ignored. `MolfileDocument` no longer
  keeps a copy of its source: `source()`, `atom_records()`, `bond_records()`,
  `property_records()`, `unsupported_records()`, and `MolfileLine` are
  removed.
- **Breaking:** `MolfileParseError`, `MolfileInterpretError`,
  `SdfParseError`, and `SdfInterpretError` report `kind()` and `location()`
  instead of `line()`. Interpretation failures that concern no source line
  have no location instead of a placeholder line. `SdfInterpretErrorKind` is
  removed; the Molfile failure is the error `source()`.
- **Breaking:** Locations inside SDF records are lines of the SDF input,
  where Molfile-level lines used to count from the record start. Atom and
  bond mappings report `location()` instead of `source_line()`, and Molfile
  interpretation warnings carry a `location`. `SdfRecord` reports `number()`
  and `location()` instead of `source_record_number()` and
  `source_start_line()`, `SdfDataField::line()` is now `location()`, and
  `SdfRecordInterpretationReport::source_start_line()` is now `location()`.
- **Breaking:** SDF records are interpreted individually, so one failing
  record no longer fails its siblings: `sdf::interpret`,
  `SdfDocument::interpret`, and `SdfInterpretation` are removed in favor of
  `SdfRecord::interpret`. Records are parsed from the input lines instead of
  re-joining each record's Molfile text. `MolfileInterpretation::into_parts`
  also returns the omitted content.
- **Breaking:** Every format namespace pairs `interpret(document)` with
  `interpret_with_options(document, options)`. `mmcif::interpret` and
  `mmcif::interpret_ensemble` now take only the document; pass options to
  the new `mmcif::interpret_with_options` and
  `mmcif::interpret_ensemble_with_options`. `mmcif::interpret_block`,
  `interpret_ensemble_block`, and `interpret_conformations_block` are removed;
  call the same-named methods on `MmcifBlock` instead.
- **Breaking:** `MolWriteError` and `MolWriteErrorKind`, returned by the
  Molfile and SMILES writers, are renamed `MoleculeWriteError` and
  `MoleculeWriteErrorKind`.
- **Breaking:** SMILES records are parsed in one pass that produces both the
  document tokens and the parsed program, and CXSMILES fields are read while
  parsing instead of during interpretation. A malformed extension field, an
  invalid radical code, or an out-of-range group number now fails parsing,
  so `interpret_base` no longer projects records whose extension is
  malformed. Unknown and unsupported fields still parse and fail
  interpretation by default. `SmilesDocumentTokenKind::Unsupported`, which
  no parsed document contained, is removed.
- **Breaking:** `SmilesParseError` and `SmilesInterpretError` report `kind()`
  and `location()` (the record's input line and the byte column) besides the
  record-relative `offset()`, and display the location instead of the offset.
- The default `MmcifParseOptions` limits are consistent with each other: the
  token budget admits the default 5,000,000 atom-site rows of a typical
  21-column table, and the input limit is 1 GiB. Previously the default token
  budget rejected atom-site tables above about 476,000 rows.
- `Topology` is one flat owner. Instances, dense orders, hierarchy, and
  properties sit beside the definitions instead of in a shared allocation,
  and layout identity is an explicit token minted at publication. Perception
  snapshots still share static storage, including properties until one
  snapshot edits them.
- **Breaking:** `Conformation` and `ConformationMut` are removed. `Model`,
  `EnsembleMember`, and `TrajectoryFrame` own their structure fields
  directly, and each holds only the fields meaningful for it: all three have
  positions, an optional cell, and structure properties; models and members
  also have occupancies and B-factors; frames also have velocities, forces,
  time, and step.
  - `Model::new`, `EnsembleMember::new`, and `TrajectoryFrame::new` take
    `Positions`; each type has `with_cell`.
  - Read structure state from the owner (`positions()`, `cell()`,
    `properties()`, `occupancies()`, ...) and change it with the owner's
    dimension-preserving setters (`set_positions`, `set_position_at`,
    `set_cell`, `set_properties`, `properties_mut`, `set_occupancy`, ...).
    The same setters are on the member and frame editors returned by
    `get_mut`, and on `ModelBuilder`. The index-based `set_position` of the
    former `ConformationMut` is now `set_position_at`, beside the existing
    ID-based `Model::set_position`.
  - `conformation()`, `conformation_mut()`, `into_conformation()`,
    `Model::replace_conformation`, `Model::into_parts`, and `ModelView::new`
    are removed. `ModelView` no longer dereferences to a conformation; it has
    its own accessors, and a view of a trajectory frame reports no
    occupancies or B-factors. `Realization` no longer requires `Deref`.
  - Trajectory frames cannot hold occupancies or B-factors, so trajectory
    writers no longer reject them, and `FrameBufferData` borrows only
    structure properties from a source frame.
- **Breaking:** `Ensemble` and `Trajectory` each own their topology,
  collection properties, and items and implement their own operations, so
  the two can diverge. The generic `Realization` trait and the
  `RealizationView`, `RealizationMut`, and `RealizationIter` types are
  removed. `EnsembleMemberView`, `EnsembleMemberMut`, `TrajectoryFrameView`,
  and `TrajectoryFrameMut` are now concrete types with the same methods, and
  `iter()` returns `EnsembleIter` or `TrajectoryIter`. Use
  `EnsembleMemberView::new` or `TrajectoryFrameView::new` in place of
  `RealizationView::new`.
- **Breaking:** "Structure" is the common term for a model, ensemble member,
  or trajectory frame, replacing "realization" and "conformation":
  `RealizationProperties` and `RealizationPropertiesMut` are now
  `StructureProperties` and `StructurePropertiesMut`,
  `ModelAtomView::realization_property` is `structure_property`, and
  `ConformationError` is `StructureError`. Error variants that wrapped it are
  renamed from `Conformation` to `Structure`.
- **Breaking:** `RealizationError` is split into `EnsembleError` and
  `TrajectoryError`, each with only the variants its collection can produce;
  `ItemCountMismatch` becomes `MemberCountMismatch` or `FrameCountMismatch`.
  The `kekule-traj` error formerly named `TrajectoryError` is now
  `TrajectoryIoError`, and its `Realization` variant is `Trajectory`.

### Fixed

- Structure writers no longer silently drop what an item is. `ModelView`
  gains `kind()`, a `StructureKind` of `Model`, `Member { weight, energy }`
  or `Frame { time, step, velocities, forces }`. Before, `sdf::write(&trajectory, …)`
  dropped each frame's time and step, and `SdfRecordSource::model(&member)`
  dropped the weight and energy that `SdfRecordSource::from(member)` wrote.
  - SDF writes a member's weight and energy as data fields however its view
    is passed, and fails on a frame's time, step, velocities or forces
    unless the new `SdfWriteOptions::with_frame_metadata(UnrepresentablePolicy::Omit)`.
  - Molfile fails on an ensemble member (its weight is relative to members
    the record does not hold, and it may have an energy) and on frame state,
    unless `MolfileWriteOptions::with_ensemble_metadata` or
    `with_frame_metadata` omits it.
  - An mmCIF model block of an ensemble member fails unless
    `MmcifWriteOptions::with_ensemble_metadata` omits its weight and energy,
    as ensemble blocks with unequal weights already did; a frame's state
    fails unless the new `with_frame_metadata` omits it.
  - The failure is `DiagnosticKind::Unsupported` with the new
    `FormatErrorDetail::UnrepresentableState(StructureField)`.
- CIP assignment no longer ranks an unperceived inferred hydrogen count as
  zero hydrogens. Unperceived `C[C@@H](O)[CH3]` was labelled R, although the
  centre carries two methyls. `assign_cip_descriptors` now fails with the new
  `CipAssignmentIssue::UnknownHydrogenCount { atom }` for every atom whose
  implicit count is inferred but not perceived, and leaves perception
  unchanged; run valence perception first. An unknown count found while
  ranking is reported the same way instead of as `UnresolvedPriority`.
- Canonical SMILES of a molecule with a ring now requires installed
  aromaticity perception and otherwise fails with
  `MoleculeWriteErrorKind::MissingPerception`. Unperceived, the two Kekule
  forms of o-xylene were written as two different "canonical" strings,
  `CC1=C(C)C=CC=C1` and `CC1=CC=CC=C1C`; perceived, both are `Cc1ccccc1C`,
  as RDKit writes them. Acyclic molecules need no aromaticity. A molecule
  with a stereo center or E/Z endpoint on an aromatic atom is still written
  in its stored Kekule form, now documented on `SmilesWriteMode::Canonical`.
- Molfile query hydrogen counts are no longer read as fixed hydrogen counts.
  The V2000 hydrogen-count field and V3000 `HCOUNT` are CTfile query
  constraints ("at least n hydrogens"), yet a V2000 code of 3 became
  `ImplicitHydrogens::Fixed(2)` without any report. A nonzero value is now
  unsupported content (`MolfileUnsupportedContent::AtomField` with field
  `hydrogen count`, or `AtomOption` with key `HCOUNT`): rejected by default,
  or reported and omitted with `UnsupportedContent::Report`, leaving the
  atom's count to its declared valence or to inference. A declared valence
  is the only Molfile hydrogen declaration. `MolfileAtom::hydrogen_count` is
  removed.
- PDB and mmCIF read deuterium and tritium. The element symbols `D` and `T`
  failed as unknown elements, so neutron structures and heavy water (`DOD`)
  could not be read; they now read as hydrogen with mass number 2 and 3, as
  component dictionary entries already did, through one symbol reader shared
  by the three formats. The mmCIF writer writes these isotopes as `D` and `T`
  instead of rejecting them, and still rejects isotopes it cannot name.
- Structure interpretation no longer drops a stated unit cell silently. PDB
  `CRYST1` records are typed (`PdbDocument::cryst1`, `PdbCryst1`) and give
  the model its periodic cell by the rule mmCIF `_cell` follows: space group
  P 1 or none, and not the 1 Å placeholder. Before, PDB entries never had a
  cell. Cell parameters in either format that are not the model's cell, being
  incomplete, under another space group, or the placeholder, are reported as
  the new `StructureInterpretIssue::CellOmitted` with a
  `CellOmissionReason` and the source location.
- mmCIF glycan links (`_pdbx_entity_branch_link`) bond the residues of every
  copy of a branched entity. Each end of a link resolved to the first matching
  residue of any copy, so when copies shared a molecule only the first got its
  glycosidic bonds, and a copy missing a residue was bonded to another copy's.
  Links now apply within each asym that `_pdbx_branch_scheme` lists, only
  within bonding distance; a longer gap is reported as
  `StructureInterpretIssue::ChainBreak`, as for polymer links.
- `MoleculeEditor` keeps represented stereo valid symmetrically and locally.
  `add_bond` never pruned, so a bond to a stereo center left an invalid
  assertion that failed `finish` with `InvalidStereo`, while every
  `delete_bond` pruned all invalid stereo in the draft, so adding the same
  bond and deleting any unrelated bond published with the center's assertion
  silently gone. `add_bond` and `delete_bond` now prune exactly the elements
  that the edited bond invalidates, among those whose focus, axis endpoints
  or carriers include its endpoints. The new
  `MoleculeEditor::pruned_stereo_elements` and `take_pruned_stereo_elements`
  report every element that structural edits pruned: adding, deleting or
  rewiring bonds, changing bond orders, deleting atoms, and oxohalogen
  normalization.
- Models read from PDB can be written as mmCIF. A PDB author chain holds its
  protein, ligands, and water, but the writer planned one entity per chain,
  so `mmcif::write` failed with `ConflictingAsymEntityClassifications` on any
  such chain (1AKE among them), and polymer residues without a
  `label_seq_id` made bonds between them ambiguous. The writer now gives each
  entity instance of a chain without a source `label_asym_id` its own label
  asym: the chain's polymer keeps its label, each ligand residue and the
  chain's water get labels of their own such as `A_2`, and every one keeps
  the chain's `auth_asym_id`. A source asym is still never split. Polymer
  residues without a `label_seq_id` are numbered from 1 in chain order. Under
  canonical classification, a metal, saccharide, or water residue bonded to a
  polymer is no longer written as a polymer residue; explicit classifications
  hold as given.
- The XTC writer encodes compressed coordinates itself, following GROMACS
  `xdr3dfcoord`, instead of through molly 0.6.1. molly packed integers wider
  than 64 bits with their axes reversed, so frames with wide spans at high
  precision could not be read back, such as a 660 nm span at 2^14 per nm or
  a 400 nm box at 10^4 per nm. Ordinary frames are written byte for byte as
  before.
- The XTC writer no longer writes a corrupt frame when the coordinates on one
  axis span more than the selected precision can encode, which it did
  silently in release builds. It now rejects such a frame with
  `TrajectoryCodecErrorKind::InvalidFrame` before writing anything, as
  GROMACS does, bounding the scaled coordinates and each axis's span by
  `i32::MAX - 2`. At a precision of 10^6 per nm, a frame from -1500 to
  1500 nm is rejected.
- Isomeric and canonical SMILES write explicitly unknown stereo unmarked,
  as RDKit does, instead of refusing the molecule. SMILES has no notation for
  an unknown configuration, so an amidine C=N read with a Molfile "either"
  mark is written `CNC(C)=N`; it reads back as unspecified.
- Canonical SMILES no longer panics on an E/Z bond whose only carrier at one
  endpoint is an explicit hydrogen, such as the N-H of `[H]/N=C/C` read from
  an RDKit-embedded SDF record. Hydrogen normalization collapsed that
  hydrogen, leaving no atom to carry the direction mark; it now stays
  explicit, and the output matches RDKit's (`[H]/N=C/C`).
- Structure writers no longer drop a model's periodic cell silently. mmCIF
  writes a cell periodic along all three axes in the standard
  crystallographic frame as `_cell` in space group P 1, and mmCIF reading now
  takes `_cell` as the model's cell when the space group is P 1 or absent
  (not for PDB's 1 Å placeholder). A cell mmCIF cannot hold (partially
  periodic, rotated, or differing between ensemble members) fails with
  `MmcifWriteError::UnrepresentableCell`. Molfile and SDF cannot hold any
  cell and fail with `UnsupportedRepresentation`. The new
  `MmcifWriteOptions::cell`, `MolfileWriteOptions::cell` and
  `SdfWriteOptions::cell` set to `UnrepresentablePolicy::Omit` write the
  structure without it.
- Writing an ensemble no longer drops its members' weights and energies
  silently. SDF writes them to each member's `kekule.weight` and
  `kekule.energy.kJ_per_mol` data fields (`sdf::WEIGHT_FIELD`,
  `sdf::ENERGY_FIELD`), and the new `sdf::ensemble` assembles interpreted
  records back into a weighted ensemble. mmCIF cannot carry them, so writing
  unequal weights or any energy fails with
  `MmcifWriteError::UnrepresentableEnsembleMember` unless
  `MmcifWriteOptions::ensemble_metadata` is `UnrepresentablePolicy::Omit`.
- V2000 atom fields are read from their CTfile columns, so a three-digit atom
  map that touches the field before it, as Kekule itself writes maps 100 to
  999, is no longer lost on reading.
- A present but nonnumeric V2000 hydrogen-count, valence, charge, atom-map,
  mass-difference or bond stereo field is a syntax error instead of reading
  as absent.
- Nonzero V2000 stereo care, H0 designator, inversion/retention and
  exact-change atom fields and ring/chain topology and reacting-center bond
  fields are reported as unsupported content, like their V3000 options,
  instead of being dropped silently.
- A V2000 atom-block mass difference sets the isotope when no `M  ISO`
  record replaces it (`MolfileAtom::mass_difference`). Like RDKit, it counts
  from the element's most abundant isotope; an element without a naturally
  abundant isotope fails interpretation instead.
- Molfile writing fails with `UnsupportedRepresentation` for an atom map of
  0, which both versions write and read as no atom map, instead of silently
  dropping it.
- Proper and improper torsions whose amplitudes are all zero contribute zero
  energy and gradient at any geometry. They used to be rejected as a
  `DegenerateDihedral` when their atoms were collinear, so `OpenFfPotential`
  could not evaluate linear groups such as acetylene, whose torsions Rosemary
  zeroes. Torsions with a nonzero amplitude still need a defined dihedral.
- A PDB `TER` record ends the residue it follows as well as the chain's
  polymer. A segment that restarted with the same chain, residue number, and
  insertion code used to join that residue, and its repeated atom names failed
  as duplicate alternate-location records; it is now a new residue, reported
  with its occurrence like a number that recurs after other residues.
- `MmcifParseOptions::max_atom_site_rows` stops parsing at the first excess
  `_atom_site` loop row instead of after buffering the whole loop, so memory
  stays bounded by the limit. The error now names the line of that row rather
  than the `loop_` keyword, and takes precedence over syntax errors later in
  the loop.
- mmCIF parsing no longer takes time quadratic in the values on one line.
  Each value searched for the end of its line from the line's start, so a
  loop written on one long line was tokenized in O(length × values).
- The XYZ reader's `max_scratch_bytes` covers its reusable line buffer as well
  as its positions, as the limit's combined-scratch contract states. The
  buffer grows only within the budget the positions leave, so a line longer
  than that fails with `ResourceLimitExceeded` before it is buffered, instead
  of growing the buffer up to the line limits regardless of the scratch limit.
- The XYZ reader's `max_frame_bytes` check counts the bytes each line takes in
  the source, terminator included. A CRLF line was counted one byte short, and
  a final line without a terminator one byte long, so the limit depended on
  line endings.
- Average-linkage dendrograms no longer overflow to infinite merge heights
  when finite distances approach `f64::MAX`, which also changed height cuts.
  Merged distances are computed as a weighted mean that never exceeds the
  larger input.
- `Ensemble::cluster_k_medoids` and the representatives of
  `Ensemble::cluster_hierarchical` choose the member with the least weighted
  distance sum when those sums exceed `f64::MAX`, instead of the lowest index
  among infinite sums. Such sums are compared on distances rescaled by a power
  of two, so ordinary inputs and their tie-breaking are unchanged.
- A 3D Molfile or SDF record no longer gives E/Z to a double bond whose
  endpoint carries two identical terminal ligands, such as a vinyl `=CH2`
  drawn with its hydrogens. The configuration was arbitrary, and canonical
  SMILES depended on atom order. Flat drawings still assert E/Z at every
  drawn double bond. Writing a conformation with such an assertion now fails
  as unsupported instead of emitting it.
- Molfile and SDF wedges, `CFG` parities, and the writer accept a nitrogen
  whose lone pair is a stereo carrier: a nitrogen in a three-membered ring or
  at a bridgehead, as in quinine. Writing a 3D model with such a center,
  including one Kekule had just read, failed with an unencodable stereo
  element, and a drawing wedged at the nitrogen failed to interpret.
- A V3000 atom `CFG` parity that contradicts its atom's 3D coordinates no
  longer overrides them silently. The coordinates decide the configuration, as
  they already did for wedges, and interpretation reports
  `AtomParityContradictsCoordinates`. Flat drawings still take configurations
  from parities.
- The Molfile writer's error for asserted stereo its coordinates contradict
  names the center or bond and the contradiction instead of reporting an
  unencodable stereo element.
- `OverwritePolicy::Forbid` path writers check that the destination's
  filesystem supports hard links, which publication needs to never replace a
  file, when the writer is created. On FAT, exFAT, and some network shares they
  now fail before any frame is written, with an error naming
  `OverwritePolicy::Replace`, instead of after writing every frame.
- `OverwritePolicy::Replace` replaces an existing destination on Windows too,
  where it used to be refused.
- On Unix, finishing a path writer synchronizes the destination directory, so
  a trajectory reported as published survives a crash.
- DCD writing with `DcdTimePolicy::HeaderDelta` accepts frame times that
  match step × `DELTA` to f32 precision, the precision `DELTA` has in the file.
  A 0.002 ps timestep used to be rejected from 2 ps on, because rounding it to
  f32 alone exceeded the old 1e-9 relative tolerance.
- XYZ, DCD, TRR, and XTC writers refuse further frames and finishing after a
  write fails partway through a frame or DCD header update, instead of
  appending to, and reporting as complete, a stream that may hold a partial
  frame. A rejected frame still writes nothing and leaves the writer usable.
- Cell-validation errors from the DCD, TRR, and XTC writers state only the
  cause, no longer the reused reader's "while attempting to read frame"
  context.
- V3000 records longer than the CTfile's 80 columns, such as atoms with
  round-trip coordinates and several properties or large stereo collections,
  continue on following lines ending in `-` instead of running past the limit.
- Molfiles set the chiral flag when they assert absolute stereo: a specified
  tetrahedral or axis center that is ungrouped or in an `ABS` group. An unset
  flag made BIOVIA-convention readers, and RDKit's
  `TranslateChiralFlagToStereoGroups`, read those centers as a racemic `AND`
  group.
- mmCIF values starting with `$`, `[`, or `]`, which CIF 1.1 reserves, are
  quoted.
- Automatic Molfile version choice promotes to V3000 only when V2000 cannot
  represent the content, and skips the V2000 attempt when the atom and bond
  counts or enhanced stereo groups already rule it out. Other failures are
  reported as they occur instead of being retried as V3000.
- Molfile and SDF output states its dimensional code in the header's program
  line, so readers such as RDKit treat a written conformation as 3D instead of
  as a 2D drawing. The code is `2D` exactly when every written z is zero,
  including a coordinate-free molecule, since RDKit takes any nonzero z as 3D
  and then ignores wedges. A flat drawing whose wedge or crossed bond carries
  stereo is therefore written on the z = 0 plane; other coordinates are written
  as given. The program name moves to the CTfile program columns 3-10, leaving
  the user-initials field blank.
- mmCIF residues whose bonds come from the block's own `_chem_comp_bond` table
  attach the hydrogens the table omits, such as N-terminal `H3` or a
  5'-terminal `HO5'`, to their nearest heavy atom, as residues bonded from
  templates do, instead of leaving each as a one-atom molecule. Charges then
  follow the attached hydrogens.
- A ring nitrogen with one hydrogen more or fewer than its template, such as
  the N3-protonated cytosine of i-motif DNA or an imidazolate histidine, has
  its charge changed by one instead of failing interpretation, whatever its
  template bond orders.
- Publishing a topology no longer takes time quadratic in the residues of one
  chain or the atom sites of one residue. A 200,000-residue chain now validates
  in about 0.1 s instead of 14 s, which also dominated mmCIF interpretation of
  large structures.
- V3000 Molfiles may omit the BOND section when their COUNTS line declares no
  bonds, as RDKit writes molecules without bonds. Declared bonds still
  require the section.
- `TopologyEditor` and `ModelEditor` can bond fragments of two occurrences
  that were both edited in the same draft, for example after deleting a bond
  in each. Merging copied the second occurrence's atom and bond annotations
  as if its draft molecule were compact, so rows kept for deleted atoms or
  bonds made `add_bond` fail with a property `LengthMismatch`. Only live rows
  are copied now, each to its new slot.

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
  `TrajectoryFrameView`/`TrajectoryFrameMut`. `subset` replaces
  `Trajectory::slice` for atom selections, and the new `select` copies items in
  any order, including ranges, strides, and repeated indices.
  `Ensemble::from_models` consumes models.
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
