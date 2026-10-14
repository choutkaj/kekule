# Changelog

All notable changes to Kekule are documented in this file.

## [Unreleased]

### Added

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
  format-neutral names and a namespace of their own, `kekule::macromolecular`;
  `kekule::mmcif` still re-exports them. `MmcifInterpretOptions`,
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
  record replaces it (`MolfileAtom::mass_difference`). It counts from the
  most abundant isotope where that agrees with the rounded standard atomic
  weight; for other elements, such as Br, interpretation fails instead of
  guessing.
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
