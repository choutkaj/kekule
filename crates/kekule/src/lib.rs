//! Pure-Rust foundations for molecular graphs, chemical perception, structure
//! I/O, structural bioinformatics, and molecular modelling.
//!
//! # Object model
//!
//! Kekule keeps chemical identity, system organization, and geometry separate:
//!
//! - [`molecule::Molecule`] is one non-empty connected molecular graph. It owns
//!   represented chemistry, derived [`perception::Perception`], and
//!   geometry-independent properties.
//! - [`topology::Topology`] is one coordinate-free system containing one or
//!   more explicit molecule instances plus an optional biological hierarchy.
//! - [`structure::Model`] is one structure of a topology: positions, an
//!   optional periodic cell, and structure-scoped properties.
//! - [`structure::Ensemble`] stores weighted samples of one shared topology,
//!   grouped into states that may carry energies.
//! - [`structure::Trajectory`] stores ordered
//!   [`structure::TrajectoryFrame`]s of one shared topology, each with
//!   positions and optional cell, velocities, forces, time, and step. File
//!   codecs and streaming analysis for trajectories live in the `kekule-traj`
//!   companion crate.
//!
//! A disconnected salt, solvent box, or protein-ligand complex is therefore a
//! [`topology::Topology`] containing several connected molecules, not one
//! disconnected [`molecule::Molecule`]. Instances may explicitly share one
//! reusable molecule definition. Geometry and numeric arrays address atoms by
//! dense topology index in the topology's stored layout order; molecule atom
//! IDs and hierarchy IDs are separate semantic identities.
//!
//! # Modules
//!
//! | Area | Modules |
//! | --- | --- |
//! | Molecules and systems | [`molecule`], [`topology`] (with [`topology::edit`]), [`structure`] (with [`structure::measure`]), [`prelude`] |
//! | Values and diagnostics | [`units`], [`properties`], [`geometry`], [`diagnostics`] |
//! | Formats | [`smiles`], [`molfile`], [`sdf`], [`mmcif`], [`pdb`], [`xyz`], [`macromolecular`], [`components`], [`native`] |
//! | Chemical perception | [`perception`], [`stereo`] ([`stereo::cip`], [`stereo::coordinates`], [`stereo::cleanup`]), [`hydrogens`], [`descriptors`] |
//! | Queries | [`query`], [`smarts`], [`substructure`], [`select`] |
//! | Structure analysis and preparation | [`alignment`], [`clustering`], [`statistics`], [`dssp`], [`solvation`] |
//!
//! # Typical workflow
//!
//! Format APIs deliberately separate parsing from chemical interpretation.
//! Interpretation publishes represented chemistry but does not run perception
//! implicitly. Geometry is supplied only when constructing a model.
//!
//! ```
//! use kekule::{
//!     geometry::Point3,
//!     smiles,
//!     structure::{Model, Positions},
//!     units::{Quantity, ANGSTROM},
//! };
//!
//! let mut molecules = smiles::to_molecules("CCO")?;
//! let mut ethanol = molecules.pop().expect("one connected component");
//! ethanol.perceive()?;
//!
//! let positions = Positions::new(Quantity::new(
//!     vec![
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.8, 0.0, 0.0),
//!     ],
//!     ANGSTROM,
//! ))?;
//! let model = Model::from_molecule(ethanol, &positions)?;
//!
//! assert_eq!(model.topology().instance_count(), 1);
//! assert_eq!(model.atom_count(), 3);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Explicit operations
//!
//! Kekule avoids hidden chemistry changes. Parsing, interpretation,
//! perception, hydrogen transforms, coordinate stereo materialization, and
//! force-field preparation are separate operations. Structural mutation uses
//! transactional builders or editors so published values retain their
//! invariants.
//!
//! | Task | Interface | Publication |
//! | --- | --- | --- |
//! | Construct or edit one connected molecule | [`molecule::MoleculeEditor`] | `finish()` |
//! | Assemble complete molecules and reusable instances | [`topology::TopologyBuilder`] | `build()` |
//! | Assemble a system with explicit geometry | [`structure::ModelBuilder`] | `build()` |
//! | Edit system atoms, bonds, components, and hierarchy | [`topology::edit::TopologyEditor`] | `finish()` |
//! | Edit a system together with its coordinates | [`structure::ModelEditor`] | `finish()` |
//! | Change geometry or structure annotations | [`structure::Model`] setters | Immediate checked update |
//!
//! Builders and editors offer non-consuming `validate()` and recoverable
//! `try_build()` / `try_finish()`. Use `to_editor()` for a detached draft or
//! `into_editor()` to move an owner into one. System editors resolve source IDs
//! to stable draft-only editing handles and publish completed values through
//! `finish()`. A bond deletion can split a system molecule;
//! a bond addition can join two occurrences. New model atoms require coordinates.
//!
//! Coordinate-dependent algorithms consume [`structure::ModelView`]. A model,
//! ensemble sample, trajectory frame, or reusable trajectory buffer can
//! therefore share kernels without copying coordinates. APIs that address atoms
//! by index accept any topology snapshot sharing their layout
//! ([`topology::Topology::shares_layout`]), so perception never detaches
//! selections, buffers, or prepared potentials. Use
//! [`topology::Topology::same_layout`] only when complete static layout equality
//! between independent publications is intended.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

macro_rules! fixed_u32_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            /// Wraps a raw value. Nothing is validated here; a lookup with a
            /// value its owner does not hold returns `None`.
            pub const fn new(raw: u32) -> Self {
                Self(raw)
            }

            /// Returns the raw value. It is this identifier's stable
            /// encoding for persistence: `new(id.raw())` restores `id`, and
            /// the encoding does not change between Kekule versions.
            pub const fn raw(self) -> u32 {
                self.0
            }

            /// Returns the raw value as a `usize` position.
            pub const fn index(self) -> usize {
                self.0 as usize
            }
        }
    };
    ($(#[$meta:meta])* $name:ident, $display:literal) => {
        fixed_u32_id!($(#[$meta])* $name);

        impl std::fmt::Display for $name {
            fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(formatter, concat!($display, " {}"), self.0)
            }
        }
    };
}

mod algorithms;
pub mod alignment;
mod chemistry;
pub mod clustering;
pub mod descriptors;
pub mod diagnostics;
mod draft;
pub mod dssp;
pub mod geometry;
mod io;
pub mod molecule;
pub mod properties;
pub mod query;
mod selection_language;
pub mod solvation;
pub mod statistics;
pub mod structure;
pub mod topology;
pub mod units;

pub use draft::DraftError;

/// Syntax-independent substructure matching algorithms.
///
/// Matching consumes [`query::QueryGraph`] and current target perception state;
/// it never invokes parsing, interpretation, or perception implicitly.
/// [`PreparedTarget`](substructure::PreparedTarget) borrows one molecule and
/// answers repeated queries against it;
/// [`PreparedTopologyTarget`](substructure::PreparedTopologyTarget) does the
/// same for every molecule instance of a topology.
///
/// ```
/// use kekule::substructure::{PreparedTarget, SubstructureMatchOptions};
/// use kekule::{smarts, smiles};
///
/// let mut toluene = smiles::to_molecules("Cc1ccccc1")?.remove(0);
/// toluene.perceive()?;
/// let aromatic_carbon = smarts::parse_str("c", Default::default())?;
/// let target = PreparedTarget::new(&toluene);
/// let matches = target.all(&aromatic_carbon, SubstructureMatchOptions::default())?;
/// assert_eq!(matches.len(), 6);
/// let nitrogen = smarts::parse_str("N", Default::default())?;
/// assert!(target.first(&nitrogen, Default::default())?.is_none());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod substructure {
    pub use crate::algorithms::{
        MatchCompletion, PreparedTarget, PreparedTopologyTarget, QueryMatch, QueryPerception,
        SubstructureMatchError, SubstructureMatchOptions, SubstructureMatchWork, TaggedMatchError,
        TaggedQuery, TopologyQueryMatch, TupleMode, MAX_SUBSTRUCTURE_QUERY_ATOMS,
    };
}

/// SMILES parsing, interpretation, and molecule/topology writing.
///
/// [`smiles::parse_str`] preserves source syntax in a
/// [`smiles::SmilesDocument`], while
/// [`SmilesDocument::interpret`](smiles::SmilesDocument::interpret) publishes
/// canonical connected molecules. [`smiles::to_molecules`] is the concise
/// parse-and-interpret path. Dot-separated components remain separate molecules
/// and perception is never run implicitly. [`smiles::write`] writes a molecule
/// or a whole topology as ordinary, isomeric, or canonical SMILES.
///
/// ```
/// use kekule::smiles::{self, SmilesWriteOptions};
///
/// let molecules = smiles::to_molecules("CCO.[Na+].[Cl-]")?;
/// assert_eq!(molecules.len(), 3);
/// assert_eq!(molecules[0].atom_count(), 3);
/// assert_eq!(smiles::write(&molecules[0], SmilesWriteOptions::ordinary())?, "CCO");
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod smiles {
    use std::sync::Arc;

    use crate::diagnostics::{DiagnosticKind, Format, FormatError};
    pub use crate::io::{
        SmilesAtomMapping, SmilesBondMapping, SmilesComponentInterpretation, SmilesCxField,
        SmilesCxFieldKind, SmilesCxStereoGroupKind, SmilesDocument, SmilesDocumentToken,
        SmilesDocumentTokenKind, SmilesInterpretOptions, SmilesInterpretation,
        SmilesInterpretationReport, SmilesParseOptions, SmilesRecords, MAX_CANONICAL_GRAPH_SLOTS,
        MAX_CANONICAL_INPUT_COMPLEXITY, MAX_CANONICAL_PENDING_LABELS,
        MAX_CANONICAL_REFINEMENT_WORK, MAX_CANONICAL_SEARCH_STATES, MAX_DIRECTIONAL_COMBINATIONS,
        MAX_DIRECTIONAL_VISITS,
    };
    use crate::molecule::Molecule;
    use crate::topology::Topology;

    /// SMILES output policy of [`write()`].
    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    #[non_exhaustive]
    pub enum SmilesWriteMode {
        /// Ordinary non-canonical SMILES without stereo. Declared and inferred
        /// hydrogen counts may be combined in a bracket atom; the original
        /// inference policy is not serialized.
        #[default]
        Ordinary,
        /// Non-canonical SMILES preserving represented stereo.
        ///
        /// SMILES has no notation for an explicitly unknown configuration, so
        /// this and canonical mode write such centers, bonds and axes
        /// unmarked, as RDKit does; they read back as unspecified.
        ///
        /// Tetrahedral absolute, AND, OR and relative groups are emitted as
        /// CXSMILES (`a`, `&`, `o`, `r`). Unsupported member geometries,
        /// quantitative racemic groups and multiple independent relative groups
        /// return an error rather than losing relationships. The global `r`
        /// flag shields other specified centers with explicit absolute
        /// membership. Directional encoding is bounded by 4,096 carrier
        /// combinations and 50,000,000 graph visits; exceeding either returns a
        /// resource-limit error.
        Isomeric,
        /// Deterministic canonical isomeric SMILES.
        ///
        /// Atom priorities, stereo refinement, branch traversal, ring closures
        /// and AND/OR representative selection follow RDKit-style conventions;
        /// byte-for-byte RDKit compatibility is not guaranteed, and canonical
        /// spelling can change when these ordering rules are improved. Output is
        /// invariant under atom numbering and preserves supported stereo,
        /// isotopes, formal charges, and atom maps. Neutral unmapped
        /// nonisotopic terminal hydrogen vertices may collapse into hydrogen
        /// counts. Supplied local stereo assertions are retained even when CIP
        /// perception finds no stereogenic unit. Group numbering and member
        /// order are canonical, and inverting every member of a non-absolute
        /// group produces the same output.
        ///
        /// Canonical output of a molecule with a ring requires installed
        /// aromaticity perception, such as [`Molecule::perceive`] installs, and
        /// fails with
        /// [`DiagnosticKind::MissingPerception`]
        /// without it.
        /// Perceived
        /// aromatic atoms and bonds are written aromatic, so every Kekule form
        /// of a molecule gives one string under the installed aromaticity
        /// model. The exception is a molecule with a stereo center or E/Z
        /// endpoint on an aromatic atom: SMILES cannot mark that configuration
        /// on aromatic syntax, so the whole molecule is written in its stored
        /// Kekule form.
        ///
        /// Complete canonical labeling is bounded by 100,000 search states,
        /// 50,000,000 atom/edge/twin visits, and 2,000,000 pending atom labels.
        /// Serialization additionally bounds the input to 2,000,000 combined
        /// atom/bond slots and a complexity score `2*n*(n+2*m)` of at most
        /// 50,000,000 for `n` live atoms and `m` live bonds. Exceeding a bound
        /// returns
        /// [`DiagnosticKind::Limit`]
        /// without a partial result.
        Canonical,
    }

    /// Options of [`write()`] and [`write_to()`]; the default writes
    /// [`SmilesWriteMode::Ordinary`] SMILES.
    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    pub struct SmilesWriteOptions {
        mode: SmilesWriteMode,
    }

    impl SmilesWriteOptions {
        /// The output policy.
        pub const fn mode(&self) -> SmilesWriteMode {
            self.mode
        }

        /// Sets [`Self::mode`].
        #[must_use]
        pub const fn with_mode(mut self, mode: SmilesWriteMode) -> Self {
            self.mode = mode;
            self
        }
    }

    impl SmilesWriteOptions {
        /// Options for [`SmilesWriteMode::Ordinary`] output.
        pub const fn ordinary() -> Self {
            Self {
                mode: SmilesWriteMode::Ordinary,
            }
        }

        /// Options for [`SmilesWriteMode::Isomeric`] output.
        pub const fn isomeric() -> Self {
            Self {
                mode: SmilesWriteMode::Isomeric,
            }
        }

        /// Options for [`SmilesWriteMode::Canonical`] output.
        pub const fn canonical() -> Self {
            Self {
                mode: SmilesWriteMode::Canonical,
            }
        }
    }

    /// Parses one record: base SMILES, an optional CX extension with its fields
    /// typed, and an optional name. `options` bounds the input and the parsed
    /// graph.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind [`DiagnosticKind::Limit`] when the
    /// input exceeds `options`, and [`DiagnosticKind::Syntax`],
    /// [`DiagnosticKind::InvalidValue`] or [`DiagnosticKind::Unsupported`]
    /// when the record is malformed or uses syntax Kekule does not read.
    pub fn parse_str(
        input: &str,
        options: SmilesParseOptions,
    ) -> Result<SmilesDocument, FormatError> {
        crate::io::parse_smiles_document(input, options)
    }

    /// Reads and parses one record from `reader`, as [`parse_str`] parses the
    /// same text. Use [`records`] for line-delimited multi-record input.
    ///
    /// # Errors
    ///
    /// As [`parse_str`]. The configured byte limit is enforced while reading;
    /// a read failure or non-UTF-8 input fails with [`DiagnosticKind::Io`].
    pub fn parse_reader(
        reader: impl std::io::Read,
        options: SmilesParseOptions,
    ) -> Result<SmilesDocument, FormatError> {
        crate::io::parse_smiles_reader(reader, options)
    }

    /// Streams the records of a multi-record SMILES input such as a `.smi`
    /// file, one record per line, each bounded by `options`.
    ///
    /// Blank lines are skipped. Each item is one record or the failure of one
    /// line, located at its input line; reading continues after a malformed
    /// line. Memory is bounded by one line.
    pub fn records<R: std::io::Read>(reader: R, options: SmilesParseOptions) -> SmilesRecords<R> {
        SmilesRecords::new(reader, options)
    }

    /// Parses and interprets one SMILES record into source-ordered connected components.
    ///
    /// This is the concise form of [`parse_str`] followed by
    /// [`SmilesDocument::interpret`], both with default options.
    /// Dot-delimited components remain separate, and no perception is run implicitly.
    ///
    /// # Errors
    ///
    /// The errors of [`parse_str`], and of interpretation: chemistry that
    /// cannot be published fails with [`DiagnosticKind::Chemistry`].
    pub fn to_molecules(input: &str) -> Result<Vec<Molecule>, FormatError> {
        let document = parse_str(input, SmilesParseOptions::default())?;
        Ok(document
            .interpret(SmilesInterpretOptions::default())?
            .into_molecules())
    }

    /// Parses and interprets one SMILES record as a shared coordinate-free
    /// topology.
    ///
    /// Every connected component becomes one explicit molecule
    /// occurrence in source order. No hierarchy or perception is fabricated.
    ///
    /// # Errors
    ///
    /// The errors of [`to_molecules`]. A topology that cannot be assembled
    /// fails with [`DiagnosticKind::Chemistry`] and the
    /// [`TopologyBuildError`](crate::topology::TopologyBuildError) as its
    /// source.
    pub fn to_topology(input: &str) -> Result<Arc<Topology>, FormatError> {
        let document = parse_str(input, SmilesParseOptions::default())?;
        let topology = document
            .interpret(SmilesInterpretOptions::default())?
            .into_topology()
            .map_err(|error| {
                FormatError::interpret(
                    Format::Smiles,
                    DiagnosticKind::Chemistry,
                    "cannot assemble topology",
                )
                .with_source(error)
            })?;
        Ok(Arc::new(topology))
    }

    /// What one SMILES record describes: one connected molecule, or every
    /// molecule occurrence of a topology joined with `.`, with an optional
    /// record name written after the SMILES as `.smi` files carry it.
    ///
    /// Convert from `&Molecule`, `&Topology`, or `&Arc<Topology>`, or build
    /// one with [`Self::molecule`] or [`Self::topology`] to add a name.
    ///
    /// ```
    /// use kekule::smiles::{self, SmilesSource};
    /// let ethanol = smiles::to_molecules("CCO")?.remove(0);
    /// let line = smiles::write(SmilesSource::molecule(&ethanol).with_name("ethanol"), Default::default())?;
    /// assert_eq!(line, "CCO ethanol");
    /// let document = smiles::parse_str(&line, Default::default())?;
    /// assert_eq!(document.name(), Some("ethanol"));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    #[derive(Debug, Clone, Copy)]
    pub struct SmilesSource<'a> {
        content: SmilesContent<'a>,
        name: Option<&'a str>,
    }

    #[derive(Debug, Clone, Copy)]
    enum SmilesContent<'a> {
        Molecule(&'a Molecule),
        /// Every molecule occurrence of a topology.
        Topology(&'a Topology),
    }

    impl<'a> SmilesSource<'a> {
        /// One connected molecule.
        pub const fn molecule(molecule: &'a Molecule) -> Self {
            Self {
                content: SmilesContent::Molecule(molecule),
                name: None,
            }
        }

        /// Every molecule occurrence of a topology, joined with `.`.
        pub const fn topology(topology: &'a Topology) -> Self {
            Self {
                content: SmilesContent::Topology(topology),
                name: None,
            }
        }

        /// Writes `name` after the SMILES, separated by one space. Writing
        /// fails for an empty name, one with surrounding whitespace, or one
        /// spanning lines, which would not read back as the same name.
        #[must_use]
        pub const fn with_name(mut self, name: &'a str) -> Self {
            self.name = Some(name);
            self
        }
    }

    impl<'a> From<&'a Molecule> for SmilesSource<'a> {
        fn from(molecule: &'a Molecule) -> Self {
            Self::molecule(molecule)
        }
    }

    impl<'a> From<&'a Topology> for SmilesSource<'a> {
        fn from(topology: &'a Topology) -> Self {
            Self::topology(topology)
        }
    }

    impl<'a> From<&'a Arc<Topology>> for SmilesSource<'a> {
        fn from(topology: &'a Arc<Topology>) -> Self {
            Self::topology(topology)
        }
    }

    /// Writes one SMILES record.
    ///
    /// [`SmilesWriteMode`] selects ordinary, isomeric, or canonical output.
    /// Canonical output of a molecule with a ring requires installed
    /// aromaticity perception.
    /// Every mode preserves isotope labels and total hydrogen counts. An atom
    /// that permits inferred hydrogens and requires bracket syntax must have
    /// hydrogen perception installed; otherwise writing returns an error. Call
    /// [`Molecule::perceive`] explicitly before exporting such atoms. Radical
    /// electrons are encoded by bracket valence, with redundant CX annotations
    /// for one, two or three electrons (`^1`, `^2`, `^5`); explicit spin
    /// multiplicities remain unsupported.
    ///
    /// A topology emits each occurrence of a reused definition once, joined
    /// with `.`; canonical mode sorts components, other modes keep instance
    /// order. Enhanced groups share one CX extension with record-global atom
    /// indices.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind [`DiagnosticKind::MissingPerception`]
    /// when required perception is not installed,
    /// [`DiagnosticKind::Unsupported`] when the mode cannot represent the
    /// chemistry without loss, and [`DiagnosticKind::Limit`] when a
    /// canonicalization or stereo-encoding bound of [`SmilesWriteMode`] is
    /// exceeded.
    ///
    /// ```
    /// use kekule::smiles::{self, SmilesWriteOptions};
    /// let canonical = |source: &str| -> Result<String, Box<dyn std::error::Error>> {
    ///     let mut molecule = smiles::to_molecules(source)?.remove(0);
    ///     molecule.perceive()?;
    ///     Ok(smiles::write(&molecule, SmilesWriteOptions::canonical())?)
    /// };
    /// assert_eq!(canonical("OCC")?, "CCO");
    /// assert_eq!(canonical("CC1=C(C)C=CC=C1")?, canonical("CC1=CC=CC=C1C")?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write<'a>(
        source: impl Into<SmilesSource<'a>>,
        options: SmilesWriteOptions,
    ) -> Result<String, FormatError> {
        let source = source.into();
        if let Some(name) = source.name {
            if name.is_empty() || name.trim() != name || name.contains(['\r', '\n']) {
                return Err(FormatError::write(
                    Format::Smiles,
                    DiagnosticKind::Unsupported,
                    "SMILES record names must be nonempty, trimmed, and on one line",
                )
                .with_detail(crate::diagnostics::FormatErrorDetail::Title));
            }
        }
        let mut record = match source.content {
            SmilesContent::Molecule(molecule) => match options.mode {
                SmilesWriteMode::Ordinary => crate::io::write_smiles(molecule),
                SmilesWriteMode::Isomeric => crate::io::write_isomeric_smiles(molecule),
                SmilesWriteMode::Canonical => crate::io::write_canonical_smiles(molecule),
            },
            SmilesContent::Topology(topology) => {
                crate::io::smiles::write_topology(topology, options)
            }
        }?;
        if let Some(name) = source.name {
            record.push(' ');
            record.push_str(name);
        }
        Ok(record)
    }

    /// Writes one SMILES record as a line terminated by `\n`, so repeated
    /// calls produce a line-delimited SMILES file. See [`write()`].
    ///
    /// # Errors
    ///
    /// The errors of [`write()`], before anything reaches `writer`; a failing
    /// `writer` gives [`DiagnosticKind::Io`].
    pub fn write_to<'a>(
        writer: &mut impl std::io::Write,
        source: impl Into<SmilesSource<'a>>,
        options: SmilesWriteOptions,
    ) -> Result<(), FormatError> {
        let mut line = write(source, options)?;
        line.push('\n');
        writer
            .write_all(line.as_bytes())
            .map_err(|error| FormatError::sink(Format::Smiles, error))
    }
}

/// The atom selection language: text such as `protein and name CA` or
/// `within 5 of resname LIG and not water`, evaluated to an
/// [`AtomSelection`](crate::topology::AtomSelection).
///
/// The keyword grammar follows VMD, PyMOL and MDAnalysis where they agree:
/// `and`, `or`, `not`, parentheses, value lists (`resname ALA GLY`), ranges
/// (`resid 10-50`), comparisons (`bfactor > 50`), quoted values, `within`
/// distances in ångströms unless a unit follows, expansion (`byres`,
/// `bychain`, `bymolecule`, `bonded to`), and macros such as `protein`,
/// `water`, `ligand`, `backbone`, and `sidechain`. Chain, residue, and atom
/// identifiers are author identifiers, falling back to label identifiers;
/// `label_*` terms match label identifiers only. The full grammar is in
/// `docs/selection-language.md`.
///
/// [`select::parse`] returns a [`select::SelectionExpr`] whose `Display` is
/// its canonical text. [`select::SelectionExpr::evaluate`] evaluates
/// topology terms against a topology;
/// [`select::SelectionExpr::evaluate_structure`] also evaluates spatial and
/// site terms against a model, ensemble sample, or trajectory frame. Chemistry
/// terms read installed perception and never perceive. Applications add
/// `#...`, `$name`, and keyword terms through a [`select::SelectionResolver`],
/// which may answer spatial references with another object's points and add
/// such points to any `within` reference, as in `within 5 of ligand`.
///
/// ```
/// use kekule::select::{self, SelectEvaluateOptions, SelectParseOptions};
/// use kekule::smiles;
///
/// let topology = smiles::to_topology("c1ccccc1O.[OH2]")?.perceived()?;
/// let expression = select::parse("aromatic or (element O and not water)", SelectParseOptions::default())?;
/// assert_eq!(expression.to_string(), "aromatic or element O and not water");
/// let atoms = expression.evaluate(topology, SelectEvaluateOptions::default())?;
/// assert_eq!(atoms.len(), 7);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod select {
    use crate::diagnostics::FormatError;
    pub use crate::selection_language::{
        HostRequest, HostRole, HostTerm, HostTermError, HostTermKind, HostValue,
        SelectEvaluateError, SelectEvaluateOptions, SelectParseOptions, SelectionExpr,
        SelectionResolver, DEFAULT_SELECT_MAX_DEPTH, DEFAULT_SELECT_MAX_INPUT_BYTES,
        DEFAULT_SELECT_MAX_SMARTS_MATCHES,
    };

    /// Parses selection text. SMARTS patterns inside it are checked now, so
    /// evaluation fails only on the evaluated structure, never on the text.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of format
    /// [`Format::Selection`](crate::diagnostics::Format::Selection) with the
    /// byte [`span`](FormatError::span) of the offending text:
    /// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
    /// for text outside the grammar,
    /// [`DiagnosticKind::InvalidValue`](crate::diagnostics::DiagnosticKind::InvalidValue)
    /// for a value a keyword cannot take, and
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when `options` is exceeded.
    pub fn parse(text: &str, options: SelectParseOptions) -> Result<SelectionExpr, FormatError> {
        SelectionExpr::parse(text, &options)
    }
}

/// SMARTS parsing into syntax-independent substructure queries.
///
/// [`smarts::parse_str`] returns a [`query::QueryGraph`], the canonical query
/// representation that [`substructure`] matching consumes. Parsing a query and
/// matching it are separate operations; matching does not perceive the target
/// implicitly. The dialect, including its explicit exclusions, is documented on
/// [`query`].
///
/// ```
/// use kekule::{smarts, smiles, substructure::PreparedTarget};
///
/// let mut succinic_acid = smiles::to_molecules("OC(=O)CCC(=O)O")?.remove(0);
/// succinic_acid.perceive()?;
/// let carboxylic_acid = smarts::parse_str("C(=O)[OX2H1]", Default::default())?;
/// let target = PreparedTarget::new(&succinic_acid);
/// let matches = target.all(&carboxylic_acid, Default::default())?;
/// assert_eq!(matches.len(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod smarts {
    pub use crate::query::smarts::{parse_str, SmartsParseOptions};
}

/// Molfile parsing, interpretation, and V2000/V3000 writing.
///
/// A parsed [`molfile::MolfileDocument`] is format state. Its interpretation can
/// project to connected molecules or to a geometry-bearing
/// [`crate::structure::Model`]. Writers reject chemistry that the selected
/// Molfile version cannot represent.
///
/// ```
/// use kekule::{molfile, smiles};
///
/// let ethanol = smiles::to_molecules("CCO")?.remove(0);
/// let text = molfile::write(&ethanol, Default::default())?;
/// let document = molfile::parse_str(&text, Default::default())?;
/// let interpretation = document.interpret(Default::default())?;
/// assert_eq!(interpretation.model().atom_count(), 3);
/// assert_eq!(interpretation.into_molecules().len(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod molfile {
    use crate::diagnostics::FormatError;
    pub use crate::io::{
        MolfileAtom, MolfileAtomMapping, MolfileAtomParity, MolfileBond, MolfileBondMapping,
        MolfileBondOrder, MolfileBondStereo, MolfileDocument, MolfileHeader,
        MolfileInterpretOptions, MolfileInterpretation, MolfileInterpretationReport,
        MolfileInterpretationWarning, MolfileParseOptions, MolfileRadical, MolfileSource,
        MolfileStereoGroup, MolfileStereoGroupKind, MolfileUnsupported, MolfileUnsupportedContent,
        MolfileVersion, MolfileWriteOptions, MolfileWriteVersion,
    };

    /// Parses one V2000 or V3000 Molfile into its typed connection table
    /// without assigning canonical meaning. `options` bounds the input and the
    /// parsed table.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when the input exceeds `options`, and `Syntax`, `InvalidValue`,
    /// `InvalidReference`, `Inconsistent` or `Unsupported` when the connection
    /// table is malformed or uses content Kekule does not read.
    pub fn parse_str(
        input: &str,
        options: MolfileParseOptions,
    ) -> Result<MolfileDocument, FormatError> {
        crate::io::parse_molfile_document(input, options)
    }

    /// Reads and parses one Molfile from `reader`.
    ///
    /// # Errors
    ///
    /// As [`parse_str`]. The configured byte limit is enforced while reading;
    /// a read failure or non-UTF-8 input fails with
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn parse_reader(
        reader: impl std::io::Read,
        options: MolfileParseOptions,
    ) -> Result<MolfileDocument, FormatError> {
        crate::io::parse_molfile_reader(reader, options)
    }

    /// Writes one CTAB from a molecule (zero coordinates) or a model.
    ///
    /// [`MolfileWriteVersion::Auto`] uses V2000 when it can represent the
    /// chemistry and V3000 otherwise. V2000 rounds coordinates to four decimal
    /// places in angstroms; V3000 preserves round-trip decimal text. Stereo is
    /// projected against the coordinates emitted by the chosen version, so
    /// specified stereo requires a geometry-bearing model.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind
    /// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
    /// when the selected version cannot represent the chemistry, or
    /// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
    /// when the source is incomplete.
    pub fn write<'a>(
        source: impl Into<MolfileSource<'a>>,
        options: MolfileWriteOptions,
    ) -> Result<String, FormatError> {
        crate::io::write_molfile(source, options)
    }

    /// Writes the CTAB [`write()`] would return. It renders the whole CTAB
    /// first, since automatic version selection learns whether V2000 can hold
    /// it only by rendering it, so nothing reaches `writer` on failure.
    ///
    /// # Errors
    ///
    /// The errors of [`write()`]; a failing `writer` gives
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn write_to<'a>(
        writer: &mut impl std::io::Write,
        source: impl Into<MolfileSource<'a>>,
        options: MolfileWriteOptions,
    ) -> Result<(), FormatError> {
        crate::io::write_molfile_to(writer, source, options)
    }
}

/// Record-oriented SDF parsing, interpretation, and writing.
///
/// An SDF input is a sequence of independent records. [`sdf::records`] streams
/// them one at a time, so a malformed record fails alone; [`sdf::parse_str`]
/// collects every record into a [`sdf::SdfDocument`]. Interpret each record
/// explicitly with [`sdf::SdfRecord::interpret`]; sibling records are not
/// merged into one model or reinterpreted as an ensemble.
///
/// ```
/// use kekule::geometry::Point3;
/// use kekule::structure::{Model, Positions};
/// use kekule::units::{Quantity, ANGSTROM};
/// use kekule::{sdf, smiles};
///
/// let water = smiles::to_molecules("O")?.remove(0);
/// let positions = Positions::new(Quantity::new(vec![Point3::new(0.0, 0.0, 0.0)], ANGSTROM))?;
/// let model = Model::from_molecule(water, &positions)?;
///
/// let text = sdf::write([&model, &model], Default::default())?;
/// let document = sdf::parse_str(&text, Default::default())?;
/// assert_eq!(document.records().len(), 2);
/// let first = document.records()[0].interpret(Default::default())?;
/// assert_eq!(first.model().atom_count(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod sdf {
    use crate::diagnostics::FormatError;
    pub use crate::io::{
        SdfDataField, SdfDocument, SdfParseOptions, SdfRecord, SdfRecordInterpretation,
        SdfRecordInterpretationReport, SdfRecordSource, SdfRecords, SdfWriteOptions,
        SDF_ENERGY_FIELD as ENERGY_FIELD, SDF_REPRESENTATIVE_FIELD as REPRESENTATIVE_FIELD,
        SDF_STATE_FIELD as STATE_FIELD, SDF_TEMPERATURE_FIELD as TEMPERATURE_FIELD,
        SDF_WEIGHT_FIELD as WEIGHT_FIELD,
    };

    /// Assembles interpreted records, in order, into one ensemble whose
    /// samples share the first record's topology, record `i` becoming sample
    /// `i`. Sibling records never form an ensemble implicitly.
    ///
    /// Each sample takes its relative weight from [`WEIGHT_FIELD`] (1.0 when
    /// absent). Records without [`STATE_FIELD`], such as a conformer set
    /// from another program, are one state each, with the energy of
    /// [`ENERGY_FIELD`] in kJ/mol (none when absent). Records with it, as
    /// [`write()`] emits for ensemble samples, are grouped by it: states are
    /// numbered in ascending order of its values, each represented by the
    /// one record marked with [`REPRESENTATIVE_FIELD`], which alone may carry
    /// the state's energy. [`TEMPERATURE_FIELD`], in kelvin, is the
    /// ensemble's temperature, the same on every record or on none.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind
    /// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
    /// when the records do not form an ensemble (none, or differing
    /// topologies), a field is repeated, some records name a state and
    /// others do not, a state has no or two representatives, an energy is
    /// on a record that does not represent its state, or temperatures
    /// differ, and
    /// [`DiagnosticKind::InvalidValue`](crate::diagnostics::DiagnosticKind::InvalidValue)
    /// when a field is not a valid number, state index, weight, energy,
    /// temperature, or `true`.
    pub fn ensemble<'a>(
        records: impl IntoIterator<Item = &'a SdfRecordInterpretation>,
    ) -> Result<crate::structure::Ensemble, FormatError> {
        crate::io::sdf_ensemble(records)
    }

    /// Streams the records of an SDF input from `reader`, each bounded by
    /// `options`.
    ///
    /// Each item is one record or the failure of one record; reading resumes at
    /// the next `$$$$` delimiter after a malformed, oversized, or non-UTF-8
    /// record. Memory is bounded by one record.
    pub fn records<R: std::io::Read>(reader: R, options: SdfParseOptions) -> SdfRecords<R> {
        SdfRecords::new(reader, options)
    }

    /// Parses every record of an SDF input, failing on the first malformed
    /// record. Use [`records`] to skip past malformed records instead.
    ///
    /// # Errors
    ///
    /// Returns the [`FormatError`] of the first record that fails: kind
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when the input or a record exceeds `options`, and the errors of
    /// [`molfile::parse_str`](crate::molfile::parse_str) for a malformed
    /// connection table or data field.
    pub fn parse_str(input: &str, options: SdfParseOptions) -> Result<SdfDocument, FormatError> {
        crate::io::parse_sdf_document(input, options)
    }

    /// Reads and parses every record of an SDF input from `reader`.
    ///
    /// # Errors
    ///
    /// As [`parse_str`]. The configured byte limit is enforced while reading;
    /// a read failure or non-UTF-8 input fails with
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn parse_reader(
        reader: impl std::io::Read,
        options: SdfParseOptions,
    ) -> Result<SdfDocument, FormatError> {
        crate::io::parse_sdf_reader(reader, options)
    }

    /// Writes independent records in input order.
    ///
    /// Records convert from models, ensemble samples and trajectory frames
    /// (`sdf::write(&ensemble, options)` writes one record per sample, which
    /// [`ensemble`] reads back with its weight, state, state energy, and
    /// temperature), or
    /// interpreted records, which keep their titles and data fields.
    ///
    /// # Errors
    ///
    /// The errors of [`molfile::write`](crate::molfile::write) for any record,
    /// and
    /// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
    /// or
    /// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
    /// for a title or data field the format cannot hold.
    pub fn write<'a, R: Into<SdfRecordSource<'a>>>(
        records: impl IntoIterator<Item = R>,
        options: SdfWriteOptions,
    ) -> Result<String, FormatError> {
        crate::io::write_sdf(records, options)
    }

    /// Streams the records [`write()`] would return, each as soon as it is
    /// rendered; wrap an unbuffered sink such as a `File` in a
    /// [`BufWriter`](std::io::BufWriter).
    ///
    /// # Errors
    ///
    /// The errors of [`write()`]; records before the failing one have already
    /// reached `writer`. A failing `writer` gives
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn write_to<'a, R: Into<SdfRecordSource<'a>>>(
        writer: &mut impl std::io::Write,
        records: impl IntoIterator<Item = R>,
        options: SdfWriteOptions,
    ) -> Result<(), FormatError> {
        crate::io::write_sdf_to(writer, records, options)
    }
}

/// Chemical component templates: residue and ligand definitions in the model of
/// the wwPDB Chemical Component Dictionary.
///
/// Structural formats name each residue's component but rarely carry its bonds.
/// [`ComponentLibrary::bundled`](components::ComponentLibrary::bundled) holds
/// the standard amino acids and nucleotides, water, common ions, cofactors, and
/// crystallization additives; read further definitions, up to the complete
/// dictionary, with
/// [`ComponentLibrary::from_mmcif`](components::ComponentLibrary::from_mmcif).
/// Structures that molecular-mechanics tools write use force-field residue and
/// atom names; [`ComponentLibrary::amber`](components::ComponentLibrary::amber)
/// and [`ComponentLibrary::charmm`](components::ComponentLibrary::charmm)
/// define them.
///
/// ```
/// use kekule::components::ComponentLibrary;
///
/// let alanine = ComponentLibrary::bundled().get("ALA").expect("bundled");
/// assert_eq!(alanine.name(), "ALANINE");
/// assert_eq!(alanine.atoms().len(), 13);
/// assert!(alanine.atom("CA").is_some());
/// assert!(ComponentLibrary::bundled().get("XYZ").is_none());
/// ```
pub mod components {
    pub use crate::io::{
        Component, ComponentAtom, ComponentBond, ComponentBondOrder, ComponentLibrary,
        ComponentLinking,
    };
}

/// Interpretation types shared by the macromolecular structure formats.
///
/// Macromolecular formats interpret atom sites, alternate locations, coordinate
/// models, and connectivity through one pipeline. Its options, reports, issues,
/// and source provenance live here, and only here; the format modules [`mmcif`]
/// and [`pdb`] take and return them.
pub mod macromolecular {
    pub use crate::io::{
        AltLocDecision, AltLocPolicy, AltLocPreference, AltLocResidue, AltLocSelection,
        AltLocSelectionReason, AtomProvenance, BondEvidence, BondProvenance, CellOmissionReason,
        ConflictingSelector, ConnectionResolutionReason, ConnectivityOptions, EntityKind,
        EntityMetadataPolicy, InstanceProvenance, ModelSelection, SourceResidueId,
        SourceResiduePosition, StructureEnsembleInterpretOptions, StructureEnsembleInterpretation,
        StructureInterpretIssue, StructureInterpretOptions, StructureInterpretation,
        StructureInterpretationReport,
    };
}

/// Structural mmCIF parsing, interpretation, and writing.
///
/// Data blocks are independent interpretation scopes. One selected coordinate
/// model naturally produces a [`crate::structure::Model`], while compatible
/// coordinate models from one block may produce an
/// [`crate::structure::Ensemble`]. Ordinary writing derives mmCIF entity kinds
/// from canonical topology classification; source reports and explicit
/// classifications remain available for faithful round trips and expert
/// overrides.
///
/// Interpretation keeps selected `_atom_site` rows in source order as dense atom
/// order, and writers emit `_atom_site` rows in dense atom order, so coordinate
/// files in the same order address the same atoms after a round trip.
///
/// ```
/// use kekule::mmcif;
///
/// let text = "\
/// data_water
/// loop_
/// _atom_site.group_PDB
/// _atom_site.id
/// _atom_site.type_symbol
/// _atom_site.label_atom_id
/// _atom_site.label_comp_id
/// _atom_site.label_asym_id
/// _atom_site.label_entity_id
/// _atom_site.label_seq_id
/// _atom_site.Cartn_x
/// _atom_site.Cartn_y
/// _atom_site.Cartn_z
/// _atom_site.auth_seq_id
/// _atom_site.auth_asym_id
/// _atom_site.pdbx_PDB_model_num
/// HETATM 1 O O   HOH A 1 . 0.000 0.000 0.000 1 A 1
/// HETATM 2 H H1  HOH A 1 . 0.960 0.000 0.000 1 A 1
/// HETATM 3 H H2  HOH A 1 . -0.240 0.930 0.000 1 A 1
/// ";
/// let document = mmcif::parse_str(text, Default::default())?;
/// let interpretation = document.interpret(Default::default())?;
/// assert_eq!(interpretation.model().atom_count(), 3);
///
/// let (written, report) = mmcif::write([interpretation.model()], Default::default())?;
/// assert!(report.omissions().is_empty());
/// let reread = mmcif::parse_str(&written, Default::default())?.interpret(Default::default())?;
/// assert_eq!(reread.model().atom_count(), 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod mmcif {
    use crate::diagnostics::FormatError;
    pub use crate::io::{
        MmcifAtomField, MmcifBlock, MmcifBlockSource, MmcifCellReason, MmcifDocument,
        MmcifEntityClassifications, MmcifEntry, MmcifItem, MmcifLoopTable, MmcifParseOptions,
        MmcifSampleField, MmcifValue, MmcifWriteIssue, MmcifWriteOptions,
    };

    /// Parses a structural mmCIF data document without assigning molecular
    /// meaning. `options` sets the syntax and resource limits.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when the input exceeds `options`,
    /// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
    /// for malformed CIF, and
    /// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
    /// for contradictory structure, such as a repeated block or item.
    pub fn parse_str(
        input: &str,
        options: MmcifParseOptions,
    ) -> Result<MmcifDocument, FormatError> {
        crate::io::parse_mmcif_str(input, options)
    }

    /// Reads and parses CIF/mmCIF text from `reader`.
    ///
    /// # Errors
    ///
    /// As [`parse_str`]. The configured byte limit is enforced while reading,
    /// before the input is held in memory; a read failure or non-UTF-8 input
    /// fails with [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn parse_reader(
        reader: impl std::io::Read,
        options: MmcifParseOptions,
    ) -> Result<MmcifDocument, FormatError> {
        crate::io::parse_mmcif_reader(reader, options)
    }

    /// Writes one data block per source and returns the text with a
    /// [`WriteReport`](crate::diagnostics::WriteReport) of what the options
    /// chose to omit, summed over the blocks.
    ///
    /// Sources convert from models (one block), ensembles (one multi-model
    /// block), and interpretations (keeping their source reports). Entity kinds
    /// derive from canonical topology classification unless a source carries
    /// explicit classifications or reports ([`MmcifBlockSource`]). Blocks are
    /// named by [`MmcifWriteOptions::block_name`], numbered when there are
    /// several, unless [`MmcifBlockSource::with_block_name`] names them.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind
    /// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
    /// when mmCIF cannot hold the content, and `Inconsistent`,
    /// `InvalidReference` or `InvalidValue` when a source's classifications or
    /// reports do not fit its topology, and `InvalidValue` when two blocks
    /// would share a name.
    ///
    /// ```no_run
    /// use kekule::mmcif::{self, MmcifBlockSource, MmcifWriteOptions};
    /// # fn example(model: &kekule::structure::Model, other: &kekule::structure::Model,
    /// #     interpretation: &kekule::macromolecular::StructureInterpretation)
    /// #     -> Result<(), kekule::diagnostics::FormatError> {
    /// let (one, _) = mmcif::write([model], MmcifWriteOptions::default())?;
    /// let (two, _) = mmcif::write([model, other], MmcifWriteOptions::default())?;
    /// let (faithful, _) = mmcif::write([interpretation], MmcifWriteOptions::default())?;
    /// # let _ = (one, two, faithful);
    /// # Ok(())
    /// # }
    /// ```
    pub fn write<'a, B: Into<MmcifBlockSource<'a>>>(
        blocks: impl IntoIterator<Item = B>,
        options: MmcifWriteOptions,
    ) -> Result<(String, crate::diagnostics::WriteReport), FormatError> {
        crate::io::write_mmcif(blocks, options)
    }

    /// Streams the blocks [`write()`] would return and returns its report.
    /// Values reach `writer` as
    /// they are formatted, so memory stays bounded however large the
    /// structure; wrap an unbuffered sink such as a `File` in a
    /// [`BufWriter`](std::io::BufWriter).
    ///
    /// # Errors
    ///
    /// The errors of [`write()`]; output before the failure has already
    /// reached `writer`. A failing `writer` gives
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn write_to<'a, B: Into<MmcifBlockSource<'a>>>(
        writer: &mut impl std::io::Write,
        blocks: impl IntoIterator<Item = B>,
        options: MmcifWriteOptions,
    ) -> Result<crate::diagnostics::WriteReport, FormatError> {
        crate::io::write_mmcif_to(writer, blocks, options)
    }
}

/// PDB-format parsing, interpretation, and writing.
///
/// A PDB file holds one entry, whose coordinate models share one
/// interpretation scope. Parsing types the unit-cell, coordinate,
/// chain-termination, model, and connection records and keeps every other
/// record verbatim;
/// interpretation runs the pipeline [`mmcif`] uses and returns the same
/// [`macromolecular`] reports.
///
/// The format has no entity or label records. A chain's polymer is the run of
/// residues its `TER` record ends, or, without one, the residues from its first
/// to its last `ATOM` record; water and single residues are never polymer.
/// Residues in such a run are adjacent in sequence. `LINK` and `SSBOND`
/// records are covalent connections; a `LINK` to a metal is coordination and
/// becomes a dative bond. `CONECT` records list connectivity without bond
/// types: within a residue they supplement its template and stand in for
/// distance inference, a partner listed two or three times has a double or
/// triple bond, and a bond between a metal and a non-metal is dative. Atoms
/// without an element symbol take the one their name's alignment implies,
/// and the report lists them. A `CRYST1` cell is the model's periodic cell
/// under the rule mmCIF `_cell` follows; otherwise the report says why not.
///
/// Interpretation keeps coordinate records in source order as dense atom order,
/// and [`write()`](pdb::write) emits records in dense atom order, so a written
/// entry interprets back to the same atoms in the same order.
///
/// ```
/// use kekule::pdb;
///
/// let text = "\
/// HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
/// HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H
/// HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H
/// END
/// ";
/// let document = pdb::parse_str(text, Default::default())?;
/// let interpretation = document.interpret(Default::default())?;
/// assert_eq!(interpretation.model().atom_count(), 3);
/// assert_eq!(interpretation.topology().instance_count(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod pdb {
    use crate::diagnostics::FormatError;
    pub use crate::io::{
        PdbAtom, PdbAtomField, PdbAtomKind, PdbCellReason, PdbConect, PdbCryst1, PdbDocument,
        PdbLink, PdbModel, PdbParseOptions, PdbPartner, PdbRecord, PdbSource, PdbSsbond, PdbTer,
        PdbWriteIssue, PdbWriteOptions,
    };

    /// Parses a PDB entry into typed records without assigning molecular
    /// meaning. `options` sets the resource limits.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of kind
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when the input exceeds `options`, and `Syntax`, `InvalidValue`,
    /// `Inconsistent` or `Unsupported` when a typed record is malformed,
    /// contradicts another, or uses content Kekule does not read.
    pub fn parse_str(input: &str, options: PdbParseOptions) -> Result<PdbDocument, FormatError> {
        crate::io::parse_pdb_str(input, options)
    }

    /// Reads and parses a PDB entry from `reader`.
    ///
    /// # Errors
    ///
    /// As [`parse_str`]. The configured byte limit is enforced while reading,
    /// before the input is held in memory; a read failure or non-UTF-8 input
    /// fails with [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn parse_reader(
        reader: impl std::io::Read,
        options: PdbParseOptions,
    ) -> Result<PdbDocument, FormatError> {
        crate::io::parse_pdb_reader(reader, options)
    }

    /// Writes one entry: a model or collection item without `MODEL` records,
    /// or an ensemble with one `MODEL` per member. Returns the text with a
    /// [`WriteReport`](crate::diagnostics::WriteReport) of what the options
    /// chose to omit.
    ///
    /// Records follow dense atom order, which interpretation keeps, so
    /// interpreting the output gives the same atoms in the same order.
    /// Identity comes from the hierarchy:
    /// - author chain, residue, and atom identifiers, falling back to label
    ///   identifiers;
    /// - residues of protein and nucleic-acid instances with an amino-acid or
    ///   nucleotide class are `ATOM` records, everything else `HETATM`;
    /// - a `TER` record, which takes the next serial number, follows each
    ///   chain's last polymer atom;
    /// - an instance without atom sites becomes one `MOL` (or `HOH`) residue
    ///   without a chain identifier, and atoms without a name are named by
    ///   their element and 1-based molecule-local number.
    ///
    /// Serial numbers beyond 99,999 and residue numbers beyond 9,999 use
    /// hybrid-36. Coordinates are written to 0.001 Å, occupancies and B
    /// factors to 0.01, and absent occupancies and B factors are blank. A
    /// representable periodic cell is a `CRYST1` record in space group P 1.
    ///
    /// `CONECT` lists every bond a reader could not derive, each from both
    /// atoms and once per step of bond order:
    /// - bonds within `HETATM` residues and residues the bundled
    ///   [`ComponentLibrary`](crate::components::ComponentLibrary) does not
    ///   cover;
    /// - bonds between residues other than the linkage of standard residues
    ///   adjacent in one chain, such as disulfides and ligand attachments.
    ///
    /// Bonds within standard residues the dictionary covers, polymer
    /// linkage, and bonds within water are derived by readers and not listed.
    ///
    /// # Errors
    ///
    /// A write-stage [`FormatError`] of kind
    /// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
    /// with a [`PdbWriteIssue`] when PDB cannot hold the content. Identifiers
    /// that do not fit their columns always fail: chain identifiers longer
    /// than one character, residue names longer than three, atom names longer
    /// than four, insertion codes longer than one, residue numbers outside
    /// -999..=2,436,111, more than 87,440,031 atom and `TER` records, and
    /// residues whose atoms are not consecutive. So do values outside their
    /// columns and formal charges outside -9..=9. Fixed implicit-hydrogen
    /// counts, other represented chemistry, ensemble or frame state, and
    /// cells PDB cannot carry fail unless [`PdbWriteOptions`] omits them,
    /// each under its own policy, and every omission is in the report.
    /// Nothing is written when validation fails.
    ///
    /// ```
    /// use kekule::pdb;
    ///
    /// let text = "\
    /// HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
    /// HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H
    /// HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H
    /// END
    /// ";
    /// let model = pdb::parse_str(text, Default::default())?
    ///     .interpret(Default::default())?
    ///     .into_model();
    /// let (written, report) = pdb::write(&model, Default::default())?;
    /// assert_eq!(written, text);
    /// assert!(report.omissions().is_empty());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write<'a>(
        source: impl Into<PdbSource<'a>>,
        options: PdbWriteOptions,
    ) -> Result<(String, crate::diagnostics::WriteReport), FormatError> {
        crate::io::write_pdb(source, options)
    }

    /// Streams the entry [`write()`] would return and returns its report.
    /// The entry is validated
    /// before the first record reaches `writer`; wrap an unbuffered sink
    /// such as a `File` in a [`BufWriter`](std::io::BufWriter).
    ///
    /// # Errors
    ///
    /// The errors of [`write()`]; a failing `writer` gives
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn write_to<'a>(
        writer: &mut impl std::io::Write,
        source: impl Into<PdbSource<'a>>,
        options: PdbWriteOptions,
    ) -> Result<crate::diagnostics::WriteReport, FormatError> {
        crate::io::write_pdb_to(writer, source, options)
    }
}

/// XYZ parsing and structure interpretation.
///
/// An XYZ file is a sequence of frames, each an atom count, a comment line,
/// and one element and three coordinates per atom. It states no bonds, so
/// interpretation infers connectivity from interatomic distances with the
/// rule macromolecular interpretation applies to residues without a
/// template: two atoms are bonded when they lie within the sum of their
/// covalent radii (Cordero et al. 2008) plus 0.45 Å, and farther apart than
/// 0.4 Å; an element without a tabulated radius is taken as 0.77 Å.
///
/// Interpretation never guesses chemistry the file does not state. Every
/// inferred bond is single, except that a bond between a metal and a
/// non-metal is dative, as for PDB `CONECT` records; bond orders, formal
/// charges, radicals, and stereo are not assigned. XYZ lists every atom, so
/// every atom has a fixed implicit hydrogen count of zero and perception
/// never adds hydrogens. Each connected component becomes one molecule
/// instance, under a synthetic hierarchy of one `UNL` residue per molecule
/// in chain `A`. The [`XyzInterpretationReport`](xyz::XyzInterpretationReport)
/// lists every inferred bond with its distance and cutoff, atoms that used
/// the fallback radius, atoms too close to bond, and omitted content. The
/// comment line is kept verbatim and not interpreted, and columns after the
/// coordinates, such as extended XYZ properties, are rejected unless
/// [`XyzInterpretOptions::with_extra_columns`](xyz::XyzInterpretOptions::with_extra_columns)
/// reports them as omitted.
///
/// One frame becomes a [`Model`](crate::structure::Model); several frames of
/// the same atoms become an equally weighted
/// [`Ensemble`](crate::structure::Ensemble), and every frame must infer the
/// same bonds. Interpretation keeps atom lines in source order as dense atom
/// order, so a trajectory reader such as `kekule-traj`'s `XyzReader` can
/// stream the frames of the same file against the interpreted topology when
/// the file spells every element as a symbol in its standard letter case,
/// the only element tokens that reader accepts; interpretation also accepts
/// other letter cases and atomic numbers.
///
/// ```
/// use kekule::molecule::BondOrder;
/// use kekule::xyz;
///
/// let text = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n";
/// let document = xyz::parse_str(text, Default::default())?;
/// let interpretation = document.interpret(Default::default())?;
/// assert_eq!(interpretation.topology().instance_count(), 1);
/// let report = interpretation.report();
/// assert_eq!(report.comment(), "water");
/// assert_eq!(report.inferred_bonds().len(), 2);
/// assert!(report
///     .inferred_bonds()
///     .iter()
///     .all(|bond| bond.order() == BondOrder::Single));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod xyz {
    use crate::diagnostics::FormatError;
    pub use crate::io::{
        XyzAtom, XyzDocument, XyzEnsembleInterpretOptions, XyzEnsembleInterpretation, XyzFrame,
        XyzFrameSelection, XyzInferredBond, XyzInterpretIssue, XyzInterpretOptions,
        XyzInterpretation, XyzInterpretationReport, XyzParseOptions, DEFAULT_XYZ_MAX_ATOM_CONTACTS,
    };

    /// Parses XYZ frames without assigning molecular meaning. `options` sets
    /// the resource limits.
    ///
    /// Blank lines may follow the last frame; any other content after them
    /// is a syntax error. Empty input has no frames.
    ///
    /// # Errors
    ///
    /// Returns a [`FormatError`] of format
    /// [`Format::Xyz`](crate::diagnostics::Format::Xyz) and kind
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when the input exceeds `options`,
    /// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
    /// for a malformed count line or a truncated frame, and
    /// [`DiagnosticKind::InvalidValue`](crate::diagnostics::DiagnosticKind::InvalidValue)
    /// for a frame of zero atoms or a coordinate that is not a finite number.
    pub fn parse_str(input: &str, options: XyzParseOptions) -> Result<XyzDocument, FormatError> {
        crate::io::parse_xyz_str(input, options)
    }

    /// Reads and parses XYZ frames from `reader`.
    ///
    /// # Errors
    ///
    /// As [`parse_str`]. The configured byte limit is enforced while reading,
    /// before the input is held in memory; a read failure or non-UTF-8 input
    /// fails with [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
    pub fn parse_reader(
        reader: impl std::io::Read,
        options: XyzParseOptions,
    ) -> Result<XyzDocument, FormatError> {
        crate::io::parse_xyz_reader(reader, options)
    }
}

/// Kekule's native binary serialization: exact, versioned persistence.
///
/// Scientific formats exchange chemistry with other tools; they are not exact
/// persistence. The native format stores one [`Molecule`](crate::molecule::Molecule),
/// [`Topology`](crate::topology::Topology), [`Model`](crate::structure::Model),
/// [`Ensemble`](crate::structure::Ensemble), or in-memory
/// [`Trajectory`](crate::structure::Trajectory) losslessly: represented
/// chemistry, definitions and instances, semantic IDs, dense order,
/// hierarchy, classifications with their explicit or inferred assignment,
/// properties at every scope with units, coordinates as exact bits, and every
/// structure and collection field. Installed perception is stored with its
/// provenance unless
/// [`NativeWriteOptions::with_perception`](native::NativeWriteOptions::with_perception)
/// omits it.
///
/// [`native::write_to`] streams to a [`std::io::Write`] and [`native::read`] reads exactly one
/// object from a [`std::io::Read`]; neither buffers the caller's I/O, so wrap
/// unbuffered files in [`std::io::BufWriter`] or [`std::io::BufReader`].
/// Decoding validates untrusted input through the owners' reconstruction rules
/// under the resource limits of [`NativeReadOptions`](native::NativeReadOptions) and fails with a
/// [`FormatError`](crate::diagnostics::FormatError), never a panic. A decoded
/// topology keeps every semantic ID but is a new publication with a fresh
/// layout identity. The byte layout and versioning policy are specified in
/// `docs/native-format.md`.
///
/// ```
/// use kekule::native::{self, NativeReadOptions, NativeWriteOptions};
/// use kekule::smiles;
///
/// let mut molecule = smiles::to_molecules("c1ccccc1O")?.remove(0);
/// molecule.perceive()?;
/// let bytes = native::write(&molecule, NativeWriteOptions::default())?;
/// let decoded = native::read(bytes.as_slice(), NativeReadOptions::default())?.into_molecule()?;
/// assert_eq!(decoded, molecule);
/// assert_eq!(decoded.perception(), molecule.perception());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod native {
    use crate::diagnostics::FormatError;
    pub use crate::io::{
        NativeKind, NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions,
        DEFAULT_NATIVE_MAX_ATOMS, DEFAULT_NATIVE_MAX_BONDS, DEFAULT_NATIVE_MAX_INPUT_BYTES,
        DEFAULT_NATIVE_MAX_ITEMS, DEFAULT_NATIVE_MAX_STRING_BYTES, NATIVE_FORMAT_VERSION,
        NATIVE_MAGIC,
    };

    /// Encodes one object into a new byte vector.
    ///
    /// # Errors
    ///
    /// As [`write_to`], except that a vector never fails as a sink.
    pub fn write<'a>(
        source: impl Into<NativeSource<'a>>,
        options: NativeWriteOptions,
    ) -> Result<Vec<u8>, FormatError> {
        crate::io::write_native(source, options)
    }

    /// Encodes one object into `writer`, staging output through a fixed
    /// 64 KiB buffer rather than holding the whole encoding.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind
    /// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io) when
    /// `writer` fails, or
    /// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
    /// when a list or string is too long for its 32-bit length field.
    pub fn write_to<'a>(
        writer: &mut impl std::io::Write,
        source: impl Into<NativeSource<'a>>,
        options: NativeWriteOptions,
    ) -> Result<(), FormatError> {
        crate::io::write_native_to(writer, source, options)
    }

    /// Decodes one object from `reader`, reading exactly its bytes, so
    /// further objects may follow in the same stream.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of format
    /// [`Format::Native`](crate::diagnostics::Format::Native) with a
    /// [`byte_offset`](FormatError::byte_offset):
    /// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
    /// for truncated input or wrong magic bytes,
    /// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
    /// with [`FormatErrorDetail::NativeVersion`](crate::diagnostics::FormatErrorDetail::NativeVersion)
    /// for another format version,
    /// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
    /// when `options` is exceeded, and `InvalidValue`, `InvalidReference`,
    /// `Inconsistent`, or `Chemistry` when stored values do not reconstruct
    /// a valid owner or the checksum does not match.
    pub fn read(
        mut reader: impl std::io::Read,
        options: NativeReadOptions,
    ) -> Result<NativeObject, FormatError> {
        crate::io::read_native(&mut reader, options)
    }
}

/// Derived valence, rings, aromaticity, conjugation and explicit resonance work.
///
/// [`crate::molecule::Molecule::perceive`] installs Kekule's default transactional perception
/// profile. The nested modules expose the individual expert algorithms.
/// Perception does not alter represented graph chemistry and is invalidated by
/// relevant graph edits.
///
/// A [`Perception`](perception::Perception) holds one optional section per
/// step (valence, rings, aromaticity, conjugation, resonance and CIP). Each
/// section records the model or policy that produced it, and an absent
/// section has not been perceived.
///
/// ```
/// use kekule::smiles;
///
/// let mut phenol = smiles::to_molecules("Oc1ccccc1")?.remove(0);
/// assert!(phenol.perception().aromaticity().is_none());
/// phenol.perceive()?;
///
/// let perception = phenol.perception();
/// let rings = perception.rings().and_then(|rings| rings.ring_set());
/// assert_eq!(rings.map(|rings| rings.len()), Some(1));
/// let aromaticity = perception.aromaticity().expect("perceived");
/// assert_eq!(aromaticity.atoms().len(), 6);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod perception {
    /// Connected conjugated groups and explicit RDKit-like contributor enumeration.
    pub mod resonance {
        pub use crate::algorithms::{
            enumerate_resonance, perceive_resonance, ResonanceContributor, ResonanceError,
            ResonanceFlags, ResonanceOptions, ResonanceStructures,
        };
        pub use crate::molecule::delocalization::{ResonanceGroup, ResonancePerception};
    }
    /// Bond conjugation, including aromatic bonds, without contributor enumeration.
    pub mod conjugation {
        pub use crate::algorithms::{perceive_conjugation, ConjugationError};
        pub use crate::molecule::delocalization::{ConjugationModel, ConjugationPerception};
    }
    pub use crate::chemistry::PerceptionError;
    pub use crate::molecule::perception::{
        DelocalizationIssue, Perception, PerceptionBuildError, PerceptionBuilder,
        PerceptionComponent, PerceptionInstallError,
    };

    /// Expert valence perception for canonical represented chemistry.
    ///
    /// The RDKit-like model derives the inferred contribution to implicit H
    /// from ordinary localized bond orders and represented atom state. It does
    /// not require installed ring or aromaticity perception. Source aromatic
    /// bonds are localized during format interpretation before this layer is
    /// reached.
    pub mod valence {
        pub use crate::algorithms::{
            perceive_valence, represented_valence, ValenceError, ValenceIssue, ValenceOptions,
        };
        pub use crate::molecule::perception::{ValenceModel, ValencePerception, ValenceValidation};
    }

    /// Ring membership and deterministic ring-basis perception.
    ///
    /// Ring perception installs derived state and never changes represented
    /// bond orders.
    pub mod rings {
        pub use crate::algorithms::{
            perceive_rings, RingPerceptionError, RingPerceptionOptions, RingPerceptionScope,
        };
        pub use crate::molecule::perception::{
            MalformedRingReason, Ring, RingBasisModel, RingBasisState, RingMembership,
            RingPerception, RingSet,
        };
    }

    /// Aromaticity perception over canonical localized bond orders.
    ///
    /// Aromatic membership is derived state; Kekule does not store an aromatic
    /// bond order in represented graph chemistry.
    pub mod aromaticity {
        pub use crate::algorithms::{perceive_aromaticity, AromaticityError, AromaticityOptions};
        pub use crate::molecule::perception::{AromaticityModel, AromaticityPerception};
    }
}

/// Represented stereochemistry, its validation, and stereo algorithms.
///
/// The root holds the stereo elements and groups a [`Molecule`](crate::molecule::Molecule)
/// represents, the issues their validation reports, and the shared
/// symmetry-search bounds. Edits, publication, and
/// [`MoleculeEditor::validate`](crate::molecule::MoleculeEditor::validate)
/// share one validator.
/// Expert algorithms live in submodules: [`cip`](crate::stereo::cip) assigns
/// opt-in CIP descriptors, [`coordinates`](crate::stereo::coordinates) infers
/// and materializes stereo from coordinates, and
/// [`cleanup`](crate::stereo::cleanup) detects stereogenic candidates and removes
/// ineligible assertions. Coordinate inference and candidate detection are
/// read-only; coordinate materialization is an explicitly named representation
/// transform.
///
/// ```
/// use kekule::smiles;
/// use kekule::stereo::{cip, StereoDescriptor};
///
/// // L-alanine.
/// let mut alanine = smiles::to_molecules("C[C@H](N)C(=O)O")?.remove(0);
/// assert_eq!(alanine.stereo_elements().count(), 1);
///
/// alanine.perceive()?;
/// cip::assign_cip_descriptors(&mut alanine, Default::default())?;
/// let (center, _) = alanine.stereo_elements().next().unwrap();
/// let cip_state = alanine.perception().cip().expect("assigned");
/// assert_eq!(cip_state.cip_descriptor(center), Some(StereoDescriptor::S));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod stereo {
    pub use crate::algorithms::{
        StereoSearchError, StereoSearchOptions, StereoValidationError, StereoValidationIssue,
    };
    pub use crate::molecule::ids::{StereoElementId, StereoGroupId};
    pub use crate::molecule::stereo::{
        AxisOrientation, AxisStereo, DoubleBondOrientation, DoubleBondStereo, StereoCarrier,
        StereoDescriptor, StereoElement, StereoElementKind, StereoGroup, StereoGroupKind,
        TetrahedralOrientation, TetrahedralStereo,
    };

    /// CIP descriptor assignment and its installed perception state.
    ///
    /// CIP descriptors remain opt-in derived state. Assignment requires known
    /// implicit hydrogen counts: run valence perception first, or it reports
    /// [`CipAssignmentIssue::UnknownHydrogenCount`](crate::stereo::cip::CipAssignmentIssue::UnknownHydrogenCount)
    /// instead of ranking with
    /// none.
    pub mod cip {
        pub use crate::algorithms::{
            assign_cip_descriptors, CipAssignment, CipAssignmentError, CipAssignmentIssue,
            CipAssignmentOptions, CipAssignmentReport, CipRankingError, CipSkipped,
            CipSkippedReason,
        };
        pub use crate::molecule::perception::CipPerception;
    }

    /// Read-only stereo inference from coordinates and its explicit
    /// materialization into represented stereo.
    pub mod coordinates {
        pub use crate::algorithms::{
            infer_coordinate_stereo, materialize_coordinate_stereo, CoordinateStereoError,
            CoordinateStereoMaterializationReport, CoordinateStereoOptions, CoordinateStereoResult,
        };
    }

    /// Stereogenic candidate detection and chemical cleanup of represented
    /// stereo assertions.
    pub mod cleanup {
        pub use crate::algorithms::{
            cleanup_stereo, detect_stereo_candidates, StereoCandidate, StereoCleanupReport,
        };
    }
}

/// Explicit small-molecule hydrogen topology transforms.
///
/// Explicit hydrogens are graph atoms; an atom's implicit hydrogens are either
/// a fixed count or wholly inferred by valence perception. Addition converts
/// implicit H to explicit atoms; removal suppresses eligible explicit atoms
/// without losing chemical information. These operations do not change
/// protonation states.
///
/// Addition requires current inference for inferred counts, unless
/// [`AddHydrogensScope::FixedOnly`](hydrogens::AddHydrogensScope::FixedOnly) is selected. Fixed counts need no perception. Removal needs current inference
/// on affected parents whose counts are inferred, and keeps a parent inferred
/// only when inference reproduces its collapsed count; otherwise it fixes it.
/// Successful topology changes invalidate perception state; recompute it before
/// reading inferred counts.
///
/// ```
/// use kekule::smiles;
///
/// let mut ethanol = smiles::to_molecules("CCO")?.remove(0);
/// ethanol.perceive()?;
/// let added = ethanol.add_hydrogens(Default::default())?;
/// assert_eq!(added.added().len(), 6);
/// assert_eq!(ethanol.atom_count(), 9);
///
/// ethanol.perceive()?;
/// let removed = ethanol.remove_hydrogens()?;
/// assert_eq!(removed.removed().len(), 6);
/// assert_eq!(ethanol.atom_count(), 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod hydrogens {
    pub use crate::algorithms::{
        AddHydrogensOptions, AddHydrogensReport, AddHydrogensScope, AddedHydrogen,
        AddedHydrogenOrigin, HydrogenCountAdjustment, HydrogenTransformError,
        RemoveHydrogensReport, RemovedHydrogen, RetainedHydrogen, RetainedHydrogenReason,
    };
}

/// Common foundational types for small examples and interactive use.
///
/// The prelude is intentionally small: the molecule, its editor and parts, the
/// topology and model that hold it, and [`Quantity`](crate::units::Quantity).
/// It is the only module that repeats items exported elsewhere. Format APIs,
/// identifiers, and algorithms remain available through their focused modules.
///
/// ```
/// use kekule::prelude::*;
/// use kekule::{geometry::Point3, units::ANGSTROM};
///
/// let mut editor = MoleculeEditor::new();
/// let carbon = editor.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
/// let oxygen = Atom::new(Element::O).with_hydrogens(ImplicitHydrogens::Fixed(1));
/// let oxygen = editor.add_atom(oxygen)?;
/// editor.add_bond(carbon, oxygen, BondOrder::Single)?;
/// let molecule: Molecule = editor.finish()?;
/// let (_, bond): (_, &Bond) = molecule.bonds().next().unwrap();
/// assert_eq!(bond.order, BondOrder::Single);
///
/// let topology = Topology::from_molecule(molecule.clone())?;
/// let positions = Positions::new(Quantity::new(
///     vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.43, 0.0, 0.0)],
///     ANGSTROM,
/// ))?;
/// let model = Model::from_molecule(molecule, &positions)?;
/// assert_eq!(model.atom_count(), topology.atom_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub mod prelude {
    pub use crate::molecule::{
        Atom, Bond, BondOrder, Element, ImplicitHydrogens, Molecule, MoleculeEditor,
    };
    pub use crate::structure::{Model, Positions};
    pub use crate::topology::Topology;
    pub use crate::units::Quantity;
}

#[cfg(test)]
mod tests;

// Compile and run the README examples so they cannot drift from the API. The
// README is the workspace root's in the repository and the package root's once
// packaged; `CARGO_PKG_README` names it relative to the manifest in both.
#[cfg(doctest)]
#[doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/", env!("CARGO_PKG_README")))]
struct ReadmeExamples;
