//! Pure-Rust foundations for molecular graphs, chemical perception, structure
//! I/O, structural bioinformatics, and molecular modelling.
//!
//! # Object model
//!
//! Kekule keeps chemical identity, system organization, and geometry separate:
//!
//! - [`core::Molecule`] is one non-empty connected molecular graph. It owns
//!   represented chemistry, derived [`core::Perception`], and
//!   geometry-independent properties.
//! - [`topology::Topology`] is one coordinate-free system containing one or
//!   more explicit molecule instances plus an optional biological hierarchy.
//! - [`structure::Model`] is one realization of a topology: positions, an
//!   optional periodic cell, and realization-scoped properties.
//! - [`structure::Ensemble`] stores several non-temporal realizations of one
//!   shared topology. Ordered trajectories live in the `kekule-traj`
//!   companion crate.
//!
//! A disconnected salt, solvent box, or protein-ligand complex is therefore a
//! [`topology::Topology`] containing several connected molecules, not one
//! disconnected [`core::Molecule`].
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
//! | Construct or edit one connected molecule | [`core::MoleculeEditor`] | `finish()` |
//! | Assemble complete molecules and reusable instances | [`topology::TopologyBuilder`] | `build()` |
//! | Assemble a system with explicit geometry | [`structure::ModelBuilder`] | `build()` |
//! | Edit system atoms, bonds, components, and hierarchy | [`topology::TopologyEditor`] | `finish()` |
//! | Edit structure while coordinating one realization | [`structure::ModelEditor`] | `finish()` |
//! | Change geometry or realization annotations | [`structure::Model`] setters | Immediate checked update |
//!
//! Builders and editors offer non-consuming `validate()` and recoverable
//! `try_build()` / `try_finish()`. Use `edit()` for a detached draft or
//! `into_editor()` to move an owner into one. System editors resolve source IDs
//! to stable draft-only editing handles and publish completed values through
//! `finish()`. A bond deletion can split a system molecule;
//! a bond addition can join two occurrences. New model atoms require coordinates.
//!
//! Coordinate-dependent algorithms consume [`structure::ModelView`]. A model,
//! ensemble member, trajectory frame, or reusable trajectory buffer can
//! therefore share kernels without copying coordinates. APIs that address atoms
//! by index accept any topology snapshot sharing their layout
//! ([`topology::Topology::shares_layout`]), so perception never detaches
//! selections, buffers, or prepared potentials. Use
//! [`topology::Topology::same_layout`] only when complete static layout equality
//! between independent publications is intended.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

macro_rules! fixed_u32_id {
    ($name:ident) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(u32);

        impl $name {
            pub const fn new(raw: u32) -> Self {
                Self(raw)
            }

            pub const fn raw(self) -> u32 {
                self.0
            }

            pub const fn index(self) -> usize {
                self.0 as usize
            }
        }
    };
    ($name:ident, $display:literal) => {
        fixed_u32_id!($name);

        impl std::fmt::Display for $name {
            fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(formatter, concat!($display, "{}"), self.0)
            }
        }
    };
}

mod algorithms;
pub mod alignment;
mod chemistry;
pub mod core;
pub mod descriptors;
pub mod dssp;
pub mod geometry;
mod io;
pub mod properties;
pub mod query;
pub mod structure;
pub mod topology;
pub mod units;

/// Syntax-independent substructure matching algorithms.
///
/// Matching consumes `query::QueryGraph` and current target perception state;
/// it never invokes parsing, interpretation, or perception implicitly.
pub mod substructure {
    pub use crate::algorithms::{
        find_match, find_matches, find_matches_with_options, find_topology_matches,
        find_topology_matches_with_options, visit_matches, visit_matches_with_options,
        MatchCompletion, PreparedTarget, PreparedTopologyTarget, QueryMatch, QueryPerception,
        SubstructureMatchError, SubstructureMatchOptions, SubstructureMatchWork, TaggedMatchError,
        TaggedQuery, TopologyQueryMatch, MAX_SUBSTRUCTURE_QUERY_ATOMS,
    };
}

/// SMILES parsing, interpretation, and molecule/topology writing.
///
/// [`smiles::parse_str`] preserves source syntax in a
/// [`smiles::SmilesDocument`], while [`smiles::interpret`] publishes canonical
/// connected molecules. [`smiles::to_molecules`] is the concise
/// parse-and-interpret path. Dot-separated components remain separate molecules
/// and perception is never run implicitly.
pub mod smiles {
    use std::fmt;
    use std::sync::Arc;

    use crate::core::Molecule;
    pub use crate::io::{
        MolWriteError, MolWriteErrorKind, SmilesAtomMapping, SmilesBondMapping,
        SmilesComponentCountError, SmilesComponentInterpretation, SmilesDocument,
        SmilesDocumentToken, SmilesDocumentTokenKind, SmilesInterpretError, SmilesInterpretation,
        SmilesInterpretationReport, SmilesParseError, SmilesParseOptions,
    };
    use crate::topology::{Topology, TopologyBuildError};

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
        /// Complete canonical labeling is bounded by 100,000 search states,
        /// 50,000,000 atom/edge/twin visits, and 2,000,000 pending atom labels.
        /// Serialization additionally bounds the input to 2,000,000 combined
        /// atom/bond slots and a complexity score `2*n*(n+2*m)` of at most
        /// 50,000,000 for `n` live atoms and `m` live bonds. Exceeding a bound
        /// returns [`MolWriteErrorKind::ResourceLimit`] without a partial result.
        Canonical,
    }

    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    pub struct SmilesWriteOptions {
        pub mode: SmilesWriteMode,
    }

    impl SmilesWriteOptions {
        pub const fn ordinary() -> Self {
            Self {
                mode: SmilesWriteMode::Ordinary,
            }
        }

        pub const fn isomeric() -> Self {
            Self {
                mode: SmilesWriteMode::Isomeric,
            }
        }

        pub const fn canonical() -> Self {
            Self {
                mode: SmilesWriteMode::Canonical,
            }
        }
    }

    /// Error produced by the concise SMILES parse-and-interpret convenience.
    #[derive(Debug, Clone, PartialEq)]
    #[non_exhaustive]
    pub enum SmilesReadError {
        Parse(SmilesParseError),
        Interpret(SmilesInterpretError),
        Topology(Box<TopologyBuildError>),
    }

    impl fmt::Display for SmilesReadError {
        fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Self::Parse(error) => write!(formatter, "{error}"),
                Self::Interpret(error) => write!(formatter, "{error}"),
                Self::Topology(error) => write!(formatter, "{error}"),
            }
        }
    }

    impl std::error::Error for SmilesReadError {
        fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
            match self {
                Self::Parse(error) => Some(error),
                Self::Interpret(error) => Some(error),
                Self::Topology(error) => Some(error.as_ref()),
            }
        }
    }

    impl From<SmilesParseError> for SmilesReadError {
        fn from(error: SmilesParseError) -> Self {
            Self::Parse(error)
        }
    }

    impl From<SmilesInterpretError> for SmilesReadError {
        fn from(error: SmilesInterpretError) -> Self {
            Self::Interpret(error)
        }
    }

    impl From<TopologyBuildError> for SmilesReadError {
        fn from(error: TopologyBuildError) -> Self {
            Self::Topology(Box::new(error))
        }
    }

    /// Parses one record, preserving base SMILES, an optional CX extension, and its name.
    /// Extension text is retained without validating its field syntax or semantics.
    pub fn parse_str(input: &str) -> Result<SmilesDocument, SmilesParseError> {
        crate::io::parse_smiles_document(input)
    }

    /// Parses source text with explicit syntax and resource-limit options.
    pub fn parse_str_with_options(
        input: &str,
        options: SmilesParseOptions,
    ) -> Result<SmilesDocument, SmilesParseError> {
        crate::io::parse_smiles_document_with_options(input, options)
    }

    /// Interprets parsed syntax as canonical represented chemistry.
    ///
    /// Bracket atoms have fixed hydrogen counts. After localizing source aromatic
    /// bonds, their valence deficits determine radical-electron occupancy using
    /// the RDKit-like octet/duet and allowed-valence conventions. This does not
    /// assert a spin multiplicity. No perception is installed implicitly.
    /// Explicit CX radical declarations override inferred occupancy and preserve
    /// specified spin. CX absolute, AND and OR groups reference source atom
    /// indices; the `r` flag relates otherwise ungrouped tetrahedral centers.
    /// Stereo relationships across disconnected molecules are not representable.
    /// Unsupported CX fields return an error instead of being discarded.
    pub fn interpret(
        document: &SmilesDocument,
    ) -> Result<SmilesInterpretation, SmilesInterpretError> {
        document.interpret()
    }

    /// Explicitly interprets only the base SMILES and records the omitted CX extension.
    /// This projection can lose chemical information, including stereochemical groups.
    pub fn interpret_base(
        document: &SmilesDocument,
    ) -> Result<SmilesInterpretation, SmilesInterpretError> {
        document.interpret_base()
    }

    /// Parses and interprets one SMILES record into source-ordered connected components.
    ///
    /// This is the concise form of [`parse_str`] followed by [`interpret`].
    /// Dot-delimited components remain separate, and no perception is run implicitly.
    pub fn to_molecules(input: &str) -> Result<Vec<Molecule>, SmilesReadError> {
        let document = parse_str(input)?;
        Ok(document.interpret()?.into_molecules())
    }

    /// Parses and interprets one SMILES record as a shared coordinate-free
    /// topology.
    ///
    /// Every connected component becomes one explicit molecule
    /// occurrence in source order. No hierarchy or perception is fabricated.
    pub fn to_topology(input: &str) -> Result<Arc<Topology>, SmilesReadError> {
        let document = parse_str(input)?;
        Ok(Arc::new(document.interpret()?.into_topology()?))
    }

    /// What one SMILES record describes: one connected molecule, or every
    /// molecule occurrence of a topology joined with `.`.
    #[derive(Debug, Clone, Copy)]
    pub enum SmilesSource<'a> {
        Molecule(&'a Molecule),
        Topology(&'a Topology),
    }

    impl<'a> From<&'a Molecule> for SmilesSource<'a> {
        fn from(molecule: &'a Molecule) -> Self {
            Self::Molecule(molecule)
        }
    }

    impl<'a> From<&'a Topology> for SmilesSource<'a> {
        fn from(topology: &'a Topology) -> Self {
            Self::Topology(topology)
        }
    }

    impl<'a> From<&'a Arc<Topology>> for SmilesSource<'a> {
        fn from(topology: &'a Arc<Topology>) -> Self {
            Self::Topology(topology)
        }
    }

    /// Writes one SMILES record.
    ///
    /// [`SmilesWriteMode`] selects ordinary, isomeric, or canonical output.
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
    /// ```
    /// use kekule::smiles::{self, SmilesWriteOptions};
    /// let molecule = smiles::to_molecules("OCC")?.remove(0);
    /// let canonical = smiles::write(&molecule, SmilesWriteOptions::canonical())?;
    /// assert_eq!(canonical, smiles::write(&smiles::to_molecules("CCO")?[0], SmilesWriteOptions::canonical())?);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn write<'a>(
        source: impl Into<SmilesSource<'a>>,
        options: SmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        match source.into() {
            SmilesSource::Molecule(molecule) => match options.mode {
                SmilesWriteMode::Ordinary => crate::io::write_smiles(molecule),
                SmilesWriteMode::Isomeric => crate::io::write_isomeric_smiles(molecule),
                SmilesWriteMode::Canonical => crate::io::write_canonical_smiles(molecule),
            },
            SmilesSource::Topology(topology) => {
                crate::io::smiles::write_topology(topology, options)
            }
        }
    }
}

/// Molfile parsing, interpretation, and V2000/V3000 writing.
///
/// A parsed [`molfile::MolfileDocument`] is format state. Its interpretation can
/// project to connected molecules or to a geometry-bearing
/// [`crate::structure::Model`]. Writers reject chemistry that the selected
/// Molfile version cannot represent.
pub mod molfile {
    pub use crate::io::{
        MolWriteError, MolWriteErrorKind, MolfileAtomMapping, MolfileBondMapping, MolfileDocument,
        MolfileHeader, MolfileInterpretError, MolfileInterpretation, MolfileInterpretationReport,
        MolfileInterpretationWarning, MolfileLine, MolfileParseError, MolfileParseOptions,
        MolfileSource, MolfileVersion, MolfileWriteOptions, MolfileWriteVersion,
    };

    /// Parses one V2000 or V3000 Molfile without assigning canonical meaning.
    pub fn parse_str(input: &str) -> Result<MolfileDocument, MolfileParseError> {
        crate::io::parse_molfile_document(input)
    }

    /// Parses one Molfile with explicit syntax and resource-limit options.
    pub fn parse_str_with_options(
        input: &str,
        options: MolfileParseOptions,
    ) -> Result<MolfileDocument, MolfileParseError> {
        crate::io::parse_molfile_document_with_options(input, options)
    }

    /// Interprets a parsed Molfile into canonical molecules, geometry, and a report.
    /// Unmarked tetrahedral configurations in 3D coordinates use the shared
    /// chemical eligibility and symmetry model. Explicit unknown assertions
    /// remain unknown; derived perception is not installed on the result.
    /// Valence outside that model is reported as a warning without discarding
    /// the represented graph. Resource exhaustion remains an error.
    pub fn interpret(
        document: &MolfileDocument,
    ) -> Result<MolfileInterpretation, MolfileInterpretError> {
        document.interpret()
    }

    /// Writes one CTAB from a molecule (zero coordinates) or a model.
    ///
    /// [`MolfileWriteVersion::Auto`] uses V2000 when it can represent the
    /// chemistry and V3000 otherwise. V2000 rounds coordinates to four decimal
    /// places in angstroms; V3000 preserves round-trip decimal text. Stereo is
    /// projected against the coordinates emitted by the chosen version, so
    /// specified stereo requires a geometry-bearing model.
    pub fn write<'a>(
        source: impl Into<MolfileSource<'a>>,
        options: MolfileWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_molfile(source, options)
    }

    pub fn write_to<'a>(
        writer: &mut impl std::io::Write,
        source: impl Into<MolfileSource<'a>>,
        options: MolfileWriteOptions,
    ) -> Result<(), MolWriteError> {
        crate::io::write_molfile_to(writer, source, options)
    }
}

/// Record-oriented SDF parsing, interpretation, and writing.
///
/// [`sdf::SdfDocument`] preserves independent records. Interpret or write those
/// records explicitly; sibling records are not merged into one model or
/// reinterpreted as an ensemble.
pub mod sdf {
    pub use crate::io::{
        MolWriteError, MolfileWriteVersion, SdfDataField, SdfDocument, SdfInterpretError,
        SdfInterpretErrorKind, SdfInterpretation, SdfParseError, SdfParseOptions, SdfRecord,
        SdfRecordInterpretation, SdfRecordInterpretationReport, SdfRecordSource, SdfWriteError,
        SdfWriteOptions,
    };

    /// Parses an SDF document while preserving independent record boundaries.
    pub fn parse_str(input: &str) -> Result<SdfDocument, SdfParseError> {
        parse_str_with_options(input, SdfParseOptions::default())
    }

    /// Parses an SDF document with explicit syntax and resource-limit options.
    pub fn parse_str_with_options(
        input: &str,
        options: SdfParseOptions,
    ) -> Result<SdfDocument, SdfParseError> {
        crate::io::parse_sdf_document(input, options)
    }

    /// Interprets every SDF record independently in source order.
    pub fn interpret(document: &SdfDocument) -> Result<SdfInterpretation, SdfInterpretError> {
        document.interpret()
    }

    /// Writes independent records in input order.
    ///
    /// Records convert from models, ensemble members and trajectory frames
    /// (`sdf::write(&ensemble, options)` writes one record per member), or
    /// interpreted records, which keep their titles and data fields.
    pub fn write<'a, R: Into<SdfRecordSource<'a>>>(
        records: impl IntoIterator<Item = R>,
        options: SdfWriteOptions,
    ) -> Result<String, SdfWriteError> {
        crate::io::write_sdf(records, options)
    }

    pub fn write_to<'a, R: Into<SdfRecordSource<'a>>>(
        writer: &mut impl std::io::Write,
        records: impl IntoIterator<Item = R>,
        options: SdfWriteOptions,
    ) -> Result<(), SdfWriteError> {
        crate::io::write_sdf_to(writer, records, options)
    }
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
pub mod mmcif {
    pub use crate::io::{
        MmcifAltLocDecision, MmcifAltLocPolicy, MmcifAltLocPreference, MmcifAltLocResidue,
        MmcifAltLocSelection, MmcifAltLocSelectionReason, MmcifAtomProvenance, MmcifBlock,
        MmcifBlockSource, MmcifConnectionResolutionReason, MmcifDocument,
        MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions, MmcifEnsembleInterpretation,
        MmcifEntityClassifications, MmcifEntityKind, MmcifEntry, MmcifInstanceProvenance,
        MmcifInterpretError, MmcifInterpretIssue, MmcifInterpretOptions, MmcifInterpretation,
        MmcifInterpretationReport, MmcifItem, MmcifLoopTable, MmcifModelSelection, MmcifParseError,
        MmcifParseOptions, MmcifResidueId, MmcifResiduePosition, MmcifValue, MmcifWriteError,
        MmcifWriteOptions,
    };

    /// Parses a structural mmCIF data document without assigning molecular meaning.
    pub fn parse_str(input: &str) -> Result<MmcifDocument, MmcifParseError> {
        parse_str_with_options(input, MmcifParseOptions::default())
    }

    /// Parses CIF/mmCIF source text with explicit syntax and resource limits.
    pub fn parse_str_with_options(
        input: &str,
        options: MmcifParseOptions,
    ) -> Result<MmcifDocument, MmcifParseError> {
        crate::io::parse_mmcif_str(input, options)
    }

    /// Interprets a document containing exactly one atom-site block.
    ///
    /// Documents with multiple atom-site blocks require explicit selection via
    /// [`interpret_block`].
    pub fn interpret(
        document: &MmcifDocument,
        options: MmcifInterpretOptions,
    ) -> Result<MmcifInterpretation, MmcifInterpretError> {
        document.interpret_with_options(options)
    }

    /// Interprets one CIF/mmCIF data block as one selected coordinate model.
    ///
    /// Coordinate-model selection applies only within `block`; sibling blocks
    /// in the source document are independent interpretation scopes.
    pub fn interpret_block(
        block: &MmcifBlock,
        options: MmcifInterpretOptions,
    ) -> Result<MmcifInterpretation, MmcifInterpretError> {
        block.interpret_with_options(options)
    }

    /// Interprets the exactly one atom-site block in a document as an ensemble.
    ///
    /// Documents with multiple atom-site blocks require explicit selection via
    /// [`interpret_ensemble_block`].
    pub fn interpret_ensemble(
        document: &MmcifDocument,
        options: MmcifEnsembleInterpretOptions,
    ) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
        document.interpret_ensemble_with_options(options)
    }

    /// Interprets coordinate models in one block as a shared-topology ensemble.
    pub fn interpret_ensemble_block(
        block: &MmcifBlock,
        options: MmcifEnsembleInterpretOptions,
    ) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
        block.interpret_ensemble_with_options(options)
    }

    /// Builds an unweighted ensemble from explicit coordinate-model and altloc
    /// selections. No alternate configurations are enumerated automatically.
    /// Selections must produce identical atom identities and compatible topology.
    pub fn interpret_conformations(
        document: &MmcifDocument,
        selections: &[MmcifInterpretOptions],
    ) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
        document.interpret_conformations(selections)
    }

    /// Builds an unweighted ensemble from explicit selections within one block.
    pub fn interpret_conformations_block(
        block: &MmcifBlock,
        selections: &[MmcifInterpretOptions],
    ) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
        block.interpret_conformations(selections)
    }

    /// Writes one data block per source.
    ///
    /// Sources convert from models (one block), ensembles (one multi-model
    /// block), and interpretations (keeping their source reports). Entity kinds
    /// derive from canonical topology classification unless a source carries
    /// explicit classifications or reports ([`MmcifBlockSource`]).
    ///
    /// ```no_run
    /// use kekule::mmcif::{self, MmcifBlockSource, MmcifWriteOptions};
    /// # fn example(model: &kekule::structure::Model, other: &kekule::structure::Model,
    /// #     interpretation: &mmcif::MmcifInterpretation) -> Result<(), mmcif::MmcifWriteError> {
    /// let one = mmcif::write([model], MmcifWriteOptions::default())?;
    /// let two = mmcif::write([model, other], MmcifWriteOptions::default())?;
    /// let faithful = mmcif::write([interpretation], MmcifWriteOptions::default())?;
    /// # let _ = (one, two, faithful);
    /// # Ok(())
    /// # }
    /// ```
    pub fn write<'a, B: Into<MmcifBlockSource<'a>>>(
        blocks: impl IntoIterator<Item = B>,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif(blocks, options)
    }

    pub fn write_to<'a, B: Into<MmcifBlockSource<'a>>>(
        writer: &mut impl std::io::Write,
        blocks: impl IntoIterator<Item = B>,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_to(writer, blocks, options)
    }
}

/// Derived valence, rings, aromaticity, conjugation and explicit resonance work.
///
/// [`crate::core::Molecule::perceive`] installs Kekule's default transactional perception
/// profile. The nested modules expose the individual expert algorithms.
/// Perception does not alter represented graph chemistry and is invalidated by
/// relevant graph edits.
pub mod perception {
    /// Connected conjugated groups and explicit RDKit-like contributor enumeration.
    pub mod resonance {
        pub use crate::algorithms::{
            enumerate_resonance, perceive_resonance, ResonanceContributor, ResonanceError,
            ResonanceFlags, ResonanceOptions, ResonanceStructures,
        };
        pub use crate::core::{ResonanceGroup, ResonancePerception};
    }
    /// Bond conjugation, including aromatic bonds, without contributor enumeration.
    pub mod conjugation {
        pub use crate::algorithms::{perceive_conjugation, ConjugationError};
        pub use crate::core::{ConjugationModel, ConjugationPerception};
    }
    pub use crate::chemistry::PerceptionError;

    /// Expert valence perception for canonical represented chemistry.
    ///
    /// The RDKit-like model derives the inferred contribution to implicit H
    /// from ordinary localized bond orders and represented atom state. It does
    /// not require installed ring or aromaticity perception. Source aromatic
    /// bonds are localized during format interpretation before this layer is
    /// reached.
    pub mod valence {
        pub use crate::algorithms::{
            perceive_valence, perceive_valence_with_options, represented_valence, ValenceError,
            ValenceIssue, ValenceOptions,
        };
        pub use crate::core::ValenceModel;
    }

    /// Ring membership and deterministic ring-basis perception.
    ///
    /// Ring perception installs derived state and never changes represented
    /// bond orders.
    pub mod rings {
        pub use crate::algorithms::{
            perceive_ring_membership, perceive_ring_set, perceive_ring_set_with_options,
            RingPerceptionError, RingPerceptionOptions,
        };
        pub use crate::core::{Ring, RingBasisModel, RingBasisState, RingMembership, RingSet};
    }

    /// Aromaticity perception over canonical localized bond orders.
    ///
    /// Aromatic membership is derived state; Kekule does not store an aromatic
    /// bond order in represented graph chemistry.
    pub mod aromaticity {
        pub use crate::algorithms::{
            perceive_aromaticity, perceive_aromaticity_with_options,
            perceive_aromaticity_with_ring_options, AromaticityError, AromaticityOptions,
        };
        pub use crate::core::AromaticityModel;
    }
}

/// Focused stereochemistry validation, inference, transforms, and CIP assignment.
///
/// Coordinate inference and candidate detection are read-only. Coordinate
/// materialization is an explicitly named representation transform, while CIP
/// descriptors remain opt-in derived state.
pub mod stereo {
    pub use crate::algorithms::{
        assign_cip_descriptors, assign_cip_descriptors_with_options, cleanup_stereo,
        detect_stereo_candidates, detect_stereo_candidates_with_options, infer_coordinate_stereo,
        infer_coordinate_stereo_with_options, materialize_coordinate_stereo,
        materialize_coordinate_stereo_with_options, validate_stereo, CipAssignment,
        CipAssignmentError, CipAssignmentIssue, CipAssignmentOptions, CipAssignmentReport,
        CipRankingError, CipSkipped, CipSkippedReason, CoordinateStereoError,
        CoordinateStereoMaterializationReport, CoordinateStereoOptions, CoordinateStereoResult,
        StereoCandidate, StereoCleanupReport, StereoPerceptionError, StereoPerceptionOptions,
        StereoValidationError, StereoValidationIssue,
    };
}

/// Canonical atom ranking for represented molecular graphs.
///
/// Ranking is a graph-derived result and does not reorder or mutate a molecule.
pub mod canon {
    pub use crate::algorithms::CanonicalAtomRanking;

    use crate::core::Molecule;

    /// Computes deterministic canonical equivalence classes without mutation.
    ///
    /// Hydrogen comparison uses declared counts plus installed implicit counts;
    /// equal totals are equivalent regardless of how the counts are stored.
    /// Run valence perception first to include inferred hydrogens. Without it,
    /// only declared hydrogens contribute. Source declarations remain unchanged.
    pub fn atom_ranking(molecule: &Molecule) -> CanonicalAtomRanking {
        crate::algorithms::canonical_atom_ranking(molecule)
    }
}

/// Read-only rotatable-bond detection.
///
/// This facade identifies configurable represented single-bond axes in
/// canonical chemistry. It does not mutate the molecule or install perception
/// state.
pub mod rotatable_bonds {
    pub use crate::algorithms::{RotatableBondOptions, RotatableBondSet};
    pub use crate::chemistry::PerceptionError;

    use crate::core::Molecule;

    /// Detects rotatable bonds using the supplied options.
    ///
    /// Resonance exclusions reuse installed RDKit-like valence and aromaticity
    /// with the default Figueras ring model, or compute the default profile on
    /// a temporary copy. Strict valence validation applies in either case.
    /// Model-neutral, incomplete, or incompatible perception is not reused.
    /// Valence and resource failures are returned without changing the source
    /// graph or its perception. Disabling resonance exclusions requires only
    /// graph ring membership and does not run chemical perception.
    pub fn detect(
        molecule: &Molecule,
        options: RotatableBondOptions,
    ) -> Result<RotatableBondSet, PerceptionError> {
        crate::algorithms::detect_rotatable_bonds(molecule, options)
    }
}

/// Explicit small-molecule hydrogen topology transforms.
///
/// Explicit hydrogens are graph atoms; implicit hydrogens are represented by
/// specified counts and/or valence inference. Addition converts implicit H to
/// explicit atoms; removal suppresses eligible explicit atoms without losing
/// chemical information. These operations do not change protonation states.
///
/// Addition requires current inference for inference-enabled atoms, unless
/// `specified_only` is selected. Fixed counts need no perception. Removal needs
/// current inference on affected parents when their counts are not fixed.
/// Successful topology changes invalidate perception state; recompute it before
/// reading inferred counts. Removal verifies its plan on a temporary copy.
pub mod hydrogens {
    pub use crate::algorithms::{
        AddHydrogensOptions, AddHydrogensReport, AddedHydrogen, AddedHydrogenOrigin,
        HydrogenCountAdjustment, HydrogenTransformError, RemoveHydrogensReport, RemovedHydrogen,
        RetainedHydrogen, RetainedHydrogenReason,
    };
}

/// Common foundational types for small examples and interactive use.
///
/// The prelude is intentionally small. Format APIs, structure types, and
/// algorithms remain available through their focused modules.
pub mod prelude {
    pub use crate::core::{
        Atom, AtomId, Bond, BondId, BondOrder, Element, HydrogenDeclaration, Molecule,
    };
    pub use crate::topology::Hierarchy;
}

#[cfg(test)]
mod tests;
