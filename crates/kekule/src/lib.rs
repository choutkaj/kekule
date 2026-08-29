//! Pure-Rust molecular graph, structure I/O, perception, bioinformatics, and
//! modelling foundations.
//!
//! The supported public surface is organized into focused facade modules.
//! Parsing produces format documents, interpretation publishes canonical
//! represented domain objects plus reports, and perception or modelling
//! preparation remains explicit.
//!
//! System structure is coordinate-free and immutable in [`topology`].
//! Dynamic structure follows three explicit relationships:
//!
//! - [`structure::Model`] = one [`topology::Topology`] plus dense positions,
//!   an optional periodic cell, and unified per-atom and per-bond properties;
//! - [`structure::Ensemble`] = one topology plus finite non-temporal members;
//!
//! Ordered multi-frame trajectories and their file codecs live in the
//! one-way `kekule-traj` companion crate.
//!
//! Coordinate-dependent kernels consume borrowed [`structure::ModelView`]
//! values, allowing the same analysis or prepared potential to operate over a
//! model, ensemble member, trajectory frame, or reusable frame buffer without
//! copying coordinates. Models, selections, buffers, and prepared systems may
//! require the same shared `Arc<Topology>` allocation; primitive dense arrays
//! carry no topology identity.
//! [`topology::Topology::same_layout`]
//! compares complete static layout, including semantic IDs and dense order;
//! general order-independent structural equivalence and isomorphism mapping
//! remain future capabilities.
//!
//! Collection-backed public identifiers are fixed-width and every insertion
//! that creates one is fallible. Exhausting an atom, bond, stereo,
//! hierarchy, definition, instance, or dense topology index space returns the
//! corresponding structured capacity error before canonical state is changed.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]
// Kekule consistently names owned conversions `to_*`, including consuming ones.
#![allow(clippy::wrong_self_convention)]

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
pub mod modeling;
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
        find_substructure_match, find_substructure_matches, find_substructure_matches_with_options,
        QueryMatch, QueryPerception, SubstructureMatchError, SubstructureMatchOptions,
        SubstructureMatchWork, MAX_SUBSTRUCTURE_QUERY_ATOMS,
    };
}

pub mod smiles {
    use std::fmt;

    use crate::core::Molecule;
    pub use crate::io::{
        MolWriteError, SmilesAtomMapping, SmilesBondMapping, SmilesComponentCountError,
        SmilesComponentInterpretation, SmilesDocument, SmilesDocumentToken,
        SmilesDocumentTokenKind, SmilesInterpretError, SmilesInterpretation,
        SmilesInterpretationReport, SmilesParseError, SmilesParseOptions,
    };
    use crate::topology::{Topology, TopologyBuildError};

    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    #[non_exhaustive]
    pub enum SmilesWriteMode {
        #[default]
        Ordinary,
        Isomeric,
        Canonical,
    }

    #[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
    pub struct SmilesWriteOptions {
        pub mode: SmilesWriteMode,
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

    pub fn parse_str(input: &str) -> Result<SmilesDocument, SmilesParseError> {
        crate::io::parse_smiles_document(input)
    }

    pub fn parse_str_with_options(
        input: &str,
        options: SmilesParseOptions,
    ) -> Result<SmilesDocument, SmilesParseError> {
        crate::io::parse_smiles_document_with_options(input, options)
    }

    pub fn interpret(
        document: &SmilesDocument,
    ) -> Result<SmilesInterpretation, SmilesInterpretError> {
        document.interpret()
    }

    /// Parses and interprets one SMILES record into source-ordered connected components.
    ///
    /// This is the concise form of [`parse_str`] followed by [`interpret`].
    /// Dot-delimited components remain separate, and no perception is run implicitly.
    pub fn to_molecules(input: &str) -> Result<Vec<Molecule>, SmilesReadError> {
        let document = parse_str(input)?;
        Ok(document.interpret()?.to_molecules())
    }

    /// Parses and interprets one SMILES record as a coordinate-free topology.
    ///
    /// Every dot-delimited connected component becomes one explicit molecule
    /// occurrence in source order. No hierarchy or perception is fabricated.
    pub fn to_topology(input: &str) -> Result<Topology, SmilesReadError> {
        let document = parse_str(input)?;
        Ok(document.interpret()?.to_topology()?)
    }

    pub fn write(molecule: &Molecule) -> Result<String, MolWriteError> {
        crate::io::write_smiles(molecule)
    }

    /// Writes one connected molecule using an explicit SMILES policy.
    pub fn write_molecule(
        molecule: &Molecule,
        options: SmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        match options.mode {
            SmilesWriteMode::Ordinary => crate::io::write_smiles(molecule),
            SmilesWriteMode::Isomeric => crate::io::write_isomeric_smiles(molecule),
            SmilesWriteMode::Canonical => crate::io::write_canonical_smiles(molecule),
        }
    }

    /// Writes every explicit topology molecule instance in authoritative order.
    ///
    /// Reused definitions are emitted once per occurrence and joined with `.`.
    pub fn write_topology(
        topology: &Topology,
        options: SmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        topology
            .molecules()
            .map(|occurrence| write_molecule(occurrence.molecule(), options))
            .collect::<Result<Vec<_>, _>>()
            .map(|components| components.join("."))
    }

    pub fn write_isomeric(molecule: &Molecule) -> Result<String, MolWriteError> {
        crate::io::write_isomeric_smiles(molecule)
    }

    pub fn write_canonical(molecule: &Molecule) -> Result<String, MolWriteError> {
        crate::io::write_canonical_smiles(molecule)
    }
}

pub mod molfile {
    pub use crate::io::{
        MolWriteError, MolWriteErrorKind, MolfileAtomMapping, MolfileBondMapping,
        MolfileComponentInterpretation, MolfileDocument, MolfileHeader, MolfileInterpretError,
        MolfileInterpretation, MolfileInterpretationReport, MolfileInterpretationWarning,
        MolfileLine, MolfileModelError, MolfileParseError, MolfileParseOptions, MolfileVersion,
        MolfileWriteOptions, MolfileWriteVersion,
    };

    use crate::core::Molecule;
    use crate::structure::Model;

    pub fn parse_str(input: &str) -> Result<MolfileDocument, MolfileParseError> {
        crate::io::parse_molfile_document(input)
    }

    pub fn parse_str_with_options(
        input: &str,
        options: MolfileParseOptions,
    ) -> Result<MolfileDocument, MolfileParseError> {
        crate::io::parse_molfile_document_with_options(input, options)
    }

    pub fn interpret(
        document: &MolfileDocument,
    ) -> Result<MolfileInterpretation, MolfileInterpretError> {
        document.interpret()
    }

    pub fn write_v2000(molecule: &Molecule) -> Result<String, MolWriteError> {
        crate::io::write_mol_v2000(molecule)
    }

    pub fn write_v3000(molecule: &Molecule) -> Result<String, MolWriteError> {
        crate::io::write_mol_v3000(molecule)
    }

    /// Writes a coordinate-free molecule using zero coordinates.
    pub fn write_molecule(
        molecule: &Molecule,
        options: MolfileWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_molfile_molecule(molecule, options)
    }

    /// Writes one geometry-bearing model as one possibly disconnected CTAB.
    pub fn write_model(
        model: &Model,
        options: MolfileWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_molfile_model(model.view(), options)
    }

    pub fn write_model_to(
        writer: &mut impl std::io::Write,
        model: &Model,
        options: MolfileWriteOptions,
    ) -> Result<(), MolWriteError> {
        crate::io::write_molfile_model_to(writer, model.view(), options)
    }

    pub fn write_model_v2000(model: &Model) -> Result<String, MolWriteError> {
        crate::io::write_model_v2000(model.view())
    }

    pub fn write_model_v3000(model: &Model) -> Result<String, MolWriteError> {
        crate::io::write_model_v3000(model.view())
    }
}

pub mod sdf {
    pub use crate::io::{
        MolWriteError, MolfileWriteVersion, SdfDataField, SdfDocument, SdfInterpretError,
        SdfInterpretErrorKind, SdfInterpretation, SdfInterpretationReport, SdfParseError,
        SdfParseOptions, SdfRecord, SdfRecordInterpretation, SdfRecordInterpretationReport,
        SdfWriteError, SdfWriteOptions,
    };

    use crate::structure::{Ensemble, Model};

    pub fn parse_str(input: &str) -> Result<SdfDocument, SdfParseError> {
        parse_str_with_options(input, SdfParseOptions::default())
    }

    pub fn parse_str_with_options(
        input: &str,
        options: SdfParseOptions,
    ) -> Result<SdfDocument, SdfParseError> {
        crate::io::parse_sdf_document(input, options)
    }

    pub fn interpret(document: &SdfDocument) -> Result<SdfInterpretation, SdfInterpretError> {
        document.interpret()
    }

    pub fn write_v2000(records: &[SdfRecordInterpretation]) -> Result<String, MolWriteError> {
        crate::io::write_sdf_v2000(records)
    }

    pub fn write_model(model: &Model, options: SdfWriteOptions) -> Result<String, SdfWriteError> {
        crate::io::write_sdf_model(model, options)
    }

    pub fn write_models(
        models: &[Model],
        options: SdfWriteOptions,
    ) -> Result<String, SdfWriteError> {
        crate::io::write_sdf_models(models, options)
    }

    pub fn write_ensemble(
        ensemble: &Ensemble,
        options: SdfWriteOptions,
    ) -> Result<String, SdfWriteError> {
        crate::io::write_sdf_ensemble(ensemble, options)
    }

    pub fn write_model_to(
        writer: &mut impl std::io::Write,
        model: &Model,
        options: SdfWriteOptions,
    ) -> Result<(), SdfWriteError> {
        crate::io::write_sdf_model_to(writer, model, options)
    }

    pub fn write_models_to(
        writer: &mut impl std::io::Write,
        models: &[Model],
        options: SdfWriteOptions,
    ) -> Result<(), SdfWriteError> {
        crate::io::write_sdf_models_to(writer, models, options)
    }

    pub fn write_ensemble_to(
        writer: &mut impl std::io::Write,
        ensemble: &Ensemble,
        options: SdfWriteOptions,
    ) -> Result<(), SdfWriteError> {
        crate::io::write_sdf_ensemble_to(writer, ensemble, options)
    }

    /// Expert/round-trip path preserving explicit SDF titles and data fields.
    pub fn write_records(
        records: &[SdfRecordInterpretation],
        options: SdfWriteOptions,
    ) -> Result<String, SdfWriteError> {
        crate::io::write_sdf_records(records, options)
    }

    pub fn write_records_to(
        writer: &mut impl std::io::Write,
        records: &[SdfRecordInterpretation],
        options: SdfWriteOptions,
    ) -> Result<(), SdfWriteError> {
        crate::io::write_sdf_records_to(writer, records, options)
    }
}

pub mod mmcif {
    pub use crate::io::{
        MmcifAltLocPolicy, MmcifAtomProvenance, MmcifBlock, MmcifConnectionResolutionReason,
        MmcifDocument, MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions,
        MmcifEnsembleInterpretation, MmcifEntityClassifications, MmcifEntityKind, MmcifEntry,
        MmcifInstanceProvenance, MmcifInterpretError, MmcifInterpretIssue, MmcifInterpretOptions,
        MmcifInterpretation, MmcifInterpretationReport, MmcifItem, MmcifLoopTable,
        MmcifModelSelection, MmcifParseError, MmcifParseOptions, MmcifValue, MmcifWriteError,
        MmcifWriteOptions,
    };

    /// Parses a structural mmCIF data document without assigning molecular meaning.
    pub fn parse_str(input: &str) -> Result<MmcifDocument, MmcifParseError> {
        parse_str_with_options(input, MmcifParseOptions::default())
    }

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

    /// Attempts to write a model without inventing mmCIF entity semantics.
    ///
    /// Generic topology state intentionally contains no mmCIF entity roles, so
    /// a non-empty model returns [`MmcifWriteError::MissingEntityClassification`].
    /// Use [`write_with_classifications`] for a generic model or
    /// [`write_with_report`] for an interpreted mmCIF model.
    pub fn write(
        model: &crate::structure::Model,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_model(model, options)
    }

    /// Attempts to write one model as one block without inventing entity semantics.
    ///
    /// This is the canonical-object spelling of [`write()`]. A generic non-empty
    /// model normally requires [`write_model_with_classifications`]; an
    /// mmCIF-derived model should use [`write_model_with_report`] to preserve
    /// source entity and asymmetry semantics.
    pub fn write_model(
        model: &crate::structure::Model,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        write(model, options)
    }

    /// Attempts to write independent models as one deterministic block per model.
    ///
    /// This strict entry point does not infer polymer, branched, non-polymer,
    /// or water classifications. Generic non-empty models normally require
    /// [`write_models_with_classifications`]; mmCIF-derived models should use
    /// [`write_models_with_reports`] to preserve source format semantics.
    pub fn write_models(
        models: &[crate::structure::Model],
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_models(models, options)
    }

    pub fn write_models_with_classifications(
        models: &[crate::structure::Model],
        classifications: &[MmcifEntityClassifications],
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_models_with_classifications(models, classifications, options)
    }

    pub fn write_models_with_reports(
        models: &[crate::structure::Model],
        reports: &[MmcifInterpretationReport],
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_models_with_reports(models, reports, options)
    }

    /// Attempts to write one shared-topology ensemble as one multi-model block.
    ///
    /// This strict entry point does not infer polymer, branched, non-polymer,
    /// or water classifications. A generic non-empty ensemble normally requires
    /// [`write_ensemble_with_classifications`]. For mmCIF-derived state, use
    /// [`write_ensemble_with_reports`] or [`write_ensemble_interpretation`] to
    /// preserve source format semantics.
    pub fn write_ensemble(
        ensemble: &crate::structure::Ensemble,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_ensemble(ensemble, options)
    }

    pub fn write_ensemble_with_classifications(
        ensemble: &crate::structure::Ensemble,
        classifications: &MmcifEntityClassifications,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_ensemble_with_classifications(ensemble, classifications, options)
    }

    pub fn write_ensemble_with_reports(
        ensemble: &crate::structure::Ensemble,
        reports: &[MmcifInterpretationReport],
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_ensemble_with_reports(ensemble, reports, options)
    }

    pub fn write_ensemble_interpretation(
        interpretation: &MmcifEnsembleInterpretation,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_ensemble_interpretation(interpretation, options)
    }

    pub fn write_models_to(
        writer: &mut impl std::io::Write,
        models: &[crate::structure::Model],
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_models_to(writer, models, options)
    }

    pub fn write_ensemble_to(
        writer: &mut impl std::io::Write,
        ensemble: &crate::structure::Ensemble,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_ensemble_to(writer, ensemble, options)
    }

    /// Writes a generic model with explicit format-specific entity semantics.
    ///
    /// Every molecule instance must occur exactly once in `classifications`.
    /// Generic writing deterministically assigns one mmCIF entity to each
    /// populated topology hierarchy chain (and one to each hierarchy-free
    /// instance). Instances touched by the same chain must therefore have the
    /// same explicit classification.
    pub fn write_with_classifications(
        model: &crate::structure::Model,
        classifications: &MmcifEntityClassifications,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_model_with_classifications(model, classifications, options)
    }

    pub fn write_model_with_classifications(
        model: &crate::structure::Model,
        classifications: &MmcifEntityClassifications,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        write_with_classifications(model, classifications, options)
    }

    /// Writes a canonical model while preserving source mmCIF entity/asymmetry semantics.
    ///
    /// Atom-level report provenance keeps one source entity and structural
    /// asymmetry consistent even when it spans multiple connected-component
    /// molecule instances. Conflicting source identity is rejected. Because the
    /// emitted foundational atom-site loop requires label identifiers, an
    /// auth-only source is deterministically normalized by copying its author
    /// atom, component, and asymmetry identifiers into the corresponding label
    /// output fields; the original report still records that those label fields
    /// were absent.
    pub fn write_with_report(
        model: &crate::structure::Model,
        report: &MmcifInterpretationReport,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_model_with_report(model, report, options)
    }

    pub fn write_model_with_report(
        model: &crate::structure::Model,
        report: &MmcifInterpretationReport,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        write_with_report(model, report, options)
    }

    pub fn write_model_to(
        writer: &mut impl std::io::Write,
        model: &crate::structure::Model,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_model_to(writer, model, options)
    }

    pub fn write_model_with_classifications_to(
        writer: &mut impl std::io::Write,
        model: &crate::structure::Model,
        classifications: &MmcifEntityClassifications,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_model_with_classifications_to(
            writer,
            model,
            classifications,
            options,
        )
    }

    pub fn write_model_with_report_to(
        writer: &mut impl std::io::Write,
        model: &crate::structure::Model,
        report: &MmcifInterpretationReport,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_model_with_report_to(writer, model, report, options)
    }

    pub fn write_models_with_classifications_to(
        writer: &mut impl std::io::Write,
        models: &[crate::structure::Model],
        classifications: &[MmcifEntityClassifications],
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_models_with_classifications_to(
            writer,
            models,
            classifications,
            options,
        )
    }

    pub fn write_models_with_reports_to(
        writer: &mut impl std::io::Write,
        models: &[crate::structure::Model],
        reports: &[MmcifInterpretationReport],
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_models_with_reports_to(writer, models, reports, options)
    }

    pub fn write_ensemble_with_classifications_to(
        writer: &mut impl std::io::Write,
        ensemble: &crate::structure::Ensemble,
        classifications: &MmcifEntityClassifications,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_ensemble_with_classifications_to(
            writer,
            ensemble,
            classifications,
            options,
        )
    }

    pub fn write_ensemble_with_reports_to(
        writer: &mut impl std::io::Write,
        ensemble: &crate::structure::Ensemble,
        reports: &[MmcifInterpretationReport],
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_ensemble_with_reports_to(writer, ensemble, reports, options)
    }

    pub fn write_ensemble_interpretation_to(
        writer: &mut impl std::io::Write,
        interpretation: &MmcifEnsembleInterpretation,
        options: MmcifWriteOptions,
    ) -> Result<(), MmcifWriteError> {
        crate::io::write_mmcif_ensemble_interpretation_to(writer, interpretation, options)
    }
}

pub mod perception {
    pub use crate::chemistry::PerceptionError;

    use crate::core::Molecule;

    /// Expert valence perception for canonical represented chemistry.
    ///
    /// The RDKit-like model derives complete implicit-hydrogen assignments
    /// from ordinary localized bond orders and represented atom state. It does
    /// not require installed ring or aromaticity perception. Source aromatic
    /// bonds are localized during format interpretation before this layer is
    /// reached.
    pub mod valence {
        pub use crate::algorithms::{
            perceive_valence, perceive_valence_with_options, ValenceError, ValenceIssue,
            ValenceOptions,
        };
        pub use crate::core::ValenceModel;
    }

    pub mod rings {
        pub use crate::algorithms::{
            perceive_ring_membership, perceive_ring_set, perceive_ring_set_with_options,
            RingPerceptionError, RingPerceptionOptions,
        };
        pub use crate::core::{Ring, RingBasisModel, RingBasisState, RingMembership, RingSet};
    }

    pub mod aromaticity {
        pub use crate::algorithms::{
            perceive_aromaticity, perceive_aromaticity_with_ring_options, AromaticityError,
        };
        pub use crate::core::AromaticityModel;
    }

    /// Install the transactional default valence, ring-set, and aromaticity profile.
    ///
    /// The represented molecule must already be canonical. This operation
    /// never rewrites represented chemistry or performs stereo or CIP work.
    pub fn perceive(molecule: &mut Molecule) -> Result<(), PerceptionError> {
        crate::chemistry::perceive_molecule(molecule)
    }
}

/// Focused stereochemistry validation, inference, transforms, and CIP assignment.
///
/// Coordinate inference and candidate detection are read-only. Coordinate
/// materialization is an explicitly named representation transform, while CIP
/// descriptors remain opt-in derived state.
pub mod stereo {
    pub use crate::algorithms::{
        assign_cip_descriptors, assign_cip_descriptors_with_options, detect_stereo_candidates,
        infer_coordinate_stereo, infer_coordinate_stereo_with_options,
        materialize_coordinate_stereo, materialize_coordinate_stereo_with_options, validate_stereo,
        CipAssignment, CipAssignmentError, CipAssignmentIssue, CipAssignmentOptions,
        CipAssignmentReport, CipSkipped, CipSkippedReason, CoordinateStereoError,
        CoordinateStereoMaterializationReport, CoordinateStereoOptions, CoordinateStereoResult,
        StereoCandidate, StereoValidationError, StereoValidationIssue,
    };
}

pub mod canon {
    pub use crate::algorithms::CanonicalAtomRanking;

    use crate::core::Molecule;

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

    use crate::core::Molecule;

    /// Detects rotatable bonds using the supplied options.
    pub fn detect(molecule: &Molecule, options: RotatableBondOptions) -> RotatableBondSet {
        crate::algorithms::detect_rotatable_bonds(molecule, options)
    }
}

/// Explicit small-molecule hydrogen topology transforms.
///
/// These functions never interpret or perceive chemistry implicitly. Addition
/// consumes current valence assignments unless `explicit_only` is selected,
/// and removal requires current valence assignments. Successful topology
/// changes invalidate perception state.
pub mod hydrogens {
    pub use crate::algorithms::{
        AddHydrogensOptions, AddHydrogensReport, AddedHydrogen, AddedHydrogenOrigin,
        HydrogenCountAdjustment, HydrogenTransformError, RemoveHydrogensReport, RemovedHydrogen,
        RetainedHydrogen, RetainedHydrogenReason,
    };

    use crate::algorithms::{add_hydrogens_to_molecule, remove_hydrogens_from_molecule};
    use crate::core::Molecule;

    /// Materialize stored explicit counts and perceived implicit hydrogens.
    pub fn add_hydrogens(
        molecule: &mut Molecule,
    ) -> Result<AddHydrogensReport, HydrogenTransformError> {
        add_hydrogens_with_options(molecule, AddHydrogensOptions::default())
    }

    /// Materialize hydrogens with an explicit growth bound and count policy.
    pub fn add_hydrogens_with_options(
        molecule: &mut Molecule,
        options: AddHydrogensOptions,
    ) -> Result<AddHydrogensReport, HydrogenTransformError> {
        add_hydrogens_to_molecule(molecule, options)
    }

    /// Collapse ordinary degree-one hydrogens without discarding protected state.
    ///
    /// Isotopic, mapped, charged, radical, property-bearing, and otherwise
    /// non-losslessly representable hydrogens remain in the graph and are
    /// described by the returned report.
    pub fn remove_hydrogens(
        molecule: &mut Molecule,
    ) -> Result<RemoveHydrogensReport, HydrogenTransformError> {
        remove_hydrogens_from_molecule(molecule)
    }
}

pub mod prelude {
    pub use crate::core::{
        Atom, AtomId, Bond, BondId, BondOrder, Element, HydrogenDeclaration, Molecule,
    };
    pub use crate::topology::Hierarchy;
}

#[cfg(test)]
mod tests;
