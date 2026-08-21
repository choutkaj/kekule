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
//! - [`structure::Model`] = one [`topology::Topology`] plus topology-bound
//!   positions, an optional periodic cell, and per-atom and per-bond data;
//! - [`structure::Ensemble`] = one topology plus finite non-temporal members;
//!
//! Ordered multi-frame trajectories and their file codecs live in the
//! one-way `kekule-traj` companion crate.
//!
//! Coordinate-dependent kernels consume borrowed [`structure::ModelView`]
//! values, allowing the same analysis or prepared potential to operate over a
//! model, ensemble member, trajectory frame, or reusable frame buffer without
//! copying coordinates. Positions, selections, buffers, and prepared systems
//! require the same shared `Arc<Topology>` allocation.
//! [`topology::Topology::same_layout`]
//! compares complete static layout, including semantic IDs and dense order;
//! general order-independent structural equivalence and isomorphism mapping
//! remain future capabilities.
//!
//! Collection-backed public identifiers are fixed-width and every insertion
//! that creates one is fallible. Exhausting an atom, bond, conformer, stereo,
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
pub mod bio;
mod chemistry;
pub mod core;
pub mod descriptors;
pub mod dssp;
pub mod geometry;
mod io;
pub mod modeling;
pub mod query;
pub mod small;
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
    pub use crate::io::{
        MolWriteError, SmilesAtomMapping, SmilesBondMapping, SmilesDocument, SmilesDocumentToken,
        SmilesDocumentTokenKind, SmilesInterpretError, SmilesInterpretation,
        SmilesInterpretationReport, SmilesParseError, SmilesParseOptions,
    };
    use crate::small::SmallMolecule;

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
        crate::io::interpret_smiles_document(document)
    }

    pub fn write(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_smiles(molecule)
    }

    pub fn write_isomeric(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_isomeric_smiles(molecule)
    }

    pub fn write_canonical(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_canonical_smiles(molecule)
    }
}

pub mod molfile {
    pub use crate::io::{
        MolWriteError, MolfileAtomMapping, MolfileBondMapping, MolfileDocument, MolfileHeader,
        MolfileInterpretError, MolfileInterpretation, MolfileInterpretationReport,
        MolfileInterpretationWarning, MolfileLine, MolfileParseError, MolfileParseOptions,
        MolfileVersion,
    };

    use crate::small::SmallMolecule;

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
        crate::io::interpret_molfile_document(document)
    }

    pub fn write_v2000(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_mol_v2000(molecule)
    }

    pub fn write_v3000(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_mol_v3000(molecule)
    }
}

pub mod sdf {
    pub use crate::io::{
        MolWriteError, SdfDataField, SdfDocument, SdfInterpretError, SdfInterpretation,
        SdfInterpretationReport, SdfParseError, SdfParseOptions, SdfRecord, SdfRecordDocument,
        SdfRecordInterpretationReport,
    };

    pub fn parse_str(input: &str, options: SdfParseOptions) -> Result<SdfDocument, SdfParseError> {
        crate::io::parse_sdf_document(input, options)
    }

    pub fn interpret(document: &SdfDocument) -> Result<SdfInterpretation, SdfInterpretError> {
        crate::io::interpret_sdf_document(document)
    }

    pub fn write_v2000(records: &[SdfRecord]) -> Result<String, MolWriteError> {
        crate::io::write_sdf_v2000(records)
    }
}

pub mod mmcif {
    pub use crate::io::{
        MmcifAltLocPolicy, MmcifAtomProvenance, MmcifConnectionResolutionReason, MmcifDataBlock,
        MmcifDocument, MmcifEnsembleInterpretError, MmcifEnsembleInterpretOptions,
        MmcifEnsembleInterpretation, MmcifEntityKind, MmcifEntry, MmcifInstanceProvenance,
        MmcifInterpretError, MmcifInterpretIssue, MmcifInterpretOptions, MmcifInterpretation,
        MmcifInterpretationReport, MmcifItem, MmcifLoopTable, MmcifModelSelection, MmcifParseError,
        MmcifParseOptions, MmcifValue, MmcifWriteError, MmcifWriteOptions,
    };

    /// Parses a structural mmCIF data document without assigning molecular meaning.
    pub fn parse_str(
        input: &str,
        options: MmcifParseOptions,
    ) -> Result<MmcifDocument, MmcifParseError> {
        crate::io::parse_mmcif_str(input, options)
    }

    /// Interprets one coordinate-containing data block as clean molecular objects.
    pub fn interpret(
        document: &MmcifDocument,
        options: MmcifInterpretOptions,
    ) -> Result<MmcifInterpretation, MmcifInterpretError> {
        crate::io::interpret_mmcif(document, options)
    }

    /// Interprets multiple coordinate models as one verified shared-topology ensemble.
    pub fn interpret_ensemble(
        document: &MmcifDocument,
        options: MmcifEnsembleInterpretOptions,
    ) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
        crate::io::interpret_mmcif_ensemble(document, options)
    }

    /// Writes one canonical molecular model as a structural mmCIF data block.
    pub fn write(
        model: &crate::structure::Model,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_model(model, options)
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
    use crate::small::SmallMolecule;

    /// Materialize stored explicit counts and perceived implicit hydrogens.
    pub fn add_hydrogens(
        molecule: &mut SmallMolecule,
    ) -> Result<AddHydrogensReport, HydrogenTransformError> {
        add_hydrogens_with_options(molecule, AddHydrogensOptions::default())
    }

    /// Materialize hydrogens with an explicit growth bound and count policy.
    pub fn add_hydrogens_with_options(
        molecule: &mut SmallMolecule,
        options: AddHydrogensOptions,
    ) -> Result<AddHydrogensReport, HydrogenTransformError> {
        add_hydrogens_to_molecule(molecule.as_molecule_mut(), options)
    }

    /// Collapse ordinary degree-one hydrogens without discarding protected state.
    ///
    /// Isotopic, mapped, charged, radical, property-bearing, and otherwise
    /// non-losslessly representable hydrogens remain in the graph and are
    /// described by the returned report.
    pub fn remove_hydrogens(
        molecule: &mut SmallMolecule,
    ) -> Result<RemoveHydrogensReport, HydrogenTransformError> {
        remove_hydrogens_from_molecule(molecule.as_molecule_mut())
    }
}

pub mod prelude {
    pub use crate::bio::{MacroMolecule, SmcraHierarchy};
    pub use crate::core::{
        Atom, AtomId, Bond, BondId, BondOrder, Conformer, Element, HydrogenDeclaration, Molecule,
    };
    pub use crate::small::SmallMolecule;
}

#[cfg(test)]
mod tests;
