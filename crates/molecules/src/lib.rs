#![forbid(unsafe_code)]

mod algorithms;
pub mod bio;
mod chemistry;
pub mod core;
pub mod dssp;
mod io;
pub mod modeling;
pub mod query;
pub mod small;
pub mod units;

/// Syntax-independent substructure matching algorithms.
///
/// Matching consumes `query::QueryGraph` and current target perception state;
/// it never invokes parsing, sanitization, or perception implicitly.
pub mod substructure {
    pub use crate::algorithms::{
        find_substructure_match, find_substructure_matches, find_substructure_matches_with_options,
        QueryMatch, QueryPerception, SubstructureMatchError, SubstructureMatchOptions,
        SubstructureMatchWork, MAX_SUBSTRUCTURE_QUERY_ATOMS,
    };
}

pub mod smiles {
    pub use crate::io::{
        CanonicalSmilesWriteOptions, IsomericSmilesWriteOptions, MolWriteError, SmilesDocument,
        SmilesDocumentToken, SmilesDocumentTokenKind, SmilesInterpretError, SmilesParseError,
        SmilesWriteOptions,
    };
    use crate::small::SmallMolecule;

    pub fn parse_str(input: &str) -> Result<SmilesDocument, SmilesParseError> {
        crate::io::parse_smiles_document(input)
    }

    pub fn interpret(document: &SmilesDocument) -> Result<SmallMolecule, SmilesInterpretError> {
        crate::io::interpret_smiles_document(document)
    }

    pub fn write(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_smiles(molecule, SmilesWriteOptions)
    }

    pub fn write_with_options(
        molecule: &SmallMolecule,
        options: SmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_smiles(molecule, options)
    }

    pub fn write_isomeric(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_isomeric_smiles(molecule, IsomericSmilesWriteOptions)
    }

    pub fn write_isomeric_with_options(
        molecule: &SmallMolecule,
        options: IsomericSmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_isomeric_smiles(molecule, options)
    }

    pub fn write_canonical(molecule: &SmallMolecule) -> Result<String, MolWriteError> {
        crate::io::write_canonical_smiles(molecule, CanonicalSmilesWriteOptions)
    }

    pub fn write_canonical_with_options(
        molecule: &SmallMolecule,
        options: CanonicalSmilesWriteOptions,
    ) -> Result<String, MolWriteError> {
        crate::io::write_canonical_smiles(molecule, options)
    }
}

pub mod molfile {
    pub use crate::io::{
        MolWriteError, MolfileDocument, MolfileHeader, MolfileInterpretError, MolfileLine,
        MolfileParseError, MolfileVersion,
    };

    use crate::small::SmallMolecule;

    pub fn parse_str(input: &str) -> Result<MolfileDocument, MolfileParseError> {
        crate::io::parse_molfile_document(input)
    }

    pub fn interpret(document: &MolfileDocument) -> Result<SmallMolecule, MolfileInterpretError> {
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
        MolWriteError, SdfDataField, SdfDocument, SdfInterpretError, SdfParseError,
        SdfParseOptions, SdfRecord, SdfRecordDocument,
    };

    pub fn parse_str(input: &str, options: SdfParseOptions) -> Result<SdfDocument, SdfParseError> {
        crate::io::parse_sdf_document(input, options)
    }

    pub fn interpret(document: &SdfDocument) -> Result<Vec<SdfRecord>, SdfInterpretError> {
        crate::io::interpret_sdf_document(document)
    }

    pub fn write_v2000(records: &[SdfRecord]) -> Result<String, MolWriteError> {
        crate::io::write_sdf_v2000(records)
    }
}

pub mod mmcif {
    pub use crate::io::{
        MmcifAltLocPolicy, MmcifDataBlock, MmcifDocument, MmcifEntityKind, MmcifEntry,
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

    /// Writes one canonical molecular model as a structural mmCIF data block.
    pub fn write(
        model: &crate::modeling::Model,
        options: MmcifWriteOptions,
    ) -> Result<String, MmcifWriteError> {
        crate::io::write_mmcif_model(model, options)
    }
}

pub mod perception {
    pub use crate::chemistry::{SanitizeError, SanitizeOptions, SanitizeReport};

    use crate::small::SmallMolecule;

    pub mod valence {
        pub use crate::algorithms::{
            perceive_valence, perceive_valence_with_options, ValenceIssue, ValenceModel,
            ValenceOptions, ValenceReport,
        };
    }

    pub mod rings {
        pub use crate::algorithms::{
            perceive_ring_membership, perceive_ring_set, perceive_ring_set_with_options, Ring,
            RingMembership, RingPerceptionError, RingPerceptionOptions, RingSet, RingWork,
        };
    }

    pub mod aromaticity {
        pub use crate::algorithms::{
            perceive_aromaticity, perceive_aromaticity_with_ring_options, AromaticityError,
            AromaticityModel,
        };
    }

    pub mod stereo {
        pub use crate::algorithms::{
            assign_cip_descriptors, assign_cip_descriptors_with_options, perceive_stereo,
            perceive_stereo_with_options, validate_stereo, validate_stereo_with_options,
            CipAssignment, CipAssignmentIssue, CipAssignmentOptions, CipAssignmentReport,
            CipSkipped, CipSkippedReason, StereoCandidate, StereoPerceptionIssue,
            StereoPerceptionOptions, StereoPerceptionReport,
        };
    }

    pub fn sanitize(molecule: &mut SmallMolecule) -> Result<SanitizeReport, SanitizeError> {
        crate::chemistry::sanitize_small_molecule(molecule, SanitizeOptions::default())
    }

    pub fn sanitize_with_options(
        molecule: &mut SmallMolecule,
        options: SanitizeOptions,
    ) -> Result<SanitizeReport, SanitizeError> {
        crate::chemistry::sanitize_small_molecule(molecule, options)
    }

    pub fn sanitize_with_ring_options(
        molecule: &mut SmallMolecule,
        options: SanitizeOptions,
        ring_options: rings::RingPerceptionOptions,
    ) -> Result<SanitizeReport, SanitizeError> {
        crate::chemistry::sanitize_small_molecule_with_ring_options(molecule, options, ring_options)
    }
}

pub mod canon {
    pub use crate::algorithms::CanonicalAtomRanking;

    use crate::core::Molecule;

    pub fn atom_ranking(molecule: &Molecule) -> CanonicalAtomRanking {
        crate::algorithms::canonical_atom_ranking(molecule)
    }
}

/// Explicit small-molecule hydrogen topology transforms.
///
/// These functions never sanitize implicitly. Addition consumes current
/// valence assignments unless `explicit_only` is selected, and removal
/// requires current valence assignments. Successful topology changes
/// invalidate perception state.
pub mod hydrogens {
    pub use crate::algorithms::{
        AddHydrogensOptions, AddHydrogensReport, AddedHydrogen, AddedHydrogenOrigin,
        HydrogenCountAdjustment, HydrogenNormalizationError, RemoveHydrogensReport,
        RemovedHydrogen, RetainedHydrogen, RetainedHydrogenReason,
    };

    use crate::algorithms::{add_hydrogens_to_molecule, remove_hydrogens_from_molecule};
    use crate::small::SmallMolecule;

    /// Materialize stored explicit counts and perceived implicit hydrogens.
    pub fn add_hydrogens(
        molecule: &mut SmallMolecule,
    ) -> Result<AddHydrogensReport, HydrogenNormalizationError> {
        add_hydrogens_with_options(molecule, AddHydrogensOptions::default())
    }

    /// Materialize hydrogens with an explicit growth bound and count policy.
    pub fn add_hydrogens_with_options(
        molecule: &mut SmallMolecule,
        options: AddHydrogensOptions,
    ) -> Result<AddHydrogensReport, HydrogenNormalizationError> {
        add_hydrogens_to_molecule(molecule.graph_mut_raw(), options)
    }

    /// Collapse ordinary degree-one hydrogens without discarding protected state.
    ///
    /// Isotopic, mapped, charged, radical, property-bearing, and otherwise
    /// non-losslessly representable hydrogens remain in the graph and are
    /// described by the returned report.
    pub fn remove_hydrogens(
        molecule: &mut SmallMolecule,
    ) -> Result<RemoveHydrogensReport, HydrogenNormalizationError> {
        remove_hydrogens_from_molecule(molecule.graph_mut_raw())
    }
}

pub mod prelude {
    pub use crate::bio::{MacroMolecule, SmcraHierarchy};
    pub use crate::core::{Atom, AtomId, Bond, BondId, BondOrder, Conformer, Element, Molecule};
    pub use crate::small::{SanitizeOptions, SanitizeReport, SmallMolecule};
    pub use crate::smiles::{
        CanonicalSmilesWriteOptions, IsomericSmilesWriteOptions, SmilesWriteOptions,
    };
}

#[cfg(test)]
mod tests;
