//! Represented stereochemistry, its validation, and stereo algorithms.
//!
//! The root holds the stereo elements and groups a [`Molecule`](crate::molecule::Molecule)
//! represents, the issues their validation reports, and the shared
//! symmetry-search bounds. Edits, publication, and
//! [`MoleculeEditor::validate`](crate::molecule::MoleculeEditor::validate)
//! share one validator.
//! Expert algorithms live in submodules: [`cip`] assigns
//! opt-in CIP descriptors, [`coordinates`] infers
//! and materializes stereo from coordinates, and
//! [`cleanup`] detects stereogenic candidates and removes
//! ineligible assertions. Coordinate inference and candidate detection are
//! read-only; coordinate materialization is an explicitly named representation
//! transform.
//!
//! ```
//! use kekule::smiles;
//! use kekule::stereo::{cip, StereoDescriptor};
//!
//! // L-alanine.
//! let mut alanine = smiles::to_molecules("C[C@H](N)C(=O)O")?.remove(0);
//! assert_eq!(alanine.stereo_elements().count(), 1);
//!
//! alanine.perceive()?;
//! cip::assign_cip_descriptors(&mut alanine, Default::default())?;
//! let (center, _) = alanine.stereo_elements().next().unwrap();
//! let cip_state = alanine.perception().cip().expect("assigned");
//! assert_eq!(cip_state.cip_descriptor(center), Some(StereoDescriptor::S));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

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
        CipAssignmentOptions, CipAssignmentReport, CipRankingError, CipSkipped, CipSkippedReason,
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
