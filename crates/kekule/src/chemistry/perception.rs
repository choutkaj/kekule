use std::fmt;

use crate::algorithms::{
    perceive_aromaticity_in_place, perceive_rings, perceive_valence, AromaticityError,
    RingPerceptionError, RingPerceptionOptions, ValenceError, ValenceOptions,
};
use crate::algorithms::{perceive_conjugation, ConjugationError};
use crate::molecule::ConjugationModel;
use crate::molecule::{AromaticityModel, Molecule, ValenceModel};

/// Failure from the default discrete chemical perception profile.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionError {
    /// Conjugation perception failed; see the source.
    Conjugation(ConjugationError),
    /// Valence perception failed; see the source.
    Valence(ValenceError),
    /// Ring perception failed; see the source.
    Rings(RingPerceptionError),
    /// Aromaticity perception failed; see the source.
    Aromaticity(AromaticityError),
}

impl fmt::Display for PerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Conjugation(_) => f.write_str("conjugation perception failed"),
            Self::Valence(_) => f.write_str("valence perception failed"),
            Self::Rings(_) => f.write_str("ring perception failed"),
            Self::Aromaticity(_) => f.write_str("aromaticity perception failed"),
        }
    }
}

impl std::error::Error for PerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Conjugation(error) => Some(error),
            Self::Valence(error) => Some(error),
            Self::Rings(error) => Some(error),
            Self::Aromaticity(error) => Some(error),
        }
    }
}

/// Install the default discrete perception state on canonical chemistry.
///
/// The fixed profile is RDKit-like valence and implicit hydrogens, followed by
/// the default deterministic ring set, RDKit-like aromaticity and conjugation. The complete
/// operation is transactional and never rewrites represented chemistry.
pub(crate) fn perceive_molecule(molecule: &mut Molecule) -> Result<(), PerceptionError> {
    let previous = molecule.perception().clone();
    if let Err(error) = perceive_molecule_in_place(molecule) {
        molecule
            .install_perception(previous)
            .expect("previous perception state must remain valid");
        return Err(error);
    }
    Ok(())
}

pub(crate) fn perceive_molecule_in_place(molecule: &mut Molecule) -> Result<(), PerceptionError> {
    perceive_valence(molecule, ValenceModel::RdkitLike, ValenceOptions::default())
        .map_err(PerceptionError::Valence)?;
    let ring_options = RingPerceptionOptions::for_molecule(molecule);
    perceive_rings(molecule, ring_options).map_err(PerceptionError::Rings)?;
    perceive_aromaticity_in_place(molecule, AromaticityModel::RdkitLike)
        .map_err(PerceptionError::Aromaticity)?;
    perceive_conjugation(molecule, ConjugationModel::RdkitLike)
        .map_err(PerceptionError::Conjugation)
}
