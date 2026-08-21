use std::fmt;

use crate::algorithms::{
    perceive_aromaticity_in_place, perceive_ring_set, perceive_valence, AromaticityError,
    RingPerceptionError, ValenceError,
};
use crate::core::{AromaticityModel, Molecule, ValenceModel};

/// Failure from the default discrete chemical perception profile.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionError {
    Valence(ValenceError),
    Rings(RingPerceptionError),
    Aromaticity(AromaticityError),
}

impl fmt::Display for PerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Valence(error) => write!(f, "{error}"),
            Self::Rings(error) => write!(f, "{error}"),
            Self::Aromaticity(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for PerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Valence(error) => Some(error),
            Self::Rings(error) => Some(error),
            Self::Aromaticity(error) => Some(error),
        }
    }
}

/// Install the default discrete perception state on canonical chemistry.
///
/// The fixed profile is RDKit-like valence and implicit hydrogens, followed by
/// the default deterministic ring set and RDKit-like aromaticity. The complete
/// operation is transactional and never rewrites represented chemistry.
pub fn perceive_molecule(molecule: &mut Molecule) -> Result<(), PerceptionError> {
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
    perceive_valence(molecule, ValenceModel::RdkitLike).map_err(PerceptionError::Valence)?;
    perceive_ring_set(molecule).map_err(PerceptionError::Rings)?;
    perceive_aromaticity_in_place(molecule, AromaticityModel::RdkitLike)
        .map_err(PerceptionError::Aromaticity)
}
