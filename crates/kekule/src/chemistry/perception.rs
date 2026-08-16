use std::fmt;

use crate::algorithms::{
    perceive_aromaticity, perceive_ring_set, perceive_valence, AromaticityError,
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

/// Install the default discrete perception state on normalized chemistry.
///
/// The fixed profile is RDKit-like valence and implicit hydrogens, followed by
/// the default deterministic ring set and RDKit-like aromaticity. The complete
/// operation is transactional and never normalizes represented chemistry.
pub fn perceive_molecule(molecule: &mut Molecule) -> Result<(), PerceptionError> {
    let mut staged = molecule.clone();
    perceive_valence(&mut staged, ValenceModel::RdkitLike).map_err(PerceptionError::Valence)?;
    perceive_ring_set(&mut staged).map_err(PerceptionError::Rings)?;
    perceive_aromaticity(&mut staged, AromaticityModel::RdkitLike)
        .map_err(PerceptionError::Aromaticity)?;
    *molecule = staged;
    Ok(())
}
