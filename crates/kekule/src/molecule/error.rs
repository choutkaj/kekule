use std::fmt;

use crate::diagnostics::{CapacityExceeded, Capitalized};
use crate::properties::PropertyError;
use crate::stereo::StereoViolation;

use super::{AtomId, BondId, StereoElementId, StereoGroupId};

/// Why a molecule or editor operation was rejected.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq)]
pub enum MoleculeError {
    /// The atom does not exist.
    InvalidAtomId(AtomId),
    /// The bond does not exist.
    InvalidBondId(BondId),
    /// The stereo element does not exist.
    InvalidStereoElementId(StereoElementId),
    /// The stereo group does not exist.
    InvalidStereoGroupId(StereoGroupId),
    /// A stereo element or group is inconsistent with the graph or with the
    /// other stereo assertions.
    InvalidStereoReference(StereoViolation),
    /// A bond would join an atom to itself.
    SelfBond(AtomId),
    /// A bond between these atoms already exists.
    DuplicateBond {
        /// One endpoint.
        a: AtomId,
        /// The other endpoint.
        b: AtomId,
    },
    /// A new atom, bond, stereo element or stereo group cannot be given a
    /// fixed-width ID.
    CapacityExceeded(CapacityExceeded),
    /// A property operation failed.
    Property(PropertyError),
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(id) => write!(f, "{} does not exist.", Capitalized(id)),
            Self::InvalidBondId(id) => write!(f, "{} does not exist.", Capitalized(id)),
            Self::InvalidStereoElementId(id) => write!(f, "{} does not exist.", Capitalized(id)),
            Self::InvalidStereoGroupId(id) => write!(f, "{} does not exist.", Capitalized(id)),
            Self::InvalidStereoReference(issue) => write!(f, "{issue}"),
            Self::SelfBond(id) => write!(f, "Cannot bond {id} to itself."),
            Self::DuplicateBond { a, b } => {
                write!(f, "{} and {b} are already bonded.", Capitalized(a))
            }
            Self::CapacityExceeded(error) => error.fmt(f),
            Self::Property(_) => f.write_str("The molecule property is invalid."),
        }
    }
}

impl std::error::Error for MoleculeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Property(error) => Some(error),
            _ => None,
        }
    }
}

/// Failure to publish checked molecule construction or editing state.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum MoleculePublicationError {
    /// The draft has no atoms.
    EmptyGraph,
    /// The draft is not one connected graph; every published molecule is.
    Disconnected {
        /// The number of graph components.
        components: usize,
    },
    /// Represented stereo is inconsistent with the graph or with itself.
    InvalidStereo(StereoViolation),
    /// A neutral chlorine, bromine, or iodine with only oxygen neighbors and
    /// explicit valence 3, 5, or 7 still has `X=O` bonds. Canonical chemistry
    /// separates them as `[X+]-[O-]`; call
    /// [`MoleculeEditor::normalize_oxohalogens`] before publishing.
    UnnormalizedOxohalogen {
        /// The halogen.
        atom: AtomId,
    },
}

impl fmt::Display for MoleculePublicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGraph => formatter.write_str("The molecule has no atoms."),
            Self::Disconnected { components } => write!(
                formatter,
                "The molecule is not connected. It has {components} components."
            ),
            Self::InvalidStereo(issue) => write!(formatter, "{issue}"),
            Self::UnnormalizedOxohalogen { atom } => write!(
                formatter,
                "Halogen {atom} has oxo double bonds. Normalize oxohalogens first."
            ),
        }
    }
}

impl std::error::Error for MoleculePublicationError {}
