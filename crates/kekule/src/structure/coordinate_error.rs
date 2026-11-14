//! [`CoordinateError`]: a failed Cartesian measurement or coordinate edit.

use std::fmt;

use crate::diagnostics::{Capitalized, TopologyMismatch};
use crate::topology::InstanceAtomId;
use crate::units::UnitError;

/// A Cartesian measurement or coordinate edit failed. A failed edit publishes
/// no coordinates.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum CoordinateError {
    /// A measurement or reference atom is not in the topology.
    InvalidAtomId(InstanceAtomId),
    /// A reference atom of an edit appears more than once.
    RepeatedAtom(InstanceAtomId),
    /// A selection, or the structure given to a prepared edit or bond
    /// dihedral, does not share the topology layout of the operation.
    TopologyMismatch(TopologyMismatch),
    /// Consecutive atoms lack a direct bond, in argument order.
    MissingBond {
        /// The earlier atom of the pair.
        a: InstanceAtomId,
        /// The later atom of the pair.
        b: InstanceAtomId,
    },
    /// Removing the traversal bond did not separate moving and fixed references.
    InseparableFragment,
    /// The last reference must be selected and the first must be excluded.
    InvalidMovingSelection,
    /// A quantity could not be converted to the canonical unit; see the source.
    Unit(UnitError),
    /// An edit target is not finite, a target distance is not positive, or a
    /// target angle is outside [0, pi].
    InvalidTarget,
    /// The distance cutoff of a spatial selection is negative or not finite.
    InvalidCutoff,
    /// A rigid-motion origin, axis, or displacement is nonfinite.
    NonFiniteInput,
    /// A direction, rotation plane, angle, or dihedral is undefined because a
    /// direction has zero length.
    DegenerateGeometry,
    /// Arithmetic overflow or loss of precision prevents the measurement or
    /// the edit target.
    NumericalFailure,
}

impl fmt::Display for CoordinateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(f, "{} is not in the topology.", Capitalized(atom)),
            Self::RepeatedAtom(atom) => {
                write!(f, "{} appears twice among the reference atoms.", Capitalized(atom))
            }
            Self::TopologyMismatch(mismatch) => mismatch.fmt(f),
            Self::MissingBond { a, b } => write!(f, "{} and {b} are not bonded.", Capitalized(a)),
            Self::InseparableFragment => {
                f.write_str("The moving fragment reaches a fixed reference atom.")
            }
            Self::InvalidMovingSelection => f.write_str(
                "The moving selection lacks the last reference atom or contains the first.",
            ),
            Self::Unit(error) => error.fmt(f),
            Self::InvalidTarget => f.write_str(
                "The target is not finite, or its distance is not positive, or its angle is outside 0 to 180 degrees.",
            ),
            Self::InvalidCutoff => f.write_str("The distance cutoff is negative or not finite."),
            Self::NonFiniteInput => f.write_str("The origin, axis or displacement is not finite."),
            Self::DegenerateGeometry => {
                f.write_str("The atoms define no direction, plane, angle or dihedral.")
            }
            Self::NumericalFailure => {
                f.write_str("The coordinates exceed the numerical range or precision.")
            }
        }
    }
}

impl std::error::Error for CoordinateError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<UnitError> for CoordinateError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
