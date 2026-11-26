//! The failures of [`Weights`](super::Weights) operations.

use std::fmt;

use crate::diagnostics::IndexOutOfRange;
use crate::properties::PropertyKey;
use crate::structure::Requirement;
use crate::units::UnitError;

use super::Unweighable;

/// Failure to build or read [`Weights`](super::Weights).
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum WeightsError {
    /// A geometry index is outside the weights.
    IndexOutOfRange(IndexOutOfRange),
    /// A log weight is not finite.
    NonFiniteLogWeight {
        /// The index of the geometry the log weight belongs to.
        geometry: usize,
    },
    /// The weights cover another number of geometries than the structure
    /// set they are used with has, or more than
    /// [`Weights::extended_to`](super::Weights::extended_to) was given.
    CoverageMismatch {
        /// The number of geometries the weights cover.
        covered: usize,
        /// The set's geometry count.
        geometries: usize,
    },
    /// Probabilities were asked of weights without samples.
    NoSamples,
    /// A population fraction is not in `(0, 1]`.
    InvalidFraction,
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// A free energy needs a temperature, and the weights state none: equal
    /// weights, or given weights without one.
    MissingTemperature,
    /// A temperature could not be converted to the canonical unit; see the
    /// source.
    Unit(UnitError),
    /// The structure set has no column under the key
    /// ([`Weights::boltzmann`](super::Weights::boltzmann)). A column without
    /// values is not stored.
    MissingColumn {
        /// The key.
        key: PropertyKey,
    },
    /// A column read as energies holds something else, or real values whose
    /// unit is not an energy.
    NonEnergyColumn {
        /// The column's key.
        key: PropertyKey,
    },
    /// The structure set lacks what
    /// [`Weights::from_bias`](super::Weights::from_bias) reads: a sampling
    /// record, or steps when the record starts after step 0.
    Missing(Requirement),
    /// No geometry of the set can be weighed
    /// ([`Weights::from_bias`](super::Weights::from_bias)). The reason is that
    /// of the first geometry that got furthest through the checks: before the
    /// first sampled step, then a missing bias value.
    NoWeighableSample(Unweighable),
    /// A geometry's bias over the thermal energy is not finite
    /// ([`Weights::from_bias`](super::Weights::from_bias)): its bias values
    /// are not physical.
    InvalidBias {
        /// The geometry's index.
        geometry: usize,
    },
}

impl fmt::Display for WeightsError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IndexOutOfRange(error) => error.fmt(formatter),
            Self::NonFiniteLogWeight { geometry } => write!(
                formatter,
                "The log weight of geometry {geometry} is not finite."
            ),
            Self::CoverageMismatch {
                covered,
                geometries,
            } => write!(
                formatter,
                "The weights cover {}, but the structure set has {geometries}.",
                crate::diagnostics::counted(*covered, "geometry", "geometries")
            ),
            Self::NoSamples => formatter.write_str("The weights have no samples."),
            Self::InvalidFraction => {
                formatter.write_str("The population fraction is not greater than 0 and at most 1.")
            }
            Self::InvalidTemperature => {
                formatter.write_str("The temperature is not finite and positive.")
            }
            Self::MissingTemperature => formatter.write_str("The weights have no temperature."),
            Self::Unit(error) => error.fmt(formatter),
            Self::MissingColumn { key } => {
                write!(formatter, "The structure set has no column \"{key}\".")
            }
            Self::NonEnergyColumn { key } => {
                write!(formatter, "The column \"{key}\" holds no energies.")
            }
            Self::Missing(requirement) => write!(formatter, "The structure set {requirement}."),
            Self::NoWeighableSample(reason) => match reason {
                Unweighable::BeforeSampling => formatter.write_str(
                    "No geometry can be weighed. Every geometry lies before the first sampled step.",
                ),
                Unweighable::MissingBias { key } => write!(
                    formatter,
                    "No geometry can be weighed. Geometries lack values in the bias column \"{key}\"."
                ),
            },
            Self::InvalidBias { geometry } => {
                write!(formatter, "The bias of geometry {geometry} is out of range.")
            }
        }
    }
}

impl std::error::Error for WeightsError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<UnitError> for WeightsError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
