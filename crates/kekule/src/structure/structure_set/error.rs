//! The failures of [`super::StructureSet`] operations.

use std::fmt;

use crate::structure::GeometryError;
use crate::topology::TopologySubsetError;
use crate::units::UnitError;

/// What an operation needs and a structure set or ensemble lacks.
///
/// Its `Display` text, such as "needs a time order", can be shown as it is,
/// for example beside an operation a host cannot offer for a set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Requirement {
    /// An [`Order::Time`](super::Order::Time).
    TimeOrder,
    /// Any [`Order`](super::Order).
    AnyOrder,
    /// A time order whose axis has times.
    Times,
    /// An energy on every state of the distribution.
    StateEnergies,
    /// A distribution with a temperature.
    Temperature,
    /// A periodic cell on every geometry.
    Cell,
    /// Every geometry held in memory.
    InMemory,
}

impl fmt::Display for Requirement {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::TimeOrder => "needs a time order",
            Self::AnyOrder => "needs an order",
            Self::Times => "needs times",
            Self::StateEnergies => "needs an energy for every state",
            Self::Temperature => "needs a temperature",
            Self::Cell => "needs a periodic cell on every geometry",
            Self::InMemory => "needs its geometries in memory",
        })
    }
}

/// Failure of a [`super::StructureSet`] or [`Ensemble`](crate::structure::Ensemble)
/// operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureSetError {
    /// The operation needs something the set or ensemble lacks, such as a
    /// time order or an energy for every state.
    Missing(Requirement),
    /// A constructor received no structures.
    EmptySource,
    /// Probabilities were asked of a distribution without samples.
    EmptyDistribution,
    /// A geometry index is outside the set.
    IndexOutOfRange {
        /// The requested geometry index.
        index: usize,
        /// The number of geometries.
        len: usize,
    },
    /// A state index is outside the distribution.
    StateOutOfRange {
        /// The requested state index.
        state: usize,
        /// The number of states.
        states: usize,
    },
    /// A geometry is not a sample of the distribution.
    NotSample {
        /// The geometry's index.
        geometry: usize,
    },
    /// A structure's topology layout differs from the others'
    /// ([`Topology::same_layout`](crate::topology::Topology::same_layout)).
    TopologyMismatch,
    /// A set-wide replacement has the wrong number of geometries.
    GeometryCountMismatch {
        /// The set's geometry count.
        expected: usize,
        /// The number of replacements supplied.
        actual: usize,
    },
    /// An order has a different number of values than the set has geometries.
    OrderLengthMismatch {
        /// The set's geometry count.
        expected: usize,
        /// The number of values.
        actual: usize,
    },
    /// A time axis has a different number of times and steps.
    TimeAxisMismatch {
        /// The number of times.
        times: usize,
        /// The number of steps.
        steps: usize,
    },
    /// An order value, a time or a coordinate, is not finite.
    NonFiniteOrderValue {
        /// The index of the geometry the value belongs to.
        geometry: usize,
    },
    /// A place is of another kind than the set's order, or a place was
    /// given to an unordered set, or none to an ordered one.
    PlaceMismatch {
        /// The index of the geometry the place belongs to.
        geometry: usize,
    },
    /// A geometry has no time on a time axis with times.
    MissingTime {
        /// The geometry's index.
        geometry: usize,
    },
    /// A geometry has a time on a time axis without times.
    UnexpectedTime {
        /// The geometry's index.
        geometry: usize,
    },
    /// A geometry has no MD step on a time axis with steps.
    MissingStep {
        /// The geometry's index.
        geometry: usize,
    },
    /// A geometry has an MD step on a time axis without steps.
    UnexpectedStep {
        /// The geometry's index.
        geometry: usize,
    },
    /// A geometry's time is earlier than the previous geometry's
    /// ([`super::StructureSet::validate_monotonic_time`]).
    NonMonotonicTime {
        /// The index of the out-of-order geometry.
        geometry: usize,
    },
    /// The samples of a distribution are not increasing geometry indices.
    UnsortedMembers,
    /// A distribution of every geometry ([`Members::All`](crate::structure::Members::All))
    /// covers another number of geometries than the set has.
    CoverageMismatch {
        /// The number of geometries the distribution covers.
        covered: usize,
        /// The set's geometry count.
        geometries: usize,
    },
    /// A partition groups another number of samples than the distribution
    /// lists ([`Distribution::new`](crate::structure::Distribution::new),
    /// [`DistributionMut::regroup`](crate::structure::DistributionMut::regroup)).
    PartitionMismatch {
        /// The number of samples.
        expected: usize,
        /// The partition's sample count.
        actual: usize,
    },
    /// Log weights were given in a different number than the distribution
    /// has samples.
    LogWeightCountMismatch {
        /// The number of samples.
        expected: usize,
        /// The number of log weights.
        actual: usize,
    },
    /// States were given in a different number than the partition has.
    StateCountMismatch {
        /// The partition's state count.
        expected: usize,
        /// The number of states.
        actual: usize,
    },
    /// A log weight is not finite.
    NonFiniteLogWeight {
        /// The index of the geometry the log weight belongs to.
        geometry: usize,
    },
    /// A new sample came without a log weight, and the distribution has
    /// weights. Pass one, or clear the weights first.
    LogWeightRequired,
    /// The representative of a state with other samples cannot be removed;
    /// choose another representative with
    /// [`DistributionMut::regroup`](crate::structure::DistributionMut::regroup) first.
    RepresentativeRemoval {
        /// The representative geometry.
        geometry: usize,
        /// Its state.
        state: usize,
    },
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// A Boltzmann population is too small for a finite log weight: the
    /// state's energy lies so far above the lowest, relative to the thermal
    /// energy of
    /// [`DistributionMut::assign_boltzmann_populations`](crate::structure::DistributionMut::assign_boltzmann_populations),
    /// that its population underflows even in log space. A higher
    /// temperature, or leaving the state out with `select_states`, avoids it.
    UnrepresentablePopulation {
        /// The index of the state.
        state: usize,
    },
    /// An energy window is not finite and nonnegative.
    InvalidEnergyWindow,
    /// A population fraction is not in `(0, 1]`.
    InvalidFraction,
    /// A geometry's structural state is invalid; see the source.
    Geometry(GeometryError),
    /// Projecting the topology onto a selection failed; see the source.
    Subset(TopologySubsetError),
    /// A quantity could not be converted to the required unit; see the
    /// source.
    Unit(Box<UnitError>),
}

impl fmt::Display for StructureSetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::counted;
        match self {
            Self::Missing(requirement) => write!(formatter, "The structure set {requirement}."),
            Self::EmptySource => formatter.write_str("There are no structures to collect."),
            Self::EmptyDistribution => formatter.write_str("The distribution has no samples."),
            Self::IndexOutOfRange { index, len } => write!(
                formatter,
                "Geometry {index} is outside the structure set of {}.",
                counted(*len, "geometry", "geometries")
            ),
            Self::StateOutOfRange { state, states } => write!(
                formatter,
                "State {state} is outside the distribution of {}.",
                counted(*states, "state", "states")
            ),
            Self::NotSample { geometry } => {
                write!(formatter, "Geometry {geometry} is not a sample of the distribution.")
            }
            Self::TopologyMismatch => formatter
                .write_str("The structure has a different topology layout from the structure set."),
            Self::GeometryCountMismatch { expected, actual } => write!(
                formatter,
                "The structure set has {}, but {} given.",
                counted(*expected, "geometry", "geometries"),
                given(*actual, "replacement", "replacements")
            ),
            Self::OrderLengthMismatch { expected, actual } => write!(
                formatter,
                "The order has {}, but the structure set has {}.",
                counted(*actual, "value", "values"),
                counted(*expected, "geometry", "geometries")
            ),
            Self::TimeAxisMismatch { times, steps } => write!(
                formatter,
                "The time axis has {} and {}.",
                counted(*times, "time", "times"),
                counted(*steps, "step", "steps")
            ),
            Self::NonFiniteOrderValue { geometry } => {
                write!(
                    formatter,
                    "The order value of geometry {geometry} is not finite."
                )
            }
            Self::PlaceMismatch { geometry } => write!(
                formatter,
                "The place of geometry {geometry} is not on the order of the structure set."
            ),
            Self::MissingTime { geometry } => write!(
                formatter,
                "Geometry {geometry} has no time, but the other geometries have one."
            ),
            Self::UnexpectedTime { geometry } => write!(
                formatter,
                "Geometry {geometry} has a time, but the other geometries have none."
            ),
            Self::MissingStep { geometry } => write!(
                formatter,
                "Geometry {geometry} has no step, but the other geometries have one."
            ),
            Self::UnexpectedStep { geometry } => write!(
                formatter,
                "Geometry {geometry} has a step, but the other geometries have none."
            ),
            Self::NonMonotonicTime { geometry } => write!(
                formatter,
                "Geometry {geometry} is earlier than the geometry before it."
            ),
            Self::UnsortedMembers => {
                formatter.write_str("The samples of the distribution are not in increasing order.")
            }
            Self::CoverageMismatch { covered, geometries } => write!(
                formatter,
                "The distribution covers {}, but the structure set has {geometries}.",
                counted(*covered, "geometry", "geometries")
            ),
            Self::PartitionMismatch { expected, actual } => write!(
                formatter,
                "The partition covers {}, but the distribution has {}.",
                counted(*actual, "sample", "samples"),
                counted(*expected, "sample", "samples")
            ),
            Self::LogWeightCountMismatch { expected, actual } => write!(
                formatter,
                "The distribution has {}, but {} given.",
                counted(*expected, "sample", "samples"),
                given(*actual, "log weight", "log weights")
            ),
            Self::StateCountMismatch { expected, actual } => write!(
                formatter,
                "The partition has {}, but {} given.",
                counted(*expected, "state", "states"),
                given(*actual, "state", "states")
            ),
            Self::NonFiniteLogWeight { geometry } => write!(
                formatter,
                "The log weight of geometry {geometry} is not finite."
            ),
            Self::LogWeightRequired => formatter.write_str(
                "The new sample has no log weight, but the other samples have log weights.",
            ),
            Self::RepresentativeRemoval { geometry, state } => write!(
                formatter,
                "Geometry {geometry} represents state {state}. Remove the state's other samples first."
            ),
            Self::InvalidTemperature => {
                formatter.write_str("The temperature is not finite and positive.")
            }
            Self::UnrepresentablePopulation { state } => write!(
                formatter,
                "The Boltzmann population of state {state} is too small to represent at this temperature. Raise the temperature."
            ),
            Self::InvalidEnergyWindow => {
                formatter.write_str("The energy window is not finite and nonnegative.")
            }
            Self::InvalidFraction => {
                formatter.write_str("The population fraction is not greater than 0 and at most 1.")
            }
            Self::Geometry(_) => formatter.write_str("A geometry is invalid."),
            Self::Subset(_) => formatter.write_str("Cannot subset the structure set."),
            Self::Unit(error) => error.fmt(formatter),
        }
    }
}

/// "1 log weight was" or "2 log weights were".
fn given(count: usize, one: &str, many: &str) -> String {
    let verb = if count == 1 { "was" } else { "were" };
    format!("{} {verb}", crate::diagnostics::counted(count, one, many))
}

impl std::error::Error for StructureSetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Geometry(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<Requirement> for StructureSetError {
    fn from(requirement: Requirement) -> Self {
        Self::Missing(requirement)
    }
}

impl From<GeometryError> for StructureSetError {
    fn from(error: GeometryError) -> Self {
        Self::Geometry(error)
    }
}

impl From<TopologySubsetError> for StructureSetError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<UnitError> for StructureSetError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
