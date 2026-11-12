//! The failures of [`super::SampleSet`] operations.

use std::fmt;

use crate::structure::StructureError;
use crate::topology::TopologySubsetError;
use crate::units::UnitError;

/// What an operation needs and a sample set lacks.
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
    /// A [`Distribution`](super::Distribution).
    Distribution,
    /// A distribution with weights.
    Weights,
    /// An energy on every state of the distribution.
    StateEnergies,
    /// A distribution with a temperature.
    Temperature,
    /// A periodic cell on every sample.
    Cell,
    /// Every sample held in memory.
    InMemory,
}

impl fmt::Display for Requirement {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::TimeOrder => "needs a time order",
            Self::AnyOrder => "needs an order",
            Self::Times => "needs times",
            Self::Distribution => "needs a distribution",
            Self::Weights => "needs weights",
            Self::StateEnergies => "needs an energy for every state",
            Self::Temperature => "needs a temperature",
            Self::Cell => "needs a periodic cell on every sample",
            Self::InMemory => "needs its samples in memory",
        })
    }
}

/// Failure of a [`super::SampleSet`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SampleSetError {
    /// The operation needs something the set lacks, such as a time order or
    /// a distribution.
    Missing(Requirement),
    /// A constructor received no models.
    EmptySource,
    /// Probabilities were asked of a distribution without members.
    EmptyDistribution,
    /// A sample index is outside the set.
    IndexOutOfRange {
        /// The requested sample index.
        index: usize,
        /// The number of samples.
        len: usize,
    },
    /// A state index is outside the distribution.
    StateOutOfRange {
        /// The requested state index.
        state: usize,
        /// The number of states.
        states: usize,
    },
    /// A sample is not a member of the distribution.
    NotMember {
        /// The sample's index.
        sample: usize,
    },
    /// A model's topology layout differs from the others'
    /// ([`Topology::same_layout`](crate::topology::Topology::same_layout)).
    TopologyMismatch,
    /// A set-wide replacement has the wrong number of samples.
    SampleCountMismatch {
        /// The set's sample count.
        expected: usize,
        /// The number of replacements supplied.
        actual: usize,
    },
    /// An order has a different number of values than the set has samples.
    OrderLengthMismatch {
        /// The set's sample count.
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
        /// The index of the sample the value belongs to.
        sample: usize,
    },
    /// A place is of another kind than the set's order, or a place was
    /// given to an unordered set, or none to an ordered one.
    PlaceMismatch {
        /// The index of the sample the place belongs to.
        sample: usize,
    },
    /// A sample has no time on a time axis with times.
    MissingTime {
        /// The sample's index.
        sample: usize,
    },
    /// A sample has a time on a time axis without times.
    UnexpectedTime {
        /// The sample's index.
        sample: usize,
    },
    /// A sample has no MD step on a time axis with steps.
    MissingStep {
        /// The sample's index.
        sample: usize,
    },
    /// A sample has an MD step on a time axis without steps.
    UnexpectedStep {
        /// The sample's index.
        sample: usize,
    },
    /// A sample's time is earlier than the previous sample's
    /// ([`super::SampleSet::validate_monotonic_time`]).
    NonMonotonicTime {
        /// The index of the out-of-order sample.
        sample: usize,
    },
    /// The members of a distribution are not increasing sample indices.
    UnsortedMembers,
    /// A distribution of every sample ([`Members::All`](super::Members::All))
    /// covers another number of samples than the set has.
    CoverageMismatch {
        /// The number of samples the distribution covers.
        covered: usize,
        /// The set's sample count.
        samples: usize,
    },
    /// A partition groups another number of members than the distribution
    /// lists ([`Distribution::new`](super::Distribution::new),
    /// [`DistributionMut::regroup`](super::DistributionMut::regroup)).
    PartitionMismatch {
        /// The number of members.
        expected: usize,
        /// The partition's sample count.
        actual: usize,
    },
    /// Weights were given in a different number than the distribution has
    /// members.
    WeightCountMismatch {
        /// The number of members.
        expected: usize,
        /// The number of weights.
        actual: usize,
    },
    /// States were given in a different number than the partition has.
    StateCountMismatch {
        /// The partition's state count.
        expected: usize,
        /// The number of states.
        actual: usize,
    },
    /// A weight is not finite and positive.
    InvalidWeight {
        /// The index of the sample the weight belongs to.
        sample: usize,
    },
    /// A new member came without a weight, and the distribution has
    /// weights. Pass one, or clear the weights first.
    WeightRequired,
    /// The representative of a state with other members cannot be removed;
    /// choose another representative with
    /// [`super::DistributionMut::regroup`] first.
    RepresentativeRemoval {
        /// The representative sample.
        sample: usize,
        /// Its state.
        state: usize,
    },
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// An energy window is not finite and nonnegative.
    InvalidEnergyWindow,
    /// A population fraction is not in `(0, 1]`.
    InvalidFraction,
    /// A sample's structural state is invalid; see the source.
    Structure(StructureError),
    /// Projecting the topology onto a selection failed; see the source.
    Subset(TopologySubsetError),
    /// A quantity could not be converted to the required unit; see the
    /// source.
    Unit(Box<UnitError>),
}

impl fmt::Display for SampleSetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::counted;
        match self {
            Self::Missing(requirement) => write!(formatter, "The sample set {requirement}."),
            Self::EmptySource => formatter.write_str("There are no models to collect."),
            Self::EmptyDistribution => formatter.write_str("The distribution has no members."),
            Self::IndexOutOfRange { index, len } => write!(
                formatter,
                "Sample {index} is outside the sample set of {}.",
                counted(*len, "sample", "samples")
            ),
            Self::StateOutOfRange { state, states } => write!(
                formatter,
                "State {state} is outside the distribution of {}.",
                counted(*states, "state", "states")
            ),
            Self::NotMember { sample } => {
                write!(formatter, "Sample {sample} is outside the distribution.")
            }
            Self::TopologyMismatch => formatter
                .write_str("The model has a different topology layout from the sample set."),
            Self::SampleCountMismatch { expected, actual } => write!(
                formatter,
                "The sample set has {}, but {} given.",
                counted(*expected, "sample", "samples"),
                given(*actual, "replacement", "replacements")
            ),
            Self::OrderLengthMismatch { expected, actual } => write!(
                formatter,
                "The order has {}, but the sample set has {}.",
                counted(*actual, "value", "values"),
                counted(*expected, "sample", "samples")
            ),
            Self::TimeAxisMismatch { times, steps } => write!(
                formatter,
                "The time axis has {} and {}.",
                counted(*times, "time", "times"),
                counted(*steps, "step", "steps")
            ),
            Self::NonFiniteOrderValue { sample } => {
                write!(
                    formatter,
                    "The order value of sample {sample} is not finite."
                )
            }
            Self::PlaceMismatch { sample } => write!(
                formatter,
                "The place of sample {sample} is not on the order of the sample set."
            ),
            Self::MissingTime { sample } => write!(
                formatter,
                "Sample {sample} has no time, but the other samples have one."
            ),
            Self::UnexpectedTime { sample } => write!(
                formatter,
                "Sample {sample} has a time, but the other samples have none."
            ),
            Self::MissingStep { sample } => write!(
                formatter,
                "Sample {sample} has no step, but the other samples have one."
            ),
            Self::UnexpectedStep { sample } => write!(
                formatter,
                "Sample {sample} has a step, but the other samples have none."
            ),
            Self::NonMonotonicTime { sample } => write!(
                formatter,
                "Sample {sample} is earlier than the sample before it."
            ),
            Self::UnsortedMembers => {
                formatter.write_str("The members of the distribution are not in increasing order.")
            }
            Self::CoverageMismatch { covered, samples } => write!(
                formatter,
                "The distribution covers {}, but the sample set has {samples}.",
                counted(*covered, "sample", "samples")
            ),
            Self::PartitionMismatch { expected, actual } => write!(
                formatter,
                "The partition covers {}, but the distribution has {}.",
                counted(*actual, "sample", "samples"),
                counted(*expected, "member", "members")
            ),
            Self::WeightCountMismatch { expected, actual } => write!(
                formatter,
                "The distribution has {}, but {} given.",
                counted(*expected, "member", "members"),
                given(*actual, "weight", "weights")
            ),
            Self::StateCountMismatch { expected, actual } => write!(
                formatter,
                "The partition has {}, but {} given.",
                counted(*expected, "state", "states"),
                given(*actual, "state", "states")
            ),
            Self::InvalidWeight { sample } => write!(
                formatter,
                "The weight of sample {sample} is not finite and positive."
            ),
            Self::WeightRequired => formatter
                .write_str("The new member has no weight, but the other members have weights."),
            Self::RepresentativeRemoval { sample, state } => write!(
                formatter,
                "Sample {sample} represents state {state}. Remove the state's other samples first."
            ),
            Self::InvalidTemperature => {
                formatter.write_str("The temperature is not finite and positive.")
            }
            Self::InvalidEnergyWindow => {
                formatter.write_str("The energy window is not finite and nonnegative.")
            }
            Self::InvalidFraction => {
                formatter.write_str("The population fraction is not greater than 0 and at most 1.")
            }
            Self::Structure(_) => formatter.write_str("A sample is invalid."),
            Self::Subset(_) => formatter.write_str("Cannot subset the sample set."),
            Self::Unit(error) => error.fmt(formatter),
        }
    }
}

/// "1 weight was" or "2 weights were".
fn given(count: usize, one: &str, many: &str) -> String {
    let verb = if count == 1 { "was" } else { "were" };
    format!("{} {verb}", crate::diagnostics::counted(count, one, many))
}

impl std::error::Error for SampleSetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<Requirement> for SampleSetError {
    fn from(requirement: Requirement) -> Self {
        Self::Missing(requirement)
    }
}

impl From<StructureError> for SampleSetError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<TopologySubsetError> for SampleSetError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<UnitError> for SampleSetError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
