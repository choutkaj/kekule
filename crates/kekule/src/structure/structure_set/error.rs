//! The failures of [`super::StructureSet`] operations.

use crate::diagnostics::{IndexOutOfRange, LengthMismatch, TopologyMismatch};
use std::fmt;

use crate::properties::{PropertyError, PropertyKey};
use crate::structure::GeometryError;
use crate::topology::TopologySubsetError;
use crate::units::UnitError;

/// What an operation needs and a structure set lacks.
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
    /// A step per geometry: a step order, or a time order with steps.
    Steps,
    /// An [`Order::Coordinate`](super::Order::Coordinate).
    CoordinateOrder,
    /// A periodic cell on every geometry.
    Cell,
    /// Every geometry held in memory.
    InMemory,
    /// A sampling record ([`Sampling`](crate::structure::Sampling)).
    Sampling,
}

impl fmt::Display for Requirement {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::TimeOrder => "needs a time order",
            Self::AnyOrder => "needs an order",
            Self::Times => "needs times",
            Self::Steps => "needs steps",
            Self::CoordinateOrder => "needs a coordinate order",
            Self::Cell => "needs a periodic cell on every geometry",
            Self::InMemory => "needs its geometries in memory",
            Self::Sampling => "needs a sampling record",
        })
    }
}

/// Failure of a [`super::StructureSet`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureSetError {
    /// The operation needs something the set lacks, such as a time order.
    Missing(Requirement),
    /// A constructor received no structures.
    EmptySource,
    /// A geometry index is outside the set.
    IndexOutOfRange(IndexOutOfRange),
    /// A structure's topology layout differs from the others'
    /// ([`Topology::same_layout`](crate::topology::Topology::same_layout),
    /// [`TopologyMismatch::Structure`]).
    TopologyMismatch(TopologyMismatch),
    /// A set-wide replacement has the wrong number of geometries.
    GeometryCountMismatch(LengthMismatch),
    /// An order has a different number of values than the set has geometries.
    OrderLengthMismatch(LengthMismatch),
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
    /// The time of a frame, which places it on a time order, is not finite
    /// ([`FrameView::new`](crate::traj::FrameView::new),
    /// [`FrameBuffer::set_time`](crate::traj::FrameBuffer::set_time)).
    NonFiniteTime,
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
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// A pressure is not finite and positive
    /// ([`ThermodynamicState::isothermal_isobaric`](crate::structure::ThermodynamicState::isothermal_isobaric)).
    InvalidPressure,
    /// A structure set description covers more geometries than it is
    /// extended to
    /// ([`StructureSetDescription::extended_to`](crate::native::StructureSetDescription::extended_to)).
    CoverageMismatch {
        /// The number of geometries described.
        covered: usize,
        /// The number of geometries it is extended to.
        geometries: usize,
    },
    /// The set has no column under the key
    /// ([`StructureSet::column_values`](super::StructureSet::column_values)).
    /// A column without values is not stored.
    MissingColumn {
        /// The key.
        key: PropertyKey,
    },
    /// A set column holds booleans or strings where numbers are asked for
    /// ([`StructureSet::column_values`](super::StructureSet::column_values)).
    NonNumericColumn {
        /// The column's key.
        key: PropertyKey,
    },
    /// An owner value of a geometry entering a set has another type than the
    /// set's column of that key, or a real value's unit does not convert to
    /// the column's. A set keeps each geometry's owner values in its columns
    /// ([`StructureSet::push`](super::StructureSet::push) and the other
    /// constructors); convert or rename the value first.
    ColumnTypeMismatch {
        /// The value's key.
        key: PropertyKey,
        /// The index of the geometry the value belongs to.
        geometry: usize,
    },
    /// A column edit through
    /// [`SetColumnsMut`](super::SetColumnsMut) is invalid for the table, such
    /// as a column of another length or a value of another type; see the
    /// source.
    Column(Box<PropertyError>),
    /// A column the set's [`Sampling`](crate::structure::Sampling) record
    /// names holds something other than real energies, or would after the
    /// edit ([`StructureSet::set_sampling`](super::StructureSet::set_sampling),
    /// [`SetColumnsMut`](super::SetColumnsMut)).
    SamplingColumn {
        /// The column's key.
        key: PropertyKey,
    },
    /// A sampling record names a column twice
    /// ([`Sampling::with_bias`](crate::structure::Sampling::with_bias)).
    RepeatedSamplingColumn {
        /// The column's key.
        key: PropertyKey,
    },
    /// A column the set's sampling record names cannot be removed
    /// ([`SetColumnsMut::remove`](super::SetColumnsMut::remove),
    /// [`SetColumnsMut::clear`](super::SetColumnsMut::clear)). Clear the
    /// record with [`StructureSet::set_sampling`](super::StructureSet::set_sampling)
    /// first.
    ColumnInUse {
        /// The column's key.
        key: PropertyKey,
    },
    /// The set already has another sampling record
    /// ([`StructureSet::join_frames`](super::StructureSet::join_frames)).
    SamplingMismatch,
    /// A geometry's structural state is invalid; see the source.
    Geometry(GeometryError),
    /// Projecting the topology onto a selection failed; see the source.
    Subset(TopologySubsetError),
    /// A quantity could not be converted to the required unit; see the
    /// source.
    Unit(UnitError),
}

impl fmt::Display for StructureSetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::counted;
        match self {
            Self::Missing(requirement) => write!(formatter, "The structure set {requirement}."),
            Self::EmptySource => formatter.write_str("There are no structures to collect."),
            Self::IndexOutOfRange(error) => error.fmt(formatter),
            Self::TopologyMismatch(mismatch) => mismatch.fmt(formatter),
            Self::GeometryCountMismatch(error) => error.fmt(formatter),
            Self::OrderLengthMismatch(error) => error.fmt(formatter),
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
            Self::NonFiniteTime => formatter.write_str("The time is not finite."),
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
            Self::InvalidTemperature => {
                formatter.write_str("The temperature is not finite and positive.")
            }
            Self::InvalidPressure => formatter.write_str("The pressure is not finite and positive."),
            Self::CoverageMismatch { covered, geometries } => write!(
                formatter,
                "The description covers {}, but the structure set has {geometries}.",
                counted(*covered, "geometry", "geometries")
            ),
            Self::MissingColumn { key } => {
                write!(formatter, "The structure set has no column \"{key}\".")
            }
            Self::NonNumericColumn { key } => {
                write!(formatter, "The column \"{key}\" holds no numbers.")
            }
            Self::ColumnTypeMismatch { key, geometry } => write!(
                formatter,
                "The value \"{key}\" of geometry {geometry} has another type or unit than the column \"{key}\"."
            ),
            Self::Column(_) => {
                formatter.write_str("Cannot change the columns of the structure set.")
            }
            Self::SamplingColumn { key } => {
                write!(formatter, "The column \"{key}\" holds no energies.")
            }
            Self::RepeatedSamplingColumn { key } => write!(
                formatter,
                "The sampling record names the column \"{key}\" twice."
            ),
            Self::ColumnInUse { key } => write!(
                formatter,
                "The sampling record names the column \"{key}\". Remove the record first."
            ),
            Self::SamplingMismatch => {
                formatter.write_str("The structure set has another sampling record.")
            }
            Self::Geometry(_) => formatter.write_str("A geometry is invalid."),
            Self::Subset(_) => formatter.write_str("Cannot subset the structure set."),
            Self::Unit(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for StructureSetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Column(error) => Some(error.as_ref()),
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
        Self::Unit(error)
    }
}
