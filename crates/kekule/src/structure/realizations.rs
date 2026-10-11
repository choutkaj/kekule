use std::fmt;

use crate::topology::transform::TopologySubsetError;
use crate::topology::TopologyBuildError;

use super::ConformationError;

/// Failure of an [`super::Ensemble`] or [`super::Trajectory`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum RealizationError {
    /// A constructor received no items.
    EmptySource,
    IndexOutOfRange {
        index: usize,
        len: usize,
    },
    /// An input belongs to a different topology layout.
    TopologyMismatch,
    /// A collection-wide replacement has the wrong number of items.
    ItemCountMismatch {
        expected: usize,
        actual: usize,
    },
    MissingTime {
        frame: usize,
    },
    NonMonotonicTime {
        frame: usize,
    },
    Conformation(ConformationError),
    Subset(TopologySubsetError),
    TopologyBuild(TopologyBuildError),
}

impl fmt::Display for RealizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("realization source is empty"),
            Self::IndexOutOfRange { index, len } => {
                write!(
                    formatter,
                    "index {index} is out of range for {len} realizations"
                )
            }
            Self::TopologyMismatch => {
                formatter.write_str("realization belongs to a different topology layout")
            }
            Self::ItemCountMismatch { expected, actual } => write!(
                formatter,
                "collection requires {expected} replacement items, but received {actual}"
            ),
            Self::MissingTime { frame } => write!(formatter, "frame {frame} has no time"),
            Self::NonMonotonicTime { frame } => {
                write!(formatter, "time decreases at frame {frame}")
            }
            Self::Conformation(error) => write!(formatter, "invalid realization: {error}"),
            Self::Subset(error) => write!(formatter, "cannot subset realizations: {error}"),
            Self::TopologyBuild(error) => write!(formatter, "cannot build topology: {error}"),
        }
    }
}

impl std::error::Error for RealizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Conformation(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::TopologyBuild(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ConformationError> for RealizationError {
    fn from(error: ConformationError) -> Self {
        Self::Conformation(error)
    }
}

impl From<TopologySubsetError> for RealizationError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<TopologyBuildError> for RealizationError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}
