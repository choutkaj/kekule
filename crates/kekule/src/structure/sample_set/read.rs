//! Reading the samples of any set: [`SampleRead`] and [`SampleReadError`].

use std::fmt;
use std::sync::Arc;

use crate::alignment::AlignmentError;
use crate::structure::StructureError;
use crate::topology::Topology;
use crate::traj::periodic::PeriodicError;
use crate::traj::TrajectoryIoError;

use super::{Distribution, Order, Sample, SampleView};
use crate::structure::ModelView;

/// Reads the samples of a set, one at a time and from any thread: a
/// [`SampleSet`](super::SampleSet) in memory, one over an indexed trajectory
/// file, a [`Processed`](crate::traj::Processed) reader, or a host's own
/// reader, such as a cache in front of another one.
///
/// Everything in Kekule that reads many samples takes `&dyn SampleRead`,
/// so the same analysis runs on every one of them. Reading does not change
/// the reader: [`Self::read`] takes `&self`, so threads of the host read
/// samples of one reader at once. Kekule starts no thread of its own.
///
/// # Contract
///
/// An implementor keeps these. Kekule checks that every sample it reads is
/// bound to the topology, and fails with [`SampleReadError::Structure`]
/// when one is not; it relies on the rest:
///
/// - [`Self::len`], [`Self::order`], and [`Self::distribution`] change only
///   through `&mut` access to the reader. The order and the distribution
///   describe [`Self::len`] samples, as a set's do.
/// - [`Self::read`] of an index below [`Self::len`] returns the same sample
///   every time, bound to [`Self::topology`]: it has one position per atom
///   and one property row per atom and bond. An index past the last fails
///   with [`SampleReadError::IndexOutOfRange`].
/// - [`Self::peek`] returns a view of that same sample when the reader holds
///   it and can lend it without reading, and `None` otherwise. It never
///   reads, so it is cheap. A view is made with [`SampleView::new`].
///
/// ```
/// use kekule::structure::{Positions, Sample, SampleRead, SampleSet};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let mut set = SampleSet::new(topology.clone());
/// set.push(Sample::new(Positions::zeros(topology.atom_count())), None)?;
/// let reader: &dyn SampleRead = &set;
/// assert_eq!(reader.len(), 1);
/// // A set in memory lends every sample.
/// assert!(reader.peek(0).is_some());
/// let sample = reader.read(0)?;
/// assert!(sample.shares_structure_with(&set.get(0).unwrap()));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait SampleRead: Send + Sync + fmt::Debug {
    /// The topology every sample is bound to.
    fn topology(&self) -> &Arc<Topology>;

    /// The number of samples.
    fn len(&self) -> usize;

    /// Whether there is no sample.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The order of the samples, or `None` when they are unordered.
    fn order(&self) -> Option<&Order>;

    /// The distribution over the samples, or `None` without one.
    fn distribution(&self) -> Option<&Distribution>;

    /// Sample `index` lent by reference, when the reader holds it: always
    /// for a set in memory, never for one in a file. `None` past the last
    /// sample.
    fn peek(&self, index: usize) -> Option<SampleView<'_, dyn SampleRead>>;

    /// Sample `index`, owned. A set in memory returns a clone, which shares
    /// the sample's storage; other readers decode or compute it.
    ///
    /// # Errors
    ///
    /// [`SampleReadError::IndexOutOfRange`] past the last sample, and the
    /// reader's own failures.
    fn read(&self, index: usize) -> Result<Sample, SampleReadError>;
}

/// Forwards every method to the reader behind a reference, a box, or an
/// [`Arc`], so a [`Processed`](crate::traj::Processed) reader can borrow or
/// share the reader below it.
macro_rules! forward_sample_read {
    ($($pointer:ty),+) => {$(
        impl<T: SampleRead + ?Sized> SampleRead for $pointer {
            fn topology(&self) -> &Arc<Topology> {
                (**self).topology()
            }

            fn len(&self) -> usize {
                (**self).len()
            }

            fn order(&self) -> Option<&Order> {
                (**self).order()
            }

            fn distribution(&self) -> Option<&Distribution> {
                (**self).distribution()
            }

            fn peek(&self, index: usize) -> Option<SampleView<'_, dyn SampleRead>> {
                (**self).peek(index)
            }

            fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
                (**self).read(index)
            }
        }
    )+};
}

forward_sample_read!(&T, Box<T>, Arc<T>);

/// Runs `visit` on sample `index` of `reader`: lent when the reader lends
/// it, read and checked against the reader's topology otherwise.
pub(crate) fn visit_sample<T, E>(
    reader: &dyn SampleRead,
    index: usize,
    read_error: impl FnOnce(SampleReadError) -> E,
    visit: impl FnOnce(ModelView<'_>) -> Result<T, E>,
) -> Result<T, E> {
    if let Some(view) = reader.peek(index) {
        return visit(view.as_model_view());
    }
    let sample = reader
        .read(index)
        .and_then(|sample| {
            SampleView::new(reader, index, &sample)?;
            Ok(sample)
        })
        .map_err(read_error)?;
    visit(SampleView::of_reader(reader, index, &sample).as_model_view())
}

/// Failure to read a sample from a [`SampleRead`].
///
/// Two errors are equal when they report the same failure; a host's own
/// failures ([`Self::Reader`]) are equal only to their clones.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum SampleReadError {
    /// No sample has this index.
    IndexOutOfRange {
        /// The requested index.
        index: usize,
        /// The number of samples.
        len: usize,
    },
    /// The trajectory file cannot be read; see the source.
    Trajectory(TrajectoryIoError),
    /// A periodic processing step failed on a sample; see the source.
    Periodic {
        /// The sample index.
        sample: usize,
        /// The failure.
        source: Box<PeriodicError>,
    },
    /// Fitting a sample onto its reference failed; see the source.
    Alignment {
        /// The sample index.
        sample: usize,
        /// The failure.
        source: Box<AlignmentError>,
    },
    /// A sample does not fit the reader's topology, or a step cannot store
    /// its result; see the source.
    Structure {
        /// The sample index.
        sample: usize,
        /// The failure.
        source: Box<StructureError>,
    },
    /// A host's reader cannot provide the sample, for its own reason, such
    /// as a cache whose source is gone; see the source. Made with
    /// [`Self::reader`].
    Reader {
        /// The sample index.
        sample: usize,
        /// The host's failure, shared so that the error stays cheap to
        /// clone.
        source: Arc<dyn std::error::Error + Send + Sync>,
    },
}

impl SampleReadError {
    /// A host reader's own failure to provide sample `sample`, such as a
    /// cache that cannot reach the reader below it.
    ///
    /// ```
    /// use kekule::structure::SampleReadError;
    ///
    /// let error = SampleReadError::reader(3, "The trajectory file is gone.");
    /// assert_eq!(error.to_string(), "Cannot read sample 3.");
    /// let source = std::error::Error::source(&error).unwrap();
    /// assert_eq!(source.to_string(), "The trajectory file is gone.");
    /// ```
    pub fn reader(
        sample: usize,
        source: impl Into<Box<dyn std::error::Error + Send + Sync>>,
    ) -> Self {
        Self::Reader {
            sample,
            source: Arc::from(source.into()),
        }
    }
}

impl PartialEq for SampleReadError {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (
                Self::IndexOutOfRange { index, len },
                Self::IndexOutOfRange {
                    index: other_index,
                    len: other_len,
                },
            ) => index == other_index && len == other_len,
            (Self::Trajectory(error), Self::Trajectory(other)) => error == other,
            (
                Self::Periodic { sample, source },
                Self::Periodic {
                    sample: other_sample,
                    source: other_source,
                },
            ) => sample == other_sample && source == other_source,
            (
                Self::Alignment { sample, source },
                Self::Alignment {
                    sample: other_sample,
                    source: other_source,
                },
            ) => sample == other_sample && source == other_source,
            (
                Self::Structure { sample, source },
                Self::Structure {
                    sample: other_sample,
                    source: other_source,
                },
            ) => sample == other_sample && source == other_source,
            (
                Self::Reader { sample, source },
                Self::Reader {
                    sample: other_sample,
                    source: other_source,
                },
            ) => sample == other_sample && Arc::ptr_eq(source, other_source),
            _ => false,
        }
    }
}

impl fmt::Display for SampleReadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IndexOutOfRange { index, len } => write!(
                formatter,
                "Sample {index} does not exist. There are {}.",
                crate::diagnostics::counted(*len, "sample", "samples")
            ),
            Self::Trajectory(error) => fmt::Display::fmt(error, formatter),
            Self::Periodic { sample, .. } | Self::Structure { sample, .. } => {
                write!(formatter, "Cannot process sample {sample}.")
            }
            Self::Alignment { sample, .. } => write!(formatter, "Cannot fit sample {sample}."),
            Self::Reader { sample, .. } => write!(formatter, "Cannot read sample {sample}."),
        }
    }
}

impl std::error::Error for SampleReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::IndexOutOfRange { .. } => None,
            Self::Trajectory(error) => std::error::Error::source(error),
            Self::Periodic { source, .. } => Some(source.as_ref()),
            Self::Alignment { source, .. } => Some(source.as_ref()),
            Self::Structure { source, .. } => Some(source.as_ref()),
            Self::Reader { source, .. } => Some(source.as_ref()),
        }
    }
}

impl From<TrajectoryIoError> for SampleReadError {
    fn from(error: TrajectoryIoError) -> Self {
        Self::Trajectory(error)
    }
}
