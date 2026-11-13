//! Reading the geometries of any set: [`StructureSetRead`] and [`StructureSetReadError`].

use std::fmt;
use std::sync::Arc;

use crate::alignment::AlignmentError;
use crate::properties::OwnerProperties;
use crate::structure::GeometryError;
use crate::topology::Topology;
use crate::traj::periodic::PeriodicError;
use crate::traj::TrajectoryIoError;

use super::{Distribution, Geometry, GeometryView, Order};
use crate::structure::StructureView;

/// Reads the geometries of a set, one at a time and from any thread: a
/// [`StructureSet`](super::StructureSet) in memory, one over an indexed trajectory
/// file, a [`Processed`](crate::traj::Processed) reader, or a host's own
/// reader, such as a cache in front of another one.
///
/// Everything in Kekule that reads many geometries takes `&dyn StructureSetRead`,
/// so the same analysis runs on every one of them. Reading does not change
/// the reader: [`Self::read`] takes `&self`, so threads of the host read
/// geometries of one reader at once. Kekule starts no thread of its own.
///
/// # Contract
///
/// An implementor keeps these. Kekule checks that every geometry it reads is
/// bound to the topology, and fails with [`StructureSetReadError::Geometry`]
/// when one is not; it relies on the rest:
///
/// - [`Self::len`], [`Self::order`], and [`Self::distribution`] change only
///   through `&mut` access to the reader. The order and the distribution
///   describe [`Self::len`] geometries, as a set's do.
/// - [`Self::read`] of an index below [`Self::len`] returns the same geometry
///   every time, bound to [`Self::topology`]: it has one position per atom
///   and one property row per atom and bond. An index past the last fails
///   with [`StructureSetReadError::IndexOutOfRange`].
/// - [`Self::peek`] returns a view of that same geometry when the reader holds
///   it and can lend it without reading, and `None` otherwise. It never
///   reads, so it is cheap. A view is made with [`GeometryView::new`].
///
/// ```
/// use kekule::structure::{Positions, Geometry, StructureSetRead, StructureSet};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let mut set = StructureSet::new(topology.clone());
/// set.push(Geometry::new(Positions::zeros(topology.atom_count())), None)?;
/// let reader: &dyn StructureSetRead = &set;
/// assert_eq!(reader.len(), 1);
/// // A set in memory lends every geometry.
/// assert!(reader.peek(0).is_some());
/// let geometry = reader.read(0)?;
/// assert!(geometry.shares_structure_with(&set.get(0).unwrap()));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub trait StructureSetRead: Send + Sync + fmt::Debug {
    /// The topology every geometry is bound to.
    fn topology(&self) -> &Arc<Topology>;

    /// The number of geometries.
    fn len(&self) -> usize;

    /// Whether there is no geometry.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The order of the geometries, or `None` when they are unordered.
    fn order(&self) -> Option<&Order>;

    /// The distribution over the geometries, or `None` without one.
    fn distribution(&self) -> Option<&Distribution>;

    /// The set-level properties, for a reader that has them; `None` by
    /// default. Structure writers refuse or report what they cannot write.
    fn properties(&self) -> Option<&OwnerProperties> {
        None
    }

    /// Geometry `index` lent by reference, when the reader holds it: always
    /// for a set in memory, never for one in a file. `None` past the last
    /// geometry.
    fn peek(&self, index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>>;

    /// Geometry `index`, owned. A set in memory returns a clone, which shares
    /// the geometry's storage; other readers decode or compute it.
    ///
    /// # Errors
    ///
    /// [`StructureSetReadError::IndexOutOfRange`] past the last geometry, and the
    /// reader's own failures.
    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError>;
}

/// Forwards every method to the reader behind a reference, a box, or an
/// [`Arc`], so a [`Processed`](crate::traj::Processed) reader can borrow or
/// share the reader below it.
macro_rules! forward_sample_read {
    ($($pointer:ty),+) => {$(
        impl<T: StructureSetRead + ?Sized> StructureSetRead for $pointer {
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

            fn properties(&self) -> Option<&OwnerProperties> {
                (**self).properties()
            }

            fn peek(&self, index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
                (**self).peek(index)
            }

            fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
                (**self).read(index)
            }
        }
    )+};
}

forward_sample_read!(&T, Box<T>, Arc<T>);

/// Runs `visit` on geometry `index` of `reader`: lent when the reader lends
/// it, read and checked against the reader's topology otherwise.
pub(crate) fn visit_geometry<T, E>(
    reader: &dyn StructureSetRead,
    index: usize,
    read_error: impl FnOnce(StructureSetReadError) -> E,
    visit: impl FnOnce(StructureView<'_>) -> Result<T, E>,
) -> Result<T, E> {
    if let Some(view) = reader.peek(index) {
        return visit(view.as_structure_view());
    }
    let geometry = reader
        .read(index)
        .and_then(|geometry| {
            GeometryView::new(reader, index, &geometry)?;
            Ok(geometry)
        })
        .map_err(read_error)?;
    visit(GeometryView::of_reader(reader, index, &geometry).as_structure_view())
}

/// Failure to read a geometry from a [`StructureSetRead`].
///
/// Two errors are equal when they report the same failure; a host's own
/// failures ([`Self::Reader`]) are equal only to their clones.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum StructureSetReadError {
    /// No geometry has this index.
    IndexOutOfRange {
        /// The requested index.
        index: usize,
        /// The number of geometries.
        len: usize,
    },
    /// A geometry cannot be read from its trajectory file; see the source.
    Trajectory {
        /// The geometry index.
        geometry: usize,
        /// The failure.
        source: TrajectoryIoError,
    },
    /// A periodic processing step failed on a geometry; see the source.
    Periodic {
        /// The geometry index.
        geometry: usize,
        /// The failure.
        source: Box<PeriodicError>,
    },
    /// Fitting a geometry onto its reference failed; see the source.
    Alignment {
        /// The geometry index.
        geometry: usize,
        /// The failure.
        source: Box<AlignmentError>,
    },
    /// A geometry does not fit the reader's topology, or a step cannot store
    /// its result; see the source.
    Geometry {
        /// The geometry index.
        geometry: usize,
        /// The failure.
        source: Box<GeometryError>,
    },
    /// A host's reader cannot provide the geometry, for its own reason, such
    /// as a cache whose source is gone; see the source. Made with
    /// [`Self::reader`].
    Reader {
        /// The geometry index.
        geometry: usize,
        /// The host's failure, shared so that the error stays cheap to
        /// clone.
        source: Arc<dyn std::error::Error + Send + Sync>,
    },
}

impl StructureSetReadError {
    /// A host reader's own failure to provide geometry `geometry`, such as a
    /// cache that cannot reach the reader below it.
    ///
    /// ```
    /// use kekule::structure::StructureSetReadError;
    ///
    /// let error = StructureSetReadError::reader(3, "The trajectory file is gone.");
    /// assert_eq!(error.to_string(), "Cannot read geometry 3.");
    /// let source = std::error::Error::source(&error).unwrap();
    /// assert_eq!(source.to_string(), "The trajectory file is gone.");
    /// ```
    pub fn reader(
        geometry: usize,
        source: impl Into<Box<dyn std::error::Error + Send + Sync>>,
    ) -> Self {
        Self::Reader {
            geometry,
            source: Arc::from(source.into()),
        }
    }
}

impl PartialEq for StructureSetReadError {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (
                Self::IndexOutOfRange { index, len },
                Self::IndexOutOfRange {
                    index: other_index,
                    len: other_len,
                },
            ) => index == other_index && len == other_len,
            (
                Self::Trajectory { geometry, source },
                Self::Trajectory {
                    geometry: other_geometry,
                    source: other_source,
                },
            ) => geometry == other_geometry && source == other_source,
            (
                Self::Periodic { geometry, source },
                Self::Periodic {
                    geometry: other_geometry,
                    source: other_source,
                },
            ) => geometry == other_geometry && source == other_source,
            (
                Self::Alignment { geometry, source },
                Self::Alignment {
                    geometry: other_geometry,
                    source: other_source,
                },
            ) => geometry == other_geometry && source == other_source,
            (
                Self::Geometry { geometry, source },
                Self::Geometry {
                    geometry: other_geometry,
                    source: other_source,
                },
            ) => geometry == other_geometry && source == other_source,
            (
                Self::Reader { geometry, source },
                Self::Reader {
                    geometry: other_geometry,
                    source: other_source,
                },
            ) => geometry == other_geometry && Arc::ptr_eq(source, other_source),
            _ => false,
        }
    }
}

impl fmt::Display for StructureSetReadError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IndexOutOfRange { index, len } => write!(
                formatter,
                "Geometry {index} does not exist. There are {}.",
                crate::diagnostics::counted(*len, "geometry", "geometries")
            ),
            Self::Periodic { geometry, .. } | Self::Geometry { geometry, .. } => {
                write!(formatter, "Cannot process geometry {geometry}.")
            }
            Self::Alignment { geometry, .. } => {
                write!(formatter, "Cannot fit geometry {geometry}.")
            }
            Self::Trajectory { geometry, .. } | Self::Reader { geometry, .. } => {
                write!(formatter, "Cannot read geometry {geometry}.")
            }
        }
    }
}

impl std::error::Error for StructureSetReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::IndexOutOfRange { .. } => None,
            Self::Trajectory { source, .. } => Some(source),
            Self::Periodic { source, .. } => Some(source.as_ref()),
            Self::Alignment { source, .. } => Some(source.as_ref()),
            Self::Geometry { source, .. } => Some(source.as_ref()),
            Self::Reader { source, .. } => Some(source.as_ref()),
        }
    }
}
