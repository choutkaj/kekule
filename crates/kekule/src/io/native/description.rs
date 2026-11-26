//! What a native structure set record states besides geometries.

use std::sync::Arc;

use crate::diagnostics::{Countable, LengthMismatch};
use crate::properties::{OwnerProperties, PropertyTable};
use crate::structure::{
    GeometryIndex, GeometryStorage, Order, Sampling, StructureSet, StructureSetError,
};
use crate::topology::Topology;

/// A structure set without its geometries: the topology, set properties,
/// order, columns, and sampling record of a set whose geometries are stored
/// elsewhere, such as in a trajectory file, and the number of them.
///
/// A native structure set record whose geometries are external stores this.
/// Writing one ([`NativeSource`](super::NativeSource) converts from it)
/// gives such a record, and reading it gives
/// [`NativeObject::StructureSetDescription`](super::NativeObject::StructureSetDescription).
/// The order and the columns describe the external geometries by index, so
/// they are checked against the geometry count as a set checks them against
/// its geometries.
///
/// ```
/// use kekule::native::{self, NativeReadOptions, StructureSetDescription};
/// use kekule::structure::{Order, TimeAxis};
/// use kekule::units::{Quantity, PICOSECOND};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let times = Quantity::new(vec![0.0, 2.0, 4.0], PICOSECOND);
/// let description = StructureSetDescription::new(topology, 3)
///     .with_order(Some(Order::Time(TimeAxis::new(Some(times), None)?)))?;
/// let bytes = native::write(&description, Default::default())?;
/// let read = native::read(bytes.as_slice(), NativeReadOptions::default())?;
/// let read = read.into_structure_set_description()?;
/// assert_eq!(read.geometry_count(), 3);
/// assert_eq!(read.order(), description.order());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct StructureSetDescription {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    geometry_count: usize,
    order: Option<Order>,
    columns: PropertyTable<GeometryIndex>,
    sampling: Option<Sampling>,
}

impl StructureSetDescription {
    /// A description of `geometry_count` geometries of `topology`, unordered,
    /// without set properties, columns, or a sampling record.
    pub fn new(topology: impl Into<Arc<Topology>>, geometry_count: usize) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            geometry_count,
            order: None,
            columns: PropertyTable::new(geometry_count),
            sampling: None,
        }
    }

    /// The description of `set`, in memory or read from a file: its
    /// topology, set properties, order, columns, sampling record, and number
    /// of geometries.
    pub fn of<S: GeometryStorage>(set: &StructureSet<S>) -> Self {
        Self {
            topology: set.shared_topology(),
            properties: set.properties().clone(),
            geometry_count: set.len(),
            order: set.order().cloned(),
            columns: set.columns().clone(),
            sampling: set.sampling().cloned(),
        }
    }

    /// Returns this description with set properties.
    #[must_use]
    pub fn with_properties(mut self, properties: OwnerProperties) -> Self {
        self.properties = properties;
        self
    }

    /// Returns this description with an order, or unordered.
    ///
    /// # Errors
    ///
    /// The errors of [`StructureSet::set_order`], against the geometry count.
    pub fn with_order(mut self, order: Option<Order>) -> Result<Self, StructureSetError> {
        if let Some(order) = &order {
            order.validate(self.geometry_count)?;
        }
        self.order = order;
        Ok(self)
    }

    /// Returns this description with columns, one row per geometry.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::GeometryCountMismatch`] when the columns have
    /// another number of rows than the description has geometries, and
    /// [`StructureSetError::SamplingColumn`] when a column the sampling
    /// record names holds other values than real energies.
    pub fn with_columns(
        mut self,
        columns: PropertyTable<GeometryIndex>,
    ) -> Result<Self, StructureSetError> {
        if columns.len() != self.geometry_count {
            return Err(StructureSetError::GeometryCountMismatch(
                LengthMismatch::new(Countable::Geometry, self.geometry_count, columns.len()),
            ));
        }
        if let Some(sampling) = &self.sampling {
            sampling.check_columns(&columns)?;
        }
        self.columns = columns;
        Ok(self)
    }

    /// Returns this description with a sampling record, or without one.
    ///
    /// # Errors
    ///
    /// The errors of [`StructureSet::set_sampling`], against the columns.
    pub fn with_sampling(mut self, sampling: Option<Sampling>) -> Result<Self, StructureSetError> {
        if let Some(sampling) = &sampling {
            sampling.check_columns(&self.columns)?;
        }
        self.sampling = sampling;
        Ok(self)
    }

    /// This description for `geometry_count` geometries whose first
    /// [`Self::geometry_count`] are the geometries it describes, such as a
    /// trajectory file that grew after the description was saved, keeping
    /// the set properties, the sampling record, and the columns, whose new
    /// cells are missing. A
    /// saved order describes the geometries it was made for, so `order`
    /// replaces it: the longer set's order, or `None` for an unordered set.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::CoverageMismatch`] for fewer geometries than the
    /// description describes, and the errors of [`Self::with_order`]
    /// against `geometry_count`.
    pub fn extended_to(
        mut self,
        geometry_count: usize,
        order: Option<Order>,
    ) -> Result<Self, StructureSetError> {
        if geometry_count < self.geometry_count {
            return Err(StructureSetError::CoverageMismatch {
                covered: self.geometry_count,
                geometries: geometry_count,
            });
        }
        if let Some(order) = &order {
            order.validate(geometry_count)?;
        }
        self.order = order;
        self.geometry_count = geometry_count;
        self.columns.resize_missing(geometry_count);
        Ok(self)
    }

    /// The topology every geometry shares.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Set-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// The number of geometries described.
    pub const fn geometry_count(&self) -> usize {
        self.geometry_count
    }

    /// The order of the geometries, or `None` when they are unordered.
    pub fn order(&self) -> Option<&Order> {
        self.order.as_ref()
    }

    /// The columns of the geometries, one row per geometry.
    pub fn columns(&self) -> &PropertyTable<GeometryIndex> {
        &self.columns
    }

    /// The sampling record of the geometries, when there is one.
    pub fn sampling(&self) -> Option<&Sampling> {
        self.sampling.as_ref()
    }
}
