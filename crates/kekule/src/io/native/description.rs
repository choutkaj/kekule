//! What a native structure set record states besides its geometries.

use std::sync::Arc;

use crate::properties::OwnerProperties;
use crate::structure::{Distribution, GeometryStorage, Order, StructureSet, StructureSetError};
use crate::topology::Topology;

/// A structure set without its geometries: the topology, set properties, order,
/// and distribution of a set whose geometries are stored elsewhere, such as
/// in a trajectory file, and the number of them.
///
/// A native structure set record whose geometries are external stores this.
/// Writing one ([`NativeSource`](super::NativeSource) converts from it)
/// gives such a record, and reading it gives
/// [`NativeObject::StructureSetDescription`](super::NativeObject::StructureSetDescription).
/// The order and distribution describe the external geometries by index, so
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
    distribution: Option<Distribution>,
}

impl StructureSetDescription {
    /// A description of `geometry_count` geometries of `topology`, unordered,
    /// without a distribution or set properties.
    pub fn new(topology: impl Into<Arc<Topology>>, geometry_count: usize) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            geometry_count,
            order: None,
            distribution: None,
        }
    }

    /// The description of `set`, in memory or read from a file: its
    /// topology, set properties, order, distribution, and number of geometries.
    pub fn of<S: GeometryStorage>(set: &StructureSet<S>) -> Self {
        Self {
            topology: set.shared_topology(),
            properties: set.properties().clone(),
            geometry_count: set.len(),
            order: set.order().cloned(),
            distribution: set.distribution().cloned(),
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

    /// Returns this description with a distribution, or without one. A
    /// subset listing every geometry becomes
    /// [`Members::All`](crate::structure::Members::All).
    ///
    /// # Errors
    ///
    /// The errors of [`StructureSet::set_distribution`], against the geometry
    /// count.
    pub fn with_distribution(
        mut self,
        distribution: Option<Distribution>,
    ) -> Result<Self, StructureSetError> {
        self.distribution = distribution
            .map(|mut distribution| {
                distribution.bind(self.geometry_count)?;
                Ok::<_, StructureSetError>(distribution)
            })
            .transpose()?;
        Ok(self)
    }

    /// This description for `geometry_count` geometries whose first
    /// [`Self::geometry_count`] are the geometries it describes, such as a
    /// trajectory file that grew after the description was saved: the
    /// distribution as [`Distribution::extended_to`] gives it, keeping its
    /// weights, partition, states, and temperature bit for bit, and the set
    /// properties. A saved order describes the geometries it was made for, so
    /// `order` replaces it: the longer set's order, or `None` for an
    /// unordered set.
    ///
    /// ```
    /// use kekule::native::StructureSetDescription;
    /// use kekule::structure::{Distribution, Members, Order, Partition};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let saved = StructureSetDescription::new(topology, 3).with_distribution(Some(
    ///     Distribution::new(Members::All, Partition::one_state(3))?,
    /// ))?;
    /// let grown = saved.extended_to(5, Some(Order::Step(vec![0, 10, 20, 30, 40])))?;
    /// assert_eq!(grown.geometry_count(), 5);
    /// assert_eq!(grown.distribution().unwrap().members(), &Members::Subset(vec![0, 1, 2]));
    /// assert_eq!(grown.order(), Some(&Order::Step(vec![0, 10, 20, 30, 40])));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
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
        self.distribution = self
            .distribution
            .map(|distribution| distribution.extended_to(geometry_count))
            .transpose()?;
        self.order = order;
        self.geometry_count = geometry_count;
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

    /// The distribution over the geometries, or `None` without one.
    pub fn distribution(&self) -> Option<&Distribution> {
        self.distribution.as_ref()
    }
}
