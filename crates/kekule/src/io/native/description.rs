//! What a native structure set, ensemble, or distribution record states
//! besides geometries.

use std::sync::Arc;

use crate::diagnostics::{Countable, LengthMismatch};
use crate::properties::{OwnerProperties, PropertyTable};
use crate::structure::{
    Distribution, Ensemble, GeometryIndex, GeometryStorage, Order, StructureSet, StructureSetError,
};
use crate::topology::Topology;

/// A structure set without its geometries: the topology, set properties,
/// order, and columns of a set whose geometries are stored elsewhere, such as
/// in a trajectory file, and the number of them.
///
/// A native structure set record whose geometries are external stores this.
/// Writing one ([`NativeSource`](super::NativeSource) converts from it)
/// gives such a record, and reading it gives
/// [`NativeObject::StructureSetDescription`](super::NativeObject::StructureSetDescription).
/// The order and the columns describe the external geometries by index, so
/// they are checked against the geometry count as a set checks them against
/// its geometries.
/// [`EnsembleDescription`] adds a distribution over them.
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
}

impl StructureSetDescription {
    /// A description of `geometry_count` geometries of `topology`, unordered,
    /// without set properties or columns.
    pub fn new(topology: impl Into<Arc<Topology>>, geometry_count: usize) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            geometry_count,
            order: None,
            columns: PropertyTable::new(geometry_count),
        }
    }

    /// The description of `set`, in memory or read from a file: its
    /// topology, set properties, order, columns, and number of geometries.
    pub fn of<S: GeometryStorage>(set: &StructureSet<S>) -> Self {
        Self {
            topology: set.shared_topology(),
            properties: set.properties().clone(),
            geometry_count: set.len(),
            order: set.order().cloned(),
            columns: set.columns().clone(),
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
    /// another number of rows than the description has geometries.
    pub fn with_columns(
        mut self,
        columns: PropertyTable<GeometryIndex>,
    ) -> Result<Self, StructureSetError> {
        if columns.len() != self.geometry_count {
            return Err(StructureSetError::GeometryCountMismatch(
                LengthMismatch::new(Countable::Geometry, self.geometry_count, columns.len()),
            ));
        }
        self.columns = columns;
        Ok(self)
    }

    /// This description for `geometry_count` geometries whose first
    /// [`Self::geometry_count`] are the geometries it describes, such as a
    /// trajectory file that grew after the description was saved, keeping
    /// the set properties and the columns, whose new cells are missing. A
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
}

/// An ensemble without its geometries: the [`StructureSetDescription`] of a
/// set whose geometries are stored elsewhere, such as in a trajectory file,
/// and a distribution over them.
///
/// A native ensemble record whose geometries are external stores this.
/// Writing one ([`NativeSource`](super::NativeSource) converts from it)
/// gives such a record, and reading it gives
/// [`NativeObject::EnsembleDescription`](super::NativeObject::EnsembleDescription).
/// The distribution describes the external geometries by index, so it is
/// checked against the geometry count as an [`Ensemble`] checks it against
/// its set.
///
/// ```
/// use kekule::native::{self, EnsembleDescription, NativeReadOptions, StructureSetDescription};
/// use kekule::structure::{Distribution, Order, Partition, Samples};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// // The last two of three frames, reweighted, as one state.
/// let distribution = Distribution::new(Samples::Subset(vec![1, 2]), Partition::one_state(2))?
///     .with_log_weights(vec![-1.5, 0.0])?;
/// let saved = EnsembleDescription::new(StructureSetDescription::new(topology, 3), distribution)?;
/// let bytes = native::write(&saved, Default::default())?;
/// let read = native::read(bytes.as_slice(), NativeReadOptions::default())?;
/// let read = read.into_ensemble_description()?;
/// assert_eq!(read.distribution(), saved.distribution());
///
/// // The file grew to five frames; the new ones are not samples.
/// let grown = read.extended_to(5, None)?;
/// assert_eq!(grown.set().geometry_count(), 5);
/// assert_eq!(grown.distribution().samples(), &Samples::Subset(vec![1, 2]));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct EnsembleDescription {
    set: StructureSetDescription,
    distribution: Distribution,
}

impl EnsembleDescription {
    /// `distribution` over the geometries `set` describes. A subset listing
    /// every geometry becomes [`Samples::All`](crate::structure::Samples::All).
    ///
    /// # Errors
    ///
    /// The errors of [`Ensemble::new`], against the geometry count.
    pub fn new(
        set: StructureSetDescription,
        mut distribution: Distribution,
    ) -> Result<Self, StructureSetError> {
        distribution.bind(set.geometry_count())?;
        Ok(Self { set, distribution })
    }

    /// `distribution`, already bound to the geometries `set` describes.
    pub(super) const fn bound(set: StructureSetDescription, distribution: Distribution) -> Self {
        Self { set, distribution }
    }

    /// The description of `ensemble`, in memory or read from a file: the
    /// [`StructureSetDescription::of`] its set and its distribution.
    pub fn of<S: GeometryStorage>(ensemble: &Ensemble<S>) -> Self {
        Self {
            set: StructureSetDescription::of(ensemble.set()),
            distribution: ensemble.distribution().clone(),
        }
    }

    /// This description for `geometry_count` geometries whose first
    /// [`StructureSetDescription::geometry_count`] are the geometries it
    /// describes, such as a trajectory file that grew after the description
    /// was saved: the set as [`StructureSetDescription::extended_to`] gives
    /// it, with `order` the longer set's, and the distribution as
    /// [`Distribution::extended_to`] gives it, keeping its log weights,
    /// partition, states, and temperature bit for bit.
    ///
    /// # Errors
    ///
    /// The errors of [`StructureSetDescription::extended_to`].
    pub fn extended_to(
        self,
        geometry_count: usize,
        order: Option<Order>,
    ) -> Result<Self, StructureSetError> {
        let set = self.set.extended_to(geometry_count, order)?;
        let distribution = self.distribution.extended_to(geometry_count)?;
        Ok(Self { set, distribution })
    }

    /// What the record states of the set besides its geometries.
    pub fn set(&self) -> &StructureSetDescription {
        &self.set
    }

    /// The distribution over the geometries.
    pub fn distribution(&self) -> &Distribution {
        &self.distribution
    }

    /// The set's description and the distribution.
    pub fn into_parts(self) -> (StructureSetDescription, Distribution) {
        (self.set, self.distribution)
    }
}

/// A distribution without its set: the [`Distribution`] over the geometries
/// of a set stored elsewhere, and the number of those geometries.
///
/// A native distribution record stores this, so a host keeps a set's
/// weights, states, energies, and temperature apart from its topology and
/// geometries. Writing one ([`NativeSource`](super::NativeSource) converts
/// from it) gives such a record, and reading it gives
/// [`NativeObject::DistributionDescription`](super::NativeObject::DistributionDescription).
/// The distribution describes the geometries by index, so it is checked
/// against the geometry count as an [`Ensemble`] checks it against its set.
///
/// ```
/// use kekule::math::Point3;
/// use kekule::native::{self, DistributionDescription, NativeReadOptions};
/// use kekule::structure::{
///     Distribution, Ensemble, Geometry, Partition, Positions, Samples, StructureSet,
/// };
/// use kekule::units::{Quantity, NANOMETER};
///
/// // Frames 1 and 2 of a three-frame set, reweighted, as one state.
/// let distribution = Distribution::new(Samples::Subset(vec![1, 2]), Partition::one_state(2))?
///     .with_log_weights(vec![-1.5, 0.0])?;
/// let saved = DistributionDescription::new(3, distribution)?;
/// let bytes = native::write(&saved, Default::default())?;
/// let read = native::read(bytes.as_slice(), NativeReadOptions::default())?
///     .into_distribution_description()?;
/// assert_eq!(read, saved);
///
/// // The set it describes, kept elsewhere, has as many geometries.
/// let geometry =
///     |x| Positions::new(Quantity::new(vec![Point3::new(x, 0.0, 0.0)], NANOMETER)).map(Geometry::new);
/// let topology = kekule::smiles::to_topology("O")?;
/// let set = StructureSet::from_geometries(topology, [geometry(0.0)?, geometry(0.1)?, geometry(0.2)?])?;
/// assert_eq!(set.len(), read.geometry_count());
/// let ensemble = Ensemble::new(set, read.into_distribution())?;
/// assert_eq!(ensemble.distribution().sample_count(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct DistributionDescription {
    geometry_count: usize,
    distribution: Distribution,
}

impl DistributionDescription {
    /// `distribution` over `geometry_count` geometries. A subset listing
    /// every geometry becomes [`Samples::All`](crate::structure::Samples::All).
    ///
    /// # Errors
    ///
    /// The errors of [`Ensemble::new`], against the geometry count.
    pub fn new(
        geometry_count: usize,
        mut distribution: Distribution,
    ) -> Result<Self, StructureSetError> {
        distribution.bind(geometry_count)?;
        Ok(Self {
            geometry_count,
            distribution,
        })
    }

    /// A distribution already bound to `geometry_count` geometries.
    pub(super) const fn bound(geometry_count: usize, distribution: Distribution) -> Self {
        Self {
            geometry_count,
            distribution,
        }
    }

    /// The distribution of `ensemble`, in memory or read from a file, over
    /// its set's geometries.
    pub fn of<S: GeometryStorage>(ensemble: &Ensemble<S>) -> Self {
        Self::bound(ensemble.set().len(), ensemble.distribution().clone())
    }

    /// This description for `geometry_count` geometries whose first
    /// [`Self::geometry_count`] are the geometries it describes, such as a
    /// trajectory file that grew after the description was saved: the
    /// distribution as [`Distribution::extended_to`] gives it, keeping its
    /// log weights, partition, states, and temperature bit for bit.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::CoverageMismatch`] for fewer geometries than the
    /// description describes.
    pub fn extended_to(self, geometry_count: usize) -> Result<Self, StructureSetError> {
        if geometry_count < self.geometry_count {
            return Err(StructureSetError::CoverageMismatch {
                covered: self.geometry_count,
                geometries: geometry_count,
            });
        }
        let distribution = self.distribution.extended_to(geometry_count)?;
        Ok(Self::bound(geometry_count, distribution))
    }

    /// The number of geometries the distribution partitions.
    pub const fn geometry_count(&self) -> usize {
        self.geometry_count
    }

    /// The distribution over the geometries.
    pub fn distribution(&self) -> &Distribution {
        &self.distribution
    }

    /// The distribution, for [`Ensemble::new`] with the set it describes.
    pub fn into_distribution(self) -> Distribution {
        self.distribution
    }
}
