//! What a native sample set record states besides its samples.

use std::sync::Arc;

use crate::properties::OwnerProperties;
use crate::structure::{Distribution, Order, SampleSet, SampleSetError};
use crate::topology::Topology;

/// A sample set without its samples: the topology, set properties, order,
/// and distribution of a set whose samples are stored elsewhere, such as
/// in a trajectory file, and the number of them.
///
/// A native sample set record whose samples are external stores this.
/// Writing one ([`NativeSource`](super::NativeSource) converts from it)
/// gives such a record, and reading it gives
/// [`NativeObject::SampleSetDescription`](super::NativeObject::SampleSetDescription).
/// The order and distribution describe the external samples by index, so
/// they are checked against the sample count as a set checks them against
/// its samples.
///
/// ```
/// use kekule::native::{self, NativeReadOptions, SampleSetDescription};
/// use kekule::structure::{Order, TimeAxis};
/// use kekule::units::{Quantity, PICOSECOND};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let times = Quantity::new(vec![0.0, 2.0, 4.0], PICOSECOND);
/// let description = SampleSetDescription::new(topology, 3)
///     .with_order(Some(Order::Time(TimeAxis::new(Some(times), None)?)))?;
/// let bytes = native::write(&description, Default::default())?;
/// let read = native::read(bytes.as_slice(), NativeReadOptions::default())?;
/// let read = read.into_sample_set_description()?;
/// assert_eq!(read.sample_count(), 3);
/// assert_eq!(read.order(), description.order());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct SampleSetDescription {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    sample_count: usize,
    order: Option<Order>,
    distribution: Option<Distribution>,
}

impl SampleSetDescription {
    /// A description of `sample_count` samples of `topology`, unordered,
    /// without a distribution or set properties.
    pub fn new(topology: impl Into<Arc<Topology>>, sample_count: usize) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            sample_count,
            order: None,
            distribution: None,
        }
    }

    /// The description of `set`: its topology, set properties, order,
    /// distribution, and number of samples.
    pub fn of(set: &SampleSet) -> Self {
        Self {
            topology: set.shared_topology(),
            properties: set.properties().clone(),
            sample_count: set.len(),
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
    /// The errors of [`SampleSet::set_order`], against the sample count.
    pub fn with_order(mut self, order: Option<Order>) -> Result<Self, SampleSetError> {
        if let Some(order) = &order {
            order.validate(self.sample_count)?;
        }
        self.order = order;
        Ok(self)
    }

    /// Returns this description with a distribution, or without one. A
    /// subset listing every sample becomes
    /// [`Members::All`](crate::structure::Members::All).
    ///
    /// # Errors
    ///
    /// The errors of [`SampleSet::set_distribution`], against the sample
    /// count.
    pub fn with_distribution(
        mut self,
        distribution: Option<Distribution>,
    ) -> Result<Self, SampleSetError> {
        self.distribution = distribution
            .map(|mut distribution| {
                distribution.bind(self.sample_count)?;
                Ok::<_, SampleSetError>(distribution)
            })
            .transpose()?;
        Ok(self)
    }

    /// The topology every sample shares.
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

    /// The number of samples described.
    pub const fn sample_count(&self) -> usize {
        self.sample_count
    }

    /// The order of the samples, or `None` when they are unordered.
    pub fn order(&self) -> Option<&Order> {
        self.order.as_ref()
    }

    /// The distribution over the samples, or `None` without one.
    pub fn distribution(&self) -> Option<&Distribution> {
        self.distribution.as_ref()
    }
}
