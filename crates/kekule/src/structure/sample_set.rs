use std::sync::Arc;

use crate::diagnostics::Cancellation;
use crate::geometry::RigidTransform;
use crate::properties::{OwnerProperties, TopologyPropertiesMut};
use crate::topology::{AtomSelection, Topology, TopologyPerceptionError};

use super::state::subset_rows;
use super::{DenseArrayError, Model, Positions, StructureError};

mod context;
mod decoded;
mod distribution;
mod error;
mod order;
mod read;
mod sample;
mod storage;
mod views;

pub use context::{Membership, SampleContext};
pub(crate) use decoded::DecodedSample;
pub use distribution::{Distribution, DistributionMut, Members, State};
pub use error::{Requirement, SampleSetError};
pub use order::{Order, Place, TimeAxis};
pub(crate) use read::visit_sample;
pub use read::{SampleRead, SampleReadError};
pub(crate) use sample::SampleStructure;
pub use sample::{Sample, SampleMut};
pub(crate) use storage::sealed;
pub use storage::{InMemory, SampleStorage};
pub(crate) use views::check_bound;
pub use views::{SampleIter, SampleView, StateIter, StateView};

/// Samples of one shared topology: many coordinate sets of one system,
/// with an optional [`Order`] and an optional [`Distribution`].
///
/// The two parts are independent and say what the set is.
///
/// - The order says what the order of the samples means: a path in time
///   ([`Order::Time`], a trajectory's frames), numbered steps
///   ([`Order::Step`], an optimization's history or a path's images), or the
///   values of a scanned coordinate ([`Order::Coordinate`]). Its values
///   belong to it, one per sample; a sample taken out of the set claims no
///   place. Without an order, sample order carries no meaning.
/// - The distribution says which samples are drawn from a distribution of
///   conformations, how much each weighs, and how they group into states
///   ([`Distribution`]). It covers every sample or a stated subset, so the
///   equilibrated part of a trajectory can be its ensemble without a copy.
///
/// A trajectory is a set with a time order and no distribution, a series one
/// with a step or coordinate order, and a conformer ensemble one without an
/// order and with a distribution over every sample; a clustered trajectory
/// has both. Operations ask for what they need and fail with
/// [`SampleSetError::Missing`] naming it ([`Requirement`]).
///
/// The storage is the type parameter ([`SampleStorage`]). A set in memory,
/// the default [`InMemory`], lends its samples by reference, and reading
/// them cannot fail. A set over an indexed trajectory file,
/// `SampleSet<`[`IndexedFile`](crate::traj::io::IndexedFile)`>`, holds no
/// sample in memory: it is read through [`SampleRead`], which returns owned
/// samples and can fail. Both share the topology, the set properties, the
/// order, the distribution, and the operations on them.
///
/// Each sample's structure is shared storage (see [`Sample`]), so cloning a
/// set copies its order and distribution but no coordinates. A set held
/// elsewhere, for example behind an [`Arc`], can therefore be regrouped or
/// reweighted through a clone ([`Arc::make_mut`]) without copying
/// coordinates; editing one sample then copies that sample alone, and the
/// other set is unaffected. [`Self::select`], [`Self::select_states`], and
/// [`Self::collapse`] share the selected structures in the same way.
///
/// Samples are validated against the shared topology on insertion, so they
/// can be read as [`super::ModelView`]s without further checks. Samples are
/// stored in a stable order: regrouping, reweighting, and editing states
/// never reorder them, so sample indices held elsewhere stay valid.
///
/// ```
/// use kekule::structure::{
///     Distribution, Partition, Members, Order, Place, Positions, Sample, SampleSet,
///     TimeAxis,
/// };
/// use kekule::units::{Quantity, PICOSECOND};
/// let topology = kekule::smiles::to_topology("O")?;
/// let atoms = topology.atom_count();
/// let mut set = SampleSet::new(topology);
/// let times = Quantity::new(Vec::new(), PICOSECOND);
/// set.set_order(Some(Order::Time(TimeAxis::new(Some(times), None)?)))?;
/// for picoseconds in [0.0, 2.0, 4.0] {
///     let place = Place::Time { time: Some(Quantity::new(picoseconds, PICOSECOND)), step: None };
///     set.push(Sample::new(Positions::zeros(atoms)), Some(place))?;
/// }
/// // The last two frames are the equilibrated ensemble, without a copy.
/// let equilibrated = Distribution::new(
///     Members::Subset(vec![1, 2]),
///     Partition::one_state(2),
/// )?;
/// set.set_distribution(Some(equilibrated))?;
/// assert_eq!(set.get(0).unwrap().weight(), None);
/// assert_eq!(set.get(2).unwrap().weight(), Some(1.0));
/// assert_eq!(
///     set.get(2).unwrap().place(),
///     Some(Place::Time { time: Some(Quantity::new(4.0, PICOSECOND)), step: None }),
/// );
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// Sample editors preserve dimensions and cannot overwrite a whole sample:
///
/// ```compile_fail,E0594
/// use kekule::structure::{Positions, Sample, SampleSet};
/// fn overwrite(set: &mut SampleSet) {
///     *set.get_mut(0).unwrap() = Sample::new(Positions::zeros(1));
/// }
/// ```
#[derive(Debug, Clone)]
pub struct SampleSet<S = InMemory> {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    storage: S,
    order: Option<Order>,
    distribution: Option<Distribution>,
}

impl<S: SampleStorage> SampleSet<S> {
    /// A set of `storage`'s samples, whose order and distribution, already
    /// checked against them, are given.
    pub(crate) fn from_storage(
        topology: Arc<Topology>,
        storage: S,
        order: Option<Order>,
        distribution: Option<Distribution>,
    ) -> Self {
        Self {
            topology,
            properties: OwnerProperties::new(),
            storage,
            order,
            distribution,
        }
    }

    /// The topology shared by every sample.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology; cloning the [`Arc`] copies no
    /// topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Where the samples are kept.
    pub fn storage(&self) -> &S {
        &self.storage
    }

    /// Set-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable set-level annotations. Sample properties of a set in memory
    /// are edited through [`SampleSet::get_mut`], state properties through
    /// [`DistributionMut::state_mut`].
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    /// Number of samples.
    pub fn len(&self) -> usize {
        self.storage.sample_count()
    }

    /// Whether the set has no samples.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The order of the samples, or `None` when the set is unordered.
    pub fn order(&self) -> Option<&Order> {
        self.order.as_ref()
    }

    /// Checks that every time is non-decreasing.
    ///
    /// # Errors
    ///
    /// `Missing(TimeOrder)` for a set without a time order,
    /// `Missing(Times)` for one whose axis has no times, and
    /// [`SampleSetError::NonMonotonicTime`] for the first sample earlier
    /// than the one before it.
    pub fn validate_monotonic_time(&self) -> Result<(), SampleSetError> {
        let Some(Order::Time(axis)) = &self.order else {
            return Err(Requirement::TimeOrder.into());
        };
        let times = axis.times().ok_or(Requirement::Times)?;
        match times.value().windows(2).position(|pair| pair[1] < pair[0]) {
            Some(previous) => Err(SampleSetError::NonMonotonicTime {
                sample: previous + 1,
            }),
            None => Ok(()),
        }
    }

    /// The distribution, or `None` when the set has none.
    pub fn distribution(&self) -> Option<&Distribution> {
        self.distribution.as_ref()
    }

    /// The distribution for reweighting, regrouping, and editing states,
    /// through a guard that keeps its members ([`DistributionMut`]).
    pub fn distribution_mut(&mut self) -> Option<DistributionMut<'_>> {
        self.distribution.as_mut().map(DistributionMut::new)
    }

    /// Each sample statistics over the set must see, with its relative
    /// weight: the distribution's members ([`Distribution::weighted_members`])
    /// or, without a distribution, every sample with weight `1.0`. Items are
    /// sample indices, in increasing order. The set's own statistics, such
    /// as [`SampleSet::rmsf`], and its clusterings, such as
    /// [`SampleSet::cluster_by_radius`], read these samples alone.
    pub fn weighted_samples(&self) -> impl ExactSizeIterator<Item = (usize, f64)> + '_ {
        let distribution = self.distribution.as_ref();
        let count = distribution.map_or(self.len(), Distribution::member_count);
        (0..count).map(move |position| {
            distribution.map_or((position, 1.0), |distribution| {
                (
                    distribution.member(position),
                    distribution.member_weight(position),
                )
            })
        })
    }

    /// Replaces or clears the distribution. A subset listing every sample
    /// becomes [`Members::All`]. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`SampleSetError::IndexOutOfRange`] for a member outside the set, and
    /// [`SampleSetError::CoverageMismatch`] when a distribution of every
    /// sample covers another number of samples.
    pub fn set_distribution(
        &mut self,
        distribution: Option<Distribution>,
    ) -> Result<(), SampleSetError> {
        let len = self.len();
        let distribution = distribution
            .map(|mut distribution| {
                distribution.bind(len)?;
                Ok::<_, SampleSetError>(distribution)
            })
            .transpose()?;
        self.distribution = distribution;
        Ok(())
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every sample, the
    /// order, the distribution, and the set properties without copying.
    /// Failure leaves the set unchanged. The new snapshot shares the
    /// original layout, so selections, buffers, and prepared potentials stay
    /// usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// [`Self::perceive`], stopping partway once `cancellation` is requested.
    ///
    /// Delegates to [`Topology::perceived_cancellable`]. A cancelled call
    /// leaves the set unchanged; a completed one installs exactly what
    /// [`Self::perceive`] would.
    ///
    /// # Errors
    ///
    /// As [`Topology::perceived_cancellable`].
    pub fn perceive_cancellable(
        &mut self,
        cancellation: &Cancellation,
    ) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived_cancellable(cancellation)?);
        Ok(())
    }

    /// Length-preserving mutable access to the shared topology's static
    /// annotations.
    ///
    /// A topology shared with other owners is first replaced by a snapshot of
    /// it ([`Arc::make_mut`]), so those owners are unaffected. Either way the
    /// layout is kept, so every sample, selection, buffer, and prepared
    /// potential remains usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// The storage, for changing what it holds.
    pub(crate) fn storage_mut(&mut self) -> &mut S {
        &mut self.storage
    }

    /// The shared topology handle, with the set's lifetime.
    pub(crate) fn topology_handle(&self) -> &Arc<Topology> {
        &self.topology
    }

    /// Takes `order`, checked against the samples, after the storage grew
    /// from `previous` samples: a distribution of every sample becomes one
    /// of the samples it covered.
    pub(crate) fn grown(&mut self, previous: usize, order: Option<Order>) {
        debug_assert!(order
            .as_ref()
            .is_none_or(|order| order.validate(self.len()).is_ok()));
        self.order = order;
        if self.len() > previous {
            if let Some(distribution) = &mut self.distribution {
                distribution.push_non_member(previous);
            }
        }
    }

    /// The same set over `storage`, which holds the same samples: the
    /// topology, the set properties, the order, and the distribution are
    /// kept.
    pub(crate) fn with_storage<T: SampleStorage>(&self, storage: T) -> SampleSet<T> {
        debug_assert_eq!(storage.sample_count(), self.len());
        SampleSet {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            storage,
            order: self.order.clone(),
            distribution: self.distribution.clone(),
        }
    }
}

impl<S: SampleStorage> SampleSet<S>
where
    Self: SampleRead,
{
    /// One state of the distribution by index, with its representative,
    /// members, and population; None without a distribution or out of
    /// range. The state of a set in memory lends its representative sample
    /// ([StateView::representative_sample]).
    pub fn state(&self, index: usize) -> Option<StateView<'_, Self>> {
        (index < self.distribution.as_ref()?.state_count()).then(|| StateView::new(self, index))
    }

    /// The states of the distribution in order; empty without one.
    pub fn states(&self) -> StateIter<'_, Self> {
        StateIter::new(self)
    }
}

impl SampleSet {
    /// An empty, unordered set of one shared topology, without a
    /// distribution.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self::from_storage(topology.into(), InMemory::default(), None, None)
    }

    /// Reads the samples `indices` of `reader` into a set in memory, in the
    /// requested order: each with its order value, and the distribution
    /// restricted to them, as [`Self::select`] restricts it. The set has no
    /// set properties; [`SampleSet::load`] keeps a file's.
    ///
    /// Ranges and strides can be written as `(start..end).step_by(stride)`,
    /// and every sample as `0..reader.len()`; duplicates and reordering are
    /// allowed. The reader can be a [`Processed`](crate::traj::Processed)
    /// one, so its steps, such as an atom subset or a fit, run on each
    /// sample as it is read. Each sample keeps its velocities and forces.
    /// The host then sets the order and the distribution it wants with
    /// [`Self::set_order`] and [`Self::set_distribution`]:
    ///
    /// ```
    /// use kekule::structure::{
    ///     Distribution, Members, Partition, Positions, Sample, SampleSet,
    /// };
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let mut frames = SampleSet::new(topology.clone());
    /// for _ in 0..10 {
    ///     frames.push(Sample::new(Positions::zeros(topology.atom_count())), None)?;
    /// }
    /// // Leave out two samples of equilibration, then take every third.
    /// let mut ensemble = SampleSet::from_reader(&frames, (2..10).step_by(3))?;
    /// assert_eq!(ensemble.len(), 3);
    /// ensemble.set_order(None)?;
    /// ensemble.set_distribution(Some(Distribution::new(
    ///     Members::All,
    ///     Partition::one_state(3),
    /// )?))?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`SampleReadError::IndexOutOfRange`] for an index past the last
    /// sample, [`SampleReadError::Structure`] for a sample the reader does
    /// not bind to its topology, and the first failure to read a sample.
    pub fn from_reader<R: SampleRead + ?Sized>(
        reader: &R,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, SampleReadError> {
        let indices = indices.into_iter().collect::<Vec<_>>();
        let samples = indices
            .iter()
            .map(|&index| {
                let sample = reader.read(index)?;
                check_bound(reader.topology(), &sample).map_err(|source| {
                    SampleReadError::Structure {
                        sample: index,
                        source: Box::new(source),
                    }
                })?;
                Ok(sample)
            })
            .collect::<Result<Vec<_>, SampleReadError>>()?;
        Ok(Self::from_storage(
            Arc::clone(reader.topology()),
            InMemory::new(samples),
            reader.order().map(|order| order.select(&indices)),
            reader
                .distribution()
                .map(|distribution| distribution.select(&indices)),
        ))
    }

    /// Validates and collects each sample in order, binding its bond
    /// property rows, into an unordered set without a distribution; give it
    /// one with [`Self::set_order`] and [`Self::set_distribution`]. Fails on
    /// the first invalid sample.
    pub fn from_samples(
        topology: impl Into<Arc<Topology>>,
        samples: impl IntoIterator<Item = Sample>,
    ) -> Result<Self, SampleSetError> {
        let mut set = Self::new(topology);
        for mut sample in samples {
            set.bind(&mut sample)?;
            set.storage.samples.push(sample);
        }
        Ok(set)
    }

    /// Collects models of one topology layout as the samples of an
    /// unordered set without a distribution, on the first model's topology
    /// snapshot. Each model's positions, cell, occupancies, B-factors, and
    /// structure properties become its sample, moving without copying.
    ///
    /// Every model's topology must have the first one's layout
    /// ([`Topology::same_layout`]), so models read from separate files are
    /// accepted as well as models sharing one published layout. Only the
    /// first topology is kept.
    ///
    /// # Errors
    ///
    /// [`SampleSetError::EmptySource`] without models and
    /// [`SampleSetError::TopologyMismatch`] for a model of another layout.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, SampleSetError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(SampleSetError::EmptySource)?;
        let mut set = Self::new(first.shared_topology());
        set.storage.samples.push(Sample::from_model(first));
        for model in models {
            if !model.topology().same_layout(&set.topology) {
                return Err(SampleSetError::TopologyMismatch);
            }
            set.storage.samples.push(Sample::from_model(model));
        }
        Ok(set)
    }

    /// One topology-bound sample by stable index.
    pub fn get(&self, index: usize) -> Option<SampleView<'_>> {
        self.storage
            .samples
            .get(index)
            .map(|sample| SampleView::of_set(self, index, sample))
    }

    /// A dimension-preserving editor for one sample. Whole-sample
    /// replacement goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<SampleMut<'_>> {
        self.storage.samples.get_mut(index).map(SampleMut::new)
    }

    /// Topology-bound samples in set order.
    pub fn iter(&self) -> SampleIter<'_> {
        SampleIter::new(self)
    }

    /// Replaces or clears the order. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`SampleSetError::OrderLengthMismatch`] when the order has another
    /// number of values than the set has samples, and
    /// [`SampleSetError::NonFiniteOrderValue`] for a coordinate that is not
    /// finite.
    pub fn set_order(&mut self, order: Option<Order>) -> Result<(), SampleSetError> {
        if let Some(order) = &order {
            order.validate(self.storage.samples.len())?;
        }
        self.order = order;
        Ok(())
    }

    /// Moves one sample to `place` in the order; see [`Self::push`] for
    /// what a place must state. Failure changes nothing.
    pub fn set_place(&mut self, index: usize, place: Option<Place>) -> Result<(), SampleSetError> {
        self.check_index(index)?;
        match &mut self.order {
            Some(order) => {
                let stored = order.check(place, index)?;
                order.set(index, stored);
                Ok(())
            }
            None if place.is_none() => Ok(()),
            None => Err(SampleSetError::PlaceMismatch { sample: index }),
        }
    }

    /// Validates and appends one sample at `place` in the order, outside
    /// the distribution, returning its index.
    ///
    /// The place must fit the order: none in an unordered set or on a time
    /// axis without columns, [`Place::Time`] with a time and step exactly
    /// when the axis has those columns, [`Place::Step`] on a step order, and
    /// [`Place::Coordinate`] on a coordinate order, converted to its unit. A
    /// distribution of every sample becomes one of the samples before this
    /// one. Failure changes nothing.
    pub fn push(
        &mut self,
        mut sample: Sample,
        place: Option<Place>,
    ) -> Result<usize, SampleSetError> {
        let index = self.storage.samples.len();
        let stored = self.check_place(place, index)?;
        self.bind(&mut sample)?;
        self.append(sample, stored);
        if let Some(distribution) = &mut self.distribution {
            distribution.push_non_member(index);
        }
        Ok(index)
    }

    /// Validates and appends one sample at `place`, as [`Self::push`] does,
    /// as a member of an existing state of the distribution, returning its
    /// index.
    ///
    /// `weight` is the new member's finite, positive weight; it is required
    /// when the distribution has weights, and in one without, a weight gives
    /// the other members `1.0`. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// `Missing(Distribution)` without a distribution,
    /// [`SampleSetError::StateOutOfRange`], [`SampleSetError::WeightRequired`],
    /// [`SampleSetError::InvalidWeight`], and the errors of [`Self::push`].
    pub fn push_member(
        &mut self,
        mut sample: Sample,
        place: Option<Place>,
        state: usize,
        weight: Option<f64>,
    ) -> Result<usize, SampleSetError> {
        let index = self.storage.samples.len();
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        distribution.check_state(state)?;
        distribution.check_new_weight(weight, index)?;
        let stored = self.check_place(place, index)?;
        self.bind(&mut sample)?;
        self.append(sample, stored);
        self.distribution
            .as_mut()
            .expect("the distribution was checked")
            .push_member(index, state, weight);
        Ok(index)
    }

    /// Validates and appends one sample at `place`, as [`Self::push`] does,
    /// as the representative of a new state of the distribution, returning
    /// the state's index. `weight` is as in [`Self::push_member`]. Failure
    /// changes nothing.
    pub fn push_state(
        &mut self,
        mut sample: Sample,
        place: Option<Place>,
        state: State,
        weight: Option<f64>,
    ) -> Result<usize, SampleSetError> {
        let index = self.storage.samples.len();
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        distribution.check_new_weight(weight, index)?;
        let stored = self.check_place(place, index)?;
        self.bind(&mut sample)?;
        self.append(sample, stored);
        Ok(self
            .distribution
            .as_mut()
            .expect("the distribution was checked")
            .push_state(index, state, weight))
    }

    /// Validates and replaces one sample, returning the previous one. The
    /// sample keeps its place and membership. The index is checked first;
    /// failure leaves the set unchanged.
    pub fn replace(&mut self, index: usize, mut sample: Sample) -> Result<Sample, SampleSetError> {
        self.check_index(index)?;
        self.bind(&mut sample)?;
        Ok(std::mem::replace(&mut self.storage.samples[index], sample))
    }

    /// Removes and returns one sample, with neither its place nor its
    /// membership, shifting later samples and members down.
    ///
    /// Removing the only member of a state removes the state, shifting later
    /// states down. The representative of a state with other members cannot
    /// be removed ([`SampleSetError::RepresentativeRemoval`]); choose another
    /// representative with [`DistributionMut::regroup`] first.
    pub fn remove(&mut self, index: usize) -> Result<Sample, SampleSetError> {
        self.check_index(index)?;
        if let Some(distribution) = &self.distribution {
            distribution.check_remove(index)?;
        }
        let sample = self.storage.samples.remove(index);
        if let Some(order) = &mut self.order {
            order.remove(index);
        }
        if let Some(distribution) = &mut self.distribution {
            distribution.remove(index, self.storage.samples.len());
        }
        Ok(sample)
    }

    /// Replaces the positions of every sample transactionally, in order.
    /// Other state is kept, and no replaced positions are copied.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), SampleSetError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.storage.samples.len() {
            return Err(SampleSetError::SampleCountMismatch {
                expected: self.storage.samples.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(StructureError::Array(DenseArrayError::CountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                })
                .into());
            }
        }
        for (sample, positions) in self.storage.samples.iter_mut().zip(positions) {
            SampleStructure::replace_positions(sample.structure_mut(), positions);
        }
        Ok(())
    }

    /// Copies samples in the requested order, sharing the topology and each
    /// sample's structure. Ranges and strides can be written as
    /// `(start..end).step_by(stride)`; duplicates and reordering are
    /// allowed. Set properties are kept.
    ///
    /// Each sample keeps its order value, which then need not increase. The
    /// selected members of the distribution stay members with their weights,
    /// and the temperature is kept. Every state with a selected member
    /// becomes a state of the copy, in state order, represented by the first
    /// copy of its representative or, when that is not selected, by its
    /// first selected member. A state whose members are all selected keeps
    /// its energy and properties; one selected only in part has neither, as
    /// a split state in [`DistributionMut::regroup`]. Use [`Self::select_states`]
    /// to keep whole states.
    pub fn select(&self, indices: impl IntoIterator<Item = usize>) -> Result<Self, SampleSetError> {
        let indices = indices
            .into_iter()
            .map(|index| self.check_index(index).map(|()| index))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(self.copy(
            &indices,
            self.order.as_ref().map(|order| order.select(&indices)),
            self.distribution
                .as_ref()
                .map(|distribution| distribution.select(&indices)),
        ))
    }

    /// Copies the members of whole states in the requested order, each
    /// state with all of its members in sample order (sharing their
    /// structures), its representative, energy, and properties. Duplicates
    /// and reordering are allowed. Weights, the temperature, the samples'
    /// order values, and set properties are kept; samples outside the
    /// distribution are left out. [`Distribution::covering_population`],
    /// [`Distribution::within_energy_window`], and
    /// [`Distribution::ranked_by_population`] give states to extract.
    ///
    /// # Errors
    ///
    /// `Missing(Distribution)` without a distribution and
    /// [`SampleSetError::StateOutOfRange`] for a state outside it.
    pub fn select_states(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, SampleSetError> {
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        let indices = indices
            .into_iter()
            .map(|index| distribution.check_state(index).map(|()| index))
            .collect::<Result<Vec<_>, _>>()?;
        let (samples, distribution) = distribution.select_states(&indices);
        Ok(self.copy(
            &samples,
            self.order.as_ref().map(|order| order.select(&samples)),
            Some(distribution),
        ))
    }

    /// A new set of one sample per state, in state order: each state's
    /// representative, sharing its structure, weighted by the state's
    /// population relative to the largest member weight, so population
    /// shares are kept.
    ///
    /// It coarse-grains, losing the spread within each state, and drops the
    /// order. Each state keeps its energy and properties; the topology, the
    /// set properties, and the temperature are kept, and the source is
    /// unchanged.
    ///
    /// # Errors
    ///
    /// `Missing(Distribution)` without a distribution, and
    /// [`SampleSetError::InvalidWeight`] when a population underflows to
    /// zero.
    pub fn collapse(&self) -> Result<Self, SampleSetError> {
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        let (samples, distribution) = distribution.collapse()?;
        Ok(self.copy(&samples, None, Some(distribution)))
    }

    /// Publishes one induced topology subset and projects every sample onto
    /// it. Set, state, and per-sample owner properties are dropped; the
    /// order and the distribution's members, weights, states, energies, and
    /// temperature still describe the samples and are kept.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, SampleSetError> {
        let subset = Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let samples = self
            .storage
            .samples
            .iter()
            .map(|sample| {
                Ok(Sample::from_structure(Arc::new(
                    sample.structure().project(&atoms, &bonds)?,
                )))
            })
            .collect::<Result<Vec<_>, StructureError>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            storage: InMemory::new(samples),
            order: self.order.clone(),
            distribution: self
                .distribution
                .as_ref()
                .map(Distribution::without_state_properties),
        })
    }

    /// Consumes the set, returning its samples in order and dropping the
    /// topology, the order, the distribution, and the set properties.
    pub fn into_samples(self) -> Vec<Sample> {
        self.storage.samples
    }

    /// Applies one rigid transform per sample transactionally: positions
    /// move; cells, velocities, and forces rotate.
    pub(crate) fn transform_samples(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), StructureError> {
        debug_assert_eq!(transforms.len(), self.storage.samples.len());
        let samples = self
            .storage
            .samples
            .iter()
            .zip(transforms)
            .map(|(sample, transform)| {
                Ok(Sample::from_structure(Arc::new(
                    sample.structure().transformed(*transform)?,
                )))
            })
            .collect::<Result<Vec<_>, StructureError>>()?;
        self.storage.samples = samples;
        Ok(())
    }

    /// A copy of the samples `indices` (checked) with an order and a
    /// distribution of that copy.
    fn copy(
        &self,
        indices: &[usize],
        order: Option<Order>,
        distribution: Option<Distribution>,
    ) -> Self {
        Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            storage: InMemory::new(
                indices
                    .iter()
                    .map(|index| self.storage.samples[*index].clone())
                    .collect(),
            ),
            order,
            distribution,
        }
    }

    /// Checks `place` for a new or moved sample `sample`.
    fn check_place(
        &self,
        place: Option<Place>,
        sample: usize,
    ) -> Result<Option<order::Stored>, SampleSetError> {
        match &self.order {
            Some(order) => order.check(place, sample).map(Some),
            None if place.is_none() => Ok(None),
            None => Err(SampleSetError::PlaceMismatch { sample }),
        }
    }

    /// Appends a bound sample at a checked place.
    fn append(&mut self, sample: Sample, place: Option<order::Stored>) {
        self.storage.samples.push(sample);
        if let (Some(order), Some(place)) = (&mut self.order, place) {
            order.push(place);
        }
    }

    fn check_index(&self, index: usize) -> Result<(), SampleSetError> {
        if index >= self.storage.samples.len() {
            return Err(SampleSetError::IndexOutOfRange {
                index,
                len: self.storage.samples.len(),
            });
        }
        Ok(())
    }

    fn bind(&self, sample: &mut Sample) -> Result<(), SampleSetError> {
        Ok(SampleStructure::bind(
            sample.structure_mut(),
            &self.topology,
        )?)
    }
}

impl<'a> IntoIterator for &'a SampleSet {
    type Item = SampleView<'a>;
    type IntoIter = SampleIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl SampleRead for SampleSet {
    fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn len(&self) -> usize {
        self.storage.samples.len()
    }

    fn order(&self) -> Option<&Order> {
        self.order.as_ref()
    }

    fn distribution(&self) -> Option<&Distribution> {
        self.distribution.as_ref()
    }

    fn peek(&self, index: usize) -> Option<SampleView<'_, dyn SampleRead>> {
        self.get(index).map(SampleView::as_dyn)
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        self.storage
            .samples
            .get(index)
            .cloned()
            .ok_or(SampleReadError::IndexOutOfRange {
                index,
                len: self.storage.samples.len(),
            })
    }
}
