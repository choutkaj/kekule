use std::sync::Arc;

use crate::diagnostics::Cancellation;
use crate::math::RigidTransform;
use crate::properties::{OwnerProperties, TopologyPropertiesMut};
use crate::topology::{AtomSelection, Topology, TopologyPerceptionError};

use super::state::subset_rows;
use super::{DenseArrayError, GeometryError, Model, Positions};

mod context;
mod decoded;
mod distribution;
mod error;
mod geometry;
mod order;
mod read;
mod storage;
mod views;

pub use context::{Membership, SetContext};
pub(crate) use decoded::DecodedGeometry;
pub use distribution::{Distribution, DistributionMut, Members, State};
pub use error::{Requirement, StructureSetError};
pub(crate) use geometry::GeometryState;
pub use geometry::{Geometry, GeometryMut};
pub use order::{Order, Place, TimeAxis};
pub(crate) use read::visit_geometry;
pub use read::{StructureSetRead, StructureSetReadError};
pub(crate) use storage::sealed;
pub use storage::{GeometryStorage, InMemory};
pub(crate) use views::check_bound;
pub use views::{GeometryIter, GeometryView, StateIter, StateView};

/// Geometries of one shared topology: many coordinate sets of one system,
/// with an optional [`Order`] and an optional [`Distribution`].
///
/// The two parts are independent and say what the set is.
///
/// - The order says what the order of the geometries means: a path in time
///   ([`Order::Time`], a trajectory's frames), numbered steps
///   ([`Order::Step`], an optimization's history or a path's images), or the
///   values of a scanned coordinate ([`Order::Coordinate`]). Its values
///   belong to it, one per geometry; a geometry taken out of the set claims no
///   place. Without an order, geometry order carries no meaning.
/// - The distribution says which geometries are drawn from a distribution of
///   conformations, how much each weighs, and how they group into states
///   ([`Distribution`]). It covers every geometry or a stated subset, so the
///   equilibrated part of a trajectory can be its ensemble without a copy.
///
/// A trajectory is a set with a time order and no distribution, a series one
/// with a step or coordinate order, and a conformer ensemble one without an
/// order and with a distribution over every geometry; a clustered trajectory
/// has both. Operations ask for what they need and fail with
/// [`StructureSetError::Missing`] naming it ([`Requirement`]).
///
/// The storage is the type parameter ([`GeometryStorage`]). A set in memory,
/// the default [`InMemory`], lends its geometries by reference, and reading
/// them cannot fail. A set over an indexed trajectory file,
/// `StructureSet<`[`IndexedFile`](crate::traj::io::IndexedFile)`>`, holds no
/// geometry in memory: it is read through [`StructureSetRead`], which returns owned
/// geometries and can fail. Both share the topology, the set properties, the
/// order, the distribution, and the operations on them.
///
/// Each geometry's state is shared storage (see [`Geometry`]), so cloning a
/// set copies its order and distribution but no coordinates. A set held
/// elsewhere, for example behind an [`Arc`], can therefore be regrouped or
/// reweighted through a clone ([`Arc::make_mut`]) without copying
/// coordinates; editing one geometry then copies that geometry alone, and the
/// other set is unaffected. [`Self::select`], [`Self::select_states`], and
/// [`Self::collapse`] share the selected structures in the same way.
///
/// Geometries are validated against the shared topology on insertion, so they
/// can be read as [`super::ModelView`]s without further checks. Geometries are
/// stored in a stable order: regrouping, reweighting, and editing states
/// never reorder them, so geometry indices held elsewhere stay valid.
///
/// ```
/// use kekule::structure::{
///     Distribution, Partition, Members, Order, Place, Positions, Geometry, StructureSet,
///     TimeAxis,
/// };
/// use kekule::units::{Quantity, PICOSECOND};
/// let topology = kekule::smiles::to_topology("O")?;
/// let atoms = topology.atom_count();
/// let mut set = StructureSet::new(topology);
/// let times = Quantity::new(Vec::new(), PICOSECOND);
/// set.set_order(Some(Order::Time(TimeAxis::new(Some(times), None)?)))?;
/// for picoseconds in [0.0, 2.0, 4.0] {
///     let place = Place::Time { time: Some(Quantity::new(picoseconds, PICOSECOND)), step: None };
///     set.push(Geometry::new(Positions::zeros(atoms)), Some(place))?;
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
/// Geometry editors preserve dimensions and cannot overwrite a whole geometry:
///
/// ```compile_fail,E0594
/// use kekule::structure::{Positions, Geometry, StructureSet};
/// fn overwrite(set: &mut StructureSet) {
///     *set.get_mut(0).unwrap() = Geometry::new(Positions::zeros(1));
/// }
/// ```
#[derive(Debug, Clone)]
pub struct StructureSet<S = InMemory> {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    storage: S,
    order: Option<Order>,
    distribution: Option<Distribution>,
}

impl<S: GeometryStorage> StructureSet<S> {
    /// A set of `storage`'s geometries, whose order and distribution, already
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

    /// The topology shared by every geometry.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology; cloning the [`Arc`] copies no
    /// topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Where the geometries are kept.
    pub fn storage(&self) -> &S {
        &self.storage
    }

    /// Set-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable set-level annotations. Geometry properties of a set in memory
    /// are edited through [`StructureSet::get_mut`], state properties through
    /// [`DistributionMut::state_mut`].
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    /// Number of geometries.
    pub fn len(&self) -> usize {
        self.storage.geometry_count()
    }

    /// Whether the set has no geometries.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The order of the geometries, or `None` when the set is unordered.
    pub fn order(&self) -> Option<&Order> {
        self.order.as_ref()
    }

    /// Checks that every time is non-decreasing.
    ///
    /// # Errors
    ///
    /// `Missing(TimeOrder)` for a set without a time order,
    /// `Missing(Times)` for one whose axis has no times, and
    /// [`StructureSetError::NonMonotonicTime`] for the first geometry earlier
    /// than the one before it.
    pub fn validate_monotonic_time(&self) -> Result<(), StructureSetError> {
        let Some(Order::Time(axis)) = &self.order else {
            return Err(Requirement::TimeOrder.into());
        };
        let times = axis.times().ok_or(Requirement::Times)?;
        match times.value().windows(2).position(|pair| pair[1] < pair[0]) {
            Some(previous) => Err(StructureSetError::NonMonotonicTime {
                geometry: previous + 1,
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

    /// Each geometry statistics over the set must see, with its relative
    /// weight: the distribution's members ([`Distribution::weighted_members`])
    /// or, without a distribution, every geometry with weight `1.0`. Items are
    /// geometry indices, in increasing order. The set's own statistics, such
    /// as [`StructureSet::rmsf`], and its clusterings, such as
    /// [`StructureSet::cluster_by_radius`], read these geometries alone.
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

    /// Replaces or clears the distribution. A subset listing every geometry
    /// becomes [`Members::All`]. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::IndexOutOfRange`] for a member outside the set, and
    /// [`StructureSetError::CoverageMismatch`] when a distribution of every
    /// geometry covers another number of geometries.
    pub fn set_distribution(
        &mut self,
        distribution: Option<Distribution>,
    ) -> Result<(), StructureSetError> {
        let len = self.len();
        let distribution = distribution
            .map(|mut distribution| {
                distribution.bind(len)?;
                Ok::<_, StructureSetError>(distribution)
            })
            .transpose()?;
        self.distribution = distribution;
        Ok(())
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every geometry, the
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
    /// layout is kept, so every geometry, selection, buffer, and prepared
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

    /// Takes `order`, checked against the geometries, after the storage grew
    /// from `previous` geometries: a distribution of every geometry becomes one
    /// of the geometries it covered.
    pub(crate) fn grown(&mut self, previous: usize, order: Option<Order>) {
        debug_assert!(order
            .as_ref()
            .is_none_or(|order| order.validate(self.len()).is_ok()));
        self.order = order;
        let len = self.len();
        if len > previous {
            if let Some(distribution) = &mut self.distribution {
                distribution
                    .extend(len)
                    .expect("the distribution describes the geometries before the new ones");
            }
        }
    }
}

impl<S: GeometryStorage> StructureSet<S>
where
    Self: StructureSetRead,
{
    /// One state of the distribution by index, with its representative,
    /// members, and population; None without a distribution or out of
    /// range. The state of a set in memory lends its representative geometry
    /// ([StateView::representative_sample]).
    pub fn state(&self, index: usize) -> Option<StateView<'_, Self>> {
        (index < self.distribution.as_ref()?.state_count()).then(|| StateView::new(self, index))
    }

    /// The states of the distribution in order; empty without one.
    pub fn states(&self) -> StateIter<'_, Self> {
        StateIter::new(self)
    }
}

impl StructureSet {
    /// An empty, unordered set of one shared topology, without a
    /// distribution.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self::from_storage(topology.into(), InMemory::default(), None, None)
    }

    /// Reads the geometries `indices` of `reader` into a set in memory, in the
    /// requested order: a faithful copy. Each picked geometry is kept whole,
    /// with its velocities, forces, occupancies, and B-factors; the set
    /// keeps the reader's order values and its distribution for the picked
    /// geometries, restricted as [`Self::select`] restricts them. Set
    /// properties are not kept; [`StructureSet::load`] keeps a file's.
    ///
    /// Ranges and strides can be written as `(start..end).step_by(stride)`,
    /// and every geometry as `0..reader.len()`. Indices may repeat and come in
    /// any order. A repeated index is read once and its geometry shared.
    /// Geometries picked out of ascending order keep their times, so
    /// [`Self::validate_monotonic_time`] refuses the set when a time
    /// decreases. The reader can be a [`Processed`](crate::traj::Processed)
    /// one, so its steps, such as an atom subset or a fit, run on each
    /// geometry as it is read. A host that needs progress or cancellation
    /// reads geometry by geometry through [`StructureSetRead`] and builds the set
    /// itself.
    ///
    /// An ensemble of trajectory frames, such as one for a structure
    /// format, drops what the frames have and an ensemble does not:
    ///
    /// ```
    /// use kekule::structure::{
    ///     Distribution, Members, Order, Partition, Place, Positions, Geometry, StructureSet,
    ///     TimeAxis,
    /// };
    /// use kekule::units::{Quantity, PICOSECOND};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let mut frames = StructureSet::new(topology.clone());
    /// frames.set_order(Some(Order::Time(TimeAxis::new(
    ///     Some(Quantity::new(Vec::new(), PICOSECOND)),
    ///     None,
    /// )?)))?;
    /// for index in 0..10 {
    ///     let time = Some(Quantity::new(f64::from(index), PICOSECOND));
    ///     let geometry = Geometry::new(Positions::zeros(topology.atom_count()));
    ///     frames.push(geometry, Some(Place::Time { time, step: None }))?;
    /// }
    /// // Leave out two geometries of equilibration, then take every third.
    /// let mut ensemble = StructureSet::from_reader(&frames, (2..10).step_by(3))?;
    /// assert_eq!(ensemble.len(), 3);
    /// ensemble.set_order(None)?;
    /// ensemble.clear_velocities_and_forces();
    /// ensemble.set_distribution(Some(Distribution::new(
    ///     Members::All,
    ///     Partition::one_state(3),
    /// )?))?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`StructureSetReadError::IndexOutOfRange`] for the first index not below
    /// `reader.len()`, before anything is read;
    /// [`StructureSetReadError::Geometry`] for a geometry the reader does not bind
    /// to its topology; and otherwise the first failure to read a geometry.
    pub fn from_reader<R: StructureSetRead + ?Sized>(
        reader: &R,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, StructureSetReadError> {
        let indices = indices.into_iter().collect::<Vec<_>>();
        let len = reader.len();
        if let Some(&index) = indices.iter().find(|&&index| index >= len) {
            return Err(StructureSetReadError::IndexOutOfRange { index, len });
        }
        let read = |index: usize| {
            let geometry = reader.read(index)?;
            check_bound(reader.topology(), &geometry).map_err(|source| {
                StructureSetReadError::Geometry {
                    geometry: index,
                    source: Box::new(source),
                }
            })?;
            Ok(geometry)
        };
        let geometries = if indices.windows(2).all(|pair| pair[0] < pair[1]) {
            indices
                .iter()
                .map(|&index| read(index))
                .collect::<Result<Vec<_>, StructureSetReadError>>()?
        } else {
            // Each index once, then a shared clone wherever it is asked for.
            let mut unique = indices.clone();
            unique.sort_unstable();
            unique.dedup();
            let picked = unique
                .iter()
                .map(|&index| read(index))
                .collect::<Result<Vec<_>, StructureSetReadError>>()?;
            indices
                .iter()
                .map(|index| {
                    let slot = unique.binary_search(index).expect("every index was read");
                    picked[slot].clone()
                })
                .collect()
        };
        Ok(Self::from_storage(
            Arc::clone(reader.topology()),
            InMemory::new(geometries),
            reader.order().map(|order| order.select(&indices)),
            reader
                .distribution()
                .map(|distribution| distribution.select(&indices)),
        ))
    }

    /// Validates and collects each geometry in order, binding its bond
    /// property rows, into an unordered set without a distribution; give it
    /// one with [`Self::set_order`] and [`Self::set_distribution`]. Fails on
    /// the first invalid geometry.
    pub fn from_geometries(
        topology: impl Into<Arc<Topology>>,
        geometries: impl IntoIterator<Item = Geometry>,
    ) -> Result<Self, StructureSetError> {
        let mut set = Self::new(topology);
        for mut geometry in geometries {
            set.bind(&mut geometry)?;
            set.storage.geometries.push(geometry);
        }
        Ok(set)
    }

    /// Collects models of one topology layout as the geometries of an
    /// unordered set without a distribution, on the first model's topology
    /// snapshot. Each model's positions, cell, occupancies, B-factors, and
    /// structure properties become its geometry, moving without copying.
    ///
    /// Every model's topology must have the first one's layout
    /// ([`Topology::same_layout`]), so models read from separate files are
    /// accepted as well as models sharing one published layout. Only the
    /// first topology is kept.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptySource`] without models and
    /// [`StructureSetError::TopologyMismatch`] for a model of another layout.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, StructureSetError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(StructureSetError::EmptySource)?;
        let mut set = Self::new(first.shared_topology());
        set.storage.geometries.push(Geometry::from_model(first));
        for model in models {
            if !model.topology().same_layout(&set.topology) {
                return Err(StructureSetError::TopologyMismatch);
            }
            set.storage.geometries.push(Geometry::from_model(model));
        }
        Ok(set)
    }

    /// One topology-bound geometry by stable index.
    pub fn get(&self, index: usize) -> Option<GeometryView<'_>> {
        self.storage
            .geometries
            .get(index)
            .map(|geometry| GeometryView::of_set(self, index, geometry))
    }

    /// A dimension-preserving editor for one geometry. Whole-geometry
    /// replacement goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<GeometryMut<'_>> {
        self.storage.geometries.get_mut(index).map(GeometryMut::new)
    }

    /// Topology-bound geometries in set order.
    pub fn iter(&self) -> GeometryIter<'_> {
        GeometryIter::new(self)
    }

    /// Replaces or clears the order. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::OrderLengthMismatch`] when the order has another
    /// number of values than the set has geometries, and
    /// [`StructureSetError::NonFiniteOrderValue`] for a coordinate that is not
    /// finite.
    pub fn set_order(&mut self, order: Option<Order>) -> Result<(), StructureSetError> {
        if let Some(order) = &order {
            order.validate(self.storage.geometries.len())?;
        }
        self.order = order;
        Ok(())
    }

    /// Moves one geometry to `place` in the order; see [`Self::push`] for
    /// what a place must state. Failure changes nothing.
    pub fn set_place(
        &mut self,
        index: usize,
        place: Option<Place>,
    ) -> Result<(), StructureSetError> {
        self.check_index(index)?;
        match &mut self.order {
            Some(order) => {
                let stored = order.check(place, index)?;
                order.set(index, stored);
                Ok(())
            }
            None if place.is_none() => Ok(()),
            None => Err(StructureSetError::PlaceMismatch { geometry: index }),
        }
    }

    /// Validates and appends one geometry at `place` in the order, outside
    /// the distribution, returning its index.
    ///
    /// The place must fit the order: none in an unordered set or on a time
    /// axis without columns, [`Place::Time`] with a time and step exactly
    /// when the axis has those columns, [`Place::Step`] on a step order, and
    /// [`Place::Coordinate`] on a coordinate order, converted to its unit. A
    /// distribution of every geometry becomes one of the geometries before this
    /// one. Failure changes nothing.
    pub fn push(
        &mut self,
        mut geometry: Geometry,
        place: Option<Place>,
    ) -> Result<usize, StructureSetError> {
        let index = self.storage.geometries.len();
        let stored = self.check_place(place, index)?;
        self.bind(&mut geometry)?;
        self.append(geometry, stored);
        if let Some(distribution) = &mut self.distribution {
            distribution
                .extend(index + 1)
                .expect("the distribution describes the geometries before the new one");
        }
        Ok(index)
    }

    /// Validates and appends one geometry at `place`, as [`Self::push`] does,
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
    /// [`StructureSetError::StateOutOfRange`], [`StructureSetError::WeightRequired`],
    /// [`StructureSetError::InvalidWeight`], and the errors of [`Self::push`].
    pub fn push_member(
        &mut self,
        mut geometry: Geometry,
        place: Option<Place>,
        state: usize,
        weight: Option<f64>,
    ) -> Result<usize, StructureSetError> {
        let index = self.storage.geometries.len();
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        distribution.check_state(state)?;
        distribution.check_new_weight(weight, index)?;
        let stored = self.check_place(place, index)?;
        self.bind(&mut geometry)?;
        self.append(geometry, stored);
        self.distribution
            .as_mut()
            .expect("the distribution was checked")
            .push_member(index, state, weight);
        Ok(index)
    }

    /// Validates and appends one geometry at `place`, as [`Self::push`] does,
    /// as the representative of a new state of the distribution, returning
    /// the state's index. `weight` is as in [`Self::push_member`]. Failure
    /// changes nothing.
    pub fn push_state(
        &mut self,
        mut geometry: Geometry,
        place: Option<Place>,
        state: State,
        weight: Option<f64>,
    ) -> Result<usize, StructureSetError> {
        let index = self.storage.geometries.len();
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        distribution.check_new_weight(weight, index)?;
        let stored = self.check_place(place, index)?;
        self.bind(&mut geometry)?;
        self.append(geometry, stored);
        Ok(self
            .distribution
            .as_mut()
            .expect("the distribution was checked")
            .push_state(index, state, weight))
    }

    /// Validates and replaces one geometry, returning the previous one. The
    /// geometry keeps its place and membership. The index is checked first;
    /// failure leaves the set unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut geometry: Geometry,
    ) -> Result<Geometry, StructureSetError> {
        self.check_index(index)?;
        self.bind(&mut geometry)?;
        Ok(std::mem::replace(
            &mut self.storage.geometries[index],
            geometry,
        ))
    }

    /// Removes and returns one geometry, with neither its place nor its
    /// membership, shifting later geometries and members down.
    ///
    /// Removing the only member of a state removes the state, shifting later
    /// states down. The representative of a state with other members cannot
    /// be removed ([`StructureSetError::RepresentativeRemoval`]); choose another
    /// representative with [`DistributionMut::regroup`] first.
    pub fn remove(&mut self, index: usize) -> Result<Geometry, StructureSetError> {
        self.check_index(index)?;
        if let Some(distribution) = &self.distribution {
            distribution.check_remove(index)?;
        }
        let geometry = self.storage.geometries.remove(index);
        if let Some(order) = &mut self.order {
            order.remove(index);
        }
        if let Some(distribution) = &mut self.distribution {
            distribution.remove(index, self.storage.geometries.len());
        }
        Ok(geometry)
    }

    /// Replaces the positions of every geometry transactionally, in order.
    /// Other state is kept, and no replaced positions are copied.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), StructureSetError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.storage.geometries.len() {
            return Err(StructureSetError::GeometryCountMismatch {
                expected: self.storage.geometries.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(GeometryError::Array(DenseArrayError::CountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                })
                .into());
            }
        }
        for (geometry, positions) in self.storage.geometries.iter_mut().zip(positions) {
            GeometryState::replace_positions(geometry.structure_mut(), positions);
        }
        Ok(())
    }

    /// Copies geometries in the requested order, sharing the topology and each
    /// geometry's structure. Ranges and strides can be written as
    /// `(start..end).step_by(stride)`; duplicates and reordering are
    /// allowed. Set properties are kept.
    ///
    /// Each geometry keeps its order value, which then need not increase. The
    /// selected members of the distribution stay members with their weights,
    /// and the temperature is kept. Every state with a selected member
    /// becomes a state of the copy, in state order, represented by the first
    /// copy of its representative or, when that is not selected, by its
    /// first selected member. A state whose members are all selected keeps
    /// its energy and properties; one selected only in part has neither, as
    /// a split state in [`DistributionMut::regroup`]. Use [`Self::select_states`]
    /// to keep whole states.
    pub fn select(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, StructureSetError> {
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
    /// state with all of its members in geometry order (sharing their
    /// structures), its representative, energy, and properties. Duplicates
    /// and reordering are allowed. Weights, the temperature, the geometries'
    /// order values, and set properties are kept; geometries outside the
    /// distribution are left out. [`Distribution::covering_population`],
    /// [`Distribution::within_energy_window`], and
    /// [`Distribution::ranked_by_population`] give states to extract.
    ///
    /// # Errors
    ///
    /// `Missing(Distribution)` without a distribution and
    /// [`StructureSetError::StateOutOfRange`] for a state outside it.
    pub fn select_states(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, StructureSetError> {
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        let indices = indices
            .into_iter()
            .map(|index| distribution.check_state(index).map(|()| index))
            .collect::<Result<Vec<_>, _>>()?;
        let (geometries, distribution) = distribution.select_states(&indices);
        Ok(self.copy(
            &geometries,
            self.order.as_ref().map(|order| order.select(&geometries)),
            Some(distribution),
        ))
    }

    /// A new set of one geometry per state, in state order: each state's
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
    /// [`StructureSetError::InvalidWeight`] when a population underflows to
    /// zero.
    pub fn collapse(&self) -> Result<Self, StructureSetError> {
        let distribution = self
            .distribution
            .as_ref()
            .ok_or(Requirement::Distribution)?;
        let (geometries, distribution) = distribution.collapse()?;
        Ok(self.copy(&geometries, None, Some(distribution)))
    }

    /// Publishes one induced topology subset and projects every geometry onto
    /// it. Set, state, and per-geometry owner properties are dropped; the
    /// order and the distribution's members, weights, states, energies, and
    /// temperature still describe the geometries and are kept.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, StructureSetError> {
        let subset = Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let geometries = self
            .storage
            .geometries
            .iter()
            .map(|geometry| {
                Ok(Geometry::from_structure(Arc::new(
                    geometry.structure().project(&atoms, &bonds)?,
                )))
            })
            .collect::<Result<Vec<_>, GeometryError>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            storage: InMemory::new(geometries),
            order: self.order.clone(),
            distribution: self
                .distribution
                .as_ref()
                .map(Distribution::without_state_properties),
        })
    }

    /// Removes the velocities and forces of every geometry. A geometry without
    /// them stays shared.
    pub fn clear_velocities_and_forces(&mut self) {
        for geometry in &mut self.storage.geometries {
            geometry.take_velocities();
            geometry.take_forces();
        }
    }

    /// Consumes the set, returning its geometries in order and dropping the
    /// topology, the order, the distribution, and the set properties.
    pub fn into_geometries(self) -> Vec<Geometry> {
        self.storage.geometries
    }

    /// Applies one rigid transform per geometry transactionally: positions
    /// move; cells, velocities, and forces rotate.
    pub(crate) fn transform_samples(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), GeometryError> {
        debug_assert_eq!(transforms.len(), self.storage.geometries.len());
        let geometries = self
            .storage
            .geometries
            .iter()
            .zip(transforms)
            .map(|(geometry, transform)| {
                Ok(Geometry::from_structure(Arc::new(
                    geometry.structure().transformed(*transform)?,
                )))
            })
            .collect::<Result<Vec<_>, GeometryError>>()?;
        self.storage.geometries = geometries;
        Ok(())
    }

    /// A copy of the geometries `indices` (checked) with an order and a
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
                    .map(|index| self.storage.geometries[*index].clone())
                    .collect(),
            ),
            order,
            distribution,
        }
    }

    /// Checks `place` for a new or moved geometry `geometry`.
    fn check_place(
        &self,
        place: Option<Place>,
        geometry: usize,
    ) -> Result<Option<order::Stored>, StructureSetError> {
        match &self.order {
            Some(order) => order.check(place, geometry).map(Some),
            None if place.is_none() => Ok(None),
            None => Err(StructureSetError::PlaceMismatch { geometry }),
        }
    }

    /// Appends a bound geometry at a checked place.
    fn append(&mut self, geometry: Geometry, place: Option<order::Stored>) {
        self.storage.geometries.push(geometry);
        if let (Some(order), Some(place)) = (&mut self.order, place) {
            order.push(place);
        }
    }

    fn check_index(&self, index: usize) -> Result<(), StructureSetError> {
        if index >= self.storage.geometries.len() {
            return Err(StructureSetError::IndexOutOfRange {
                index,
                len: self.storage.geometries.len(),
            });
        }
        Ok(())
    }

    fn bind(&self, geometry: &mut Geometry) -> Result<(), StructureSetError> {
        Ok(GeometryState::bind(
            geometry.structure_mut(),
            &self.topology,
        )?)
    }
}

impl<'a> IntoIterator for &'a StructureSet {
    type Item = GeometryView<'a>;
    type IntoIter = GeometryIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl StructureSetRead for StructureSet {
    fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn len(&self) -> usize {
        self.storage.geometries.len()
    }

    fn order(&self) -> Option<&Order> {
        self.order.as_ref()
    }

    fn distribution(&self) -> Option<&Distribution> {
        self.distribution.as_ref()
    }

    fn properties(&self) -> Option<&OwnerProperties> {
        Some(&self.properties)
    }

    fn peek(&self, index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        self.get(index).map(GeometryView::as_dyn)
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        self.storage
            .geometries
            .get(index)
            .cloned()
            .ok_or(StructureSetReadError::IndexOutOfRange {
                index,
                len: self.storage.geometries.len(),
            })
    }
}
