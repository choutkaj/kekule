use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch, TopologyMismatch};
use std::sync::Arc;

use crate::math::RigidTransform;
use crate::parallel::Run;
use crate::properties::{OwnerProperties, PropertyTable, TopologyPropertiesMut};
use crate::topology::{AtomSelection, Topology, TopologyPerceptionError};

use super::subset_rows;
use super::{
    DenseArrayError, Geometry, GeometryData, GeometryError, Positions, Sampling, Structure,
};

mod columns;
mod context;
mod error;
mod order;
mod read;
mod storage;
mod views;

pub(crate) use columns::move_owner_values;
pub use columns::{ColumnRef, GeometryIndex, RowView, SetColumnsMut};
pub use context::SetContext;
pub use error::{Requirement, StructureSetError};
pub use order::{Order, Place, TimeAxis};
pub(crate) use read::visit_geometry;
pub use read::{StructureSetRead, StructureSetReadError};
pub(crate) use storage::sealed;
pub use storage::{GeometryStorage, InMemory};
pub(crate) use views::check_bound;
pub use views::{GeometryIter, GeometryView, SetGeometryMut};

/// Geometries of one shared topology: many coordinate sets of one system,
/// with an optional [`Order`].
///
/// The order says what the order of the geometries means: a path in time
/// ([`Order::Time`], a trajectory's frames), numbered steps ([`Order::Step`],
/// an optimization's history or a path's images), or the values of a scanned
/// coordinate ([`Order::Coordinate`]). Its values belong to it, one per
/// geometry; a geometry taken out of the set claims no place. Without an
/// order, geometry order carries no meaning. Per-geometry scalars, such as
/// an energy or a collective variable per frame, are the set's
/// [`columns`](Self::columns), one row per geometry, kept with the
/// geometries as the order is. They have no other home: a geometry entering
/// a set in memory moves its owner values into its row, and a geometry taken
/// out ([`GeometryView::to_structure`], [`Self::remove`],
/// [`Self::into_geometries`]) carries its row's cells as owner values again.
/// An order is recorded by the
/// source, so it belongs to the set; a distribution over the geometries is
/// inferred by an analysis, so it belongs to an
/// [`Ensemble`](crate::structure::Ensemble) that reads the set. Operations
/// ask for what they need and fail with [`StructureSetError::Missing`]
/// naming it ([`Requirement`]).
///
/// The storage is the type parameter ([`GeometryStorage`]). A set in memory,
/// the default [`InMemory`], lends its geometries by reference, and reading
/// them cannot fail. A set over an indexed trajectory file,
/// `StructureSet<`[`IndexedFile`](crate::traj::io::IndexedFile)`>`, holds no
/// geometry in memory: it is read through [`StructureSetRead`], which returns owned
/// geometries and can fail. Both share the topology, the set properties, the
/// order, the columns, and the operations on them. A geometry read from a
/// file keeps its codec's per-frame values, such as a TRR frame's lambda, as
/// its own owner values; they are not columns until
/// [`Self::from_reader`] moves them there.
///
/// Each geometry's state is shared storage (see [`Geometry`]), so cloning a
/// set copies its order but no coordinates; editing one geometry then copies
/// that geometry alone, and the other set is unaffected. [`Self::select`]
/// shares the selected geometries in the same way.
///
/// Geometries are validated against the shared topology on insertion, so they
/// can be read as [`super::StructureView`]s without further checks. Geometries are
/// stored in a stable order: editing them never reorders them, so geometry
/// indices held elsewhere stay valid.
///
/// ```
/// use kekule::structure::{Order, Place, Positions, Geometry, StructureSet, TimeAxis};
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
    columns: PropertyTable<GeometryIndex>,
    sampling: Option<Sampling>,
}

impl<S: GeometryStorage> StructureSet<S> {
    /// A set of `storage`'s geometries, whose order, already checked against
    /// them, is given.
    pub(crate) fn from_storage(topology: Arc<Topology>, storage: S, order: Option<Order>) -> Self {
        let columns = PropertyTable::new(storage.geometry_count());
        Self {
            topology,
            properties: OwnerProperties::new(),
            storage,
            order,
            columns,
            sampling: None,
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
    /// are edited through [`StructureSet::get_mut`].
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

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every geometry, the
    /// order, and the set properties without copying.
    /// Failure leaves the set unchanged. The new snapshot shares the
    /// original layout, so selections, buffers, and prepared potentials stay
    /// usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// [`Self::perceive`] under a [`Run`], through
    /// [`Topology::perceived_with`].
    ///
    /// # Errors
    ///
    /// As [`Topology::perceived_with`].
    pub fn perceive_with(&mut self, run: Run<'_>) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived_with(run)?);
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

    /// Takes `columns`, one row per geometry, which keep the sampling record's
    /// columns real energies.
    pub(crate) fn set_columns(&mut self, columns: PropertyTable<GeometryIndex>) {
        debug_assert_eq!(columns.len(), self.len());
        debug_assert!(self
            .sampling
            .as_ref()
            .is_none_or(|sampling| sampling.check_columns(&columns).is_ok()));
        self.columns = columns;
    }

    /// The sampling record: the state the geometries sample, from which step,
    /// and which columns hold their potential energy and bias. `None` when
    /// nothing records it.
    ///
    /// The record belongs to the set, as the columns do: selecting, loading
    /// through a reader, an atom subset, a [`Processed`](crate::traj::Processed)
    /// reader, and refreshing a growing file keep it.
    pub fn sampling(&self) -> Option<&Sampling> {
        self.sampling.as_ref()
    }

    /// Replaces or clears the sampling record. Failure changes nothing.
    ///
    /// A column the record names must be a real column of energies, or
    /// absent: a column without values, which reads as missing energies.
    ///
    /// ```
    /// use kekule::properties::{PropertyColumn, PropertyKey};
    /// use kekule::structure::{Geometry, Positions, Sampling, StructureSet, ThermodynamicState};
    /// use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let geometry = || Geometry::new(Positions::zeros(topology.atom_count()));
    /// let mut set = StructureSet::from_geometries(topology.clone(), [geometry(), geometry()])?;
    /// let bias = PropertyKey::new("bias")?;
    /// set.columns_mut().insert(
    ///     bias.clone(),
    ///     PropertyColumn::Real { unit: KILOJOULE_PER_MOLE, values: vec![Some(1.5), Some(0.25)] },
    /// )?;
    /// let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN))?;
    /// set.set_sampling(Some(Sampling::new(state).with_bias([bias.clone()])?))?;
    /// assert_eq!(set.sampling().unwrap().bias(), [bias.clone()]);
    /// // The record's columns stay.
    /// assert!(set.columns_mut().remove(&bias).is_err());
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`StructureSetError::SamplingColumn`] for a named column that holds
    /// other values than real energies.
    pub fn set_sampling(&mut self, sampling: Option<Sampling>) -> Result<(), StructureSetError> {
        if let Some(sampling) = &sampling {
            sampling.check_columns(&self.columns)?;
        }
        self.sampling = sampling;
        Ok(())
    }

    /// Takes `order`, checked against the geometries, after the storage
    /// grew; the new geometries' cells are missing.
    pub(crate) fn grown(&mut self, order: Option<Order>) {
        debug_assert!(order
            .as_ref()
            .is_none_or(|order| order.validate(self.len()).is_ok()));
        self.order = order;
        self.columns.resize_missing(self.len());
    }
}

impl StructureSet {
    /// An empty, unordered set of one shared topology.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self::from_storage(topology.into(), InMemory::default(), None)
    }

    /// Reads the geometries `indices` of `reader` into a set in memory, in the
    /// requested order: a faithful copy. Each picked geometry is kept whole,
    /// with its velocities, forces, occupancies, and B-factors; the set
    /// keeps the reader's order values and column rows for the picked
    /// geometries, and the reader's set properties and sampling record, when
    /// it has them. A geometry's own owner values, such as a trajectory
    /// codec's per-frame values, move into its row as [`Self::push`] moves
    /// them.
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
    ///     Ensemble, Order, Place, Positions, Geometry, StructureSet, TimeAxis,
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
    /// let mut samples = StructureSet::from_reader(&frames, (2..10).step_by(3))?;
    /// assert_eq!(samples.len(), 3);
    /// samples.set_order(None)?;
    /// samples.clear_velocities_and_forces();
    /// let ensemble = Ensemble::uniform(samples);
    /// assert_eq!(ensemble.distribution().sample_count(), 3);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`StructureSetReadError::IndexOutOfRange`] for the first index not below
    /// `reader.len()`, before anything is read;
    /// [`StructureSetReadError::ColumnRows`] for a reader whose columns do not
    /// have one row per geometry, also before anything is read;
    /// [`StructureSetReadError::Geometry`] for a geometry the reader does not bind
    /// to its topology; [`StructureSetReadError::ColumnTypeMismatch`] for an
    /// owner value that does not fit its column; and otherwise the first
    /// failure to read a geometry.
    pub fn from_reader<R: StructureSetRead + ?Sized>(
        reader: &R,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, StructureSetReadError> {
        let indices = indices.into_iter().collect::<Vec<_>>();
        let len = reader.len();
        if let Some(&index) = indices.iter().find(|&&index| index >= len) {
            return Err(StructureSetReadError::IndexOutOfRange(
                IndexOutOfRange::new(Countable::Geometry, index, len),
            ));
        }
        if let Some(columns) = reader.columns().filter(|columns| columns.len() != len) {
            return Err(StructureSetReadError::ColumnRows {
                rows: columns.len(),
                geometries: len,
            });
        }
        // Each geometry with its owner values taken out, for its row.
        let read = |index: usize| {
            let mut geometry = reader.read(index)?;
            check_bound(reader.topology(), &geometry).map_err(|source| {
                StructureSetReadError::Geometry {
                    geometry: index,
                    source,
                }
            })?;
            let owner = geometry.take_owner_values();
            Ok((geometry, owner))
        };
        let picked = if indices.windows(2).all(|pair| pair[0] < pair[1]) {
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
        let mut columns = match reader.columns() {
            Some(columns) => columns
                .select_rows(&indices)
                .expect("the columns have a row per geometry"),
            None => PropertyTable::new(indices.len()),
        };
        let mut geometries = Vec::with_capacity(picked.len());
        for (row, ((geometry, owner), &index)) in picked.into_iter().zip(&indices).enumerate() {
            move_owner_values(&mut columns, row, &owner).map_err(|error| match error {
                StructureSetError::ColumnTypeMismatch { key, .. } => {
                    StructureSetReadError::ColumnTypeMismatch {
                        key,
                        geometry: index,
                    }
                }
                StructureSetError::Column(error) => StructureSetReadError::Geometry {
                    geometry: index,
                    source: GeometryError::Property(error),
                },
                _ => unreachable!("moving owner values fails on a column alone"),
            })?;
            geometries.push(geometry);
        }
        let mut set = Self::from_storage(
            Arc::clone(reader.topology()),
            InMemory::new(geometries),
            reader.order().map(|order| order.select(&indices)),
        );
        if let Some(properties) = reader.properties() {
            set.properties_mut().clone_from(properties);
        }
        set.columns = columns;
        set.sampling = reader.sampling().cloned();
        Ok(set)
    }

    /// Validates and collects each geometry in order, binding its bond
    /// property rows, into an unordered set; give it an order with
    /// [`Self::set_order`]. Each geometry's owner values move into its row
    /// of the columns, as [`Self::push`] moves them. Fails on the first
    /// invalid geometry.
    pub fn from_geometries(
        topology: impl Into<Arc<Topology>>,
        geometries: impl IntoIterator<Item = Geometry>,
    ) -> Result<Self, StructureSetError> {
        let mut set = Self::new(topology);
        for geometry in geometries {
            set.push(geometry, None)?;
        }
        Ok(set)
    }

    /// Collects structures of one topology layout as the geometries of an
    /// unordered set, on the first structure's topology snapshot. Each
    /// structure's positions, cell, occupancies, B-factors, and atom and
    /// bond properties become its geometry, moving without copying; its
    /// owner values move into its row of the columns, as [`Self::push`]
    /// moves them.
    ///
    /// Every structure's topology must have the first one's layout
    /// ([`Topology::same_layout`]), so structures read from separate files are
    /// accepted as well as structures sharing one published layout. Only the
    /// first topology is kept.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::EmptySource`] without structures,
    /// [`StructureSetError::TopologyMismatch`] for a structure of another
    /// layout, and [`StructureSetError::ColumnTypeMismatch`] for an owner
    /// value that does not fit the column an earlier structure made.
    pub fn from_structures(
        structures: impl IntoIterator<Item = Structure>,
    ) -> Result<Self, StructureSetError> {
        let mut structures = structures.into_iter();
        let first = structures.next().ok_or(StructureSetError::EmptySource)?;
        let mut set = Self::new(first.shared_topology());
        set.take(Geometry::from_structure(first))?;
        for structure in structures {
            if !structure.topology().same_layout(&set.topology) {
                return Err(StructureSetError::TopologyMismatch(
                    TopologyMismatch::Structure,
                ));
            }
            set.take(Geometry::from_structure(structure))?;
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

    /// A dimension-preserving editor for one geometry. It edits everything
    /// but owner values, which are the geometry's row of the
    /// [`columns`](Self::columns_mut). Whole-geometry replacement goes
    /// through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<SetGeometryMut<'_>> {
        self.storage
            .geometries
            .get_mut(index)
            .map(SetGeometryMut::new)
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

    /// Validates and appends one geometry at `place` in the order, returning
    /// its index.
    ///
    /// The geometry's owner values move into its row of the set's columns:
    /// a new key adds a column whose earlier cells are missing, a real value
    /// converts to its column's unit, and the row's other cells are missing.
    /// The geometry in the set has no owner values.
    ///
    /// The place must fit the order: none in an unordered set or on a time
    /// axis without columns, [`Place::Time`] with a time and step exactly
    /// when the axis has those columns, [`Place::Step`] on a step order, and
    /// [`Place::Coordinate`] on a coordinate order, converted to its unit.
    /// Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::ColumnTypeMismatch`] for an owner value of
    /// another type or unit than its column, besides the errors of a place
    /// or geometry that does not fit.
    pub fn push(
        &mut self,
        geometry: Geometry,
        place: Option<Place>,
    ) -> Result<usize, StructureSetError> {
        self.push_with_row(geometry, place, None)
    }

    /// Appends `geometry` at `place` as [`Self::push`] does, `row`'s cells
    /// joining its owner values in its new row and replacing values of the
    /// same key. The geometry stays shared, as a frame's payload is.
    pub(crate) fn push_with_row(
        &mut self,
        mut geometry: Geometry,
        place: Option<Place>,
        row: Option<RowView<'_>>,
    ) -> Result<usize, StructureSetError> {
        let index = self.storage.geometries.len();
        let stored = self.check_place(place, index)?;
        self.bind(&mut geometry)?;
        self.append_with_row(geometry, stored, row)?;
        Ok(index)
    }

    /// Validates and replaces one geometry, returning the previous one. The
    /// geometry keeps its place and its column row; its owner values move
    /// into the row as [`Self::push`] moves them, replacing the cells of
    /// their keys; a cell whose key the geometry lacks stays. The previous
    /// geometry carries the row's cells from before as owner values. The index is checked first; failure leaves the set
    /// unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut geometry: Geometry,
    ) -> Result<Geometry, StructureSetError> {
        self.check_index(index)?;
        self.bind(&mut geometry)?;
        let row = self.row_properties(index);
        let owner = geometry.take_owner_values();
        move_owner_values(&mut self.columns, index, &owner)?;
        let previous = std::mem::replace(&mut self.storage.geometries[index], geometry);
        Ok(previous.with_owner_values(row))
    }

    /// Removes and returns one geometry with its place and column row,
    /// shifting later geometries down. The geometry carries its row's cells
    /// as owner values.
    pub fn remove(&mut self, index: usize) -> Result<Geometry, StructureSetError> {
        self.check_index(index)?;
        let row = self.row_properties(index);
        let geometry = self.storage.geometries.remove(index);
        if let Some(order) = &mut self.order {
            order.remove(index);
        }
        self.columns.remove_row(index);
        Ok(geometry.with_owner_values(row))
    }

    /// Replaces the positions of every geometry transactionally, in order.
    /// Other state is kept, and no replaced positions are copied.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), StructureSetError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.storage.geometries.len() {
            return Err(StructureSetError::GeometryCountMismatch(
                LengthMismatch::new(
                    Countable::Geometry,
                    self.storage.geometries.len(),
                    positions.len(),
                ),
            ));
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(
                    GeometryError::Array(Box::new(DenseArrayError::CountMismatch(
                        LengthMismatch::new(
                            Countable::Row,
                            self.topology.atom_count(),
                            values.len(),
                        ),
                    )))
                    .into(),
                );
            }
        }
        for (geometry, positions) in self.storage.geometries.iter_mut().zip(positions) {
            GeometryData::replace_positions(geometry.data_mut(), positions);
        }
        Ok(())
    }

    /// Copies geometries in the requested order, sharing the topology and each
    /// geometry's state. Ranges and strides can be written as
    /// `(start..end).step_by(stride)`; duplicates and reordering are
    /// allowed. Set properties and the sampling record are kept.
    ///
    /// Each geometry keeps its order value, which then need not increase,
    /// and its column row.
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
        ))
    }

    /// Publishes one induced topology subset and projects every geometry onto
    /// it. Set properties are dropped; the order, the columns, and the
    /// sampling record still describe the geometries and are kept.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, StructureSetError> {
        let subset = Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let geometries = self
            .storage
            .geometries
            .iter()
            .map(|geometry| {
                Ok(Geometry::from_data(
                    geometry.data().project(&atoms, &bonds)?,
                ))
            })
            .collect::<Result<Vec<_>, GeometryError>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            storage: InMemory::new(geometries),
            order: self.order.clone(),
            columns: self.columns.clone(),
            sampling: self.sampling.clone(),
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
    /// topology, the order, the sampling record, and the set properties.
    /// Each geometry carries its row's cells as owner values.
    pub fn into_geometries(self) -> Vec<Geometry> {
        let columns = self.columns;
        self.storage
            .geometries
            .into_iter()
            .enumerate()
            .map(|(index, geometry)| geometry.with_row(RowView::new(&columns, index)))
            .collect()
    }

    /// Applies one rigid transform per geometry transactionally: positions
    /// move; cells, velocities, and forces rotate.
    pub(crate) fn transform_geometries(
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
                Ok(Geometry::from_data(
                    geometry.data().transformed(*transform)?,
                ))
            })
            .collect::<Result<Vec<_>, GeometryError>>()?;
        self.storage.geometries = geometries;
        Ok(())
    }

    /// A copy of the geometries `indices` (checked) with an order of that
    /// copy.
    fn copy(&self, indices: &[usize], order: Option<Order>) -> Self {
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
            columns: self
                .columns
                .select_rows(indices)
                .expect("the indices are checked"),
            sampling: self.sampling.clone(),
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

    /// Appends a bound geometry at a checked place, moving its owner values
    /// into its new row. Failure changes nothing.
    fn append(
        &mut self,
        geometry: Geometry,
        place: Option<order::Stored>,
    ) -> Result<(), StructureSetError> {
        self.append_with_row(geometry, place, None)
    }

    fn append_with_row(
        &mut self,
        mut geometry: Geometry,
        place: Option<order::Stored>,
        row: Option<RowView<'_>>,
    ) -> Result<(), StructureSetError> {
        let index = self.storage.geometries.len();
        let mut owner = geometry.take_owner_values();
        for (key, value) in row.into_iter().flat_map(RowView::iter) {
            owner
                .insert(key.clone(), value.to_value())
                .expect("stored cells are finite");
        }
        self.columns.resize_missing(index + 1);
        if let Err(error) = move_owner_values(&mut self.columns, index, &owner) {
            self.columns.remove_row(index);
            return Err(error);
        }
        self.storage.geometries.push(geometry);
        if let (Some(order), Some(place)) = (&mut self.order, place) {
            order.push(place);
        }
        Ok(())
    }

    /// Appends a bound geometry to an unordered set.
    fn take(&mut self, geometry: Geometry) -> Result<(), StructureSetError> {
        debug_assert!(self.order.is_none());
        self.append(geometry, None)
    }

    /// The cells of geometry `index`'s row, as owner values.
    fn row_properties(&self, index: usize) -> OwnerProperties {
        RowView::new(&self.columns, index).map_or_else(OwnerProperties::new, RowView::to_properties)
    }

    fn check_index(&self, index: usize) -> Result<(), StructureSetError> {
        if index >= self.storage.geometries.len() {
            return Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
                Countable::Geometry,
                index,
                self.storage.geometries.len(),
            )));
        }
        Ok(())
    }

    fn bind(&self, geometry: &mut Geometry) -> Result<(), StructureSetError> {
        Ok(GeometryData::bind(geometry.data_mut(), &self.topology)?)
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

    fn properties(&self) -> Option<&OwnerProperties> {
        Some(&self.properties)
    }

    fn columns(&self) -> Option<&PropertyTable<GeometryIndex>> {
        Some(&self.columns)
    }

    fn sampling(&self) -> Option<&Sampling> {
        self.sampling.as_ref()
    }

    fn peek(&self, index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        self.get(index).map(GeometryView::as_dyn)
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        self.storage
            .geometries
            .get(index)
            .cloned()
            .ok_or(StructureSetReadError::IndexOutOfRange(
                IndexOutOfRange::new(Countable::Geometry, index, self.storage.geometries.len()),
            ))
    }
}
