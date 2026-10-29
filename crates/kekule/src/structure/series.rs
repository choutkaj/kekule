use crate::diagnostics::Cancellation;
use crate::structure::DenseArrayError;
use std::fmt;
use std::sync::Arc;

use crate::geometry::{PeriodicCell, RigidTransform};
use crate::properties::{OwnerProperties, StructureProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{AtomSelection, Topology, TopologyPerceptionError};
use crate::units::{Quantity, Unit, UnitError};

use super::scalars::ABSENT;
use super::state::{
    structure_mutators, structure_readers, subset_rows, transform_cell, transform_positions,
    validate_binding,
};
use super::{AsModelView, Model, ModelView, Positions, StructureError, StructureKind};

/// What orders the items of a [`Series`]: steps, or the values of one
/// coordinate.
///
/// Every item of a series lies on its axis: at a step of a
/// [`SeriesAxis::Step`] series, at a coordinate value of a
/// [`SeriesAxis::Coordinate`] series.
///
/// ```
/// use kekule::structure::SeriesAxis;
/// use kekule::units::DEGREE;
/// // A geometry optimization's accepted iterations, and a dihedral scan.
/// let history = SeriesAxis::Step;
/// let scan = SeriesAxis::Coordinate(DEGREE);
/// assert_ne!(history, scan);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum SeriesAxis {
    /// Items numbered by steps: the accepted iterations of a geometry
    /// optimization, or the images of a path. Each item has a step.
    Step,
    /// Items placed by the value of one coordinate, stored in this unit:
    /// the scanned distance, angle, or dihedral of a relaxed scan, or a
    /// reaction coordinate. Each item has a coordinate, which may be
    /// dimensionless, as a normalized path length is.
    Coordinate(Unit),
}

impl fmt::Display for SeriesAxis {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Step => formatter.write_str("steps"),
            Self::Coordinate(unit) => write!(formatter, "a coordinate in {unit}"),
        }
    }
}

/// An ordered series of structures of one shared topology that is neither a
/// path in time nor a distribution: the steps of a geometry optimization,
/// the points of a relaxed scan, or the images of a reaction path.
///
/// Item order is the series order, and every item lies on the series'
/// [`SeriesAxis`]: at a step, or at a value of one coordinate in the axis
/// unit. The values need not increase, as a scan swept back and forth shows,
/// and selection keeps each item's own. An item holds positions, an optional
/// periodic cell, and structure properties: scalar values such as an
/// energy or a gradient norm are its owner properties, keyed and united as
/// the host chooses, and per-atom values such as partial charges are atom
/// columns. A series carries no time, which belongs to a
/// [`super::Trajectory`], and no weights, which belong to an
/// [`super::Ensemble`]; items carry no occupancies or B factors.
///
/// Items are validated against the shared topology on insertion, so they can
/// be read as [`ModelView`]s without further checks; a model converts to an
/// item through [`Self::from_models`] and an item to a model through
/// [`SeriesItemView::to_model`]. Superposition and RMSD are inherent methods;
/// see [`crate::alignment`].
///
/// ```
/// use kekule::properties::{PropertyKey, PropertyValue};
/// use kekule::structure::{Positions, Series, SeriesAxis, SeriesItem};
/// use kekule::units::KILOJOULE_PER_MOLE;
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let atoms = topology.atom_count();
/// let energy = PropertyKey::new("energy")?;
/// let mut history = Series::new(topology, SeriesAxis::Step);
/// for (step, value) in [(0, -10.0), (1, -12.5), (2, -12.75)] {
///     let mut item = SeriesItem::at_step(Positions::zeros(atoms), step);
///     item.properties_mut()
///         .owner_mut()
///         .insert(energy.clone(), PropertyValue::real(value, KILOJOULE_PER_MOLE)?)?;
///     history.push(item)?;
/// }
/// let last = history.get(2).unwrap();
/// assert_eq!(last.step(), Some(2));
/// assert!(last.properties().owner().get(&energy).is_some());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct Series {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    axis: SeriesAxis,
    items: Vec<SeriesItem>,
}

impl Series {
    /// An empty series of one shared topology on `axis`.
    pub fn new(topology: impl Into<Arc<Topology>>, axis: SeriesAxis) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            axis,
            items: Vec::new(),
        }
    }

    /// Validates and collects each item in order, as [`Self::push`] does.
    /// Fails on the first invalid item.
    pub fn from_items(
        topology: impl Into<Arc<Topology>>,
        axis: SeriesAxis,
        items: impl IntoIterator<Item = SeriesItem>,
    ) -> Result<Self, SeriesError> {
        let mut series = Self::new(topology, axis);
        for item in items {
            series.push(item)?;
        }
        Ok(series)
    }

    /// Collects models of one topology layout, in order, as the items of a
    /// [`SeriesAxis::Step`] series, model `i` at step `i`, on the first
    /// model's topology snapshot.
    ///
    /// Each model's positions, cell, and structure properties become its
    /// item, moving without copying; its occupancies and B factors, which
    /// series items do not carry, are dropped. Every model's topology must
    /// have the first one's layout ([`Topology::same_layout`]), so models
    /// read from separate files are accepted, and only the first topology is
    /// kept. Give the items other steps with [`Self::get_mut`], or build
    /// coordinate items with [`SeriesItem::at_coordinate`].
    ///
    /// # Errors
    ///
    /// [`SeriesError::EmptySource`] without models and
    /// [`SeriesError::TopologyMismatch`] for a model of another layout.
    pub fn from_models(models: impl IntoIterator<Item = Model>) -> Result<Self, SeriesError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(SeriesError::EmptySource)?;
        let mut series = Self::new(first.shared_topology(), SeriesAxis::Step);
        series.items.push(SeriesItem::from_model(first, 0));
        for (step, model) in (1..).zip(models) {
            if !model.topology().same_layout(&series.topology) {
                return Err(SeriesError::TopologyMismatch);
            }
            series.items.push(SeriesItem::from_model(model, step));
        }
        Ok(series)
    }

    /// The topology shared by every item.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology; cloning the [`Arc`] copies no topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// The axis the items lie on.
    pub const fn axis(&self) -> SeriesAxis {
        self.axis
    }

    /// Series-level annotations, such as what the axis coordinate measures
    /// or the method that computed the items.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable series-level annotations. Item properties are edited through
    /// [`Self::get_mut`].
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    /// Number of items.
    pub fn len(&self) -> usize {
        self.items.len()
    }

    /// Whether the series has no items.
    pub fn is_empty(&self) -> bool {
        self.items.is_empty()
    }

    /// One topology-bound item by index.
    pub fn get(&self, index: usize) -> Option<SeriesItemView<'_>> {
        (index < self.items.len()).then_some(SeriesItemView {
            series: self,
            index,
        })
    }

    /// A dimension-preserving editor for one item. Whole-item replacement
    /// goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<SeriesItemMut<'_>> {
        let axis = self.axis;
        self.items
            .get_mut(index)
            .map(|item| SeriesItemMut { item, axis })
    }

    /// Topology-bound items in series order.
    pub fn iter(&self) -> SeriesIter<'_> {
        SeriesIter {
            series: self,
            indices: 0..self.items.len(),
        }
    }

    /// Validates and appends one item, binding its bond property rows.
    ///
    /// The item must lie on the series' axis: a step item on a step axis, a
    /// coordinate item on a coordinate axis, its coordinate converted to the
    /// axis unit.
    ///
    /// # Errors
    ///
    /// [`SeriesError::AxisMismatch`] for an item of the other axis,
    /// [`SeriesError::Unit`] for a coordinate of another dimension than the
    /// axis unit, and [`SeriesError::Structure`] for positions or property
    /// rows that do not fit the topology. Failure changes nothing.
    pub fn push(&mut self, mut item: SeriesItem) -> Result<(), SeriesError> {
        self.bind(&mut item)?;
        self.items.push(item);
        Ok(())
    }

    /// Validates and replaces one item, returning the previous one. The
    /// index is checked first; failure leaves the series unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut item: SeriesItem,
    ) -> Result<SeriesItem, SeriesError> {
        self.check_index(index)?;
        self.bind(&mut item)?;
        Ok(std::mem::replace(&mut self.items[index], item))
    }

    /// Removes and returns one item, shifting later ones down.
    pub fn remove(&mut self, index: usize) -> Result<SeriesItem, SeriesError> {
        self.check_index(index)?;
        Ok(self.items.remove(index))
    }

    /// Replaces the positions of every item transactionally, in order.
    /// Other state is kept.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), SeriesError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.items.len() {
            return Err(SeriesError::ItemCountMismatch {
                expected: self.items.len(),
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
        for (item, positions) in self.items.iter_mut().zip(positions) {
            item.positions = positions;
        }
        Ok(())
    }

    /// Copies items in the requested order, sharing the topology. Ranges and
    /// strides can be written as `(start..end).step_by(stride)`; duplicates
    /// and reordering are allowed, and each item keeps its step or
    /// coordinate. The axis and series properties are kept.
    pub fn select(&self, indices: impl IntoIterator<Item = usize>) -> Result<Self, SeriesError> {
        let items = indices
            .into_iter()
            .map(|index| {
                self.check_index(index)?;
                Ok(self.items[index].clone())
            })
            .collect::<Result<Vec<_>, SeriesError>>()?;
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            axis: self.axis,
            items,
        })
    }

    /// Publishes one induced topology subset and projects every item onto
    /// it. Each item keeps its step or coordinate; series and per-item owner
    /// properties, which describe the whole structure, are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, SeriesError> {
        let subset = crate::topology::Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let items = self
            .items
            .iter()
            .map(|item| item.project(&atoms, &bonds))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            axis: self.axis,
            items,
        })
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every item and the
    /// series properties without copying. Failure leaves the series
    /// unchanged. The new snapshot shares the original layout, so
    /// selections, buffers, and prepared potentials stay usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// [`Self::perceive`], stopping partway once `cancellation` is requested.
    ///
    /// Delegates to [`Topology::perceived_cancellable`]. A cancelled call
    /// leaves the series unchanged; a completed one installs exactly what
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
    /// layout is kept, so every item, selection, buffer, and prepared
    /// potential remains usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// Consumes the series, transferring the topology, series properties,
    /// axis, and every item without copying.
    pub fn into_parts(self) -> (Arc<Topology>, OwnerProperties, SeriesAxis, Vec<SeriesItem>) {
        (self.topology, self.properties, self.axis, self.items)
    }

    /// Consumes the series, returning its items in order and dropping the
    /// topology, axis, and series properties.
    pub fn into_items(self) -> Vec<SeriesItem> {
        self.items
    }

    /// Applies one rigid transform per item transactionally.
    pub(crate) fn transform_items(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), StructureError> {
        debug_assert_eq!(transforms.len(), self.items.len());
        let items = self
            .items
            .iter()
            .zip(transforms)
            .map(|(item, transform)| item.transformed(*transform))
            .collect::<Result<Vec<_>, _>>()?;
        self.items = items;
        Ok(())
    }

    fn check_index(&self, index: usize) -> Result<(), SeriesError> {
        if index >= self.items.len() {
            return Err(SeriesError::IndexOutOfRange {
                index,
                len: self.items.len(),
            });
        }
        Ok(())
    }

    /// Checks `item` against the topology and the axis, then converts its
    /// coordinate to the axis unit and binds its bond rows.
    fn bind(&self, item: &mut SeriesItem) -> Result<(), SeriesError> {
        validate_binding(item.atom_count(), &item.properties, &self.topology)?;
        let place = self.axis.place(item.place)?;
        item.place = place;
        item.properties.bind_bond_rows(self.topology.bond_count());
        Ok(())
    }
}

impl SeriesAxis {
    /// `place` on this axis: its coordinate in the axis unit.
    fn place(self, place: Place) -> Result<Place, SeriesError> {
        match (self, place) {
            (Self::Step, Place::Step(step)) => Ok(Place::Step(step)),
            (Self::Coordinate(unit), Place::Coordinate(coordinate)) => {
                Ok(Place::Coordinate(checked_coordinate(coordinate, unit)?))
            }
            _ => Err(SeriesError::AxisMismatch { axis: self }),
        }
    }
}

impl<'a> IntoIterator for &'a Series {
    type Item = SeriesItemView<'a>;
    type IntoIter = SeriesIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Where an item lies on its series' axis.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Place {
    Step(u64),
    /// Finite, in the axis unit once bound.
    Coordinate(Quantity<f64>),
}

/// A finite coordinate, converted to `unit`.
fn checked_coordinate(coordinate: Quantity<f64>, unit: Unit) -> Result<Quantity<f64>, SeriesError> {
    let coordinate = coordinate.into_unit(unit)?;
    if !coordinate.value().is_finite() {
        return Err(StructureError::NonFiniteCoordinate.into());
    }
    Ok(coordinate)
}

/// One series item: positions in dense topology atom order, an optional
/// periodic cell, structure properties, and its place on the series' axis,
/// a step or a coordinate value.
///
/// Items carry no time, velocities, forces, weights, occupancies, or
/// B-factors. The item carries no topology; it is bound when inserted into a
/// [`Series`], which checks that it lies on the series' axis.
#[derive(Debug, Clone, PartialEq)]
pub struct SeriesItem {
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: StructureProperties,
    place: Place,
}

impl SeriesItem {
    /// A non-periodic item without properties at `step`, for a
    /// [`SeriesAxis::Step`] series.
    pub fn at_step(positions: Positions, step: u64) -> Self {
        let properties = StructureProperties::new(positions.len());
        Self {
            positions,
            cell: None,
            properties,
            place: Place::Step(step),
        }
    }

    /// A non-periodic item without properties at a finite `coordinate`, for
    /// a [`SeriesAxis::Coordinate`] series, whose unit it converts to when
    /// inserted.
    ///
    /// # Errors
    ///
    /// [`StructureError::NonFiniteCoordinate`] when the coordinate is not
    /// finite.
    pub fn at_coordinate(
        positions: Positions,
        coordinate: Quantity<f64>,
    ) -> Result<Self, StructureError> {
        if !coordinate.value().is_finite() {
            return Err(StructureError::NonFiniteCoordinate);
        }
        let properties = StructureProperties::new(positions.len());
        Ok(Self {
            positions,
            cell: None,
            properties,
            place: Place::Coordinate(coordinate),
        })
    }

    /// An item holding a model's positions, cell, and properties at `step`.
    fn from_model(model: Model, step: u64) -> Self {
        Self {
            positions: model.positions,
            cell: model.cell,
            properties: model.properties,
            place: Place::Step(step),
        }
    }

    /// Returns this item with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    /// Number of atoms, the length of every dense per-atom array of the item.
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// The item's step, or `None` for an item placed by a coordinate.
    pub const fn step(&self) -> Option<u64> {
        match self.place {
            Place::Step(step) => Some(step),
            Place::Coordinate(_) => None,
        }
    }

    /// The item's coordinate, in the axis unit once it is in a series, or
    /// `None` for an item placed by a step.
    pub const fn coordinate(&self) -> Option<Quantity<f64>> {
        match self.place {
            Place::Step(_) => None,
            Place::Coordinate(coordinate) => Some(coordinate),
        }
    }

    /// Places the item at `step`, for a [`SeriesAxis::Step`] series.
    pub fn set_step(&mut self, step: u64) {
        self.place = Place::Step(step);
    }

    /// Places the item at a finite `coordinate`, for a
    /// [`SeriesAxis::Coordinate`] series.
    ///
    /// # Errors
    ///
    /// [`StructureError::NonFiniteCoordinate`] when the coordinate is not
    /// finite; the item is unchanged then.
    pub fn set_coordinate(&mut self, coordinate: Quantity<f64>) -> Result<(), StructureError> {
        if !coordinate.value().is_finite() {
            return Err(StructureError::NonFiniteCoordinate);
        }
        self.place = Place::Coordinate(coordinate);
        Ok(())
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }

    fn model_view<'a>(&'a self, topology: &'a Arc<Topology>) -> ModelView<'a> {
        ModelView::bound(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &ABSENT,
            &ABSENT,
            &self.properties,
            StructureKind::Item {
                step: self.step(),
                coordinate: self.coordinate(),
            },
        )
    }

    /// Projects every dense row through a subset correspondence, keeping the
    /// place; owner properties are dropped.
    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, StructureError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            properties: self.properties.project(atoms, bonds)?,
            place: self.place,
        })
    }

    /// Positions move and the cell rotates.
    fn transformed(&self, transform: RigidTransform) -> Result<Self, StructureError> {
        Ok(Self {
            cell: transform_cell(self.cell, transform)?,
            positions: transform_positions(&self.positions, transform)?,
            properties: self.properties.clone(),
            place: self.place,
        })
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    structure_readers!();
    structure_mutators!();
}

/// One topology-bound [`SeriesItem`] borrowed from a [`Series`].
///
/// Reads item state through `Deref` and converts to a [`ModelView`] through
/// [`AsModelView`] without copying.
#[derive(Debug, Clone, Copy)]
pub struct SeriesItemView<'a> {
    series: &'a Series,
    index: usize,
}

impl std::ops::Deref for SeriesItemView<'_> {
    type Target = SeriesItem;

    fn deref(&self) -> &Self::Target {
        &self.series.items[self.index]
    }
}

impl<'a> SeriesItemView<'a> {
    /// The item's index in its series.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The topology the item is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        &self.series.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        self.series.shared_topology()
    }

    /// The borrowed item with the view's full lifetime.
    pub fn payload(self) -> &'a SeriesItem {
        &self.series.items[self.index]
    }

    /// Borrows the item as the common kernel input, with the view's full
    /// lifetime. Its [`StructureKind::Item`] carries the item's step or
    /// coordinate.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.payload().model_view(&self.series.topology)
    }

    /// Materializes the item as an owned model, dropping its step or
    /// coordinate.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl super::model::sealed::Sealed for SeriesItemView<'_> {}

impl AsModelView for SeriesItemView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        SeriesItemView::as_model_view(*self)
    }
}

/// Dimension-preserving mutable access to one item of a [`Series`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// item through [`Series::replace`].
#[derive(Debug)]
pub struct SeriesItemMut<'a> {
    item: &'a mut SeriesItem,
    axis: SeriesAxis,
}

impl std::ops::Deref for SeriesItemMut<'_> {
    type Target = SeriesItem;

    fn deref(&self) -> &Self::Target {
        self.item
    }
}

impl SeriesItemMut<'_> {
    /// Moves the item to `step`.
    ///
    /// # Errors
    ///
    /// [`SeriesError::AxisMismatch`] unless the series is a
    /// [`SeriesAxis::Step`] series; the item is unchanged then.
    pub fn set_step(&mut self, step: u64) -> Result<(), SeriesError> {
        self.item.place = self.axis.place(Place::Step(step))?;
        Ok(())
    }

    /// Moves the item to a finite `coordinate`, converted to the axis unit.
    ///
    /// # Errors
    ///
    /// [`SeriesError::AxisMismatch`] unless the series is a
    /// [`SeriesAxis::Coordinate`] series, [`SeriesError::Unit`] for a
    /// coordinate of another dimension than the axis unit, and
    /// [`SeriesError::Structure`] for one that is not finite; the item is
    /// unchanged then.
    pub fn set_coordinate(&mut self, coordinate: Quantity<f64>) -> Result<(), SeriesError> {
        self.item.place = self.axis.place(Place::Coordinate(coordinate))?;
        Ok(())
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        self.item.apply_rigid_transform(transform)
    }

    fn storage(&mut self) -> &mut SeriesItem {
        self.item
    }

    structure_mutators!();
}

/// Iterator over the topology-bound items of a [`Series`], in series order.
#[derive(Debug, Clone)]
pub struct SeriesIter<'a> {
    series: &'a Series,
    indices: std::ops::Range<usize>,
}

impl<'a> Iterator for SeriesIter<'a> {
    type Item = SeriesItemView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let series = self.series;
        self.indices
            .next()
            .map(|index| SeriesItemView { series, index })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for SeriesIter<'_> {}

/// Failure of a [`Series`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum SeriesError {
    /// A constructor received no models.
    EmptySource,
    /// An item index is outside the series.
    IndexOutOfRange {
        /// The requested item index.
        index: usize,
        /// The number of items.
        len: usize,
    },
    /// A series-wide replacement has the wrong number of items.
    ItemCountMismatch {
        /// The series' item count.
        expected: usize,
        /// The number of replacements supplied.
        actual: usize,
    },
    /// A model's topology layout differs from the others'
    /// ([`Topology::same_layout`]).
    TopologyMismatch,
    /// An item placed by a step was given to a coordinate series, or the
    /// other way round.
    AxisMismatch {
        /// The series' axis.
        axis: SeriesAxis,
    },
    /// An item's structural state is invalid; see the source.
    Structure(StructureError),
    /// Projecting the topology onto a selection failed; see the source.
    Subset(TopologySubsetError),
    /// A coordinate could not be converted to the axis unit; see the source.
    Unit(Box<UnitError>),
}

impl fmt::Display for SeriesError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("series source is empty"),
            Self::IndexOutOfRange { index, len } => {
                write!(formatter, "index {index} is out of range for {len} items")
            }
            Self::ItemCountMismatch { expected, actual } => write!(
                formatter,
                "series requires {expected} replacement items, but received {actual}"
            ),
            Self::TopologyMismatch => {
                formatter.write_str("model belongs to a different topology layout")
            }
            Self::AxisMismatch { axis } => {
                write!(formatter, "item does not lie on a series axis of {axis}")
            }
            Self::Structure(_) => formatter.write_str("invalid series item"),
            Self::Subset(_) => formatter.write_str("cannot subset series"),
            Self::Unit(_) => formatter.write_str("invalid coordinate unit"),
        }
    }
}

impl std::error::Error for SeriesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::Unit(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<StructureError> for SeriesError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<TopologySubsetError> for SeriesError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<UnitError> for SeriesError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
