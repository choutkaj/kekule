//! Borrowed geometries of a [`StructureSet`] or another [`StructureSetRead`].

use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::math::{PeriodicCell, Point3, RigidTransform};
use crate::properties::PropertyTableMut;
use crate::structure::geometry_validation::{is_bound, validate_binding};
use crate::structure::view::sealed::Sealed;
use crate::structure::{
    AsStructureView, Forces, GeometryError, GeometryMut, Structure, StructureKind, StructureView,
    Velocities,
};
use crate::topology::{Topology, TopologyAtomIndex, TopologyBondIndex};
use crate::units::Quantity;

use super::{
    Geometry, Place, RowView, SetContext, StructureSet, StructureSetRead, StructureSetReadError,
};

/// One topology-bound [`Geometry`] borrowed from a [`StructureSet`] or another
/// [`StructureSetRead`], with its place in the reader's order.
///
/// Reads geometry state through `Deref` and converts to a [`StructureView`]
/// through [`AsStructureView`] without copying; the view's
/// [`StructureKind::InSet`] carries the geometry's [`SetContext`].
///
/// The type parameter says where the view comes from. `GeometryView<'_>` is a
/// geometry of a set in memory ([`StructureSet::get`]).
/// `GeometryView<'_, dyn StructureSetRead>` is a geometry of any reader, as
/// [`StructureSetRead::peek`] lends it; [`Self::as_dyn`] turns the first into the
/// second.
pub struct GeometryView<'a, R: ?Sized = StructureSet> {
    pub(super) reader: &'a dyn StructureSetRead,
    pub(super) geometry: &'a Geometry,
    pub(super) index: usize,
    origin: PhantomData<&'a R>,
}

impl<R: ?Sized> Clone for GeometryView<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R: ?Sized> Copy for GeometryView<'_, R> {}

impl<R: ?Sized> fmt::Debug for GeometryView<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("GeometryView")
            .field("index", &self.index)
            .field("geometry", self.geometry)
            .finish_non_exhaustive()
    }
}

impl<R: ?Sized> std::ops::Deref for GeometryView<'_, R> {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        self.geometry
    }
}

impl<'a> GeometryView<'a, dyn StructureSetRead> {
    /// A view of `geometry` as geometry `index` of `reader`, for a
    /// [`StructureSetRead::peek`] that lends a geometry the reader holds, such as a
    /// host's cache of geometries read from another reader.
    ///
    /// # Errors
    ///
    /// [`StructureSetReadError::IndexOutOfRange`] past the reader's last geometry,
    /// and [`StructureSetReadError::Geometry`] for a geometry not bound to the
    /// reader's topology: one without a position per atom and a property
    /// row per atom and bond.
    pub fn new(
        reader: &'a dyn StructureSetRead,
        index: usize,
        geometry: &'a Geometry,
    ) -> Result<Self, StructureSetReadError> {
        let len = reader.len();
        if index >= len {
            return Err(StructureSetReadError::IndexOutOfRange(
                IndexOutOfRange::new(Countable::Geometry, index, len),
            ));
        }
        check_bound(reader.topology(), geometry).map_err(|source| {
            StructureSetReadError::Geometry {
                geometry: index,
                source,
            }
        })?;
        Ok(Self::of_reader(reader, index, geometry))
    }

    /// A view of `geometry`, geometry `index` of `reader` and bound to its
    /// topology.
    pub(crate) fn of_reader(
        reader: &'a dyn StructureSetRead,
        index: usize,
        geometry: &'a Geometry,
    ) -> Self {
        Self {
            reader,
            geometry,
            index,
            origin: PhantomData,
        }
    }
}

impl<'a> GeometryView<'a> {
    /// A view of geometry `index` of an in-memory `set`, which holds `geometry`.
    pub(crate) fn of_set(set: &'a StructureSet, index: usize, geometry: &'a Geometry) -> Self {
        Self {
            reader: set,
            geometry,
            index,
            origin: PhantomData,
        }
    }
}

impl<'a, R: ?Sized> GeometryView<'a, R> {
    /// The same view as one of any reader, as [`StructureSetRead::peek`] returns
    /// it.
    pub fn as_dyn(self) -> GeometryView<'a, dyn StructureSetRead> {
        GeometryView::of_reader(self.reader, self.index, self.geometry)
    }

    /// The geometry's index in its set.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The geometry's place in the set's order, or `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        self.reader.order().map(|order| order.place(self.index))
    }

    /// The geometry's row of the set's columns: its own scalars. `None` for
    /// a reader without columns.
    pub fn row(self) -> Option<RowView<'a>> {
        self.reader
            .columns()
            .and_then(|columns| RowView::new(columns, self.index))
    }

    /// The topology the geometry is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        self.reader.topology()
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.reader.topology())
    }

    /// The bound topology's handle, with the view's full lifetime.
    pub(crate) fn topology_handle(self) -> &'a Arc<Topology> {
        self.reader.topology()
    }

    /// The borrowed geometry with the view's full lifetime.
    pub fn payload(self) -> &'a Geometry {
        self.geometry
    }

    /// Borrows the geometry as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_structure_view(self) -> StructureView<'a> {
        self.geometry.structure_view(
            self.reader.topology(),
            StructureKind::InSet(SetContext::of_view(self.as_dyn())),
        )
    }

    /// Materializes the geometry as an owned structure, dropping its place
    /// and keeping its velocities and forces. The cells of its row become
    /// the structure's owner values, so the structure is whole; a row with
    /// cells copies the geometry once.
    pub fn to_structure(self) -> Structure {
        self.as_structure_view().to_structure()
    }
}

/// Dimension-preserving mutable access to a geometry of a set in memory
/// ([`StructureSet::get_mut`]): a [`GeometryMut`] without owner values.
///
/// Reads use `Deref`. Positions, the cell, site values, velocities, forces,
/// and the per-atom and per-bond property tables are editable. The
/// geometry's own scalars are its row of the set's columns, edited through
/// [`StructureSet::columns_mut`], so this guard offers no owner values:
///
/// ```compile_fail,E0599
/// use kekule::structure::StructureSet;
/// fn label(set: &mut StructureSet) {
///     set.get_mut(0).unwrap().properties_mut();
/// }
/// ```
#[derive(Debug)]
pub struct SetGeometryMut<'a> {
    geometry: GeometryMut<'a>,
}

impl<'a> SetGeometryMut<'a> {
    pub(crate) fn new(geometry: &'a mut Geometry) -> Self {
        Self {
            geometry: GeometryMut::new(geometry),
        }
    }
}

impl std::ops::Deref for SetGeometryMut<'_> {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        &self.geometry
    }
}

impl SetGeometryMut<'_> {
    /// Sets or clears velocities; see [`Geometry::set_velocities`].
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), GeometryError> {
        self.geometry.set_velocities(velocities)
    }

    /// Sets or clears forces; see [`Geometry::set_forces`].
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), GeometryError> {
        self.geometry.set_forces(forces)
    }

    /// Removes velocities, returning them.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.geometry.take_velocities()
    }

    /// Removes forces, returning them.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.geometry.take_forces()
    }

    /// Applies a rigid transform transactionally: positions move; the cell,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), GeometryError> {
        self.geometry.apply_rigid_transform(transform)
    }

    /// Replaces every position; see [`Geometry::set_positions`].
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), GeometryError> {
        self.geometry.set_positions(positions)
    }

    /// Sets one position; see [`Geometry::set_position_at`].
    pub fn set_position_at(
        &mut self,
        atom: TopologyAtomIndex,
        position: Quantity<Point3>,
    ) -> Result<(), GeometryError> {
        self.geometry.set_position_at(atom, position)
    }

    /// Replaces the periodic cell; see [`Geometry::set_cell`].
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.geometry.set_cell(cell);
    }

    /// Mutable topology-atom columns of the geometry's properties. A shared
    /// geometry is copied when the guard is taken.
    pub fn atom_properties_mut(&mut self) -> PropertyTableMut<'_, TopologyAtomIndex> {
        self.geometry.atom_properties_mut()
    }

    /// Mutable topology-bond columns of the geometry's properties. A shared
    /// geometry is copied when the guard is taken.
    pub fn bond_properties_mut(&mut self) -> PropertyTableMut<'_, TopologyBondIndex> {
        self.geometry.bond_properties_mut()
    }

    /// Sets or clears one occupancy; see [`Geometry::set_occupancy`].
    pub fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), GeometryError> {
        self.geometry.set_occupancy(atom, value)
    }

    /// Replaces or clears every occupancy; see [`Geometry::set_occupancies`].
    pub fn set_occupancies(
        &mut self,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), GeometryError> {
        self.geometry.set_occupancies(values)
    }

    /// Sets or clears one B-factor; see [`Geometry::set_b_factor`].
    pub fn set_b_factor(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), GeometryError> {
        self.geometry.set_b_factor(atom, value)
    }

    /// Replaces or clears every B-factor; see [`Geometry::set_b_factors`].
    pub fn set_b_factors(
        &mut self,
        values: Option<Quantity<Vec<Option<f64>>>>,
    ) -> Result<(), GeometryError> {
        self.geometry.set_b_factors(values)
    }
}

/// Fails unless `geometry` is bound to `topology`.
pub(crate) fn check_bound(topology: &Topology, geometry: &Geometry) -> Result<(), GeometryError> {
    let atoms = geometry.atom_count();
    let properties = geometry.properties();
    validate_binding(atoms, properties, topology)?;
    if !is_bound(atoms, properties, topology) {
        return Err(GeometryError::BondCountMismatch(LengthMismatch::new(
            Countable::Row,
            topology.bond_count(),
            properties.bonds().len(),
        )));
    }
    Ok(())
}

impl<R: ?Sized> Sealed for GeometryView<'_, R> {}

impl<R: ?Sized> AsStructureView for GeometryView<'_, R> {
    fn as_structure_view(&self) -> StructureView<'_> {
        GeometryView::as_structure_view(*self)
    }
}

/// Iterator over the topology-bound geometries of a [`StructureSet`], in set
/// order.
#[derive(Debug, Clone)]
pub struct GeometryIter<'a> {
    set: &'a StructureSet,
    indices: std::ops::Range<usize>,
}

impl<'a> GeometryIter<'a> {
    pub(crate) fn new(set: &'a StructureSet) -> Self {
        Self {
            set,
            indices: 0..set.len(),
        }
    }
}

impl<'a> Iterator for GeometryIter<'a> {
    type Item = GeometryView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let set = self.set;
        self.indices.next().and_then(|index| set.get(index))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for GeometryIter<'_> {}
