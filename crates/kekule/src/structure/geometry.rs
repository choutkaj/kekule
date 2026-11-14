//! [`Geometry`], the one owned coordinate payload of structures, structure
//! sets, and frame buffers, and its guard [`GeometryMut`].

mod data;
mod error;

pub(crate) use data::{subset_rows, transform_cell, GeometryData};
pub use error::GeometryError;

use std::sync::Arc;

use crate::math::{PeriodicCell, Point3, RigidTransform, Vector3};
use crate::properties::{GeometryProperties, GeometryPropertiesMut};
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::Quantity;

use super::geometry_validation::check_property_rows;
use super::scalars::{
    b_factor_quantities, b_factor_quantity, canonical_b_factor, canonical_b_factors,
};
use super::{
    DenseArrayError, Forces, Positions, Structure, StructureKind, StructureView, Velocities,
};

/// One arrangement of a topology's atoms: positions in dense topology atom
/// order, an optional periodic cell, geometry properties, and, only when
/// stated, occupancies, B-factors, velocities, and forces.
///
/// A geometry is the one owned coordinate payload. A [`Structure`] is a
/// topology and one geometry, and reads it through `Deref`; a
/// [`StructureSet`](super::StructureSet) is a topology and many, and a
/// [`FrameBuffer`](crate::traj::FrameBuffer) holds one for a stream. Owners
/// that keep a geometry bound hand out a [`GeometryMut`] to edit it.
///
/// A geometry has no place and no weight: an order's values belong to the
/// set and a distribution's weights to an
/// [`Ensemble`](crate::structure::Ensemble), so a geometry taken out of its
/// set claims neither. The geometry carries no topology. It is *bound* to one when
/// it has a property row per atom and bond of that topology: a set binds the
/// geometries it takes, and [`Self::bind`] binds a geometry outside a set, as a
/// [`FrameView`](crate::traj::FrameView) of it needs.
///
/// Its state is shared storage: cloning a geometry, and therefore a set,
/// copies no coordinates. Editing a geometry whose storage is shared copies
/// that one geometry (copy on write), so the other holders are unaffected. An
/// edit that fails, and one that changes nothing such as clearing absent
/// velocities or setting the same cell, keeps the storage shared.
///
/// ```
/// use kekule::structure::{Positions, Geometry};
/// let geometry = Geometry::new(Positions::zeros(3));
/// let mut moved = geometry.clone();
/// assert!(moved.shares_storage_with(&geometry));
/// moved.set_occupancies(Some(vec![Some(0.5); 3]))?;
/// assert!(!moved.shares_storage_with(&geometry));
/// assert_eq!(geometry.occupancies(), None);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    data: Arc<GeometryData>,
}

impl Geometry {
    /// A non-periodic geometry without properties, occupancies, B-factors,
    /// velocities, or forces.
    ///
    /// The geometry is bound to no topology: it has no bond rows. A set binds
    /// it when it takes it; bind it with [`Self::bind`] to use it outside a
    /// set, such as in a [`FrameView`](crate::traj::FrameView) of a topology
    /// with bonds.
    pub fn new(positions: Positions) -> Self {
        Self::from_data(GeometryData::new(positions))
    }

    /// A structure's geometry, moved without copying; the structure's
    /// topology is dropped. See [`Structure::into_parts`].
    pub fn from_structure(structure: Structure) -> Self {
        structure.into_parts().1
    }

    /// Binds the geometry to `topology`, as a set binds the geometries it takes:
    /// checks that it has a position per atom, and gives it a property row
    /// per bond unless it has them already. Binding never copies
    /// positions of a geometry whose storage only it holds, and a bound
    /// geometry stays as it is, sharing its storage.
    ///
    /// ```
    /// use kekule::structure::{Positions, Geometry};
    /// use kekule::traj::FrameView;
    ///
    /// let topology = kekule::smiles::to_topology("CO")?;
    /// let mut geometry = Geometry::new(Positions::zeros(topology.atom_count()));
    /// assert!(FrameView::new(&topology, &geometry, None, None).is_err());
    /// geometry.bind(&topology)?;
    /// let frame = FrameView::new(&topology, &geometry, None, Some(0))?;
    /// assert_eq!(frame.step(), Some(0));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// [`GeometryError::Array`] for another atom count, and
    /// [`GeometryError::BondCountMismatch`] for bond rows of another
    /// topology. A failed binding changes nothing.
    pub fn bind(&mut self, topology: &Topology) -> Result<(), GeometryError> {
        GeometryData::bind(&mut self.data, topology)
    }

    /// This geometry projected onto the dense atom and bond rows of a
    /// topology subset.
    pub(crate) fn projected(
        &self,
        atoms: &[usize],
        bonds: &[usize],
    ) -> Result<Self, GeometryError> {
        Ok(Self::from_data(self.data.project(atoms, bonds)?))
    }

    /// The positions for writing in place: a state shared with a clone
    /// is copied first, once. The caller keeps every value finite.
    pub(crate) fn positions_in_place(&mut self) -> &mut [Point3] {
        Arc::make_mut(&mut self.data).positions.as_mut_slice()
    }

    /// Moves every position by `shift` in place, copying a shared state
    /// once; the cell, velocities, and forces do not change under a
    /// translation. On failure some positions have moved, so the caller
    /// drops the geometry.
    pub(crate) fn translate_in_place(&mut self, shift: Vector3) -> Result<(), GeometryError> {
        for (index, point) in self.positions_in_place().iter_mut().enumerate() {
            *point = *point + shift;
            if !point.is_finite() {
                return Err(DenseArrayError::NonFinite { index }.into());
            }
        }
        Ok(())
    }

    /// Applies a rigid transform in place, as
    /// [`Self::apply_rigid_transform`] does, copying a shared state
    /// once. On failure the geometry is partly transformed, so the caller
    /// drops it.
    pub(crate) fn transform_in_place(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), GeometryError> {
        let cell = transform_cell(self.data.cell, *transform)?;
        let data = Arc::make_mut(&mut self.data);
        data.cell = cell;
        for (index, point) in data.positions.as_mut_slice().iter_mut().enumerate() {
            *point = transform.transform_point(*point);
            if !point.is_finite() {
                return Err(DenseArrayError::NonFinite { index }.into());
            }
        }
        let rotate = |vector| transform.transform_vector(vector);
        if let Some(velocities) = &mut data.velocities {
            velocities.rotate_in_place(rotate)?;
        }
        if let Some(forces) = &mut data.forces {
            forces.rotate_in_place(rotate)?;
        }
        Ok(())
    }

    /// Returns this geometry with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.set_cell(Some(cell));
        self
    }

    /// Number of atoms, the length of every dense per-atom array of the
    /// geometry.
    pub fn atom_count(&self) -> usize {
        self.data.positions.len()
    }

    /// Whether this geometry and `other` share their storage, so
    /// that neither has copied the other's coordinates. Clones share until
    /// one of them is edited.
    pub fn shares_storage_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }

    /// Positions in dense topology atom order.
    pub fn positions(&self) -> &Positions {
        &self.data.positions
    }

    /// The periodic cell, or `None` for a non-periodic geometry.
    pub fn cell(&self) -> Option<&PeriodicCell> {
        self.data.cell.as_ref()
    }

    /// Geometry annotations: owner values plus one row per dense atom and
    /// bond.
    pub fn properties(&self) -> &GeometryProperties {
        &self.data.properties
    }

    /// Occupancy of one atom; `None` when it was never recorded or `atom` is
    /// out of range.
    pub fn occupancy(&self, atom: TopologyAtomIndex) -> Option<f64> {
        self.data.occupancies.get(atom)
    }

    /// Every occupancy in dense atom order, or `None` when none was recorded.
    pub fn occupancies(&self) -> Option<&[Option<f64>]> {
        self.data.occupancies.values()
    }

    /// Isotropic B-factor of one atom; `None` when it was never recorded or
    /// `atom` is out of range.
    pub fn b_factor(&self, atom: TopologyAtomIndex) -> Option<Quantity<f64>> {
        b_factor_quantity(self.data.b_factors.get(atom))
    }

    /// Every B-factor in dense atom order and square nanometres, or `None`
    /// when none was recorded.
    pub fn b_factors(&self) -> Option<Quantity<&[Option<f64>]>> {
        b_factor_quantities(&self.data.b_factors)
    }

    /// Borrows this geometry, bound to `topology`, as a structure view.
    pub(crate) fn structure_view<'a>(
        &'a self,
        topology: &'a Arc<Topology>,
        kind: StructureKind<'a>,
    ) -> StructureView<'a> {
        StructureView::bound(topology, self, kind)
    }

    /// Per-atom velocities (nm/ps), or `None` when the geometry has none.
    pub fn velocities(&self) -> Option<&Velocities> {
        self.data.velocities.as_ref()
    }

    /// Per-atom forces (kJ/(mol·nm)), or `None` when the geometry has none.
    pub fn forces(&self) -> Option<&Forces> {
        self.data.forces.as_ref()
    }

    /// Sets or clears velocities; they must have one vector per atom.
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), GeometryError> {
        match &velocities {
            Some(values) => self.data.check_len(values.len())?,
            None if self.data.velocities.is_none() => return Ok(()),
            None => {}
        }
        Arc::make_mut(&mut self.data).velocities = velocities;
        Ok(())
    }

    /// Sets or clears forces; they must have one vector per atom.
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), GeometryError> {
        match &forces {
            Some(values) => self.data.check_len(values.len())?,
            None if self.data.forces.is_none() => return Ok(()),
            None => {}
        }
        Arc::make_mut(&mut self.data).forces = forces;
        Ok(())
    }

    /// Removes velocities, returning them.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.data.velocities.as_ref()?;
        Arc::make_mut(&mut self.data).velocities.take()
    }

    /// Removes forces, returning them.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.data.forces.as_ref()?;
        Arc::make_mut(&mut self.data).forces.take()
    }

    /// Replaces every position without changing the atom count. Units,
    /// dimensions, and finite values are checked before mutation.
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), GeometryError> {
        // Borrow the caller's slice once, so that the slice validated is the
        // slice copied, and no shared array is copied for a failing edit.
        let positions = Quantity::new(positions.value().as_ref(), positions.unit());
        self.data.positions.validate_all(&positions)?;
        self.edit(|data| Ok(data.positions.set_all(positions)?))
    }

    /// Sets one position by dense atom index.
    pub fn set_position_at(
        &mut self,
        atom: TopologyAtomIndex,
        position: Quantity<Point3>,
    ) -> Result<(), GeometryError> {
        self.edit(|data| Ok(data.positions.set(atom.index(), position)?))
    }

    /// Replaces the periodic cell; `None` makes the geometry non-periodic.
    /// Positions are not wrapped or rescaled.
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        if self.data.cell != cell {
            Arc::make_mut(&mut self.data).cell = cell;
        }
    }

    /// Replaces every geometry annotation with rows of the same
    /// dimensions, for example ones copied from another geometry of the
    /// same topology. A geometry not yet in a set accepts no bond rows.
    pub fn set_properties(&mut self, properties: GeometryProperties) -> Result<(), GeometryError> {
        check_property_rows(
            self.data.positions.len(),
            &self.data.properties,
            &properties,
        )?;
        Arc::make_mut(&mut self.data).properties = properties;
        Ok(())
    }

    /// Length-preserving mutable access to geometry annotations. A shared
    /// geometry is copied when the guard is taken.
    pub fn properties_mut(&mut self) -> GeometryPropertiesMut<'_> {
        GeometryPropertiesMut::new(&mut Arc::make_mut(&mut self.data).properties)
    }

    /// Sets or clears one occupancy. Values must be finite.
    pub fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), GeometryError> {
        self.edit(|data| {
            let len = data.positions.len();
            data.occupancies.set(atom, len, value)
        })
    }

    /// Replaces every occupancy in dense atom order, or clears them all.
    pub fn set_occupancies(
        &mut self,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), GeometryError> {
        self.edit(|data| {
            let len = data.positions.len();
            data.occupancies.replace(len, values)
        })
    }

    /// Sets or clears one isotropic B-factor, converting it to square
    /// nanometres.
    pub fn set_b_factor(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), GeometryError> {
        let value = canonical_b_factor(value)?;
        self.edit(|data| {
            let len = data.positions.len();
            data.b_factors.set(atom, len, value)
        })
    }

    /// Replaces every B-factor in dense atom order, or clears them all.
    pub fn set_b_factors(
        &mut self,
        values: Option<Quantity<Vec<Option<f64>>>>,
    ) -> Result<(), GeometryError> {
        let values = canonical_b_factors(values)?;
        self.edit(|data| {
            let len = data.positions.len();
            data.b_factors.replace(len, values)
        })
    }

    /// Applies a rigid transform transactionally: positions move; the cell,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), GeometryError> {
        self.data = Arc::new(self.data.transformed(*transform)?);
        Ok(())
    }

    /// Runs an edit that changes nothing when it fails. A shared state is
    /// edited as a copy that replaces it only on success, so a failed edit
    /// keeps the state shared.
    fn edit<R>(
        &mut self,
        edit: impl FnOnce(&mut GeometryData) -> Result<R, GeometryError>,
    ) -> Result<R, GeometryError> {
        if let Some(owned) = Arc::get_mut(&mut self.data) {
            return edit(owned);
        }
        let mut copy = GeometryData::clone(&self.data);
        let result = edit(&mut copy)?;
        self.data = Arc::new(copy);
        Ok(result)
    }
}

/// Dimension-preserving mutable access to a geometry its owner keeps bound
/// to a topology: the geometry of a [`Structure`]
/// ([`Structure::geometry_mut`]), of a [`super::StructureSet`]
/// ([`super::StructureSet::get_mut`]), or of a
/// [`FrameBuffer`](crate::traj::FrameBuffer).
///
/// Reads use `Deref`. There is deliberately no `DerefMut`, so no edit changes
/// the atom count or the bond rows; replace a whole geometry through its
/// owner, such as [`super::StructureSet::replace`] or
/// [`Structure::from_parts`]. A set geometry's place belongs to the set
/// ([`super::StructureSet::set_place`]) and its weight to an ensemble
/// ([`DistributionMut::set_log_weight`](crate::structure::DistributionMut::set_log_weight)).
#[derive(Debug)]
pub struct GeometryMut<'a> {
    geometry: &'a mut Geometry,
}

impl<'a> GeometryMut<'a> {
    pub(crate) fn new(geometry: &'a mut Geometry) -> Self {
        Self { geometry }
    }
}

impl std::ops::Deref for GeometryMut<'_> {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        self.geometry
    }
}

impl GeometryMut<'_> {
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

    /// Replaces every geometry annotation; see [`Geometry::set_properties`].
    pub fn set_properties(&mut self, properties: GeometryProperties) -> Result<(), GeometryError> {
        self.geometry.set_properties(properties)
    }

    /// Length-preserving mutable access to geometry annotations; see
    /// [`Geometry::properties_mut`].
    pub fn properties_mut(&mut self) -> GeometryPropertiesMut<'_> {
        self.geometry.properties_mut()
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
