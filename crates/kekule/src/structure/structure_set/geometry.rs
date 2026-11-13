//! One geometry of a [`super::StructureSet`] and the structure it shares.

use std::sync::Arc;

use crate::math::{PeriodicCell, Point3, RigidTransform, Vector3};
use crate::properties::{GeometryProperties, GeometryPropertiesMut};
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::Quantity;

use crate::structure::model::ModelView;
use crate::structure::scalars::{canonical_b_factor, canonical_b_factors, AtomScalars};
use crate::structure::state::{
    check_property_rows, is_bound, site_readers, structure_readers, transform_cell,
    transform_positions, validate_binding,
};
use crate::structure::{
    DenseArrayError, Forces, GeometryError, Model, Positions, StructureKind, Velocities,
};

/// The structure state of a geometry, shared between clones and copied on
/// write.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeometryState {
    pub(crate) positions: Positions,
    pub(crate) cell: Option<PeriodicCell>,
    pub(crate) occupancies: AtomScalars,
    pub(crate) b_factors: AtomScalars,
    pub(crate) properties: GeometryProperties,
    pub(crate) velocities: Option<Velocities>,
    pub(crate) forces: Option<Forces>,
}

impl GeometryState {
    /// Positions alone: no cell, site fields, properties, or dynamics.
    pub(crate) fn new(positions: Positions) -> Self {
        let properties = GeometryProperties::new(positions.len());
        Self::from_frame_state(positions, None, properties)
    }

    /// Positions, a cell, and properties, without site fields or dynamics.
    pub(crate) fn from_frame_state(
        positions: Positions,
        cell: Option<PeriodicCell>,
        properties: GeometryProperties,
    ) -> Self {
        Self {
            positions,
            cell,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            velocities: None,
            forces: None,
        }
    }

    /// A model's structure state; the model's topology is dropped.
    pub(crate) fn from_model(model: Model) -> Self {
        Self {
            positions: model.positions,
            cell: model.cell,
            occupancies: model.occupancies,
            b_factors: model.b_factors,
            properties: model.properties,
            velocities: None,
            forces: None,
        }
    }

    /// Projects every dense row through a subset correspondence.
    pub(crate) fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, GeometryError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(atoms),
            b_factors: self.b_factors.select(atoms),
            properties: self.properties.project(atoms, bonds)?,
            velocities: self.velocities.as_ref().map(|values| values.select(atoms)),
            forces: self.forces.as_ref().map(|values| values.select(atoms)),
        })
    }

    /// Positions move; the cell, velocities, and forces rotate.
    pub(crate) fn transformed(&self, transform: RigidTransform) -> Result<Self, GeometryError> {
        let rotate = |vector| transform.transform_vector(vector);
        Ok(Self {
            cell: transform_cell(self.cell, transform)?,
            positions: transform_positions(&self.positions, transform)?,
            occupancies: self.occupancies.clone(),
            b_factors: self.b_factors.clone(),
            properties: self.properties.clone(),
            velocities: self
                .velocities
                .as_ref()
                .map(|values| values.rotated(rotate))
                .transpose()?,
            forces: self
                .forces
                .as_ref()
                .map(|values| values.rotated(rotate))
                .transpose()?,
        })
    }

    /// Replaces the positions of a possibly shared structure, copying the
    /// rest of it when it is shared but never the positions being replaced.
    pub(crate) fn replace_positions(structure: &mut Arc<Self>, positions: Positions) {
        match Arc::get_mut(structure) {
            Some(owned) => owned.positions = positions,
            None => {
                let shared = &**structure;
                *structure = Arc::new(Self {
                    positions,
                    cell: shared.cell,
                    occupancies: shared.occupancies.clone(),
                    b_factors: shared.b_factors.clone(),
                    properties: shared.properties.clone(),
                    velocities: shared.velocities.clone(),
                    forces: shared.forces.clone(),
                });
            }
        }
    }

    /// Checks that a per-atom array has one row per atom.
    fn check_len(&self, actual: usize) -> Result<(), GeometryError> {
        if actual != self.positions.len() {
            return Err(GeometryError::Array(DenseArrayError::CountMismatch {
                expected: self.positions.len(),
                actual,
            }));
        }
        Ok(())
    }

    /// Binds a possibly shared structure to `topology`, allocating its bond
    /// rows. A bound structure, such as one taken from another set of this
    /// topology, stays shared.
    pub(crate) fn bind(
        structure: &mut Arc<Self>,
        topology: &Topology,
    ) -> Result<(), GeometryError> {
        let atoms = structure.positions.len();
        validate_binding(atoms, &structure.properties, topology)?;
        if !is_bound(atoms, &structure.properties, topology) {
            Arc::make_mut(structure)
                .properties
                .bind_bond_rows(topology.bond_count());
        }
        Ok(())
    }

    /// Borrows the structure bound to `topology`.
    pub(crate) fn model_view<'a>(
        &'a self,
        topology: &'a Arc<Topology>,
        kind: StructureKind<'a>,
    ) -> ModelView<'a> {
        ModelView::bound(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.occupancies,
            &self.b_factors,
            &self.properties,
            kind,
        )
    }
}

/// One geometry of a [`super::StructureSet`]: positions in dense topology atom
/// order, an optional periodic cell, structure properties, and, only when
/// stated, occupancies, B-factors, velocities, and forces.
///
/// A geometry has no place and no weight: an order's values and a
/// distribution's weights belong to the set, so a geometry taken out of its set
/// claims neither. The geometry carries no topology. It is *bound* to one when
/// it has a property row per atom and bond of that topology: a set binds the
/// geometries it takes, and [`Self::bind`] binds a geometry outside a set, as a
/// [`FrameView`](crate::traj::FrameView) of it needs.
///
/// The structure is shared storage: cloning a geometry, and therefore a set,
/// copies no coordinates. Editing a geometry whose storage is shared copies
/// that one geometry (copy on write), so the other holders are unaffected. An
/// edit that fails, and one that changes nothing such as clearing absent
/// velocities or setting the same cell, keeps the storage shared.
///
/// ```
/// use kekule::structure::{Positions, Geometry};
/// let geometry = Geometry::new(Positions::zeros(3));
/// let mut moved = geometry.clone();
/// assert!(moved.shares_structure_with(&geometry));
/// moved.set_occupancies(Some(vec![Some(0.5); 3]))?;
/// assert!(!moved.shares_structure_with(&geometry));
/// assert_eq!(geometry.occupancies(), None);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Geometry {
    structure: Arc<GeometryState>,
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
        Self::from_structure(Arc::new(GeometryState::new(positions)))
    }

    /// A geometry holding a model's positions, cell, occupancies, B-factors,
    /// and structure properties, moving them without copying; the model's
    /// topology is dropped.
    pub fn from_model(model: Model) -> Self {
        Self::from_structure(Arc::new(GeometryState::from_model(model)))
    }

    pub(crate) const fn from_structure(structure: Arc<GeometryState>) -> Self {
        Self { structure }
    }

    pub(crate) fn structure(&self) -> &Arc<GeometryState> {
        &self.structure
    }

    pub(crate) fn structure_mut(&mut self) -> &mut Arc<GeometryState> {
        &mut self.structure
    }

    /// The structure for editing in place: one shared with a clone is
    /// copied first, once. The caller keeps it bound and valid.
    pub(crate) fn structure_in_place(&mut self) -> &mut GeometryState {
        Arc::make_mut(&mut self.structure)
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
        GeometryState::bind(&mut self.structure, topology)
    }

    /// This geometry projected onto the dense atom and bond rows of a
    /// topology subset.
    pub(crate) fn projected(
        &self,
        atoms: &[usize],
        bonds: &[usize],
    ) -> Result<Self, GeometryError> {
        Ok(Self::from_structure(Arc::new(
            self.structure.project(atoms, bonds)?,
        )))
    }

    /// The positions for writing in place: a structure shared with a clone
    /// is copied first, once. The caller keeps every value finite.
    pub(crate) fn positions_in_place(&mut self) -> &mut [Point3] {
        Arc::make_mut(&mut self.structure).positions.as_mut_slice()
    }

    /// Moves every position by `shift` in place, copying a shared structure
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
    /// [`Self::apply_rigid_transform`] does, copying a shared structure
    /// once. On failure the geometry is partly transformed, so the caller
    /// drops it.
    pub(crate) fn transform_in_place(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), GeometryError> {
        let cell = transform_cell(self.structure.cell, *transform)?;
        let structure = Arc::make_mut(&mut self.structure);
        structure.cell = cell;
        for (index, point) in structure.positions.as_mut_slice().iter_mut().enumerate() {
            *point = transform.transform_point(*point);
            if !point.is_finite() {
                return Err(DenseArrayError::NonFinite { index }.into());
            }
        }
        let rotate = |vector| transform.transform_vector(vector);
        if let Some(velocities) = &mut structure.velocities {
            velocities.rotate_in_place(rotate)?;
        }
        if let Some(forces) = &mut structure.forces {
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
        self.structure.positions.len()
    }

    /// Whether this geometry and `other` share one structure's storage, so
    /// that neither has copied the other's coordinates. Clones share until
    /// one of them is edited.
    pub fn shares_structure_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.structure, &other.structure)
    }

    /// Per-atom velocities (nm/ps), or `None` when the geometry has none.
    pub fn velocities(&self) -> Option<&Velocities> {
        self.structure.velocities.as_ref()
    }

    /// Per-atom forces (kJ/(mol·nm)), or `None` when the geometry has none.
    pub fn forces(&self) -> Option<&Forces> {
        self.structure.forces.as_ref()
    }

    /// Sets or clears velocities; they must have one vector per atom.
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), GeometryError> {
        match &velocities {
            Some(values) => self.structure.check_len(values.len())?,
            None if self.structure.velocities.is_none() => return Ok(()),
            None => {}
        }
        Arc::make_mut(&mut self.structure).velocities = velocities;
        Ok(())
    }

    /// Sets or clears forces; they must have one vector per atom.
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), GeometryError> {
        match &forces {
            Some(values) => self.structure.check_len(values.len())?,
            None if self.structure.forces.is_none() => return Ok(()),
            None => {}
        }
        Arc::make_mut(&mut self.structure).forces = forces;
        Ok(())
    }

    /// Removes velocities, returning them.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.structure.velocities.as_ref()?;
        Arc::make_mut(&mut self.structure).velocities.take()
    }

    /// Removes forces, returning them.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.structure.forces.as_ref()?;
        Arc::make_mut(&mut self.structure).forces.take()
    }

    /// Replaces every position without changing the atom count. Units,
    /// dimensions, and finite values are checked before mutation.
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), GeometryError> {
        self.structure.positions.validate_all(&positions)?;
        self.edit(|structure| Ok(structure.positions.set_all(positions)?))
    }

    /// Sets one position by dense atom index.
    pub fn set_position_at(
        &mut self,
        atom: TopologyAtomIndex,
        position: Quantity<Point3>,
    ) -> Result<(), GeometryError> {
        self.edit(|structure| Ok(structure.positions.set(atom.index(), position)?))
    }

    /// Replaces the periodic cell; `None` makes the geometry non-periodic.
    /// Positions are not wrapped or rescaled.
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        if self.structure.cell != cell {
            Arc::make_mut(&mut self.structure).cell = cell;
        }
    }

    /// Replaces every structure annotation with rows of the same
    /// dimensions, for example ones copied from another structure of the
    /// same topology. A geometry not yet in a set accepts no bond rows.
    pub fn set_properties(&mut self, properties: GeometryProperties) -> Result<(), GeometryError> {
        check_property_rows(
            self.structure.positions.len(),
            &self.structure.properties,
            &properties,
        )?;
        Arc::make_mut(&mut self.structure).properties = properties;
        Ok(())
    }

    /// Length-preserving mutable access to structure annotations. A shared
    /// structure is copied when the guard is taken.
    pub fn properties_mut(&mut self) -> GeometryPropertiesMut<'_> {
        GeometryPropertiesMut::new(&mut Arc::make_mut(&mut self.structure).properties)
    }

    /// Sets or clears one occupancy. Values must be finite.
    pub fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), GeometryError> {
        self.edit(|structure| {
            let len = structure.positions.len();
            structure.occupancies.set(atom, len, value)
        })
    }

    /// Replaces every occupancy in dense atom order, or clears them all.
    pub fn set_occupancies(
        &mut self,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), GeometryError> {
        self.edit(|structure| {
            let len = structure.positions.len();
            structure.occupancies.replace(len, values)
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
        self.edit(|structure| {
            let len = structure.positions.len();
            structure.b_factors.set(atom, len, value)
        })
    }

    /// Replaces every B-factor in dense atom order, or clears them all.
    pub fn set_b_factors(
        &mut self,
        values: Option<Quantity<Vec<Option<f64>>>>,
    ) -> Result<(), GeometryError> {
        let values = canonical_b_factors(values)?;
        self.edit(|structure| {
            let len = structure.positions.len();
            structure.b_factors.replace(len, values)
        })
    }

    /// Applies a rigid transform transactionally: positions move; the cell,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), GeometryError> {
        self.structure = Arc::new(self.structure.transformed(*transform)?);
        Ok(())
    }

    fn storage_ref(&self) -> &GeometryState {
        &self.structure
    }

    /// Runs an edit that changes nothing when it fails. A shared structure
    /// is edited as a copy that replaces it only on success, so a failed
    /// edit keeps the structure shared.
    fn edit<R>(
        &mut self,
        edit: impl FnOnce(&mut GeometryState) -> Result<R, GeometryError>,
    ) -> Result<R, GeometryError> {
        if let Some(owned) = Arc::get_mut(&mut self.structure) {
            return edit(owned);
        }
        let mut copy = GeometryState::clone(&self.structure);
        let result = edit(&mut copy)?;
        self.structure = Arc::new(copy);
        Ok(result)
    }

    structure_readers!();
    site_readers!();
}

/// Dimension-preserving mutable access to one geometry of a
/// [`super::StructureSet`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// geometry through [`super::StructureSet::replace`]. The geometry's place and
/// weight belong to the set ([`super::StructureSet::set_place`],
/// [`super::DistributionMut::set_weight`]).
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

    /// Replaces every structure annotation; see [`Geometry::set_properties`].
    pub fn set_properties(&mut self, properties: GeometryProperties) -> Result<(), GeometryError> {
        self.geometry.set_properties(properties)
    }

    /// Length-preserving mutable access to structure annotations; see
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
