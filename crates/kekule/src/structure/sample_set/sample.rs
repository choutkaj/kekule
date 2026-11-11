//! One sample of a [`super::SampleSet`] and the structure it shares.

use std::sync::Arc;

use crate::geometry::{PeriodicCell, Point3, RigidTransform};
use crate::properties::{StructureProperties, StructurePropertiesMut};
use crate::topology::{Topology, TopologyAtomIndex};
use crate::units::Quantity;

use crate::structure::model::ModelView;
use crate::structure::scalars::{canonical_b_factor, canonical_b_factors, AtomScalars};
use crate::structure::state::{
    check_property_rows, is_bound, site_readers, structure_readers, transform_cell,
    transform_positions, validate_binding,
};
use crate::structure::{
    DenseArrayError, Forces, Model, Positions, StructureError, StructureKind, Velocities,
};

/// The structure state of a sample, shared between clones and copied on
/// write. Ensemble samples hold the same storage, without velocities or
/// forces, so converting between the two copies no coordinates.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SampleStructure {
    pub(crate) positions: Positions,
    pub(crate) cell: Option<PeriodicCell>,
    pub(crate) occupancies: AtomScalars,
    pub(crate) b_factors: AtomScalars,
    pub(crate) properties: StructureProperties,
    pub(crate) velocities: Option<Velocities>,
    pub(crate) forces: Option<Forces>,
}

impl SampleStructure {
    /// Positions alone: no cell, site fields, properties, or dynamics.
    pub(crate) fn new(positions: Positions) -> Self {
        let properties = StructureProperties::new(positions.len());
        Self::from_frame_state(positions, None, properties)
    }

    /// Positions, a cell, and properties, without site fields or dynamics.
    pub(crate) fn from_frame_state(
        positions: Positions,
        cell: Option<PeriodicCell>,
        properties: StructureProperties,
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

    /// Whether any occupancy is recorded.
    pub(crate) fn has_occupancies(&self) -> bool {
        recorded(&self.occupancies)
    }

    /// Whether any B-factor is recorded.
    pub(crate) fn has_b_factors(&self) -> bool {
        recorded(&self.b_factors)
    }

    /// Projects every dense row through a subset correspondence.
    pub(crate) fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, StructureError> {
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
    pub(crate) fn transformed(&self, transform: RigidTransform) -> Result<Self, StructureError> {
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
    fn check_len(&self, actual: usize) -> Result<(), StructureError> {
        if actual != self.positions.len() {
            return Err(StructureError::Array(DenseArrayError::CountMismatch {
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
    ) -> Result<(), StructureError> {
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

/// Whether a per-atom column records any value.
fn recorded(column: &AtomScalars) -> bool {
    column
        .values()
        .is_some_and(|values| values.iter().any(Option::is_some))
}

/// One sample of a [`super::SampleSet`]: positions in dense topology atom
/// order, an optional periodic cell, structure properties, and, only when
/// stated, occupancies, B-factors, velocities, and forces.
///
/// A sample has no place and no weight: an order's values and a
/// distribution's weights belong to the set, so a sample taken out of its set
/// claims neither. The sample carries no topology; it is bound when inserted.
///
/// The structure is shared storage: cloning a sample, and therefore a set,
/// copies no coordinates. Editing a sample whose storage is shared copies
/// that one sample (copy on write), so the other holders are unaffected. An
/// edit that fails, and one that changes nothing such as clearing absent
/// velocities or setting the same cell, keeps the storage shared.
///
/// ```
/// use kekule::structure::{Positions, Sample};
/// let sample = Sample::new(Positions::zeros(3));
/// let mut moved = sample.clone();
/// assert!(moved.shares_structure_with(&sample));
/// moved.set_occupancies(Some(vec![Some(0.5); 3]))?;
/// assert!(!moved.shares_structure_with(&sample));
/// assert_eq!(sample.occupancies(), None);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct Sample {
    structure: Arc<SampleStructure>,
}

impl Sample {
    /// A non-periodic sample without properties, occupancies, B-factors,
    /// velocities, or forces.
    pub fn new(positions: Positions) -> Self {
        Self::from_structure(Arc::new(SampleStructure::new(positions)))
    }

    /// A sample holding a model's positions, cell, occupancies, B-factors,
    /// and structure properties, moving them without copying; the model's
    /// topology is dropped.
    pub fn from_model(model: Model) -> Self {
        Self::from_structure(Arc::new(SampleStructure::from_model(model)))
    }

    pub(crate) const fn from_structure(structure: Arc<SampleStructure>) -> Self {
        Self { structure }
    }

    pub(crate) fn structure(&self) -> &Arc<SampleStructure> {
        &self.structure
    }

    pub(crate) fn structure_mut(&mut self) -> &mut Arc<SampleStructure> {
        &mut self.structure
    }

    pub(crate) fn into_structure(self) -> Arc<SampleStructure> {
        self.structure
    }

    /// This sample projected onto the dense atom and bond rows of a
    /// topology subset.
    pub(crate) fn projected(
        &self,
        atoms: &[usize],
        bonds: &[usize],
    ) -> Result<Self, StructureError> {
        Ok(Self::from_structure(Arc::new(
            self.structure.project(atoms, bonds)?,
        )))
    }

    /// Replaces the positions without checking them again, copying the rest
    /// of a shared structure but never the positions.
    pub(crate) fn replace_checked_positions(&mut self, positions: Positions) {
        SampleStructure::replace_positions(&mut self.structure, positions);
    }

    /// Returns this sample with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.set_cell(Some(cell));
        self
    }

    /// Number of atoms, the length of every dense per-atom array of the
    /// sample.
    pub fn atom_count(&self) -> usize {
        self.structure.positions.len()
    }

    /// Whether this sample and `other` share one structure's storage, so
    /// that neither has copied the other's coordinates. Clones share until
    /// one of them is edited.
    pub fn shares_structure_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.structure, &other.structure)
    }

    /// Per-atom velocities (nm/ps), or `None` when the sample has none.
    pub fn velocities(&self) -> Option<&Velocities> {
        self.structure.velocities.as_ref()
    }

    /// Per-atom forces (kJ/(mol·nm)), or `None` when the sample has none.
    pub fn forces(&self) -> Option<&Forces> {
        self.structure.forces.as_ref()
    }

    /// Sets or clears velocities; they must have one vector per atom.
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), StructureError> {
        match &velocities {
            Some(values) => self.structure.check_len(values.len())?,
            None if self.structure.velocities.is_none() => return Ok(()),
            None => {}
        }
        Arc::make_mut(&mut self.structure).velocities = velocities;
        Ok(())
    }

    /// Sets or clears forces; they must have one vector per atom.
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), StructureError> {
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
    ) -> Result<(), StructureError> {
        self.structure.positions.validate_all(&positions)?;
        self.edit(|structure| Ok(structure.positions.set_all(positions)?))
    }

    /// Sets one position by dense atom index.
    pub fn set_position_at(
        &mut self,
        atom: TopologyAtomIndex,
        position: Quantity<Point3>,
    ) -> Result<(), StructureError> {
        self.edit(|structure| Ok(structure.positions.set(atom.index(), position)?))
    }

    /// Replaces the periodic cell; `None` makes the sample non-periodic.
    /// Positions are not wrapped or rescaled.
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        if self.structure.cell != cell {
            Arc::make_mut(&mut self.structure).cell = cell;
        }
    }

    /// Replaces every structure annotation with rows of the same
    /// dimensions, for example ones copied from another structure of the
    /// same topology. A sample not yet in a set accepts no bond rows.
    pub fn set_properties(
        &mut self,
        properties: StructureProperties,
    ) -> Result<(), StructureError> {
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
    pub fn properties_mut(&mut self) -> StructurePropertiesMut<'_> {
        StructurePropertiesMut::new(&mut Arc::make_mut(&mut self.structure).properties)
    }

    /// Sets or clears one occupancy. Values must be finite.
    pub fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), StructureError> {
        self.edit(|structure| {
            let len = structure.positions.len();
            structure.occupancies.set(atom, len, value)
        })
    }

    /// Replaces every occupancy in dense atom order, or clears them all.
    pub fn set_occupancies(
        &mut self,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), StructureError> {
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
    ) -> Result<(), StructureError> {
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
    ) -> Result<(), StructureError> {
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
    ) -> Result<(), StructureError> {
        self.structure = Arc::new(self.structure.transformed(*transform)?);
        Ok(())
    }

    fn storage_ref(&self) -> &SampleStructure {
        &self.structure
    }

    /// Runs an edit that changes nothing when it fails. A shared structure
    /// is edited as a copy that replaces it only on success, so a failed
    /// edit keeps the structure shared.
    fn edit<R>(
        &mut self,
        edit: impl FnOnce(&mut SampleStructure) -> Result<R, StructureError>,
    ) -> Result<R, StructureError> {
        if let Some(owned) = Arc::get_mut(&mut self.structure) {
            return edit(owned);
        }
        let mut copy = SampleStructure::clone(&self.structure);
        let result = edit(&mut copy)?;
        self.structure = Arc::new(copy);
        Ok(result)
    }

    structure_readers!();
    site_readers!();
}

/// Dimension-preserving mutable access to one sample of a
/// [`super::SampleSet`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// sample through [`super::SampleSet::replace`]. The sample's place and
/// weight belong to the set ([`super::SampleSet::set_place`],
/// [`super::DistributionMut::set_weight`]).
#[derive(Debug)]
pub struct SampleMut<'a> {
    sample: &'a mut Sample,
}

impl<'a> SampleMut<'a> {
    pub(crate) fn new(sample: &'a mut Sample) -> Self {
        Self { sample }
    }
}

impl std::ops::Deref for SampleMut<'_> {
    type Target = Sample;

    fn deref(&self) -> &Self::Target {
        self.sample
    }
}

impl SampleMut<'_> {
    /// Sets or clears velocities; see [`Sample::set_velocities`].
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), StructureError> {
        self.sample.set_velocities(velocities)
    }

    /// Sets or clears forces; see [`Sample::set_forces`].
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), StructureError> {
        self.sample.set_forces(forces)
    }

    /// Removes velocities, returning them.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.sample.take_velocities()
    }

    /// Removes forces, returning them.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.sample.take_forces()
    }

    /// Applies a rigid transform transactionally: positions move; the cell,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        self.sample.apply_rigid_transform(transform)
    }

    /// Replaces every position; see [`Sample::set_positions`].
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), StructureError> {
        self.sample.set_positions(positions)
    }

    /// Sets one position; see [`Sample::set_position_at`].
    pub fn set_position_at(
        &mut self,
        atom: TopologyAtomIndex,
        position: Quantity<Point3>,
    ) -> Result<(), StructureError> {
        self.sample.set_position_at(atom, position)
    }

    /// Replaces the periodic cell; see [`Sample::set_cell`].
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.sample.set_cell(cell);
    }

    /// Replaces every structure annotation; see [`Sample::set_properties`].
    pub fn set_properties(
        &mut self,
        properties: StructureProperties,
    ) -> Result<(), StructureError> {
        self.sample.set_properties(properties)
    }

    /// Length-preserving mutable access to structure annotations; see
    /// [`Sample::properties_mut`].
    pub fn properties_mut(&mut self) -> StructurePropertiesMut<'_> {
        self.sample.properties_mut()
    }

    /// Sets or clears one occupancy; see [`Sample::set_occupancy`].
    pub fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), StructureError> {
        self.sample.set_occupancy(atom, value)
    }

    /// Replaces or clears every occupancy; see [`Sample::set_occupancies`].
    pub fn set_occupancies(
        &mut self,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), StructureError> {
        self.sample.set_occupancies(values)
    }

    /// Sets or clears one B-factor; see [`Sample::set_b_factor`].
    pub fn set_b_factor(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), StructureError> {
        self.sample.set_b_factor(atom, value)
    }

    /// Replaces or clears every B-factor; see [`Sample::set_b_factors`].
    pub fn set_b_factors(
        &mut self,
        values: Option<Quantity<Vec<Option<f64>>>>,
    ) -> Result<(), StructureError> {
        self.sample.set_b_factors(values)
    }
}
