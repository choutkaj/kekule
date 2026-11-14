//! [`GeometryData`]: the state a [`super::Geometry`] shares, and the row
//! operations of drafts that stage it.

use crate::diagnostics::Countable;
use crate::diagnostics::LengthMismatch;
use std::sync::Arc;

use crate::math::{PeriodicCell, RigidTransform};
use crate::properties::GeometryProperties;
use crate::structure::geometry_validation::{is_bound, validate_binding};
use crate::structure::scalars::AtomScalars;
use crate::structure::{DenseArrayError, Forces, GeometryError, Positions, Velocities};
use crate::topology::{Topology, TopologySubset};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::Geometry;

/// The state of a geometry: shared between clones of a [`Geometry`] and
/// copied on write.
///
/// Builders and editors stage one in their own row order, which for an
/// editor includes rows of deleted slots; publication projects it into
/// dense order with [`Self::project`].
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct GeometryData {
    pub(crate) positions: Positions,
    pub(crate) cell: Option<PeriodicCell>,
    pub(crate) occupancies: AtomScalars,
    pub(crate) b_factors: AtomScalars,
    pub(crate) properties: GeometryProperties,
    pub(crate) velocities: Option<Velocities>,
    pub(crate) forces: Option<Forces>,
}

impl Default for GeometryData {
    fn default() -> Self {
        Self::new(Positions::default())
    }
}

impl GeometryData {
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

    /// Projects every row through a correspondence: `atoms[new]` and
    /// `bonds[new]` are source rows.
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

    /// Replaces the positions of a possibly shared state, copying the
    /// rest of it when it is shared but never the positions being replaced.
    pub(crate) fn replace_positions(data: &mut Arc<Self>, positions: Positions) {
        match Arc::get_mut(data) {
            Some(owned) => owned.positions = positions,
            None => {
                let shared = &**data;
                *data = Arc::new(Self {
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
    pub(crate) fn check_len(&self, actual: usize) -> Result<(), GeometryError> {
        if actual != self.positions.len() {
            return Err(GeometryError::Array(DenseArrayError::CountMismatch(
                LengthMismatch::new(Countable::Row, self.positions.len(), actual),
            )));
        }
        Ok(())
    }

    /// Binds a possibly shared state to `topology`, allocating its bond
    /// rows. A bound state, such as one taken from another set of this
    /// topology, stays shared.
    pub(crate) fn bind(data: &mut Arc<Self>, topology: &Topology) -> Result<(), GeometryError> {
        let atoms = data.positions.len();
        validate_binding(atoms, &data.properties, topology)?;
        if !is_bound(atoms, &data.properties, topology) {
            Arc::make_mut(data)
                .properties
                .bind_bond_rows(topology.bond_count());
        }
        Ok(())
    }

    pub(crate) fn try_reserve_atoms(&mut self, additional: usize) -> Result<(), GeometryError> {
        self.positions
            .try_reserve(additional)
            .map_err(|_| GeometryError::CapacityOverflow)
    }

    /// Appends atoms and bond rows for staged construction. Owner
    /// properties are kept; callers apply their own structural-edit policy.
    /// Velocities and forces cover the earlier atoms only, so they are
    /// cleared.
    pub(crate) fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        self.positions.extend_canonical(positions);
        let atoms = self.positions.len();
        self.occupancies.resize(atoms);
        self.b_factors.resize(atoms);
        let bonds = self.properties.bonds().len() + added_bonds;
        self.properties.resize(atoms, bonds);
        self.clear_dynamics_unless(positions.is_empty());
    }

    /// Grows private editor slot rows; new atom slots start at the origin
    /// and are positioned by the caller. New atom slots clear velocities and
    /// forces, which cover the earlier atoms only.
    pub(crate) fn resize_slots(&mut self, atoms: usize, bonds: usize) {
        debug_assert!(atoms >= self.positions.len());
        let grown = atoms != self.positions.len();
        self.positions.resize_canonical(atoms);
        self.occupancies.resize(atoms);
        self.b_factors.resize(atoms);
        self.properties.resize(atoms, bonds);
        self.clear_dynamics_unless(!grown);
    }

    fn clear_dynamics_unless(&mut self, keep: bool) {
        if !keep {
            self.velocities = None;
            self.forces = None;
        }
    }

    /// Clears annotations of the given dead editor atom slots, at a cost
    /// proportional to their number.
    pub(crate) fn clear_atom_rows(&mut self, slots: &[usize]) {
        for &slot in slots {
            self.occupancies.clear(slot);
            self.b_factors.clear(slot);
            self.properties.atoms_raw_mut().clear_index(slot);
        }
    }

    /// Clears annotations of the given dead editor bond slots.
    pub(crate) fn clear_bond_rows(&mut self, slots: &[usize]) {
        for &slot in slots {
            self.properties.bonds_raw_mut().clear_index(slot);
        }
    }

    /// Copies positions, occupancies, and B-factors of every atom of
    /// `source` into the given destination rows. Property rows are copied
    /// by the caller.
    pub(crate) fn copy_atom_state_from(&mut self, source: &Self, rows: &[usize]) {
        debug_assert_eq!(rows.len(), source.positions.len());
        let staged = rows
            .iter()
            .copied()
            .zip(source.positions.values().value().iter().copied())
            .collect::<Vec<_>>();
        self.positions
            .set_canonical_batch(&staged)
            .expect("source positions are finite and rows are allocated");
        let len = self.positions.len();
        self.occupancies
            .copy_rows_from(&source.occupancies, rows, len);
        self.b_factors.copy_rows_from(&source.b_factors, rows, len);
    }

    /// Reorders atom rows; `previous[new]` is the old row.
    pub(crate) fn reorder_atoms(&mut self, previous: &[usize]) {
        self.positions = self
            .positions
            .select_indices(previous)
            .expect("staged positions follow staged atoms");
        self.occupancies = self.occupancies.select(previous);
        self.b_factors = self.b_factors.select(previous);
        self.velocities = self
            .velocities
            .as_ref()
            .map(|values| values.select(previous));
        self.forces = self.forces.as_ref().map(|values| values.select(previous));
        self.properties
            .reorder_atoms(previous)
            .expect("staged structure atom rows follow staged atoms");
    }
}

impl Geometry {
    pub(crate) fn from_data(data: GeometryData) -> Self {
        Self {
            data: Arc::new(data),
        }
    }

    /// The data, moved out when only this geometry holds it.
    pub(crate) fn into_data(self) -> GeometryData {
        Arc::unwrap_or_clone(self.data)
    }

    pub(crate) fn data(&self) -> &Arc<GeometryData> {
        &self.data
    }

    pub(crate) fn data_mut(&mut self) -> &mut Arc<GeometryData> {
        &mut self.data
    }

    /// The data for editing in place: data shared with a clone is copied
    /// first, once. The caller keeps it bound and valid.
    pub(crate) fn data_in_place(&mut self) -> &mut GeometryData {
        Arc::make_mut(&mut self.data)
    }
}

/// Source dense atom and bond rows of a subset, in subset dense order.
pub(crate) fn subset_rows(subset: &TopologySubset) -> (Vec<usize>, Vec<usize>) {
    let correspondence = subset.correspondence();
    let atoms = correspondence
        .source_atom_indices()
        .iter()
        .map(|index| index.index())
        .collect();
    let bonds = correspondence
        .source_bond_indices()
        .iter()
        .map(|index| index.index())
        .collect();
    (atoms, bonds)
}

/// Applies a rigid transform to every position.
fn transform_positions(
    positions: &Positions,
    transform: RigidTransform,
) -> Result<Positions, GeometryError> {
    let points = positions
        .values()
        .value()
        .iter()
        .map(|point| transform.transform_point(*point))
        .collect::<Vec<_>>();
    Ok(Positions::new(Quantity::new(
        points,
        CANONICAL_LENGTH_UNIT,
    ))?)
}

/// Rotates cell vectors; a rigid transform's translation does not move a cell.
pub(crate) fn transform_cell(
    cell: Option<PeriodicCell>,
    transform: RigidTransform,
) -> Result<Option<PeriodicCell>, GeometryError> {
    cell.map(|cell| {
        let vectors = cell
            .vectors()
            .map(|vectors| vectors.map(|vector| transform.transform_vector(vector)));
        Ok(PeriodicCell::new(vectors, cell.periodic_axes())?)
    })
    .transpose()
}
