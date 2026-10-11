use std::collections::BTreeSet;

use crate::geometry::{PeriodicCell, Point3};
use crate::properties::{RealizationProperties, RealizationPropertiesMut};
use crate::topology::TopologyAtomIndex;
use crate::units::Quantity;

use super::sites::SiteParameters;
use super::{ConformationError, ModelView, Positions};

/// Realization state staged by [`super::ModelBuilder`] and
/// [`super::ModelEditor`] before it becomes a [`super::Model`].
///
/// Rows follow the draft's staged dense order or private editor slot order and
/// may include rows of deleted editor slots. Publication projects and moves
/// these fields into the model's own.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ModelState {
    pub(crate) positions: Positions,
    pub(crate) cell: Option<PeriodicCell>,
    pub(crate) sites: SiteParameters,
    pub(crate) properties: RealizationProperties,
}

impl Default for ModelState {
    fn default() -> Self {
        Self {
            positions: Positions::default(),
            cell: None,
            sites: SiteParameters::default(),
            properties: RealizationProperties::new(0),
        }
    }
}

impl ModelState {
    pub(crate) fn set_position(
        &mut self,
        atom: TopologyAtomIndex,
        position: Quantity<Point3>,
    ) -> Result<(), ConformationError> {
        Ok(self.positions.set_position_at(atom.index(), position)?)
    }

    pub(crate) fn properties_mut(&mut self) -> RealizationPropertiesMut<'_> {
        RealizationPropertiesMut::new(&mut self.properties)
    }

    pub(crate) fn occupancy(
        &self,
        atom: TopologyAtomIndex,
    ) -> Result<Option<f64>, ConformationError> {
        self.sites.occupancy(atom, self.positions.len())
    }

    pub(crate) fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<f64>,
    ) -> Result<(), ConformationError> {
        let len = self.positions.len();
        self.sites.set_occupancy(atom, len, value)
    }

    pub(crate) fn b_factor(
        &self,
        atom: TopologyAtomIndex,
    ) -> Result<Option<Quantity<f64>>, ConformationError> {
        self.sites.b_factor(atom, self.positions.len())
    }

    pub(crate) fn set_b_factor(
        &mut self,
        atom: TopologyAtomIndex,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ConformationError> {
        let len = self.positions.len();
        self.sites.set_b_factor(atom, len, value)
    }

    pub(crate) fn try_reserve_atoms(&mut self, additional: usize) -> Result<(), ConformationError> {
        self.positions
            .try_reserve(additional)
            .map_err(|_| ConformationError::CapacityOverflow)
    }

    /// Appends atoms and bond rows for staged construction. Owner properties
    /// are kept; callers apply their own structural-edit policy.
    pub(crate) fn extend(&mut self, positions: &Positions, added_bonds: usize) {
        self.positions.extend_canonical(positions);
        let atoms = self.positions.len();
        self.sites.resize(atoms);
        let bonds = self.properties.bonds().len() + added_bonds;
        self.properties.resize(atoms, bonds);
    }

    /// Grows private editor slot rows; new atom slots start at the origin and
    /// are positioned by the caller.
    pub(crate) fn resize_slots(&mut self, atoms: usize, bonds: usize) {
        debug_assert!(atoms >= self.positions.len());
        self.positions.resize_canonical(atoms);
        self.sites.resize(atoms);
        self.properties.resize(atoms, bonds);
    }

    /// Clears annotations of editor atom slots that are no longer live.
    pub(crate) fn clear_dead_atom_rows(&mut self, live: &[usize]) {
        let live = live.iter().copied().collect::<BTreeSet<_>>();
        for slot in (0..self.positions.len()).filter(|slot| !live.contains(slot)) {
            self.sites.clear_row(slot);
            self.properties.atoms_raw_mut().clear_index(slot);
        }
    }

    /// Clears annotations of editor bond slots that are no longer live.
    pub(crate) fn clear_dead_bond_rows(&mut self, live: &[usize]) {
        let live = live.iter().copied().collect::<BTreeSet<_>>();
        for slot in (0..self.properties.bonds().len()).filter(|slot| !live.contains(slot)) {
            self.properties.bonds_raw_mut().clear_index(slot);
        }
    }

    /// Copies positions, occupancies, and B-factors of every source atom into
    /// the given destination rows. Property rows are copied by the caller.
    pub(crate) fn copy_atom_state_from(&mut self, source: ModelView<'_>, rows: &[usize]) {
        debug_assert_eq!(rows.len(), source.atom_count());
        let staged = rows
            .iter()
            .copied()
            .zip(source.positions().values().value().iter().copied())
            .collect::<Vec<_>>();
        self.positions
            .set_canonical_batch(&staged)
            .expect("source positions are finite and rows are allocated");
        if let Some(sites) = source.sites() {
            let len = self.positions.len();
            self.sites.copy_rows_from(sites, rows, len);
        }
    }

    /// Reorders atom rows; `previous[new]` is the old dense index.
    pub(crate) fn reorder_atoms(&mut self, previous: &[usize]) {
        self.positions = self
            .positions
            .select_indices(previous)
            .expect("staged positions follow staged atoms");
        self.sites = self.sites.select(previous);
        self.properties
            .reorder_atoms(previous)
            .expect("staged realization atom rows follow staged atoms");
    }

    /// Projects atom and bond rows into published dense order. Owner
    /// properties are dropped because the projection is a new owner.
    pub(crate) fn project(
        &self,
        atoms: &[usize],
        bonds: &[usize],
    ) -> Result<Self, ConformationError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            sites: self.sites.select(atoms),
            properties: self.properties.project(atoms, bonds)?,
        })
    }
}
