//! One-realization coordination over the coordinate-free structural editor.
use super::{Model, ModelError, PositionError, Positions};
use crate::core::{Atom, BondOrder, Molecule};
use crate::geometry::{PeriodicCell, Point3};
use crate::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use crate::topology::{
    AtomSiteId, AtomSiteMetadata, ChainId, EditAtomId, EditAtomSite, EditAtomSiteId, EditBond,
    EditBondId, EditChain, EditChainId, EditMolecule, EditResidue, EditResidueId, InstanceAtomId,
    InstanceBondId, MoleculeClass, MoleculeInstanceId, ResidueClass, ResidueId,
    TopologyEditCorrespondence, TopologyEditError, TopologyEditor,
};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};
use std::fmt;

/// Detached structural editing for one geometry-bearing molecular system.
///
/// Positions accompany atom insertion and follow stable editing handles through
/// deletion, splitting and merging. The coordinate-free editor is exposed only
/// for inspection. Use the coordinated methods here for structural changes.
/// Property methods refer to this realization; explicitly named `topology_*`
/// methods edit static annotations. Generic coordinate changes preserve stored
/// annotations without asserting that arbitrary derived values remain valid.
///
/// ```
/// use kekule::{core::{Atom, Element, BondOrder}, geometry::Point3,
///     structure::ModelEditor, units::{Quantity, ANGSTROM}};
/// let mut editor = ModelEditor::new();
/// let c = editor.add_atom(Atom::new(Element::from_symbol("C").unwrap()),
///     Quantity::new(Point3::origin(), ANGSTROM))?;
/// let o = editor.add_atom(Atom::new(Element::from_symbol("O").unwrap()),
///     Quantity::new(Point3::new(1.4, 0.0, 0.0), ANGSTROM))?;
/// editor.add_bond(c, o, BondOrder::Single)?;
/// let result = editor.finish_with_correspondence()?;
/// assert_eq!(result.model().topology().instance_count(), 1);
/// assert!(result.correspondence().atom(o).is_some());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct ModelEditor {
    topology: TopologyEditor,
    positions: Vec<Point3>,
    cell: Option<PeriodicCell>,
    properties: Properties,
}
impl Model {
    pub fn edit(&self) -> ModelEditor {
        ModelEditor::from_model(self)
    }
    /// Moves realization arrays into a draft while retaining the shared topology.
    pub fn into_editor(self) -> ModelEditor {
        ModelEditor {
            topology: TopologyEditor::from_topology(self.topology),
            positions: self.positions.into_canonical_values(),
            cell: self.cell,
            properties: self.properties,
        }
    }
}

impl ModelEditor {
    pub fn new() -> Self {
        Self::default()
    }
    /// Clears structure, geometry, hierarchy, cell and annotations in the draft.
    pub fn clear(&mut self) {
        self.topology.clear();
        self.positions.clear();
        self.cell = None;
        self.properties = Properties::new();
    }
    pub fn from_model(model: &Model) -> Self {
        model.clone().into_editor()
    }
    pub fn topology_editor(&self) -> &TopologyEditor {
        &self.topology
    }
    pub fn atom_count(&self) -> usize {
        self.topology.atom_count()
    }
    pub fn bond_count(&self) -> usize {
        self.topology.bond_count()
    }
    pub fn is_empty(&self) -> bool {
        self.topology.is_empty()
    }
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = EditAtomId> + '_ {
        self.topology.atom_ids()
    }
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = EditBondId> + '_ {
        self.topology.bond_ids()
    }
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (EditAtomId, &Atom)> {
        self.topology.atoms()
    }
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (EditBondId, EditBond)> + '_ {
        self.topology.bonds()
    }
    pub fn atom(&self, id: EditAtomId) -> Result<&Atom, ModelEditError> {
        Ok(self.topology.atom(id)?)
    }
    pub fn bond(&self, id: EditBondId) -> Result<EditBond, ModelEditError> {
        Ok(self.topology.bond(id)?)
    }
    pub fn atom_handle(&self, source: InstanceAtomId) -> Result<EditAtomId, ModelEditError> {
        Ok(self.topology.atom_handle(source)?)
    }
    pub fn bond_handle(&self, source: InstanceBondId) -> Result<EditBondId, ModelEditError> {
        Ok(self.topology.bond_handle(source)?)
    }
    pub fn neighbors(
        &self,
        id: EditAtomId,
    ) -> Result<impl Iterator<Item = EditAtomId> + '_, ModelEditError> {
        Ok(self.topology.neighbors(id)?)
    }
    pub fn incident_bonds(
        &self,
        id: EditAtomId,
    ) -> Result<impl Iterator<Item = EditBondId> + '_, ModelEditError> {
        Ok(self.topology.incident_bonds(id)?)
    }
    pub fn bond_between(
        &self,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<Option<EditBondId>, ModelEditError> {
        Ok(self.topology.bond_between(a, b)?)
    }
    pub fn connected_components(&self) -> Vec<Vec<EditAtomId>> {
        self.topology.connected_components()
    }

    pub fn add_atom(
        &mut self,
        atom: Atom,
        position: Quantity<Point3>,
    ) -> Result<EditAtomId, ModelEditError> {
        let point = checked_point(position)?;
        self.positions
            .try_reserve(1)
            .map_err(|_| ModelEditError::CapacityOverflow)?;
        let id = self.structural(|topology| topology.add_atom(atom))?;
        debug_assert_eq!(self.topology.atom_slot(id)?, self.positions.len());
        self.positions.push(point);
        Ok(id)
    }
    pub fn add_molecule(
        &mut self,
        molecule: &Molecule,
        positions: &Positions,
    ) -> Result<EditMolecule, ModelEditError> {
        if positions.len() != molecule.atom_count() {
            return Err(PositionError::PositionCountMismatch {
                expected: molecule.atom_count(),
                actual: positions.len(),
            }
            .into());
        }
        self.positions
            .try_reserve(positions.len())
            .map_err(|_| ModelEditError::CapacityOverflow)?;
        let added = self.structural(|topology| topology.add_molecule(molecule))?;
        self.positions.extend_from_slice(positions.values().value());
        Ok(added)
    }
    pub fn replace_atom(&mut self, id: EditAtomId, atom: Atom) -> Result<Atom, ModelEditError> {
        self.structural(|t| t.replace_atom(id, atom))
    }
    pub fn delete_atom(&mut self, id: EditAtomId) -> Result<Atom, ModelEditError> {
        self.structural(|t| t.delete_atom(id))
    }
    pub fn delete_atoms(
        &mut self,
        ids: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<Vec<(EditAtomId, Atom)>, ModelEditError> {
        self.structural(|t| t.delete_atoms(ids))
    }
    pub fn retain_atoms(
        &mut self,
        ids: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<Vec<(EditAtomId, Atom)>, ModelEditError> {
        self.structural(|t| t.retain_atoms(ids))
    }
    /// Deletes surviving atoms of a source occurrence and their geometry.
    /// See [`TopologyEditor::delete_instance`] for semantics after splits/merges.
    pub fn delete_instance(
        &mut self,
        id: MoleculeInstanceId,
    ) -> Result<Vec<(EditAtomId, Atom)>, ModelEditError> {
        self.structural(|t| t.delete_instance(id))
    }
    pub fn add_bond(
        &mut self,
        a: EditAtomId,
        b: EditAtomId,
        order: BondOrder,
    ) -> Result<EditBondId, ModelEditError> {
        self.structural(|t| t.add_bond(a, b, order))
    }
    pub fn delete_bond(&mut self, id: EditBondId) -> Result<EditBond, ModelEditError> {
        self.structural(|t| t.delete_bond(id))
    }
    pub fn delete_bonds(
        &mut self,
        ids: impl IntoIterator<Item = EditBondId>,
    ) -> Result<Vec<(EditBondId, EditBond)>, ModelEditError> {
        self.structural(|t| t.delete_bonds(ids))
    }
    pub fn set_bond_order(
        &mut self,
        id: EditBondId,
        order: BondOrder,
    ) -> Result<(), ModelEditError> {
        self.structural(|t| t.set_bond_order(id, order))
    }
    pub fn set_bond_endpoints(
        &mut self,
        id: EditBondId,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<(), ModelEditError> {
        self.structural(|t| t.set_bond_endpoints(id, a, b))
    }
    pub fn replace_bond(
        &mut self,
        id: EditBondId,
        replacement: EditBond,
    ) -> Result<EditBond, ModelEditError> {
        self.structural(|t| t.replace_bond(id, replacement))
    }
    pub fn set_molecule_class(
        &mut self,
        atom: EditAtomId,
        class: MoleculeClass,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_molecule_class(atom, class)?)
    }

    /// Copies dense positions in current live atom-handle order.
    pub fn positions(&self) -> Positions {
        Positions::from_canonical_values(
            self.atom_ids()
                .map(|id| self.positions[self.topology.atom_slot(id).unwrap()])
                .collect(),
        )
    }
    pub fn position(&self, id: EditAtomId) -> Result<Quantity<Point3>, ModelEditError> {
        Ok(Quantity::new(
            self.positions[self.topology.atom_slot(id)?],
            CANONICAL_LENGTH_UNIT,
        ))
    }
    pub fn set_position(
        &mut self,
        id: EditAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        let point = checked_point(position)?;
        self.positions[slot] = point;
        Ok(())
    }
    /// Replaces all live positions transactionally in [`Self::atom_ids`] order.
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), ModelEditError> {
        let positions = Positions::new(positions)?;
        if positions.len() != self.atom_count() {
            return Err(PositionError::PositionCountMismatch {
                expected: self.atom_count(),
                actual: positions.len(),
            }
            .into());
        }
        let slots = self
            .atom_ids()
            .map(|id| self.topology.atom_slot(id).unwrap())
            .collect::<Vec<_>>();
        for (slot, &point) in slots.into_iter().zip(positions.values().value().iter()) {
            self.positions[slot] = point;
        }
        Ok(())
    }
    /// Checked sparse coordinate batch; repeated handles are applied in input order.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (EditAtomId, Quantity<Point3>)>,
    ) -> Result<(), ModelEditError> {
        let updates = values
            .into_iter()
            .map(|(id, value)| Ok((self.topology.atom_slot(id)?, checked_point(value)?)))
            .collect::<Result<Vec<_>, ModelEditError>>()?;
        for (slot, point) in updates {
            self.positions[slot] = point;
        }
        Ok(())
    }
    pub fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }
    /// Realization properties in private stable-slot order. Column helpers use live order.
    pub fn properties(&self) -> &Properties {
        &self.properties
    }
    pub fn atom_properties(&self) -> &PropertyTable {
        self.properties.atoms()
    }
    pub fn bond_properties(&self) -> &PropertyTable {
        self.properties.bonds()
    }
    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, ModelEditError> {
        Ok(self.properties.insert(key, value)?)
    }
    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.properties.remove(key)
    }
    pub fn clear_properties(&mut self) {
        self.properties.clear_owner();
    }
    pub fn atom_property(
        &self,
        id: EditAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelEditError> {
        Ok(self
            .properties
            .atoms()
            .value(key, self.topology.atom_slot(id)?)?)
    }
    pub fn bond_property(
        &self,
        id: EditBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelEditError> {
        Ok(self
            .properties
            .bonds()
            .value(key, self.topology.bond_slot(id)?)?)
    }
    pub fn set_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self
            .properties
            .set_realization_atom_value(key, slot, value)?)
    }
    pub fn set_bond_property(
        &mut self,
        id: EditBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.bond_slot(id)?;
        Ok(self
            .properties
            .set_realization_bond_value(key, slot, value)?)
    }
    pub fn set_atom_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditAtomId, Option<PropertyValue>)>,
    ) -> Result<(), ModelEditError> {
        let mut staged = self.properties.clone();
        for (id, value) in values {
            staged.set_realization_atom_value(key.clone(), self.topology.atom_slot(id)?, value)?;
        }
        self.properties = staged;
        Ok(())
    }
    pub fn set_bond_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditBondId, Option<PropertyValue>)>,
    ) -> Result<(), ModelEditError> {
        let mut staged = self.properties.clone();
        for (id, value) in values {
            staged.set_realization_bond_value(key.clone(), self.topology.bond_slot(id)?, value)?;
        }
        self.properties = staged;
        Ok(())
    }
    pub fn atom_property_column(
        &self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        Ok(self
            .atom_properties()
            .select_indices(&self.atom_slots())?
            .remove(key))
    }
    pub fn bond_property_column(
        &self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        Ok(self
            .bond_properties()
            .select_indices(&self.bond_slots())?
            .remove(key))
    }
    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        let previous = self.atom_property_column(&key)?;
        let column =
            column.into_editor_slots(&self.atom_slots(), self.topology.atom_slot_count())?;
        self.properties
            .insert_realization_atom_column(key, column)?;
        Ok(previous)
    }
    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        let previous = self.bond_property_column(&key)?;
        let column =
            column.into_editor_slots(&self.bond_slots(), self.topology.bond_slot_count())?;
        self.properties
            .insert_realization_bond_column(key, column)?;
        Ok(previous)
    }
    pub fn remove_atom_property_column(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        let previous = self.atom_property_column(key)?;
        self.properties.remove_realization_atom_column(key)?;
        Ok(previous)
    }
    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.bond_property_column(key).expect("live property slots");
        self.properties.remove_realization_bond_column(key);
        previous
    }
    pub fn occupancy(&self, id: EditAtomId) -> Result<Option<f64>, ModelEditError> {
        Ok(self.properties.occupancy_at(self.topology.atom_slot(id)?)?)
    }
    pub fn set_occupancy(
        &mut self,
        id: EditAtomId,
        value: Option<f64>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.properties.set_occupancy_at(slot, value)?)
    }
    pub fn b_factor(&self, id: EditAtomId) -> Result<Option<Quantity<f64>>, ModelEditError> {
        Ok(self.properties.b_factor_at(self.topology.atom_slot(id)?)?)
    }
    pub fn set_b_factor(
        &mut self,
        id: EditAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.properties.set_b_factor_at(slot, value)?)
    }

    pub fn chains(&self) -> impl ExactSizeIterator<Item = (EditChainId, &EditChain)> {
        self.topology.chains()
    }
    pub fn residues(&self) -> impl ExactSizeIterator<Item = (EditResidueId, &EditResidue)> {
        self.topology.residues()
    }
    pub fn atom_sites(&self) -> impl ExactSizeIterator<Item = (EditAtomSiteId, &EditAtomSite)> {
        self.topology.atom_sites()
    }
    pub fn chain(&self, id: EditChainId) -> Result<&EditChain, ModelEditError> {
        Ok(self.topology.chain(id)?)
    }
    pub fn residue(&self, id: EditResidueId) -> Result<&EditResidue, ModelEditError> {
        Ok(self.topology.residue(id)?)
    }
    pub fn atom_site(&self, id: EditAtomSiteId) -> Result<&EditAtomSite, ModelEditError> {
        Ok(self.topology.atom_site(id)?)
    }
    pub fn atom_site_for_atom(
        &self,
        atom: EditAtomId,
    ) -> Result<Option<EditAtomSiteId>, ModelEditError> {
        Ok(self.topology.atom_site_for_atom(atom)?)
    }
    pub fn chain_handle(&self, id: ChainId) -> Result<EditChainId, ModelEditError> {
        Ok(self.topology.chain_handle(id)?)
    }
    pub fn residue_handle(&self, id: ResidueId) -> Result<EditResidueId, ModelEditError> {
        Ok(self.topology.residue_handle(id)?)
    }
    pub fn atom_site_handle(&self, id: AtomSiteId) -> Result<EditAtomSiteId, ModelEditError> {
        Ok(self.topology.atom_site_handle(id)?)
    }
    pub fn add_chain(
        &mut self,
        label: impl Into<String>,
        author: Option<String>,
    ) -> Result<EditChainId, ModelEditError> {
        self.structural(|t| t.add_chain(label, author))
    }
    pub fn add_residue(
        &mut self,
        chain: EditChainId,
        name: impl Into<String>,
        label_seq: Option<i32>,
        author_seq: Option<String>,
        insertion: Option<String>,
    ) -> Result<EditResidueId, ModelEditError> {
        self.structural(|t| t.add_residue(chain, name, label_seq, author_seq, insertion))
    }
    pub fn add_atom_site(
        &mut self,
        residue: EditResidueId,
        atom: EditAtomId,
        metadata: AtomSiteMetadata,
    ) -> Result<EditAtomSiteId, ModelEditError> {
        self.structural(|t| t.add_atom_site(residue, atom, metadata))
    }
    pub fn set_residue_class(
        &mut self,
        id: EditResidueId,
        class: ResidueClass,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_residue_class(id, class)?)
    }
    pub fn set_residue_component_ids(
        &mut self,
        id: EditResidueId,
        label: Option<String>,
        author: Option<String>,
    ) -> Result<(), ModelEditError> {
        self.structural(|t| t.set_residue_component_ids(id, label, author))
    }
    pub fn set_chain_identifiers(
        &mut self,
        id: EditChainId,
        label: impl Into<String>,
        author: Option<String>,
    ) -> Result<(), ModelEditError> {
        self.structural(|t| t.set_chain_identifiers(id, label, author))
    }
    pub fn set_atom_site_metadata(
        &mut self,
        id: EditAtomSiteId,
        metadata: AtomSiteMetadata,
    ) -> Result<(), ModelEditError> {
        self.structural(|t| t.set_atom_site_metadata(id, metadata))
    }
    pub fn set_atom_site_residue(
        &mut self,
        id: EditAtomSiteId,
        residue: EditResidueId,
    ) -> Result<(), ModelEditError> {
        self.structural(|t| t.set_atom_site_residue(id, residue))
    }
    pub fn delete_chain(&mut self, id: EditChainId) -> Result<EditChain, ModelEditError> {
        self.structural(|t| t.delete_chain(id))
    }
    pub fn delete_residue(&mut self, id: EditResidueId) -> Result<EditResidue, ModelEditError> {
        self.structural(|t| t.delete_residue(id))
    }
    pub fn delete_atom_site(&mut self, id: EditAtomSiteId) -> Result<EditAtomSite, ModelEditError> {
        self.structural(|t| t.delete_atom_site(id))
    }
    pub fn insert_topology_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, ModelEditError> {
        Ok(self.topology.insert_property(key, value)?)
    }
    pub fn set_topology_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_atom_property(id, key, value)?)
    }
    pub fn remove_topology_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.topology.remove_property(key)
    }
    pub fn clear_topology_properties(&mut self) {
        self.topology.clear_properties();
    }
    pub fn set_topology_bond_property(
        &mut self,
        id: EditBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_bond_property(id, key, value)?)
    }
    pub fn set_chain_property(
        &mut self,
        id: EditChainId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_chain_property(id, key, value)?)
    }
    pub fn set_residue_property(
        &mut self,
        id: EditResidueId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_residue_property(id, key, value)?)
    }
    pub fn set_atom_site_property(
        &mut self,
        id: EditAtomSiteId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_atom_site_property(id, key, value)?)
    }
    pub fn set_definition_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_definition_atom_property(id, key, value)?)
    }
    pub fn set_definition_bond_property(
        &mut self,
        id: EditBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        Ok(self.topology.set_definition_bond_property(id, key, value)?)
    }

    pub fn validate(&self) -> Result<(), ModelEditError> {
        self.clone().finish().map(|_| ())
    }
    pub fn finish(self) -> Result<Model, ModelEditError> {
        self.finish_with_correspondence().map(|r| r.model)
    }
    pub fn finish_with_correspondence(self) -> Result<ModelEdit, ModelEditError> {
        let (topology, correspondence) = self.topology.finish_with_correspondence()?.into_parts();
        let positions = Positions::from_canonical_values(
            correspondence
                .atom_slots
                .iter()
                .map(|&slot| self.positions[slot])
                .collect(),
        );
        let mut properties = self
            .properties
            .project_realization(&correspondence.atom_slots, &correspondence.bond_slots)?;
        for (key, value) in self.properties.iter() {
            properties.insert(key.clone(), value.clone())?;
        }
        let model = Model::with_properties(topology, positions, self.cell, properties)?;
        Ok(ModelEdit {
            model,
            correspondence,
        })
    }
    pub fn try_finish(self) -> Result<Model, ModelFinishError> {
        self.try_finish_with_correspondence().map(|r| r.model)
    }
    pub fn try_finish_with_correspondence(self) -> Result<ModelEdit, ModelFinishError> {
        let snapshot = self.clone();
        self.finish_with_correspondence()
            .map_err(|error| ModelFinishError {
                error: Box::new(error),
                editor: Box::new(snapshot),
            })
    }

    fn structural<T>(
        &mut self,
        edit: impl FnOnce(&mut TopologyEditor) -> Result<T, TopologyEditError>,
    ) -> Result<T, ModelEditError> {
        let revision = self.topology.structural_revision;
        let atom_count = self.atom_count();
        let bond_count = self.bond_count();
        let result = edit(&mut self.topology)?;
        self.properties
            .resize_atoms(self.topology.atom_slot_count());
        self.properties
            .resize_bonds(self.topology.bond_slot_count());
        if self.topology.structural_revision != revision {
            self.properties.clear_owner();
        }
        if self.atom_count() < atom_count && self.properties.atoms().has_data() {
            let live = self
                .atom_slots()
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>();
            for slot in 0..self.properties.atoms().len() {
                if !live.contains(&slot) {
                    self.properties.atoms_mut().clear_index(slot);
                }
            }
        }
        if self.bond_count() < bond_count && self.properties.bonds().has_data() {
            let live = self
                .bond_slots()
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>();
            for slot in 0..self.properties.bonds().len() {
                if !live.contains(&slot) {
                    self.properties.bonds_mut().clear_index(slot);
                }
            }
        }
        Ok(result)
    }
    fn atom_slots(&self) -> Vec<usize> {
        self.atom_ids()
            .map(|id| self.topology.atom_slot(id).unwrap())
            .collect()
    }
    fn bond_slots(&self) -> Vec<usize> {
        self.bond_ids()
            .map(|id| self.topology.bond_slot(id).unwrap())
            .collect()
    }
}

fn checked_point(position: Quantity<Point3>) -> Result<Point3, PositionError> {
    let point = position.to_unit(CANONICAL_LENGTH_UNIT)?.to_value();
    if !point.is_finite() {
        return Err(PositionError::NonFinitePosition { index: 0 });
    }
    Ok(point)
}

#[derive(Debug, Clone)]
pub struct ModelEdit {
    model: Model,
    correspondence: TopologyEditCorrespondence,
}
impl ModelEdit {
    pub fn model(&self) -> &Model {
        &self.model
    }
    pub fn correspondence(&self) -> &TopologyEditCorrespondence {
        &self.correspondence
    }
    pub fn into_parts(self) -> (Model, TopologyEditCorrespondence) {
        (self.model, self.correspondence)
    }
}
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelEditError {
    Topology(TopologyEditError),
    Position(PositionError),
    Property(PropertyError),
    Model(Box<ModelError>),
    CapacityOverflow,
}
impl fmt::Display for ModelEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Topology(e) => e.fmt(f),
            Self::Position(e) => e.fmt(f),
            Self::Property(e) => e.fmt(f),
            Self::Model(e) => e.fmt(f),
            Self::CapacityOverflow => f.write_str("model editing exceeds coordinate capacity"),
        }
    }
}
impl std::error::Error for ModelEditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(e) => Some(e),
            Self::Position(e) => Some(e),
            Self::Property(e) => Some(e),
            Self::Model(e) => Some(e.as_ref()),
            Self::CapacityOverflow => None,
        }
    }
}
macro_rules! convert {
    ($ty:ty, $variant:ident) => {
        impl From<$ty> for ModelEditError {
            fn from(e: $ty) -> Self {
                Self::$variant(e)
            }
        }
    };
}
convert!(TopologyEditError, Topology);
convert!(PositionError, Position);
convert!(PropertyError, Property);
impl From<ModelError> for ModelEditError {
    fn from(error: ModelError) -> Self {
        Self::Model(Box::new(error))
    }
}
#[derive(Debug)]
pub struct ModelFinishError {
    error: Box<ModelEditError>,
    editor: Box<ModelEditor>,
}
impl ModelFinishError {
    pub fn error(&self) -> &ModelEditError {
        &self.error
    }
    pub fn editor(&self) -> &ModelEditor {
        &self.editor
    }
    pub fn into_editor(self) -> ModelEditor {
        *self.editor
    }
}
impl fmt::Display for ModelFinishError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.error.fmt(f)
    }
}
impl std::error::Error for ModelFinishError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(self.error.as_ref())
    }
}
