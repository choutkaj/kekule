//! One-realization coordination over the coordinate-free structural editor.
mod append;
pub use append::*;

use super::{Conformation, ConformationError, Model, ModelError, PositionError, Positions};
use crate::core::{Atom, BondOrder, Molecule};
use crate::geometry::{PeriodicCell, Point3};
use crate::properties::{
    OwnerProperties, PropertyColumn, PropertyError, PropertyKey, PropertyValue, RawPropertyTable,
};
use crate::topology::{
    AtomSiteMetadata, EditAtomId, EditAtomSite, EditAtomSiteId, EditBond, EditBondId, EditChain,
    EditChainId, EditMolecule, EditResidue, EditResidueId, MoleculeClass, MoleculeInstanceId,
    ResidueClass, TopologyAtomIndex, TopologyBondIndex, TopologyEditError, TopologyEditor,
};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};
use std::fmt;

/// Detached structural editing for one geometry-bearing molecular system.
///
/// Positions accompany atom insertion and follow stable editing handles through
/// deletion, splitting and merging. The coordinate-free [`TopologyEditor`] is
/// readable through `Deref` (`editor.atoms()`, `editor.neighbors(id)`, ...);
/// structural changes go through the coordinated methods here so geometry stays
/// in step. Realization annotations are addressed by editing handle; explicitly
/// named `topology_*` methods edit static annotations. Generic coordinate
/// changes preserve stored annotations without asserting that arbitrary derived
/// values remain valid.
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
/// let model = editor.finish()?;
/// assert_eq!(model.topology().instance_count(), 1);
/// assert_eq!(model.atom_count(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct ModelEditor {
    topology: TopologyEditor,
    /// Realization state in private editor slot order.
    slots: Conformation,
}

impl Default for ModelEditor {
    fn default() -> Self {
        Self {
            topology: TopologyEditor::default(),
            slots: Conformation::new(Positions::default()),
        }
    }
}

impl std::ops::Deref for ModelEditor {
    type Target = TopologyEditor;

    fn deref(&self) -> &Self::Target {
        &self.topology
    }
}

impl Model {
    pub fn edit(&self) -> ModelEditor {
        self.clone().into_editor()
    }
    /// Moves the conformation into a draft while retaining the shared topology.
    /// Source slots are source dense indices.
    pub fn into_editor(self) -> ModelEditor {
        ModelEditor {
            topology: TopologyEditor::from_topology(self.topology),
            slots: self.conformation,
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
        self.slots = Conformation::new(Positions::default());
    }
    /// The coordinate-free structural draft, also reachable through `Deref`.
    pub fn topology(&self) -> &TopologyEditor {
        &self.topology
    }
    /// Replaces represented stereo for an intact source occurrence; see
    /// [`TopologyEditor::replace_source_instance_stereo`] for identity semantics.
    pub fn replace_source_instance_stereo(
        &mut self,
        instance: MoleculeInstanceId,
        elements: &[crate::core::StereoElement],
    ) -> Result<(), ModelEditError> {
        self.structural(|topology| topology.replace_source_instance_stereo(instance, elements))
    }

    pub fn add_atom(
        &mut self,
        atom: Atom,
        position: Quantity<Point3>,
    ) -> Result<EditAtomId, ModelEditError> {
        let point = checked_point(position)?;
        self.slots.try_reserve_atoms(1)?;
        let id = self.structural(|topology| topology.add_atom(atom))?;
        let slot = self.topology.atom_slot(id)?;
        self.slots
            .set_position(atom_row(slot), Quantity::new(point, CANONICAL_LENGTH_UNIT))?;
        Ok(id)
    }
    pub fn add_molecule(
        &mut self,
        molecule: Molecule,
        positions: &Positions,
    ) -> Result<EditMolecule, ModelEditError> {
        if positions.len() != molecule.atom_count() {
            return Err(PositionError::PositionCountMismatch {
                expected: molecule.atom_count(),
                actual: positions.len(),
            }
            .into());
        }
        self.slots.try_reserve_atoms(positions.len())?;
        let first = self.topology.atom_slot_count();
        let added = self.structural(|topology| topology.add_molecule(molecule))?;
        for (offset, point) in positions.values().value().iter().enumerate() {
            self.slots.set_position(
                atom_row(first + offset),
                Quantity::new(*point, CANONICAL_LENGTH_UNIT),
            )?;
        }
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

    /// Copies positions in live atom-handle order, which is the dense atom order
    /// of the model this draft publishes.
    pub fn positions(&self) -> Positions {
        self.slots
            .positions()
            .select_indices(&self.atom_slots())
            .expect("live atoms have position slots")
    }
    pub fn position(&self, id: EditAtomId) -> Result<Quantity<Point3>, ModelEditError> {
        Ok(self
            .slots
            .positions()
            .position_at(self.topology.atom_slot(id)?)?)
    }
    pub fn set_position(
        &mut self,
        id: EditAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.set_position(atom_row(slot), position)?)
    }
    /// Replaces all live positions transactionally in [`TopologyEditor::atom_ids`] order.
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
        let staged = self
            .atom_slots()
            .into_iter()
            .zip(positions.values().value().iter().copied())
            .collect::<Vec<_>>();
        Ok(self.slots.positions_mut().set_canonical_batch(&staged)?)
    }
    /// Checked sparse coordinate batch; repeated handles are applied in input order.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (EditAtomId, Quantity<Point3>)>,
    ) -> Result<(), ModelEditError> {
        let staged = values
            .into_iter()
            .map(|(id, value)| Ok((self.topology.atom_slot(id)?, checked_point(value)?)))
            .collect::<Result<Vec<_>, ModelEditError>>()?;
        Ok(self.slots.positions_mut().set_canonical_batch(&staged)?)
    }
    pub fn cell(&self) -> Option<&PeriodicCell> {
        self.slots.cell()
    }
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.slots.set_cell(cell);
    }
    /// Realization-level owner annotations. Structural edits clear them.
    pub fn owner_properties(&self) -> &OwnerProperties {
        self.slots.properties().owner()
    }
    pub fn owner_properties_mut(&mut self) -> &mut OwnerProperties {
        self.slots.properties_storage_mut().owner_mut()
    }
    pub fn atom_property(
        &self,
        id: EditAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.properties().atoms().value(key, atom_row(slot))?)
    }
    pub fn bond_property(
        &self,
        id: EditBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, ModelEditError> {
        let slot = self.topology.bond_slot(id)?;
        Ok(self.slots.properties().bonds().value(key, bond_row(slot))?)
    }
    pub fn set_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self
            .slots
            .properties_mut()
            .atoms_mut()
            .set_value(key, atom_row(slot), value)?)
    }
    pub fn set_bond_property(
        &mut self,
        id: EditBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.bond_slot(id)?;
        Ok(self
            .slots
            .properties_mut()
            .bonds_mut()
            .set_value(key, bond_row(slot), value)?)
    }
    /// Applies one atom-column batch transactionally; repeated handles use the last value.
    pub fn set_atom_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditAtomId, Option<PropertyValue>)>,
    ) -> Result<(), ModelEditError> {
        let values = values
            .into_iter()
            .map(|(id, value)| Ok((atom_row(self.topology.atom_slot(id)?), value)))
            .collect::<Result<Vec<_>, ModelEditError>>()?;
        Ok(self
            .slots
            .properties_mut()
            .atoms_mut()
            .set_values(key, values)?)
    }
    /// Applies one bond-column batch transactionally; repeated handles use the last value.
    pub fn set_bond_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditBondId, Option<PropertyValue>)>,
    ) -> Result<(), ModelEditError> {
        let values = values
            .into_iter()
            .map(|(id, value)| Ok((bond_row(self.topology.bond_slot(id)?), value)))
            .collect::<Result<Vec<_>, ModelEditError>>()?;
        Ok(self
            .slots
            .properties_mut()
            .bonds_mut()
            .set_values(key, values)?)
    }
    /// One realization atom column in live atom-handle order.
    pub fn atom_property_column(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        live_column(
            self.slots.properties().atoms().raw(),
            &self.atom_slots(),
            key,
        )
    }
    /// One realization bond column in live bond-handle order.
    pub fn bond_property_column(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        live_column(
            self.slots.properties().bonds().raw(),
            &self.bond_slots(),
            key,
        )
    }
    /// Inserts a realization atom column given in live atom-handle order.
    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        let previous = self.atom_property_column(&key);
        let column =
            column.into_editor_slots(&self.atom_slots(), self.topology.atom_slot_count())?;
        self.slots
            .properties_mut()
            .atoms_mut()
            .insert(key, column)?;
        Ok(previous)
    }
    /// Inserts a realization bond column given in live bond-handle order.
    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, ModelEditError> {
        let previous = self.bond_property_column(&key);
        let column =
            column.into_editor_slots(&self.bond_slots(), self.topology.bond_slot_count())?;
        self.slots
            .properties_mut()
            .bonds_mut()
            .insert(key, column)?;
        Ok(previous)
    }
    pub fn remove_atom_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.atom_property_column(key);
        self.slots.properties_mut().atoms_mut().remove(key);
        previous
    }
    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.bond_property_column(key);
        self.slots.properties_mut().bonds_mut().remove(key);
        previous
    }
    pub fn occupancy(&self, id: EditAtomId) -> Result<Option<f64>, ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.occupancy(atom_row(slot))?)
    }
    pub fn set_occupancy(
        &mut self,
        id: EditAtomId,
        value: Option<f64>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.set_occupancy(atom_row(slot), value)?)
    }
    pub fn b_factor(&self, id: EditAtomId) -> Result<Option<Quantity<f64>>, ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.b_factor(atom_row(slot))?)
    }
    pub fn set_b_factor(
        &mut self,
        id: EditAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ModelEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.set_b_factor(atom_row(slot), value)?)
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
    /// Publishes the completed model, preserving coordinates and entity properties.
    /// Editing handles are draft-only; inspect the returned model for its final IDs.
    pub fn finish(self) -> Result<Model, ModelEditError> {
        self.finish_with_correspondence().map(|(model, _)| model)
    }
    /// Publishes the model and maps surviving draft atom/bond handles to its layout.
    pub fn finish_with_correspondence(
        self,
    ) -> Result<(Model, crate::topology::EditCorrespondence), ModelEditError> {
        let (published, correspondence) = self.topology.into_mapped_publication()?;
        let mut conformation = self
            .slots
            .project(&published.atom_slots, &published.bond_slots)?;
        *conformation.properties_storage_mut().owner_mut() =
            self.slots.properties().owner().clone();
        Ok((
            Model::new(published.topology, conformation)?,
            correspondence,
        ))
    }
    /// Publishes a model, returning the draft with any publication error.
    pub fn try_finish(self) -> Result<Model, ModelFinishError> {
        let snapshot = self.clone();
        self.finish().map_err(|error| ModelFinishError {
            error: Box::new(error),
            editor: Box::new(snapshot),
        })
    }

    /// Runs one topology edit and keeps slot rows in step: new slots get rows,
    /// dead slots lose their annotations, and a structural change clears
    /// realization owner properties.
    fn structural<T>(
        &mut self,
        edit: impl FnOnce(&mut TopologyEditor) -> Result<T, TopologyEditError>,
    ) -> Result<T, ModelEditError> {
        let revision = self.topology.structural_revision;
        let atom_count = self.atom_count();
        let bond_count = self.bond_count();
        let result = edit(&mut self.topology)?;
        self.slots.resize_slots(
            self.topology.atom_slot_count(),
            self.topology.bond_slot_count(),
        );
        if self.topology.structural_revision != revision {
            self.slots.properties_storage_mut().owner_mut().clear();
        }
        if self.atom_count() < atom_count {
            let live = self.atom_slots();
            self.slots.clear_dead_atom_rows(&live);
        }
        if self.bond_count() < bond_count {
            let live = self.bond_slots();
            self.slots.clear_dead_bond_rows(&live);
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

fn atom_row(slot: usize) -> TopologyAtomIndex {
    TopologyAtomIndex::new(u32::try_from(slot).expect("editor slots fit topology indices"))
}

fn bond_row(slot: usize) -> TopologyBondIndex {
    TopologyBondIndex::new(u32::try_from(slot).expect("editor slots fit topology indices"))
}

fn live_column(
    table: &RawPropertyTable,
    live: &[usize],
    key: &PropertyKey,
) -> Option<PropertyColumn> {
    table.get(key)?;
    table
        .select_indices(live)
        .expect("live slots are allocated rows")
        .remove(key)
}

fn checked_point(position: Quantity<Point3>) -> Result<Point3, PositionError> {
    let point = position.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
    if !point.is_finite() {
        return Err(PositionError::NonFinitePosition { index: 0 });
    }
    Ok(point)
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ModelEditError {
    Topology(TopologyEditError),
    Position(PositionError),
    Property(PropertyError),
    Conformation(ConformationError),
    Model(Box<ModelError>),
    CapacityOverflow,
    /// A periodic source cannot be imported under the destination's cell.
    IncompatibleAppendCell,
    /// An imported realization property cannot be combined with its destination table.
    AppendProperty {
        domain: &'static str,
        error: Box<PropertyError>,
    },
}
impl fmt::Display for ModelEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Topology(e) => e.fmt(f),
            Self::Position(e) => e.fmt(f),
            Self::Property(e) => e.fmt(f),
            Self::Conformation(e) => e.fmt(f),
            Self::Model(e) => e.fmt(f),
            Self::CapacityOverflow => f.write_str("model editing exceeds coordinate capacity"),
            Self::IncompatibleAppendCell => f.write_str(
                "appended model has a different periodic cell; set the intended destination cell or explicitly remove the source cell before appending",
            ),
            Self::AppendProperty { domain, error } => write!(f, "cannot append {domain} properties: {error}"),
        }
    }
}
impl std::error::Error for ModelEditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(e) => Some(e),
            Self::Position(e) => Some(e),
            Self::Property(e) => Some(e),
            Self::Conformation(e) => Some(e),
            Self::Model(e) => Some(e.as_ref()),
            Self::AppendProperty { error, .. } => Some(error.as_ref()),
            Self::CapacityOverflow | Self::IncompatibleAppendCell => None,
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
convert!(ConformationError, Conformation);
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

#[cfg(test)]
mod allocation_tests {
    use super::*;

    fn key(name: &str) -> PropertyKey {
        PropertyKey::new(name).unwrap()
    }

    fn int_column(table: &RawPropertyTable, key: &PropertyKey) -> *const Option<i64> {
        match table.get(key).unwrap() {
            PropertyColumn::Int(values) => values.as_ptr(),
            _ => panic!("expected an integer column"),
        }
    }

    // Realization batches stage only their own column, as topology batches do.
    #[test]
    fn realization_batches_keep_unrelated_column_allocations() {
        let topology = crate::smiles::to_topology("CCC").unwrap();
        let atoms = topology.atom_count();
        let mut editor = Model::new(topology, Positions::zeros(atoms))
            .unwrap()
            .into_editor();
        let atoms = editor.atom_ids().collect::<Vec<_>>();
        let bonds = editor.bond_ids().collect::<Vec<_>>();
        let untouched = key("untouched");
        editor
            .set_atom_property(atoms[0], untouched.clone(), Some(PropertyValue::Int(1)))
            .unwrap();
        editor
            .set_bond_property(bonds[0], untouched.clone(), Some(PropertyValue::Int(2)))
            .unwrap();
        let atom_column = int_column(editor.slots.properties().atoms().raw(), &untouched);
        let bond_column = int_column(editor.slots.properties().bonds().raw(), &untouched);
        editor
            .set_atom_properties(
                key("edited"),
                [
                    (atoms[0], Some(PropertyValue::Int(3))),
                    (atoms[1], Some(PropertyValue::Int(4))),
                ],
            )
            .unwrap();
        assert!(editor
            .set_atom_properties(
                key("edited"),
                [(atoms[1], Some(PropertyValue::String("bad type".into())))],
            )
            .is_err());
        editor
            .set_bond_properties(key("edited"), [(bonds[1], Some(PropertyValue::Int(5)))])
            .unwrap();
        editor
            .set_bond_properties(key("edited"), [(bonds[1], None)])
            .unwrap();
        assert_eq!(
            int_column(editor.slots.properties().atoms().raw(), &untouched),
            atom_column
        );
        assert_eq!(
            int_column(editor.slots.properties().bonds().raw(), &untouched),
            bond_column
        );
    }
}
