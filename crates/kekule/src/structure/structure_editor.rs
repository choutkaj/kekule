//! One-structure coordination over the coordinate-free structural editor.
mod append;
mod placement;
pub use append::{StructureAppend, StructureAppendReport};

use super::draft_geometry::DraftGeometry;
use super::{DenseArrayError, GeometryError, Positions, Structure, StructureError};
use crate::math::{PeriodicCell, Point3};
use crate::molecule::{Atom, BondOrder, Molecule};
use crate::properties::{OwnerProperties, PropertyError};
use crate::topology::{
    edit::{
        EditAtomId, EditAtomSite, EditAtomSiteId, EditBond, EditBondId, EditChain, EditChainId,
        EditMolecule, EditPropertiesMut, EditPropertyTable, EditPropertyTableMut, EditResidue,
        EditResidueId, TopologyEditError, TopologyEditor,
    },
    AtomSiteMetadata, MoleculeClass, MoleculeInstanceId, ResidueClass, TopologyAtomIndex,
};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};
use crate::DraftError;
use std::fmt;

/// Detached structural editing for one geometry-bearing molecular system.
///
/// Positions accompany atom insertion and follow stable editing handles through
/// deletion, splitting and merging. The coordinate-free [`TopologyEditor`] is
/// readable through [`Self::topology`] (`editor.topology().atoms()`, ...);
/// structural changes go through the coordinated methods here so geometry stays
/// in step, and [`Self::topology_mut`] reaches the topology edits that need no
/// coordination: static annotations and class assignments. Structure
/// annotations are addressed by editing handle through [`Self::properties_mut`].
/// Generic coordinate changes preserve stored annotations without asserting
/// that arbitrary derived values remain valid.
///
/// ```
/// use kekule::{molecule::{Atom, Element, BondOrder}, math::Point3,
///     structure::StructureEditor, units::{Quantity, ANGSTROM}};
/// let mut editor = StructureEditor::new();
/// let c = editor.add_atom(Atom::new(Element::from_symbol("C").unwrap()),
///     Quantity::new(Point3::origin(), ANGSTROM))?;
/// let o = editor.add_atom(Atom::new(Element::from_symbol("O").unwrap()),
///     Quantity::new(Point3::new(1.4, 0.0, 0.0), ANGSTROM))?;
/// editor.add_bond(c, o, BondOrder::Single)?;
/// let structure = editor.finish()?;
/// assert_eq!(structure.topology().instance_count(), 1);
/// assert_eq!(structure.atom_count(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct StructureEditor {
    topology: TopologyEditor,
    /// Geometry state in private editor slot order.
    slots: DraftGeometry,
}

impl Structure {
    /// Starts a detached draft from a copy of this structure.
    pub fn to_editor(&self) -> StructureEditor {
        self.clone().into_editor()
    }
    /// Moves the geometry state into a draft while retaining the shared
    /// topology. Source slots are source dense indices.
    pub fn into_editor(self) -> StructureEditor {
        let (topology, slots) = self.into_state();
        StructureEditor {
            topology: TopologyEditor::from_topology(topology),
            slots,
        }
    }
}

impl StructureEditor {
    /// An empty draft with no atoms, positions, cell or annotations.
    pub fn new() -> Self {
        Self::default()
    }
    /// Clears structure, geometry, hierarchy, cell and annotations in the draft.
    pub fn clear(&mut self) {
        self.topology.clear();
        self.slots = DraftGeometry::default();
    }
    /// The coordinate-free structural draft.
    pub fn topology(&self) -> &TopologyEditor {
        &self.topology
    }
    /// The topology edits that need no coordination with geometry: static
    /// annotations and class assignments.
    pub fn topology_mut(&mut self) -> StructureEditorTopologyMut<'_> {
        StructureEditorTopologyMut {
            topology: &mut self.topology,
        }
    }
    /// Structure annotations, addressed by editing handle.
    pub fn properties(&self) -> StructureEditorProperties<'_> {
        StructureEditorProperties { editor: self }
    }
    /// Mutable structure annotations, addressed by editing handle. Structural
    /// edits clear the owner annotations.
    pub fn properties_mut(&mut self) -> StructureEditorPropertiesMut<'_> {
        StructureEditorPropertiesMut { editor: self }
    }
    /// Replaces represented stereo for an intact source occurrence; see
    /// [`TopologyEditor::replace_source_instance_stereo`] for identity semantics.
    pub fn replace_source_instance_stereo(
        &mut self,
        instance: MoleculeInstanceId,
        elements: &[crate::stereo::StereoElement],
    ) -> Result<(), StructureEditError> {
        self.structural(|topology| topology.replace_source_instance_stereo(instance, elements))
    }

    /// Adds an isolated atom at `position`, accepting any length unit.
    ///
    /// Bond it to an existing atom or publish it as a new molecule instance; see
    /// [`TopologyEditor::add_atom`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Position`] when the unit is not a length,
    /// [`StructureEditError::NonFinitePosition`] when the converted position is not
    /// finite, and [`StructureEditError::Geometry`] when capacity is exhausted. The
    /// draft is unchanged on error.
    pub fn add_atom(
        &mut self,
        atom: Atom,
        position: Quantity<Point3>,
    ) -> Result<EditAtomId, StructureEditError> {
        let point = checked_point(position, None)?;
        self.slots.try_reserve_atoms(1)?;
        let id = self.structural(|topology| topology.add_atom(atom))?;
        let slot = self.topology.atom_slot(id)?;
        self.slots
            .set_position(atom_row(slot), Quantity::new(point, CANONICAL_LENGTH_UNIT))?;
        Ok(id)
    }
    /// Adds a complete molecule occurrence with positions in the molecule's atom
    /// order; see [`TopologyEditor::add_molecule`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Position`] when the position count differs from
    /// the molecule's atom count, [`StructureEditError::Geometry`] when capacity is
    /// exhausted, and [`StructureEditError::Topology`] when the molecule is rejected.
    pub fn add_molecule(
        &mut self,
        molecule: Molecule,
        positions: &Positions,
    ) -> Result<EditMolecule, StructureEditError> {
        if positions.len() != molecule.atom_count() {
            return Err(DenseArrayError::CountMismatch {
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
    /// Replaces the chemistry of one atom, keeping its position and annotations,
    /// and returns the previous atom.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live atom.
    pub fn replace_atom(&mut self, id: EditAtomId, atom: Atom) -> Result<Atom, StructureEditError> {
        self.structural(|t| t.replace_atom(id, atom))
    }
    /// Deletes one atom, its incident bonds, hierarchy site, position and
    /// annotations, and returns the deleted atom.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live atom.
    pub fn delete_atom(&mut self, id: EditAtomId) -> Result<Atom, StructureEditError> {
        self.structural(|t| t.delete_atom(id))
    }
    /// Deletes several atoms as [`Self::delete_atom`] does, returning them in
    /// handle order. Every handle is validated before any atom is deleted;
    /// duplicates are ignored.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when a handle names no live atom.
    pub fn delete_atoms(
        &mut self,
        ids: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<Vec<(EditAtomId, Atom)>, StructureEditError> {
        self.structural(|t| t.delete_atoms(ids))
    }
    /// Deletes every live atom not in `ids`, returning the deleted atoms in
    /// handle order.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when a handle names no live atom;
    /// nothing is deleted then.
    pub fn retain_atoms(
        &mut self,
        ids: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<Vec<(EditAtomId, Atom)>, StructureEditError> {
        self.structural(|t| t.retain_atoms(ids))
    }
    /// Deletes surviving atoms of a source occurrence and their geometry.
    /// See [`TopologyEditor::delete_instance`] for semantics after splits/merges.
    pub fn delete_instance(
        &mut self,
        id: MoleculeInstanceId,
    ) -> Result<Vec<(EditAtomId, Atom)>, StructureEditError> {
        self.structural(|t| t.delete_instance(id))
    }
    /// Adds an asserted bond between two live atoms; atoms of different
    /// occurrences are merged. See [`TopologyEditor::add_bond`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or a bond the
    /// molecule rejects.
    pub fn add_bond(
        &mut self,
        a: EditAtomId,
        b: EditAtomId,
        order: BondOrder,
    ) -> Result<EditBondId, StructureEditError> {
        self.structural(|t| t.add_bond(a, b, order))
    }
    /// Deletes one bond and its annotations, returning it. Positions are kept.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live bond.
    pub fn delete_bond(&mut self, id: EditBondId) -> Result<EditBond, StructureEditError> {
        self.structural(|t| t.delete_bond(id))
    }
    /// Deletes several bonds, returning them in handle order. Every handle is
    /// validated first; duplicates are ignored.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when a handle names no live bond.
    pub fn delete_bonds(
        &mut self,
        ids: impl IntoIterator<Item = EditBondId>,
    ) -> Result<Vec<(EditBondId, EditBond)>, StructureEditError> {
        self.structural(|t| t.delete_bonds(ids))
    }
    /// Changes the represented order of one bond; see
    /// [`TopologyEditor::set_bond_order`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or an order the
    /// molecule rejects.
    pub fn set_bond_order(
        &mut self,
        id: EditBondId,
        order: BondOrder,
    ) -> Result<(), StructureEditError> {
        self.structural(|t| t.set_bond_order(id, order))
    }
    /// Rewires one bond, keeping its handle and annotations; see
    /// [`TopologyEditor::set_bond_endpoints`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or a rejected
    /// bond.
    pub fn set_bond_endpoints(
        &mut self,
        id: EditBondId,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<(), StructureEditError> {
        self.structural(|t| t.set_bond_endpoints(id, a, b))
    }
    /// Replaces the endpoints and order of one bond transactionally, returning
    /// the previous bond.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or a rejected
    /// bond; the draft is unchanged then.
    pub fn replace_bond(
        &mut self,
        id: EditBondId,
        replacement: EditBond,
    ) -> Result<EditBond, StructureEditError> {
        self.structural(|t| t.replace_bond(id, replacement))
    }
    /// Copies positions in live atom-handle order, which is the dense atom order
    /// of the structure this draft publishes.
    pub fn positions(&self) -> Positions {
        self.slots
            .positions
            .select_indices(&self.atom_slots())
            .expect("live atoms have position slots")
    }
    /// The draft position of one atom, or `None` when `id` names no live atom.
    pub fn position(&self, id: EditAtomId) -> Option<Quantity<Point3>> {
        self.slots.positions.get(self.topology.atom_slot(id).ok()?)
    }
    /// Replaces the draft position of one atom, accepting any length unit.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live atom and
    /// [`StructureEditError::Geometry`] when the unit is not a length or the converted
    /// position is not finite.
    pub fn set_position(
        &mut self,
        id: EditAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), StructureEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.set_position(atom_row(slot), position)?)
    }
    /// Replaces all live positions transactionally in [`TopologyEditor::atom_ids`] order.
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), StructureEditError> {
        let positions = Positions::new(positions)?;
        if positions.len() != self.topology.atom_count() {
            return Err(DenseArrayError::CountMismatch {
                expected: self.topology.atom_count(),
                actual: positions.len(),
            }
            .into());
        }
        let staged = self
            .atom_slots()
            .into_iter()
            .zip(positions.values().value().iter().copied())
            .collect::<Vec<_>>();
        Ok(self.slots.positions.set_canonical_batch(&staged)?)
    }
    /// Checked sparse coordinate batch; repeated handles are applied in input order.
    pub fn set_atom_positions(
        &mut self,
        values: impl IntoIterator<Item = (EditAtomId, Quantity<Point3>)>,
    ) -> Result<(), StructureEditError> {
        let staged = values
            .into_iter()
            .map(|(id, value)| {
                Ok((
                    self.topology.atom_slot(id)?,
                    checked_point(value, Some(id))?,
                ))
            })
            .collect::<Result<Vec<_>, StructureEditError>>()?;
        Ok(self.slots.positions.set_canonical_batch(&staged)?)
    }
    /// The draft periodic cell, or `None` for a non-periodic structure.
    pub fn cell(&self) -> Option<&PeriodicCell> {
        self.slots.cell.as_ref()
    }
    /// Replaces the draft periodic cell; positions are not wrapped or rescaled.
    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.slots.cell = cell;
    }
    /// The draft occupancy of one atom; `None` when it is unset or `id` names
    /// no live atom.
    pub fn occupancy(&self, id: EditAtomId) -> Option<f64> {
        self.slots
            .occupancy(atom_row(self.topology.atom_slot(id).ok()?))
    }
    /// Sets or, with `None`, clears the draft occupancy of one atom.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live atom and
    /// [`StructureEditError::Geometry`] when the occupancy is not finite.
    pub fn set_occupancy(
        &mut self,
        id: EditAtomId,
        value: Option<f64>,
    ) -> Result<(), StructureEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.set_occupancy(atom_row(slot), value)?)
    }
    /// The draft B-factor of one atom; `None` when it is unset or `id` names
    /// no live atom.
    pub fn b_factor(&self, id: EditAtomId) -> Option<Quantity<f64>> {
        self.slots
            .b_factor(atom_row(self.topology.atom_slot(id).ok()?))
    }
    /// Sets or, with `None`, clears the draft isotropic B-factor of one atom,
    /// accepting any area unit.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live atom and
    /// [`StructureEditError::Geometry`] when the unit is not an area or the value is
    /// not finite.
    pub fn set_b_factor(
        &mut self,
        id: EditAtomId,
        value: Option<Quantity<f64>>,
    ) -> Result<(), StructureEditError> {
        let slot = self.topology.atom_slot(id)?;
        Ok(self.slots.set_b_factor(atom_row(slot), value)?)
    }

    /// Adds a hierarchy chain with a label and optional author identifier; see
    /// [`TopologyEditor::add_chain`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when the hierarchy rejects the chain.
    pub fn add_chain(
        &mut self,
        label: impl Into<String>,
        author: Option<String>,
    ) -> Result<EditChainId, StructureEditError> {
        self.structural(|t| t.add_chain(label, author))
    }
    /// Adds a residue to a chain; see [`TopologyEditor::add_residue`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `chain` names no live chain.
    pub fn add_residue(
        &mut self,
        chain: EditChainId,
        name: impl Into<String>,
        label_seq: Option<i32>,
        author_seq: Option<String>,
        insertion: Option<String>,
    ) -> Result<EditResidueId, StructureEditError> {
        self.structural(|t| t.add_residue(chain, name, label_seq, author_seq, insertion))
    }
    /// Places one atom in a residue as an atom site; see
    /// [`TopologyEditor::add_atom_site`].
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle, an atom that
    /// already has a site, or a site the hierarchy rejects.
    pub fn add_atom_site(
        &mut self,
        residue: EditResidueId,
        atom: EditAtomId,
        metadata: AtomSiteMetadata,
    ) -> Result<EditAtomSiteId, StructureEditError> {
        self.structural(|t| t.add_atom_site(residue, atom, metadata))
    }
    /// See [`TopologyEditor::set_residue_component_ids`].
    pub fn set_residue_component_ids(
        &mut self,
        id: EditResidueId,
        label: Option<String>,
        author: Option<String>,
    ) -> Result<(), StructureEditError> {
        self.structural(|t| t.set_residue_component_ids(id, label, author))
    }
    /// Replaces the label and author identifiers of one chain.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or rejected
    /// identifiers.
    pub fn set_chain_identifiers(
        &mut self,
        id: EditChainId,
        label: impl Into<String>,
        author: Option<String>,
    ) -> Result<(), StructureEditError> {
        self.structural(|t| t.set_chain_identifiers(id, label, author))
    }
    /// Replaces the metadata of one atom site.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or rejected
    /// metadata.
    pub fn set_atom_site_metadata(
        &mut self,
        id: EditAtomSiteId,
        metadata: AtomSiteMetadata,
    ) -> Result<(), StructureEditError> {
        self.structural(|t| t.set_atom_site_metadata(id, metadata))
    }
    /// Moves one atom site to another residue.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] for an invalid handle or a move the
    /// hierarchy rejects.
    pub fn set_atom_site_residue(
        &mut self,
        id: EditAtomSiteId,
        residue: EditResidueId,
    ) -> Result<(), StructureEditError> {
        self.structural(|t| t.set_atom_site_residue(id, residue))
    }
    /// Deletes one hierarchy chain, returning it; see
    /// [`TopologyEditor::delete_chain`]. Atoms and positions are kept.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live chain.
    pub fn delete_chain(&mut self, id: EditChainId) -> Result<EditChain, StructureEditError> {
        self.structural(|t| t.delete_chain(id))
    }
    /// Deletes one hierarchy residue, returning it; see
    /// [`TopologyEditor::delete_residue`]. Atoms and positions are kept.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live residue.
    pub fn delete_residue(&mut self, id: EditResidueId) -> Result<EditResidue, StructureEditError> {
        self.structural(|t| t.delete_residue(id))
    }
    /// Removes one atom site from the hierarchy, returning it. The atom and its
    /// position are kept.
    ///
    /// # Errors
    ///
    /// Returns [`StructureEditError::Topology`] when `id` names no live atom site.
    pub fn delete_atom_site(
        &mut self,
        id: EditAtomSiteId,
    ) -> Result<EditAtomSite, StructureEditError> {
        self.structural(|t| t.delete_atom_site(id))
    }
    /// Checks that [`Self::finish`] would succeed, without consuming the draft.
    ///
    /// # Errors
    ///
    /// The error [`Self::finish`] would return.
    pub fn validate(&self) -> Result<(), StructureEditError> {
        self.clone().finish().map(|_| ())
    }
    /// Publishes the completed structure, preserving coordinates and entity properties.
    /// Editing handles are draft-only; inspect the returned structure for its final IDs.
    pub fn finish(self) -> Result<Structure, StructureEditError> {
        self.finish_with_correspondence()
            .map(|(structure, _)| structure)
    }
    /// Publishes the structure and maps surviving draft atom/bond handles to its layout.
    pub fn finish_with_correspondence(
        self,
    ) -> Result<(Structure, crate::topology::edit::EditCorrespondence), StructureEditError> {
        let (published, correspondence) = self.topology.into_mapped_publication()?;
        let mut state = self
            .slots
            .project(&published.atom_slots, &published.bond_slots)?;
        *state.properties.owner_mut() = self.slots.properties.owner().clone();
        Ok((
            Structure::from_state(published.topology, state)?,
            correspondence,
        ))
    }
    /// Publishes a structure, returning the draft with any publication error.
    pub fn try_finish(self) -> Result<Structure, DraftError<Self, StructureEditError>> {
        let snapshot = self.clone();
        self.finish()
            .map_err(|error| DraftError::new(error, snapshot))
    }

    /// Runs one topology edit and keeps slot rows in step: new slots get rows,
    /// dead slots lose their annotations, and a structural change clears
    /// structure owner properties.
    fn structural<T>(
        &mut self,
        edit: impl FnOnce(&mut TopologyEditor) -> Result<T, TopologyEditError>,
    ) -> Result<T, StructureEditError> {
        let revision = self.topology.structural_revision;
        let mark = self.topology.deletion_mark();
        let result = edit(&mut self.topology)?;
        self.slots.resize_slots(
            self.topology.atom_slot_count(),
            self.topology.bond_slot_count(),
        );
        if self.topology.structural_revision != revision {
            self.slots.properties.owner_mut().clear();
        }
        // Only the rows this edit freed are cleared; earlier ones already were.
        let (atoms, bonds) = self.topology.deleted_slots_since(mark);
        self.slots.clear_atom_rows(atoms);
        self.slots.clear_bond_rows(bonds);
        Ok(result)
    }
    fn atom_slots(&self) -> Vec<usize> {
        self.topology.live_atom_slots()
    }
}

fn atom_row(slot: usize) -> TopologyAtomIndex {
    TopologyAtomIndex::new(u32::try_from(slot).expect("editor slots fit topology indices"))
}

/// The topology edits of a [`StructureEditor`] that need no coordination with
/// geometry: static annotations and class assignments.
#[derive(Debug)]
pub struct StructureEditorTopologyMut<'a> {
    topology: &'a mut TopologyEditor,
}

impl StructureEditorTopologyMut<'_> {
    /// Static topology annotations; see [`TopologyEditor::properties_mut`].
    pub fn properties_mut(&mut self) -> EditPropertiesMut<'_> {
        self.topology.properties_mut()
    }
    /// See [`TopologyEditor::set_molecule_class`].
    pub fn set_molecule_class(
        &mut self,
        atom: EditAtomId,
        class: MoleculeClass,
    ) -> Result<(), TopologyEditError> {
        self.topology.set_molecule_class(atom, class)
    }
    /// See [`TopologyEditor::set_residue_class`].
    pub fn set_residue_class(
        &mut self,
        residue: EditResidueId,
        class: ResidueClass,
    ) -> Result<(), TopologyEditError> {
        self.topology.set_residue_class(residue, class)
    }
    /// See [`TopologyEditor::clear_residue_class`].
    pub fn clear_residue_class(&mut self, residue: EditResidueId) -> Result<(), TopologyEditError> {
        self.topology.clear_residue_class(residue)
    }
}

/// Structure annotations of a [`StructureEditor`], addressed by editing handle.
#[derive(Debug, Clone, Copy)]
pub struct StructureEditorProperties<'a> {
    editor: &'a StructureEditor,
}

impl<'a> StructureEditorProperties<'a> {
    /// Structure-level owner annotations.
    pub fn owner(&self) -> &'a OwnerProperties {
        self.editor.slots.properties.owner()
    }
    /// Atom annotations, addressed by editing handle.
    pub fn atoms(&self) -> EditPropertyTable<'a, EditAtomId> {
        EditPropertyTable::new(
            self.editor.slots.properties.atoms().raw(),
            self.editor.topology.atom_rows(),
        )
    }
    /// Bond annotations, addressed by editing handle.
    pub fn bonds(&self) -> EditPropertyTable<'a, EditBondId> {
        EditPropertyTable::new(
            self.editor.slots.properties.bonds().raw(),
            self.editor.topology.bond_rows(),
        )
    }
}

/// Mutable structure annotations of a [`StructureEditor`], addressed by editing
/// handle; see [`EditPropertyTableMut`] for the table methods.
#[derive(Debug)]
pub struct StructureEditorPropertiesMut<'a> {
    editor: &'a mut StructureEditor,
}

impl StructureEditorPropertiesMut<'_> {
    /// Structure-level owner annotations.
    pub fn owner(&self) -> &OwnerProperties {
        self.editor.slots.properties.owner()
    }
    /// Structure-level owner annotations. Structural edits clear them.
    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        self.editor.slots.properties.owner_mut()
    }
    /// Atom annotations, addressed by editing handle.
    pub fn atoms_mut(&mut self) -> EditPropertyTableMut<'_, EditAtomId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.slots.properties.atoms_raw_mut(),
            editor.topology.atom_rows(),
            None,
        )
    }
    /// Bond annotations, addressed by editing handle.
    pub fn bonds_mut(&mut self) -> EditPropertyTableMut<'_, EditBondId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.slots.properties.bonds_raw_mut(),
            editor.topology.bond_rows(),
            None,
        )
    }
}

/// A canonical finite point for editor atom `atom`, or for an atom being
/// added when `atom` is `None`.
fn checked_point(
    position: Quantity<Point3>,
    atom: Option<EditAtomId>,
) -> Result<Point3, StructureEditError> {
    let point = position
        .into_unit(CANONICAL_LENGTH_UNIT)
        .map_err(DenseArrayError::from)?
        .into_value();
    if !point.is_finite() {
        return Err(StructureEditError::NonFinitePosition { atom });
    }
    Ok(point)
}

/// A structure edit was rejected.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum StructureEditError {
    /// Editing the topology failed; see the source.
    Topology(TopologyEditError),
    /// Dense per-atom state is invalid; see the source.
    Position(DenseArrayError),
    /// A position is not finite: of editor atom `atom`, or of the atom being
    /// added when `atom` is `None`.
    NonFinitePosition {
        /// The editor atom.
        atom: Option<EditAtomId>,
    },
    /// A property edit failed; see the source.
    Property(PropertyError),
    /// Geometry state is invalid; see the source.
    Geometry(GeometryError),
    /// The edited structure is invalid; see the source.
    Structure(Box<StructureError>),
    /// The edit exceeds coordinate capacity.
    CapacityOverflow,
    /// A periodic source cannot be imported under the destination's cell.
    IncompatibleAppendCell,
    /// An imported geometry property cannot be combined with its destination table.
    AppendProperty {
        /// The table whose rows could not be appended.
        domain: crate::properties::PropertyDomain,
        /// The cause, also the source.
        error: Box<PropertyError>,
    },
}
impl fmt::Display for StructureEditError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Topology(error) => error.fmt(f),
            Self::Position(error) => error.fmt(f),
            Self::NonFinitePosition { atom: Some(atom) } => {
                write!(f, "The position of {atom} is not finite.")
            }
            Self::NonFinitePosition { atom: None } => {
                f.write_str("The position of the added atom is not finite.")
            }
            Self::Property(error) => error.fmt(f),
            Self::Geometry(error) => error.fmt(f),
            Self::Structure(error) => error.fmt(f),
            Self::CapacityOverflow => {
                f.write_str("The edit exceeds the addressable number of atoms.")
            }
            Self::IncompatibleAppendCell => f.write_str(
                "The appended structure has a different periodic cell. Set the destination cell or remove the source cell first.",
            ),
            Self::AppendProperty { domain, .. } => {
                write!(f, "Cannot append the {domain} properties.")
            }
        }
    }
}
impl std::error::Error for StructureEditError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Topology(e) => e.source(),
            Self::Position(e) => e.source(),
            Self::Property(e) => e.source(),
            Self::Geometry(e) => e.source(),
            Self::Structure(e) => e.source(),
            Self::AppendProperty { error, .. } => Some(error.as_ref()),
            Self::CapacityOverflow
            | Self::IncompatibleAppendCell
            | Self::NonFinitePosition { .. } => None,
        }
    }
}
macro_rules! convert {
    ($ty:ty, $variant:ident) => {
        impl From<$ty> for StructureEditError {
            fn from(e: $ty) -> Self {
                Self::$variant(e)
            }
        }
    };
}
convert!(TopologyEditError, Topology);
convert!(DenseArrayError, Position);
convert!(PropertyError, Property);
convert!(GeometryError, Geometry);
impl From<StructureError> for StructureEditError {
    fn from(error: StructureError) -> Self {
        Self::Structure(Box::new(error))
    }
}

#[cfg(test)]
mod allocation_tests {
    use super::*;
    use crate::properties::{PropertyColumn, PropertyKey, PropertyValue, RawPropertyTable};

    fn key(name: &str) -> PropertyKey {
        PropertyKey::new(name).unwrap()
    }

    fn int_column(table: &RawPropertyTable, key: &PropertyKey) -> *const Option<i64> {
        match table.get(key).unwrap() {
            PropertyColumn::Int(values) => values.as_ptr(),
            _ => panic!("expected an integer column"),
        }
    }

    // Structure batches stage only their own column, as topology batches do.
    #[test]
    fn structure_batches_keep_unrelated_column_allocations() {
        let topology = crate::smiles::to_topology("CCC").unwrap();
        let atoms = topology.atom_count();
        let mut editor = Structure::new(topology, Positions::zeros(atoms))
            .unwrap()
            .into_editor();
        let atoms = editor.topology().atom_ids().collect::<Vec<_>>();
        let bonds = editor.topology().bond_ids().collect::<Vec<_>>();
        let untouched = key("untouched");
        let mut properties = editor.properties_mut();
        properties
            .atoms_mut()
            .set_value(untouched.clone(), atoms[0], Some(PropertyValue::Int(1)))
            .unwrap();
        properties
            .bonds_mut()
            .set_value(untouched.clone(), bonds[0], Some(PropertyValue::Int(2)))
            .unwrap();
        let atom_column = int_column(editor.slots.properties.atoms().raw(), &untouched);
        let bond_column = int_column(editor.slots.properties.bonds().raw(), &untouched);
        let mut properties = editor.properties_mut();
        properties
            .atoms_mut()
            .set_values(
                key("edited"),
                [
                    (atoms[0], Some(PropertyValue::Int(3))),
                    (atoms[1], Some(PropertyValue::Int(4))),
                ],
            )
            .unwrap();
        assert!(properties
            .atoms_mut()
            .set_values(
                key("edited"),
                [(atoms[1], Some(PropertyValue::String("bad type".into())))],
            )
            .is_err());
        properties
            .bonds_mut()
            .set_values(key("edited"), [(bonds[1], Some(PropertyValue::Int(5)))])
            .unwrap();
        properties
            .bonds_mut()
            .set_values(key("edited"), [(bonds[1], None)])
            .unwrap();
        assert_eq!(
            int_column(editor.slots.properties.atoms().raw(), &untouched),
            atom_column
        );
        assert_eq!(
            int_column(editor.slots.properties.bonds().raw(), &untouched),
            bond_column
        );
    }
}
