//! [`MoleculeEditor`], the draft that owns every structural change to a
//! molecule.
//!
//! The draft graph may be empty, disconnected, or hold deleted slots.
//! Structural edits invalidate perception and clear owner properties while
//! keeping annotations on surviving atoms and bonds; property edits leave
//! perception intact. Graph edits prune the represented stereo they
//! invalidate and record it ([`MoleculeEditor::pruned_stereo_elements`]).
//! Publication validates the draft and renumbers it densely.
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::{Deref, DerefMut};

use crate::diagnostics::Countable;
use crate::fixed_u32_id::checked_id;
use crate::properties::{LiveRows, MoleculeProperties, MoleculePropertiesMut};

use super::*;

mod publication;
mod stereo;

pub(crate) use publication::canonicalize_represented_chemistry;
pub(crate) use stereo::validate_represented_stereo;

/// Transactional construction and structural editing state.
///
/// The working graph may be empty, disconnected, or chemically incomplete.
/// [`Self::finish`] validates graph references, connectedness, stereochemistry,
/// and canonical represented chemistry before publishing a [`Molecule`]. A
/// failed finish returns structured errors and publishes nothing.
/// Use [`Self::try_finish`] when a failed draft must remain available for repair.
///
/// | Task | Operations |
/// | --- | --- |
/// | Start | [`Self::new`], [`Molecule::to_editor`], [`Molecule::into_editor`] |
/// | Inspect | [`DraftGraph`] reads through `Deref` ([`Graph::atoms`], [`Graph::bonds`], [`Graph::neighbors`]), [`Self::connected_components`] |
/// | Change graph | [`Self::add_atom`], [`Self::add_bond`], [`Self::replace_atom`], [`Self::replace_bond`], [`Self::delete_atoms`], [`Self::retain_atoms`] |
/// | Combine fragments | [`Self::append_molecule`] with returned ID mappings |
/// | Annotate | [`Self::properties_mut`] (rows are draft atom and bond slots) |
/// | Edit stereo | [`Self::add_stereo_element`], [`Self::replace_stereo_element`], [`Self::replace_stereo_group`] |
///
/// Draft atom and bond IDs remain stable until [`Self::clear`]; iteration
/// skips deleted slots. Publication renumbers every ID space densely, so a
/// published [`Molecule`] never contains deleted slots;
/// [`Self::finish_with_correspondence`] maps draft IDs to published IDs.
///
/// An unfinished editor cannot be passed to algorithms requiring a published molecule:
/// ```compile_fail
/// use kekule::molecule::{Molecule, MoleculeEditor};
/// fn published(editor: &MoleculeEditor) -> &Molecule { editor }
/// ```
///
/// # Example
///
/// ```
/// use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
///
/// let mut editor = MoleculeEditor::new();
/// let carbon = editor.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
/// let oxygen = editor.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
/// editor.add_bond(carbon, oxygen, BondOrder::Single)?;
/// let molecule = editor.finish()?;
///
/// assert_eq!(molecule.atom_count(), 2);
/// assert_eq!(molecule.bond_count(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeEditor {
    graph: DraftGraph,
    perception: Perception,
    properties: MoleculeProperties,
    /// Stereo elements that structural edits removed, in removal order.
    pruned_stereo: Vec<StereoElementId>,
}

impl Deref for MoleculeEditor {
    type Target = DraftGraph;

    fn deref(&self) -> &DraftGraph {
        &self.graph
    }
}

impl Default for MoleculeEditor {
    fn default() -> Self {
        Self {
            graph: DraftGraph::default(),
            perception: Perception::default(),
            properties: MoleculeProperties::new(0, 0),
            pruned_stereo: Vec::new(),
        }
    }
}

impl Molecule {
    /// Starts a detached transaction from this published molecule.
    pub fn to_editor(&self) -> MoleculeEditor {
        self.clone().into_editor()
    }

    /// Moves this molecule into an editor whose draft IDs equal its IDs,
    /// without cloning its atoms, bonds, or properties.
    pub fn into_editor(self) -> MoleculeEditor {
        MoleculeEditor {
            graph: self.graph.into_draft(),
            perception: self.perception,
            properties: self.properties,
            pruned_stereo: Vec::new(),
        }
    }
}

/// Mutable access to one live atom of a draft, from
/// [`MoleculeEditor::atom_mut`].
///
/// Read the atom through `Deref`. Change one field with a setter, or replace
/// the whole atom through `DerefMut` (`*guard = atom`). Each setter call and
/// each `DerefMut` borrow clears installed perception and owner properties
/// at once, even when the value is unchanged; borrowing or reading the guard
/// clears nothing. Use [`MoleculeEditor::replace_atom`] to keep cached state
/// for an identical replacement.
pub struct AtomMut<'a> {
    editor: &'a mut MoleculeEditor,
    id: AtomId,
}

impl fmt::Debug for AtomMut<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("AtomMut")
            .field("id", &self.id)
            .field("atom", &**self)
            .finish_non_exhaustive()
    }
}

impl Deref for AtomMut<'_> {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        self.editor
            .graph
            .atom(self.id)
            .expect("validated atom must remain live while borrowed")
    }
}

impl DerefMut for AtomMut<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.editor.chemistry_changed();
        self.editor
            .graph
            .atom_mut(self.id)
            .expect("validated atom must remain live while borrowed")
    }
}

impl AtomMut<'_> {
    /// Changes the element.
    pub fn set_element(&mut self, element: Element) {
        self.element = element;
    }

    /// Changes the mass number, or natural abundance for `None`.
    pub fn set_isotope(&mut self, isotope: Option<NonZeroU16>) {
        self.isotope = isotope;
    }

    /// Changes the formal charge.
    pub fn set_formal_charge(&mut self, formal_charge: i8) {
        self.formal_charge = formal_charge;
    }

    /// Changes the radical electrons, or removes them for `None`.
    pub fn set_radical(&mut self, radical: Option<AtomRadical>) {
        self.radical = radical;
    }

    /// Changes the implicit-hydrogen declaration.
    pub fn set_hydrogens(&mut self, hydrogens: ImplicitHydrogens) {
        self.hydrogens = hydrogens;
    }

    /// Changes the reaction atom map, or removes it for `None`.
    pub fn set_atom_map(&mut self, atom_map: Option<NonZeroU32>) {
        self.atom_map = atom_map;
    }
}

/// Order-only mutable access to one live bond of a draft, from
/// [`MoleculeEditor::bond_mut`].
///
/// Endpoints stay fixed so adjacency remains consistent; read the bond
/// through `Deref` and change its order with [`BondMut::set_order`].
pub struct BondMut<'a> {
    editor: &'a mut MoleculeEditor,
    id: BondId,
}

impl fmt::Debug for BondMut<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("BondMut")
            .field("id", &self.id)
            .field("bond", &**self)
            .finish_non_exhaustive()
    }
}

impl Deref for BondMut<'_> {
    type Target = Bond;

    fn deref(&self) -> &Self::Target {
        self.editor
            .graph
            .bond(self.id)
            .expect("validated bond must remain live while borrowed")
    }
}

impl BondMut<'_> {
    /// Changes represented order without exposing mutable bond endpoints.
    /// Use [`MoleculeEditor::set_bond_endpoints`] to maintain adjacency when rewiring.
    /// A changed order removes stereo assertions focused on this bond, including
    /// their relation-group memberships, and records them as pruned. Assigning
    /// the current order is a no-op.
    pub fn set_order(&mut self, order: BondOrder) {
        if self.order == order {
            return;
        }
        self.editor.chemistry_changed();
        self.editor
            .graph
            .bond_mut(self.id)
            .expect("validated bond must remain live while borrowed")
            .order = order;
        self.editor.prune_stereo_for_bond(self.id);
    }
}

impl MoleculeEditor {
    /// An empty draft with no atoms, perception, or properties.
    pub fn new() -> Self {
        Self::default()
    }

    /// The draft graph, also reachable through `Deref`.
    pub fn graph(&self) -> &DraftGraph {
        &self.graph
    }

    /// Draft annotations; rows are draft atom and bond slots.
    pub fn properties(&self) -> &MoleculeProperties {
        &self.properties
    }

    /// Length-preserving mutable access to draft annotations. Rows are draft
    /// slots, so they follow draft IDs; rows of deleted atoms and bonds accept
    /// only missing values and are dropped at publication.
    pub fn properties_mut(&mut self) -> MoleculePropertiesMut<'_> {
        MoleculePropertiesMut::new(
            &mut self.properties,
            LiveRows::Slots(&self.graph.atoms),
            LiveRows::Slots(&self.graph.bonds),
        )
    }

    /// Whole-table access for interpreters and topology edits that assemble
    /// draft annotations from another owner's rows.
    pub(crate) fn properties_storage_mut(&mut self) -> &mut MoleculeProperties {
        &mut self.properties
    }

    /// The draft's installed perception. Chemical edits clear it.
    pub fn perception(&self) -> &Perception {
        &self.perception
    }

    /// See [`Molecule::explicit_hydrogens`].
    pub fn explicit_hydrogens(&self, atom: AtomId) -> Option<usize> {
        explicit_hydrogen_count(&self.graph, atom)
    }

    /// See [`Molecule::implicit_hydrogens`]. Chemical edits invalidate inference.
    pub fn implicit_hydrogens(&self, atom: AtomId) -> Option<usize> {
        implicit_hydrogen_count(&self.graph, &self.perception, atom)
    }

    /// See [`Molecule::total_hydrogens`].
    pub fn total_hydrogens(&self, atom: AtomId) -> Option<usize> {
        Some(self.explicit_hydrogens(atom)? + self.implicit_hydrogens(atom)?)
    }

    /// Removes all installed perception without changing represented chemistry.
    pub fn clear_perception(&mut self) {
        self.perception = Perception::default();
    }

    fn invalidate_cip(&mut self) {
        self.perception.cip = None;
    }

    /// Records a chemical edit: perception and owner properties no longer
    /// describe the draft.
    fn chemistry_changed(&mut self) {
        self.clear_perception();
        self.properties.owner_mut().clear();
    }

    /// Raw draft storage for interpreters that rewrite bond orders or atoms
    /// in bulk. Obtaining it clears perception and owner properties; writes
    /// skip stereo pruning, which publication's validation backs up.
    pub(crate) fn graph_mut(&mut self) -> &mut DraftGraph {
        self.chemistry_changed();
        &mut self.graph
    }

    /// Whether the draft has no live atoms.
    pub fn is_empty(&self) -> bool {
        self.graph.atom_count() == 0
    }

    /// Live connected components in stable atom order, including isolated atoms.
    pub fn connected_components(&self) -> Vec<Vec<AtomId>> {
        self.graph.connected_components()
    }

    /// Empty editing state is not a connected molecule.
    pub fn is_connected(&self) -> bool {
        self.graph.connected_components().len() == 1
    }

    /// Mutable represented atom state. Mutable access through the guard
    /// invalidates perception and owner properties immediately, even if the
    /// guard is forgotten. Use [`Self::replace_atom`] when an identical
    /// replacement should preserve cached state.
    ///
    /// Returns `None` when `atom` names no live atom of this draft.
    pub fn atom_mut(&mut self, atom: AtomId) -> Option<AtomMut<'_>> {
        self.graph.atom(atom)?;
        Some(AtomMut {
            editor: self,
            id: atom,
        })
    }

    pub(crate) fn require_atom_mut(&mut self, atom: AtomId) -> Result<AtomMut<'_>, MoleculeError> {
        self.atom_mut(atom)
            .ok_or(MoleculeError::InvalidAtomId(atom))
    }

    /// Borrows a bond with checked order editing through [`BondMut::set_order`].
    /// Change connectivity with [`Self::set_bond_endpoints`] or [`Self::replace_bond`].
    ///
    /// ```compile_fail
    /// use kekule::molecule::{AtomId, Bond, BondId, BondOrder, MoleculeEditor};
    /// fn bypass(editor: &mut MoleculeEditor, id: BondId) {
    ///     *editor.bond_mut(id).unwrap() = Bond::new(AtomId::new(0), AtomId::new(99), BondOrder::Single);
    /// }
    /// ```
    ///
    /// Returns `None` when `bond` names no live bond of this draft.
    pub fn bond_mut(&mut self, bond: BondId) -> Option<BondMut<'_>> {
        self.graph.bond(bond)?;
        Some(BondMut {
            editor: self,
            id: bond,
        })
    }

    pub(crate) fn require_bond_mut(&mut self, bond: BondId) -> Result<BondMut<'_>, MoleculeError> {
        self.bond_mut(bond)
            .ok_or(MoleculeError::InvalidBondId(bond))
    }

    /// Adds an unbonded atom and returns its draft ID. Perception and owner
    /// properties are cleared.
    ///
    /// # Errors
    ///
    /// [`MoleculeError::CapacityExceeded`] when the draft already
    /// holds the maximum number of atom slots.
    pub fn add_atom(&mut self, atom: Atom) -> Result<AtomId, MoleculeError> {
        self.add_atom_at_slot(atom, self.graph.atom_slot_count())
    }

    /// [`Self::add_atom`] into slot `slot`, which must be the next one; tests
    /// pass a larger slot to reach the capacity bound without allocating it.
    fn add_atom_at_slot(&mut self, atom: Atom, slot: usize) -> Result<AtomId, MoleculeError> {
        checked_id::<AtomId>(slot, Countable::Atom).map_err(MoleculeError::CapacityExceeded)?;
        debug_assert_eq!(slot, self.graph.atom_slot_count());
        let id = self.graph.push_atom(atom);
        self.properties
            .atoms_mut()
            .resize_missing(self.graph.atom_slot_count());
        self.chemistry_changed();
        Ok(id)
    }

    /// Deletes an atom and its bonds, removing every stereo element that
    /// references the atom and those its deleted bonds invalidate. See
    /// [`Self::pruned_stereo_elements`].
    pub fn delete_atom(&mut self, atom: AtomId) -> Result<Atom, MoleculeError> {
        self.graph.require_atom(atom)?;
        let incident = self.graph.adjacency[atom.index()].clone();
        for bond in incident {
            self.delete_bond(bond)?;
        }
        let removed = self
            .graph
            .take_atom(atom)
            .ok_or(MoleculeError::InvalidAtomId(atom))?;
        self.properties.atoms_mut().clear_index(atom.index());
        self.prune_stereo_for_atom(atom);
        self.chemistry_changed();
        Ok(removed)
    }

    /// Adds a bond. Stereo elements it invalidates are removed, such as a
    /// tetrahedral center that gains a neighbor its carriers do not list, or
    /// a double-bond endpoint that gains a third substituent; elements away
    /// from the bond's endpoints are untouched. See
    /// [`Self::pruned_stereo_elements`].
    pub fn add_bond(
        &mut self,
        a: AtomId,
        b: AtomId,
        order: BondOrder,
    ) -> Result<BondId, MoleculeError> {
        let bond = self.add_bond_unpruned(a, b, order)?;
        self.prune_invalid_stereo_at([a, b], bond);
        Ok(bond)
    }

    /// [`Self::add_bond`] leaving existing stereo as it is, for transforms
    /// that rewrite the carriers of the new bond's endpoints themselves.
    pub(crate) fn add_bond_unpruned(
        &mut self,
        a: AtomId,
        b: AtomId,
        order: BondOrder,
    ) -> Result<BondId, MoleculeError> {
        self.graph.require_atom(a)?;
        self.graph.require_atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.graph.bond_between(a, b).is_some() {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        checked_id::<BondId>(self.graph.bond_slot_count(), Countable::Bond)
            .map_err(MoleculeError::CapacityExceeded)?;
        let id = self.graph.push_bond(Bond::new(a, b, order));
        self.properties
            .bonds_mut()
            .resize_missing(self.graph.bond_slot_count());
        self.chemistry_changed();
        Ok(id)
    }

    /// Deletes a bond. As for [`Self::add_bond`], stereo elements it
    /// invalidates are removed, such as an assertion whose carrier it bonded
    /// to the focus, and every other element is untouched. See
    /// [`Self::pruned_stereo_elements`].
    pub fn delete_bond(&mut self, bond: BondId) -> Result<Bond, MoleculeError> {
        let removed = self
            .graph
            .take_bond(bond)
            .ok_or(MoleculeError::InvalidBondId(bond))?;
        // A deleted bond may connect a focus to a carrier even when it is not
        // itself a stereo focus. Ring opening can leave both atoms alive and
        // connected through another path, so live IDs alone are insufficient.
        self.prune_invalid_stereo_at([removed.a, removed.b], bond);
        self.properties.bonds_mut().clear_index(bond.index());
        self.chemistry_changed();
        Ok(removed)
    }

    /// Stereo elements that structural edits of this draft removed because
    /// they no longer described valid stereo, in removal order: those pruned
    /// by adding, deleting, or rewiring bonds, by changing a bond order, by
    /// deleting atoms, and by oxohalogen normalization. Explicit
    /// [`Self::remove_stereo_element`] calls are not listed. IDs are draft
    /// IDs, which are never reused within a draft.
    pub fn pruned_stereo_elements(&self) -> &[StereoElementId] {
        &self.pruned_stereo
    }

    /// Returns [`Self::pruned_stereo_elements`] and clears the record, so a
    /// caller can check one edit at a time.
    pub fn take_pruned_stereo_elements(&mut self) -> Vec<StereoElementId> {
        std::mem::take(&mut self.pruned_stereo)
    }

    /// Replaces represented atom state while retaining its ID and properties.
    pub fn replace_atom(&mut self, id: AtomId, atom: Atom) -> Result<Atom, MoleculeError> {
        if self.graph.require_atom(id)? == &atom {
            return Ok(atom);
        }
        Ok(std::mem::replace(&mut *self.require_atom_mut(id)?, atom))
    }

    /// Replaces a bond's endpoints and order, retaining its ID and properties.
    /// Invalid endpoints, self-bonds, and duplicate bonds leave state unchanged.
    /// Rewiring removes stereo assertions referencing the bond or changed
    /// endpoints, whose chemical neighborhoods have changed.
    /// Changing only the order removes assertions focused on this bond.
    pub fn replace_bond(&mut self, id: BondId, replacement: Bond) -> Result<Bond, MoleculeError> {
        let previous = self.graph.require_bond(id)?.clone();
        let (a, b) = replacement.endpoints();
        self.graph.require_atom(a)?;
        self.graph.require_atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self
            .graph
            .require_bond_between(a, b)?
            .is_some_and(|other| other != id)
        {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        if previous == replacement {
            return Ok(previous);
        }
        if previous.connects(a, b) {
            self.require_bond_mut(id)?.set_order(replacement.order);
            return Ok(previous);
        }
        self.graph.rewire_bond(id, replacement);
        let invalid = self
            .graph
            .stereo_elements()
            .filter_map(|(element_id, element)| {
                (element.references_bond(id)
                    || [previous.a(), previous.b(), a, b]
                        .into_iter()
                        .any(|atom| element.references_atom(atom)))
                .then_some(element_id)
            })
            .collect::<Vec<_>>();
        for &element in &invalid {
            self.remove_stereo_element(element)?;
        }
        self.pruned_stereo.extend(invalid);
        self.chemistry_changed();
        Ok(previous)
    }

    /// Changes represented order, removing stereo assertions focused on this bond.
    /// Assigning the current order leaves chemistry, stereo, and annotations intact.
    pub fn set_bond_order(&mut self, id: BondId, order: BondOrder) -> Result<(), MoleculeError> {
        self.require_bond_mut(id)?.set_order(order);
        Ok(())
    }

    /// Rewires a bond to new endpoints, keeping its order; see
    /// [`Self::replace_bond`] for the stereo it removes.
    ///
    /// # Errors
    ///
    /// As [`Self::replace_bond`]: [`MoleculeError::InvalidBondId`],
    /// [`MoleculeError::InvalidAtomId`], [`MoleculeError::SelfBond`], or
    /// [`MoleculeError::DuplicateBond`], leaving the draft unchanged.
    pub fn set_bond_endpoints(
        &mut self,
        id: BondId,
        a: AtomId,
        b: AtomId,
    ) -> Result<(), MoleculeError> {
        self.replace_bond(id, Bond::new(a, b, self.graph.require_bond(id)?.order))?;
        Ok(())
    }

    /// Deletes a set of live atoms and their incident bonds. IDs are validated
    /// before mutation; duplicates are ignored. Surviving IDs are never renumbered.
    pub fn delete_atoms(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Result<Vec<(AtomId, Atom)>, MoleculeError> {
        let atoms = atoms.into_iter().collect::<BTreeSet<_>>();
        for &id in &atoms {
            self.graph.require_atom(id)?;
        }
        atoms
            .into_iter()
            .map(|id| self.delete_atom(id).map(|atom| (id, atom)))
            .collect()
    }

    /// Deletes a set of live bonds after validating every ID. Duplicates are ignored.
    pub fn delete_bonds(
        &mut self,
        bonds: impl IntoIterator<Item = BondId>,
    ) -> Result<Vec<(BondId, Bond)>, MoleculeError> {
        let bonds = bonds.into_iter().collect::<BTreeSet<_>>();
        for &id in &bonds {
            self.graph.require_bond(id)?;
        }
        bonds
            .into_iter()
            .map(|id| self.delete_bond(id).map(|bond| (id, bond)))
            .collect()
    }

    /// Retains the induced graph on these atoms; an empty set leaves an empty editor.
    pub fn retain_atoms(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Result<Vec<(AtomId, Atom)>, MoleculeError> {
        let retained = atoms.into_iter().collect::<BTreeSet<_>>();
        for &id in &retained {
            self.graph.require_atom(id)?;
        }
        let removed = self
            .graph
            .atom_ids()
            .filter(|id| !retained.contains(id))
            .collect::<Vec<_>>();
        self.delete_atoms(removed)
    }

    /// Clears graph, properties, and perception, restarting the local ID space.
    pub fn clear(&mut self) {
        *self = Self::new();
    }

    /// Appends a published fragment transactionally and returns semantic ID maps.
    /// Copies live atom/bond properties and represented stereo, including groups.
    /// Source owner properties and perception are not transferred. Successful
    /// append clears target owner properties and perception as a structural edit.
    /// A temporary disconnected result is allowed; connect it before finishing.
    /// This compound operation stages a clone of the target for rollback on error.
    pub fn append_molecule(
        &mut self,
        source: &Molecule,
    ) -> Result<MoleculeCorrespondence, MoleculeError> {
        self.append_graph(&source.graph, &source.properties)
    }

    /// Appends another draft, as topology edits do when they join two
    /// instances that are both being edited.
    pub(crate) fn append_draft(
        &mut self,
        source: &MoleculeEditor,
    ) -> Result<MoleculeCorrespondence, MoleculeError> {
        self.append_graph(&source.graph, &source.properties)
    }

    fn append_graph<L: SlotLiveness>(
        &mut self,
        source: &Graph<L>,
        properties: &MoleculeProperties,
    ) -> Result<MoleculeCorrespondence, MoleculeError> {
        let mut staged = self.clone();
        let mut map = MoleculeCorrespondence::default();
        for (id, atom) in source.atoms() {
            map.atoms.insert(id, staged.add_atom(atom.clone())?);
        }
        for (id, bond) in source.bonds() {
            map.bonds.insert(
                id,
                staged.add_bond(map.atoms[&bond.a()], map.atoms[&bond.b()], bond.order)?,
            );
        }
        // Property rows are slots. A draft source keeps rows for deleted slots,
        // so copy only live rows.
        let (live, atom_rows): (Vec<_>, Vec<_>) = map
            .atoms
            .iter()
            .map(|(source, draft)| (source.index(), draft.index()))
            .unzip();
        let atoms = properties
            .atoms()
            .raw()
            .select_indices(&live)
            .map_err(MoleculeError::Property)?;
        staged
            .properties
            .atoms_mut()
            .copy_rows_from(&atoms, &atom_rows)
            .map_err(MoleculeError::Property)?;
        let (live, bond_rows): (Vec<_>, Vec<_>) = map
            .bonds
            .iter()
            .map(|(source, draft)| (source.index(), draft.index()))
            .unzip();
        let bonds = properties
            .bonds()
            .raw()
            .select_indices(&live)
            .map_err(MoleculeError::Property)?;
        staged
            .properties
            .bonds_mut()
            .copy_rows_from(&bonds, &bond_rows)
            .map_err(MoleculeError::Property)?;
        for (id, element) in source.stereo_elements() {
            let mut element = StereoElement::new(element.kind.clone());
            element.map_references(|atom| map.atoms[&atom], |bond| map.bonds[&bond]);
            map.stereo_elements
                .insert(id, staged.add_stereo_element(element)?);
        }
        for (id, group) in source.stereo_groups() {
            let mut group = group.clone();
            for member in &mut group.members {
                *member = map.stereo_elements[member];
            }
            map.stereo_groups
                .insert(id, staged.add_stereo_group(group)?);
        }
        *self = staged;
        Ok(map)
    }
}

/// Identifier correspondence from one molecular graph to another.
///
/// Every source atom, bond, stereo element, and stereo group with a
/// counterpart has one entry; deleted sources have none. It is returned by
/// [`MoleculeEditor::append_molecule`] (fragment IDs to draft IDs),
/// [`MoleculeEditor::finish_with_correspondence`] (draft IDs to published IDs),
/// and hydrogen removal (input IDs to transformed IDs).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MoleculeCorrespondence {
    atoms: BTreeMap<AtomId, AtomId>,
    bonds: BTreeMap<BondId, BondId>,
    stereo_elements: BTreeMap<StereoElementId, StereoElementId>,
    stereo_groups: BTreeMap<StereoGroupId, StereoGroupId>,
}

impl MoleculeCorrespondence {
    pub(crate) fn from_compaction(compaction: &SlotCompaction) -> Self {
        fn entries<Id: Copy + Ord>(
            slots: &[Option<Id>],
            source: impl Fn(u32) -> Id,
        ) -> BTreeMap<Id, Id> {
            (0..=u32::MAX)
                .zip(slots)
                .filter_map(|(raw, target)| target.map(|target| (source(raw), target)))
                .collect()
        }
        Self {
            atoms: entries(&compaction.atoms, AtomId::new),
            bonds: entries(&compaction.bonds, BondId::new),
            stereo_elements: entries(&compaction.stereo_elements, StereoElementId::new),
            stereo_groups: entries(&compaction.stereo_groups, StereoGroupId::new),
        }
    }

    /// Maps every entity of `molecule` to itself.
    pub(crate) fn identity(molecule: &Molecule) -> Self {
        Self {
            atoms: molecule.atom_ids().map(|id| (id, id)).collect(),
            bonds: molecule.bond_ids().map(|id| (id, id)).collect(),
            stereo_elements: molecule.stereo_element_ids().map(|id| (id, id)).collect(),
            stereo_groups: molecule.stereo_groups().map(|(id, _)| (id, id)).collect(),
        }
    }

    /// The counterpart of a source atom, or `None` when it has none.
    pub fn atom(&self, source: AtomId) -> Option<AtomId> {
        self.atoms.get(&source).copied()
    }
    /// The counterpart of a source bond, or `None` when it has none.
    pub fn bond(&self, source: BondId) -> Option<BondId> {
        self.bonds.get(&source).copied()
    }
    /// The counterpart of a source stereo element, or `None` when it has none.
    pub fn stereo_element(&self, source: StereoElementId) -> Option<StereoElementId> {
        self.stereo_elements.get(&source).copied()
    }
    /// The counterpart of a source stereo group, or `None` when it has none.
    pub fn stereo_group(&self, source: StereoGroupId) -> Option<StereoGroupId> {
        self.stereo_groups.get(&source).copied()
    }
    /// Every atom entry, keyed by source ID.
    pub fn atoms(&self) -> &BTreeMap<AtomId, AtomId> {
        &self.atoms
    }
    /// Every bond entry, keyed by source ID.
    pub fn bonds(&self) -> &BTreeMap<BondId, BondId> {
        &self.bonds
    }
    /// Every stereo element entry, keyed by source ID.
    pub fn stereo_elements(&self) -> &BTreeMap<StereoElementId, StereoElementId> {
        &self.stereo_elements
    }
    /// Every stereo group entry, keyed by source ID.
    pub fn stereo_groups(&self) -> &BTreeMap<StereoGroupId, StereoGroupId> {
        &self.stereo_groups
    }
}

#[cfg(test)]
mod tests {
    use super::publication::{validate_graph, GraphValidationError};
    use super::*;
    use crate::stereo::StereoViolation;
    use crate::stereo::{
        DoubleBondOrientation, DoubleBondStereo, StereoCarrier, StereoElement, StereoElementKind,
    };

    fn carbon() -> Atom {
        Atom::new(Element::from_atomic_number(6).expect("carbon"))
    }

    fn atom(symbol: &str) -> Atom {
        Atom::new(Element::from_symbol(symbol).expect("fixture element"))
    }

    #[test]
    fn graph_validation_checks_property_domains() {
        for bonds in [false, true] {
            let mut editor = MoleculeEditor::new();
            editor.add_atom(carbon()).unwrap();
            if bonds {
                editor.properties.bonds_mut().resize_missing(1);
            } else {
                editor.properties.atoms_mut().resize_missing(0);
            }
            let expected = if bonds {
                GraphValidationError::BondPropertySlotCount {
                    expected: 0,
                    actual: 1,
                }
            } else {
                GraphValidationError::AtomPropertySlotCount {
                    expected: 1,
                    actual: 0,
                }
            };
            assert_eq!(validate_graph(&editor), Err(expected));
        }
    }

    #[test]
    fn stereo_publication_issues_have_diagnostic_display_messages() {
        let error =
            MoleculePublicationError::InvalidStereo(StereoViolation::DuplicateStereoGroupMember {
                group: StereoGroupId::new(2),
                element: StereoElementId::new(7),
            });
        assert_eq!(
            error.to_string(),
            "Stereo group 2 lists stereo element 7 twice."
        );
    }

    #[test]
    fn builder_rejects_disconnected_final_graph() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        builder.add_atom(carbon()).expect("first atom");
        builder.add_atom(carbon()).expect("second atom");
        assert_eq!(
            builder.finish(),
            Err(MoleculePublicationError::Disconnected { components: 2 })
        );
    }

    #[test]
    fn builder_accepts_connected_graph() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        builder
            .add_bond(a, b, BondOrder::Single)
            .expect("connecting bond");
        assert_eq!(
            builder
                .finish()
                .expect("connected molecule")
                .connected_components()
                .len(),
            1
        );
    }

    #[test]
    fn publication_rejects_a_stereo_focus_changed_by_internal_graph_code() {
        let molecule = crate::smiles::to_molecules("F/C=C/F")
            .unwrap()
            .pop()
            .unwrap();
        let (element, stereo) = molecule.stereo_elements().next().unwrap();
        let StereoElementKind::DoubleBond(stereo) = &stereo.kind else {
            panic!("expected double-bond stereo");
        };
        let bond = stereo.bond;
        let mut editor = molecule.into_editor();
        // Internal algorithms can touch graph storage directly; publication
        // must still reject an assertion invalidated outside the public setter.
        editor.graph.bond_mut(bond).unwrap().order = BondOrder::Single;
        let error = editor.try_finish().unwrap_err();
        assert_eq!(
            error.error(),
            &MoleculePublicationError::InvalidStereo(StereoViolation::InvalidDoubleBondOrder {
                element,
                bond,
                order: BondOrder::Single,
            })
        );
        let mut repaired = error.into_draft();
        repaired.remove_stereo_element(element).unwrap();
        assert_eq!(repaired.finish().unwrap().stereo_elements().count(), 0);
    }

    #[test]
    fn publication_rejects_nonadjacent_stereo_carriers_from_internal_graph_code() {
        let molecule = crate::smiles::to_molecules("F[C@]1(Cl)CCNC1")
            .unwrap()
            .pop()
            .unwrap();
        let (element, original) = molecule.stereo_elements().next().unwrap();
        let mut invalid = original.clone();
        let StereoElementKind::Tetrahedral(stereo) = &mut invalid.kind else {
            panic!("expected tetrahedral stereo");
        };
        let center = stereo.center;
        stereo.carriers[0] = StereoCarrier::Atom(center);
        let mut editor = molecule.into_editor();
        *editor.graph.stereo_element_mut(element).unwrap() = invalid;
        let failure = editor.try_finish().unwrap_err();
        assert_eq!(
            failure.error(),
            &MoleculePublicationError::InvalidStereo(
                StereoViolation::TetrahedralCarrierNotAdjacent {
                    element,
                    center,
                    carrier: StereoCarrier::Atom(center),
                }
            )
        );
        let mut repaired = failure.into_draft();
        repaired.remove_stereo_element(element).unwrap();
        repaired.finish().unwrap();
    }

    #[test]
    fn canonical_bond_order_rewriting_prunes_obsolete_stereo_focus() {
        let mut editor = MoleculeEditor::new();
        let chlorine = editor.add_atom(atom("Cl")).unwrap();
        let oxo = editor.add_atom(atom("O")).unwrap();
        let hydroxyl = editor.add_atom(atom("O")).unwrap();
        let bond = editor.add_bond(chlorine, oxo, BondOrder::Double).unwrap();
        editor
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .unwrap();
        editor
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                DoubleBondStereo {
                    bond,
                    left: chlorine,
                    right: oxo,
                    left_carrier: StereoCarrier::Atom(hydroxyl),
                    right_carrier: StereoCarrier::ImplicitHydrogen,
                    orientation: None,
                },
            )))
            .unwrap();
        assert_eq!(editor.normalize_oxohalogens(), [chlorine]);
        let molecule = editor.finish().unwrap();
        assert_eq!(molecule.bond(bond).unwrap().order, BondOrder::Single);
        assert_eq!(molecule.stereo_elements().count(), 0);
    }

    #[test]
    fn explicit_normalization_publishes_canonical_hypervalent_representation() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
        let oxo = builder.add_atom(atom("O")).expect("oxo oxygen");
        let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
        builder
            .add_bond(chlorine, oxo, BondOrder::Double)
            .expect("source-convention oxo bond");
        builder
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .expect("hydroxyl bond");

        let before = builder.clone();
        assert_eq!(
            builder.validate(),
            Err(MoleculePublicationError::UnnormalizedOxohalogen { atom: chlorine })
        );
        assert_eq!(builder, before, "validation never rewrites chemistry");
        let error = builder.try_finish().unwrap_err();
        assert_eq!(
            error.error(),
            &MoleculePublicationError::UnnormalizedOxohalogen { atom: chlorine }
        );
        let mut builder = error.into_draft();
        assert_eq!(builder, before);
        assert_eq!(builder.normalize_oxohalogens(), [chlorine]);
        assert!(
            builder.normalize_oxohalogens().is_empty(),
            "normalization is idempotent"
        );
        let molecule = builder.finish().expect("canonical publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap();
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
        assert_eq!(molecule.perception(), &Perception::default());
    }

    #[test]
    fn failed_editor_commit_leaves_original_unchanged() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        let bond = builder
            .add_bond(a, b, BondOrder::Single)
            .expect("connecting bond");
        let molecule = builder.finish().expect("connected molecule");
        let before = molecule.clone();

        let mut editor = molecule.to_editor();
        editor.delete_bond(bond).expect("delete in working copy");
        assert_eq!(
            editor.finish(),
            Err(MoleculePublicationError::Disconnected { components: 2 })
        );
        assert_eq!(molecule, before);
    }

    #[test]
    fn editor_can_pass_through_temporarily_disconnected_state() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let a = builder.add_atom(carbon()).expect("first atom");
        let b = builder.add_atom(carbon()).expect("second atom");
        let c = builder.add_atom(carbon()).expect("third atom");
        let ab = builder
            .add_bond(a, b, BondOrder::Single)
            .expect("first bond");
        builder
            .add_bond(b, c, BondOrder::Single)
            .expect("second bond");
        let molecule = builder.finish().expect("connected molecule");

        let mut editor = molecule.to_editor();
        editor.delete_bond(ab).expect("temporary disconnect");
        editor
            .add_bond(a, c, BondOrder::Single)
            .expect("reconnect through another edge");
        let molecule = editor.finish().expect("connected final graph");
        assert_eq!(molecule.connected_components().len(), 1);
    }

    #[test]
    fn editor_rejects_unnormalized_oxohalogen_edits_until_normalized() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
        let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
        builder
            .add_bond(chlorine, hydroxyl, BondOrder::Single)
            .expect("initial bond");
        let molecule = builder.finish().expect("initial molecule");

        let mut editor = molecule.to_editor();
        let oxo = editor.add_atom(atom("O")).expect("oxo oxygen");
        editor
            .add_bond(chlorine, oxo, BondOrder::Double)
            .expect("source-convention oxo bond");
        assert!(matches!(
            editor.validate(),
            Err(MoleculePublicationError::UnnormalizedOxohalogen { .. })
        ));
        editor.normalize_oxohalogens();
        let molecule = editor.finish().expect("canonical edit publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap();
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    }

    #[test]
    fn explicit_normalization_covers_direct_bond_order_changes() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let anchor = builder.add_atom(atom("O")).unwrap();
        let oxo = builder.add_atom(atom("O")).unwrap();
        builder
            .add_bond(chlorine, anchor, BondOrder::Single)
            .unwrap();
        let oxo_bond = builder.add_bond(chlorine, oxo, BondOrder::Single).unwrap();
        let molecule = builder.finish().unwrap();

        let mut editor = molecule.to_editor();
        editor
            .bond_mut(oxo_bond)
            .unwrap()
            .set_order(BondOrder::Double);
        assert_eq!(editor.normalize_oxohalogens(), [chlorine]);
        let molecule = editor.finish().expect("canonical edit publication");

        assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
        assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
        assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    }

    #[test]
    fn editor_recanonicalizes_stereo_references_after_topology_changes() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let left = builder.add_atom(carbon()).unwrap();
        let right = builder.add_atom(carbon()).unwrap();
        let hydrogen = builder.add_atom(atom("H")).unwrap();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let bromine = builder.add_atom(atom("Br")).unwrap();
        let double_bond = builder.add_bond(left, right, BondOrder::Double).unwrap();
        builder.add_bond(left, hydrogen, BondOrder::Single).unwrap();
        builder
            .add_bond(right, chlorine, BondOrder::Single)
            .unwrap();
        builder.add_bond(right, bromine, BondOrder::Single).unwrap();
        let molecule = builder.finish().unwrap();
        let mut editor = molecule.to_editor();
        let element = editor
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                DoubleBondStereo {
                    bond: double_bond,
                    left,
                    right,
                    left_carrier: StereoCarrier::Atom(hydrogen),
                    right_carrier: StereoCarrier::Atom(chlorine),
                    orientation: Some(DoubleBondOrientation::Together),
                },
            )))
            .unwrap();
        let molecule = editor.finish().unwrap();

        let mut editor = molecule.to_editor();
        let fluorine = editor.add_atom(atom("F")).unwrap();
        editor.add_bond(left, fluorine, BondOrder::Single).unwrap();
        let molecule = editor.finish().expect("canonical edit publication");

        let StereoElementKind::DoubleBond(stereo) = &molecule.stereo_element(element).unwrap().kind
        else {
            panic!("expected double-bond stereo");
        };
        assert_eq!(stereo.left_carrier, StereoCarrier::Atom(fluorine));
        assert_eq!(stereo.right_carrier, StereoCarrier::Atom(chlorine));
        assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
    }

    #[test]
    fn editor_publication_leaves_out_of_model_oxohalogens_for_valence_validation() {
        let mut builder = crate::molecule::MoleculeEditor::new();
        let chlorine = builder.add_atom(atom("Cl")).unwrap();
        let anchor = builder.add_atom(atom("O")).unwrap();
        builder
            .add_bond(chlorine, anchor, BondOrder::Single)
            .unwrap();
        let mut oxo_bonds = Vec::new();
        for _ in 0..128 {
            let oxygen = builder.add_atom(atom("O")).unwrap();
            oxo_bonds.push(
                builder
                    .add_bond(chlorine, oxygen, BondOrder::Single)
                    .unwrap(),
            );
        }
        let molecule = builder.finish().unwrap();
        let before = molecule.clone();

        let mut editor = molecule.to_editor();
        for &bond in &oxo_bonds {
            editor.bond_mut(bond).unwrap().set_order(BondOrder::Double);
        }
        let mut edited = editor
            .finish()
            .expect("publication does not validate valence");
        assert_eq!(edited.atom(chlorine).unwrap().formal_charge, 0);
        assert!(oxo_bonds
            .iter()
            .all(|&bond| edited.bond(bond).unwrap().order == BondOrder::Double));
        assert!(edited.perceive().is_err(), "valence must reject this graph");
        assert_eq!(molecule, before);
    }
}

#[cfg(all(test, target_pointer_width = "64"))]
mod capacity_tests {
    use super::*;
    use crate::diagnostics::CapacityExceeded;

    fn first_unsupported_slot() -> usize {
        usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize")
    }

    #[test]
    fn capacity_rejection_does_not_mutate_molecule_state() {
        let mut molecule = crate::molecule::MoleculeEditor::new();
        let before = molecule.clone();
        assert_eq!(
            molecule.add_atom_at_slot(
                Atom::new(Element::from_atomic_number(6).expect("carbon")),
                first_unsupported_slot(),
            ),
            Err(MoleculeError::CapacityExceeded(CapacityExceeded::new(
                Countable::Atom
            )))
        );
        assert_eq!(molecule, before);
    }
}
