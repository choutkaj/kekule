//! Transactional system edits, with molecular chemistry staged per affected occurrence.
mod append;
mod error;
mod hierarchy;
mod identity;
mod overlay;
mod placement;
mod properties;
mod publication;
pub use publication::{DenseCorrespondence, EditCorrespondence};

use super::dense_map::{DenseKey, DenseMap};
use super::shared_rows::SharedRow;
use super::{
    AtomSiteId, AtomSiteMetadata, ChainId, HierarchyError, InstanceAtomId, InstanceBondId,
    MoleculeClass, MoleculeDefinitionId, MoleculeInstanceId, ResidueClass, ResidueId,
    ResidueLabels, Topology, TopologyBuildError, TopologyBuilder,
};
use crate::molecule::{
    Atom, AtomId, BondId, BondOrder, Molecule, MoleculeEditor, MoleculeError,
    MoleculePublicationError,
};
use crate::properties::{
    OwnerProperties, PropertyColumn, PropertyError, PropertyKey, PropertyValue, RawPropertyTable,
    TopologyProperties,
};
pub(crate) use append::AppendMapping;
pub use error::TopologyEditError;
pub use hierarchy::{EditAtomSite, EditChain, EditResidue};
use hierarchy::{EditHierarchy, SourceHierarchy};
pub use identity::{EditAtomId, EditAtomSiteId, EditBondId, EditChainId, EditResidueId};
use overlay::{HandleBlock, Overlay};
pub(crate) use placement::{AttachedInOccurrence, DraftHydrogenError};
pub use properties::{EditProperties, EditPropertiesMut, EditPropertyTable, EditPropertyTableMut};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Arc;

#[derive(Debug, Clone)]
enum GroupChemistry {
    Source(MoleculeDefinitionId),
    /// An added molecule, and whether it holds exactly the default
    /// perception, as a definition imported from a perceived snapshot does.
    Added(Arc<Molecule>, bool),
    Draft(Box<MoleculeEditor>),
}
/// One staging group: an untouched source occurrence read from the source,
/// a stored group, or a group merged into another.
#[derive(Debug, Clone)]
enum GroupSlot {
    /// The occurrence of the source instance with this group's index, which
    /// no edit has touched.
    Source,
    Live(Box<Group>),
    Merged,
}

#[derive(Debug, Clone)]
struct Group {
    chemistry: GroupChemistry,
    atoms: DenseMap<AtomId, EditAtomId>,
    bonds: DenseMap<BondId, EditBondId>,
    instance_slot: Option<usize>,
    changed: bool,
    class: Option<MoleculeClass>,
    class_explicit: bool,
}
#[derive(Debug, Clone, Copy)]
pub(crate) struct Location<Id> {
    group: usize,
    local: Id,
    slot: usize,
}

/// Editing handles and storage slots for one registered group, listed in the
/// group molecule's local atom and bond ID order.
///
/// Slots are draft rows and define publication order: surviving atoms keep
/// their source dense order and added atoms follow in creation order. Handles
/// are allocated in slot order, so handle order matches slot order too.
struct GroupIdentity {
    atoms: Vec<(EditAtomId, usize)>,
    bonds: Vec<(EditBondId, usize)>,
}

impl GroupIdentity {
    /// Fresh handles and slots appended in local ID order.
    fn appended(molecule: &Molecule, atom_start: usize, bond_start: usize) -> Self {
        Self {
            atoms: (atom_start..)
                .zip(molecule.atom_ids())
                .map(|(slot, _)| (EditAtomId::new(), slot))
                .collect(),
            bonds: (bond_start..)
                .zip(molecule.bond_ids())
                .map(|(slot, _)| (EditBondId::new(), slot))
                .collect(),
        }
    }

    /// One instance of a dense source layout, with handles taken from the
    /// given blocks in source dense order and slots offset by the given bases.
    fn dense(
        source: &Topology,
        instance: MoleculeInstanceId,
        handles: (HandleBlock, HandleBlock),
        bases: (usize, usize),
    ) -> Self {
        let molecule = source
            .instance(instance)
            .expect("registered source instance exists")
            .molecule();
        let atom = |local| {
            let dense = source
                .atom_index(InstanceAtomId::new(instance, local))
                .expect("source instance atom has a dense index")
                .index();
            (handles.0.handle(dense), bases.0 + dense)
        };
        let bond = |local| {
            let dense = source
                .bond_index(InstanceBondId::new(instance, local))
                .expect("source instance bond has a dense index")
                .index();
            (handles.1.handle(dense), bases.1 + dense)
        };
        Self {
            atoms: molecule.atom_ids().map(atom).collect(),
            bonds: molecule.bond_ids().map(bond).collect(),
        }
    }

    /// One consecutive handle block per dense source atom and bond, in dense
    /// order.
    fn source_handles(source: &Topology) -> (HandleBlock, HandleBlock) {
        (
            HandleBlock::reserve(source.atom_count()),
            HandleBlock::reserve(source.bond_count()),
        )
    }
}

/// Detached coordinate-free structural editing state.
///
/// Chemical edits affect individual occurrences, even when definitions are reused.
/// Deleting bonds can split molecules; adding bonds can merge them. Publication
/// constructs valid connected definitions and one immutable topology snapshot.
/// Its dense atom order keeps surviving atoms in source order, followed by added
/// atoms in creation order; [`Self::atom_ids`] already iterates in that order.
/// Stable opaque handles survive these changes within the draft. Use `*_handle`
/// to resolve source IDs. [`Self::finish`] returns the completed topology; editing
/// handles do not identify entities in the published result.
///
/// Graph edits clear changed owner annotations. Surviving entity annotations are
/// transferred explicitly; newly added entities have missing property values.
/// Untouched definitions retain perception and classification. Hierarchy changes
/// affecting their atoms invalidate inferred classifications; explicit assignments
/// survive metadata changes. Changes to represented chemistry or residue
/// composition require a fresh override. Unchanged inferred classes, including
/// classes retained by a prior subset, are preserved without reclassification.
/// No-op publication
/// retains the exact input `Arc<Topology>` supplied to [`Self::from_topology`].
///
/// ```
/// use kekule::{smiles, topology::{Topology, edit::TopologyEditor}};
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Topology::from_molecule(molecule)?;
/// let source_bond = topology.bond_ids()[0];
/// let mut editor = TopologyEditor::from_topology(topology);
/// editor.delete_bond(editor.bond_handle(source_bond).expect("source bond"))?;
/// let edited = editor.finish()?;
/// assert_eq!(edited.instance_count(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct TopologyEditor {
    source: Option<Arc<Topology>>,
    // The group of source instance `i` is group `i`, read from the source
    // until an edit touches it.
    groups: Vec<GroupSlot>,
    // Source atoms and bonds are read from the source until an edit moves or
    // deletes them: a source handle's block offset is its source dense index,
    // which is also its slot.
    atoms: Overlay<EditAtomId, Location<AtomId>>,
    bonds: Overlay<EditBondId, Location<BondId>>,
    molecule_classes: BTreeMap<EditAtomId, MoleculeClass>,
    hierarchy: EditHierarchy,
    properties: TopologyProperties,
    revision: u64,
    pub(crate) structural_revision: u64,
    // Storage slots of atoms deleted from this draft, so a slot's liveness is
    // known without walking every live atom.
    dead_atom_slots: BTreeSet<usize>,
    // Whether the draft was cleared, so no slot names a source entity.
    cleared: bool,
    // Slots of deleted atoms and bonds in deletion order, so coordinated
    // state clears exactly the rows one edit freed.
    deleted_slots: (Vec<usize>, Vec<usize>),
}

/// The unchanged location of the source atom at dense index `index`.
fn source_atom(source: Option<&Topology>, index: usize) -> Location<AtomId> {
    let id = source.expect("source handles have a source").atom_ids()[index];
    Location {
        group: id.instance().index(),
        local: id.atom(),
        slot: index,
    }
}

/// The unchanged location of the source bond at dense index `index`.
fn source_bond(source: Option<&Topology>, index: usize) -> Location<BondId> {
    let id = source.expect("source handles have a source").bond_ids()[index];
    Location {
        group: id.instance().index(),
        local: id.bond(),
        slot: index,
    }
}

impl Topology {
    /// Starts a detached draft from a copy of this topology.
    ///
    /// To edit a shared snapshot without copying it, pass the `Arc` to
    /// [`TopologyEditor::from_topology`]: molecular data is then copied only
    /// when edited, and a no-op publication returns that same allocation. Use
    /// [`Self::into_editor`] to transfer an owned topology.
    pub fn to_editor(&self) -> TopologyEditor {
        TopologyEditor::from_topology(self.clone())
    }
    /// Transfers an owned topology into a structural draft. Successful consuming
    /// publication can move its unchanged molecular definitions into the result.
    pub fn into_editor(self) -> TopologyEditor {
        TopologyEditor::from_topology(self)
    }
}

impl TopologyEditor {
    /// An empty draft with no source topology.
    ///
    /// Publishing it fails until atoms are added.
    pub fn new() -> Self {
        Self::default()
    }

    /// Clears the draft while retaining source correspondence for deleted entities.
    /// Previously allocated handles cannot be reused for newly added entities.
    pub fn clear(&mut self) {
        self.groups.clear();
        self.atoms.clear();
        self.bonds.clear();
        self.molecule_classes.clear();
        self.hierarchy.clear();
        self.properties = TopologyProperties::default();
        self.dead_atom_slots.clear();
        self.deleted_slots = Default::default();
        self.cleared = true;
        self.changed();
    }

    /// Retains the supplied immutable source snapshot without cloning chemistry.
    ///
    /// The draft reads the source until an edit changes something, so
    /// starting it costs nothing per atom: each edit copies only the
    /// occurrences, atoms and hierarchy nodes it touches.
    pub fn from_topology(topology: impl Into<Arc<Topology>>) -> Self {
        let source = topology.into();
        // Source slots are source dense indices, so an unedited draft
        // republishes the source order exactly. Source handles are one block
        // in dense order, so a handle's offset is its source index.
        let (atoms, bonds) = GroupIdentity::source_handles(&source);
        let mut editor = Self {
            properties: source.properties().clone(),
            groups: vec![GroupSlot::Source; source.instance_count()],
            atoms: Overlay::new(atoms),
            bonds: Overlay::new(bonds),
            source: Some(Arc::clone(&source)),
            ..Self::default()
        };
        editor.import_hierarchy(&source);
        editor
    }

    /// The source snapshot, or `None` for a draft started with [`Self::new`].
    pub fn source_topology(&self) -> Option<&Topology> {
        self.source.as_deref()
    }
    /// The number of live draft atoms.
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
    /// The number of live draft bonds.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }
    /// Whether the draft has no live atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
    /// Live atom handles in publication order.
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = EditAtomId> + '_ {
        self.atoms.keys()
    }
    /// Live bond handles in handle order.
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = EditBondId> + '_ {
        self.bonds.keys()
    }
    /// Live atoms with their handles, in publication order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (EditAtomId, &Atom)> + '_ {
        self.atom_locations().map(|(id, loc)| {
            (
                id,
                self.molecule(loc.group)
                    .atom(loc.local)
                    .expect("live draft atom"),
            )
        })
    }
    /// Inspects bonds using editing handles rather than molecule-local endpoints.
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (EditBondId, EditBond)> + '_ {
        self.bond_ids()
            .map(|id| (id, self.bond(id).expect("live draft bond")))
    }
    /// The draft handle of a source atom, or `None` when the atom is not in the
    /// source or was deleted.
    pub fn atom_handle(&self, source: InstanceAtomId) -> Option<EditAtomId> {
        let index = self.source.as_ref()?.atom_index(source)?.index();
        Some(self.atoms.block().handle(index)).filter(|id| self.atoms.contains_key(id))
    }
    /// The draft handle of a source bond, or `None` when the bond is not in the
    /// source or was deleted.
    pub fn bond_handle(&self, source: InstanceBondId) -> Option<EditBondId> {
        let index = self.source.as_ref()?.bond_index(source)?.index();
        Some(self.bonds.block().handle(index)).filter(|id| self.bonds.contains_key(id))
    }
    /// One draft atom, or `None` when `id` names no live atom of this draft.
    pub fn atom(&self, id: EditAtomId) -> Option<&Atom> {
        let loc = self.atom_location(id).ok()?;
        self.molecule(loc.group).atom(loc.local)
    }
    /// One draft bond with handle endpoints, or `None` when `id` names no live
    /// bond of this draft.
    pub fn bond(&self, id: EditBondId) -> Option<EditBond> {
        let loc = self.bond_location(id).ok()?;
        let bond = self.molecule(loc.group).bond(loc.local)?;
        Some(EditBond {
            a: self.group_atom(loc.group, bond.a()),
            b: self.group_atom(loc.group, bond.b()),
            order: bond.order,
        })
    }
    pub(crate) fn require_atom_handle(
        &self,
        source: InstanceAtomId,
    ) -> Result<EditAtomId, TopologyEditError> {
        self.atom_handle(source)
            .ok_or(TopologyEditError::InvalidSourceAtom(source))
    }
    pub(crate) fn require_atom(&self, id: EditAtomId) -> Result<&Atom, TopologyEditError> {
        self.atom(id).ok_or(TopologyEditError::InvalidAtom(id))
    }
    pub(crate) fn require_bond(&self, id: EditBondId) -> Result<EditBond, TopologyEditError> {
        self.bond(id).ok_or(TopologyEditError::InvalidBond(id))
    }
    /// Handles of the atoms bonded to `id`.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when `id` is deleted or foreign.
    pub fn neighbors(
        &self,
        id: EditAtomId,
    ) -> Result<impl Iterator<Item = EditAtomId> + '_, TopologyEditError> {
        let loc = self.atom_location(id)?;
        Ok(self
            .molecule(loc.group)
            .require_neighbors(loc.local)?
            .map(move |id| self.group_atom(loc.group, id)))
    }
    /// Handles of the bonds incident to `id`.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when `id` is deleted or foreign.
    pub fn incident_bonds(
        &self,
        id: EditAtomId,
    ) -> Result<impl Iterator<Item = EditBondId> + '_, TopologyEditError> {
        let loc = self.atom_location(id)?;
        Ok(self
            .molecule(loc.group)
            .require_incident_bonds(loc.local)?
            .map(move |(id, _)| self.group_bond(loc.group, id)))
    }
    /// The bond joining `a` and `b`, or `None` when they are not bonded.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when either handle is deleted or foreign.
    pub fn bond_between(
        &self,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<Option<EditBondId>, TopologyEditError> {
        let left = self.atom_location(a)?;
        let right = self.atom_location(b)?;
        if left.group != right.group {
            return Ok(None);
        }
        Ok(self
            .molecule(left.group)
            .require_bond_between(left.local, right.local)?
            .map(|id| self.group_bond(left.group, id)))
    }
    /// The current connected components as atom-handle lists.
    ///
    /// Each becomes one molecule instance at publication.
    pub fn connected_components(&self) -> Vec<Vec<EditAtomId>> {
        (0..self.groups.len())
            .filter(|&group| self.group_is_live(group))
            .flat_map(|group| {
                self.molecule(group)
                    .connected_components()
                    .into_iter()
                    .map(|atoms| {
                        atoms
                            .into_iter()
                            .map(|atom| self.group_atom(group, atom))
                            .collect()
                    })
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    /// Adds an isolated atom. Bond it to an existing atom or publish a new occurrence.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::Molecule`] when the molecular draft rejects the atom.
    pub fn add_atom(&mut self, atom: Atom) -> Result<EditAtomId, TopologyEditError> {
        let mut draft = MoleculeEditor::new();
        let local = draft.add_atom(atom)?;
        let id = EditAtomId::new();
        let slot = self.properties.atoms().len();
        let group = self.groups.len();
        self.atoms.insert(id, Location { group, local, slot });
        self.groups.push(GroupSlot::Live(Box::new(Group {
            chemistry: GroupChemistry::Draft(Box::new(draft)),
            atoms: DenseMap::from_iter([(local, id)]),
            bonds: DenseMap::new(),
            instance_slot: None,
            changed: true,
            class: None,
            class_explicit: false,
        })));
        self.properties.resize_atoms(slot + 1);
        self.changed();
        Ok(id)
    }

    /// Adds a complete occurrence, retaining represented stereo and definition properties.
    ///
    /// The returned [`EditMolecule`] maps the molecule's local IDs to new handles.
    ///
    /// # Errors
    ///
    /// Currently never fails; the `Result` is reserved for checked additions.
    pub fn add_molecule(&mut self, molecule: Molecule) -> Result<EditMolecule, TopologyEditError> {
        let owned = Arc::new(molecule);
        let (atom_start, bond_start) = (self.atom_slot_count(), self.bond_slot_count());
        let identity = GroupIdentity::appended(&owned, atom_start, bond_start);
        let group = self.register_group(
            GroupChemistry::Added(Arc::clone(&owned), false),
            &owned,
            None,
            false,
            identity,
        );
        let result = self.group_handles(group);
        self.properties
            .resize_atoms(atom_start + owned.atom_count());
        self.properties
            .resize_bonds(bond_start + owned.bond_count());
        self.changed();
        Ok(result)
    }

    /// Replaces one atom's chemistry and returns the previous atom.
    ///
    /// Replacing an atom with an equal one leaves the draft unchanged. A change
    /// invalidates the class of the atom's component and residue.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when `id` is deleted or foreign.
    /// [`TopologyEditError::Molecule`] when the molecular draft rejects the atom.
    pub fn replace_atom(&mut self, id: EditAtomId, atom: Atom) -> Result<Atom, TopologyEditError> {
        let loc = self.atom_location(id)?;
        if self.require_atom(id)? == &atom {
            return Ok(atom);
        }
        let previous = self.draft_mut(loc.group).replace_atom(loc.local, atom)?;
        self.mark_group_changed(loc.group);
        self.invalidate_component_classes([id]);
        self.invalidate_residue_for_atom(id);
        self.changed();
        Ok(previous)
    }

    /// Deletes one atom with its incident bonds and atom site, returning the atom.
    ///
    /// Deleting a bridging atom splits its molecule at publication.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when `id` is deleted or foreign.
    pub fn delete_atom(&mut self, id: EditAtomId) -> Result<Atom, TopologyEditError> {
        let loc = self.atom_location(id)?;
        let bonds = self.incident_bonds(id)?.collect::<Vec<_>>();
        let neighbors = self.neighbors(id)?.collect::<Vec<_>>();
        let previous = self.draft_mut(loc.group).delete_atom(loc.local)?;
        for bond in bonds {
            let source = self.source.as_deref();
            let bond_loc = self
                .bonds
                .remove(&bond, |index| source_bond(source, index))
                .unwrap();
            self.group_mut(loc.group).bonds.remove(&bond_loc.local);
            self.properties.bonds_mut().clear_index(bond_loc.slot);
            self.deleted_slots.1.push(bond_loc.slot);
        }
        let source = self.source.as_deref();
        self.atoms.remove(&id, |index| source_atom(source, index));
        self.dead_atom_slots.insert(loc.slot);
        self.deleted_slots.0.push(loc.slot);
        self.group_mut(loc.group).atoms.remove(&loc.local);
        self.properties.atoms_mut().clear_index(loc.slot);
        self.remove_site_for_atom(id);
        self.mark_group_changed(loc.group);
        self.molecule_classes.remove(&id);
        self.invalidate_component_classes(neighbors);
        self.changed();
        Ok(previous)
    }

    /// Deletes several atoms, returning them in handle order.
    ///
    /// Validates every handle before deleting any atom; duplicates are ignored.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when any handle is deleted or foreign;
    /// nothing is deleted then.
    pub fn delete_atoms(
        &mut self,
        ids: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<Vec<(EditAtomId, Atom)>, TopologyEditError> {
        let ids = ids.into_iter().collect::<BTreeSet<_>>();
        for &id in &ids {
            self.atom_location(id)?;
        }
        ids.into_iter()
            .map(|id| self.delete_atom(id).map(|a| (id, a)))
            .collect()
    }
    /// Deletes every live atom not in `ids`, returning the deleted atoms in
    /// handle order.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when any retained handle is deleted or
    /// foreign; nothing is deleted then.
    pub fn retain_atoms(
        &mut self,
        ids: impl IntoIterator<Item = EditAtomId>,
    ) -> Result<Vec<(EditAtomId, Atom)>, TopologyEditError> {
        let retained = ids.into_iter().collect::<BTreeSet<_>>();
        for &id in &retained {
            self.atom_location(id)?;
        }
        let removed = self
            .atom_ids()
            .filter(|id| !retained.contains(id))
            .collect::<Vec<_>>();
        self.delete_atoms(removed)
    }
    /// Deletes surviving atoms that belonged to this source occurrence.
    /// After merges, atoms from other source occurrences remain. To delete a
    /// current component, pass its handles from [`Self::connected_components`]
    /// to [`Self::delete_atoms`].
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidSourceInstance`] when `source` is not an
    /// instance of the source topology.
    pub fn delete_instance(
        &mut self,
        source: MoleculeInstanceId,
    ) -> Result<Vec<(EditAtomId, Atom)>, TopologyEditError> {
        let Some(instance) = self.source.as_ref().and_then(|t| t.instance(source)) else {
            return Err(TopologyEditError::InvalidSourceInstance(source));
        };
        let ids = instance
            .atoms()
            .filter_map(|atom| self.atom_handle(atom.id()))
            .collect::<Vec<_>>();
        self.delete_atoms(ids)
    }

    /// Adds an asserted bond. Different occurrences are merged transactionally.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when either handle is deleted or
    /// foreign, or [`TopologyEditError::Molecule`] for a self bond, an
    /// existing bond, or conflicting fragment properties. The draft is
    /// unchanged on error.
    pub fn add_bond(
        &mut self,
        a: EditAtomId,
        b: EditAtomId,
        order: BondOrder,
    ) -> Result<EditBondId, TopologyEditError> {
        let left = self.atom_location(a)?;
        let right = self.atom_location(b)?;
        if left.group == right.group {
            if a == b {
                return Err(MoleculeError::SelfBond(left.local).into());
            }
            if self.bond_between(a, b)?.is_some() {
                return Err(MoleculeError::DuplicateBond {
                    a: left.local,
                    b: right.local,
                }
                .into());
            }
            let local = self
                .draft_mut(left.group)
                .add_bond(left.local, right.local, order)?;
            return Ok(self.register_bond(left.group, local, a, b));
        }
        // Stage only the two affected chemical groups. Conflicting fragment
        // properties fail before any handle, hierarchy, or coordinate changes.
        let target = left.group.min(right.group);
        let other = left.group.max(right.group);
        let mut draft = self.molecule(target).to_editor();
        let mapping = draft.append_working(self.molecule(other))?;
        let local_a = if left.group == target {
            left.local
        } else {
            mapping.atoms()[&left.local]
        };
        let local_b = if right.group == target {
            right.local
        } else {
            mapping.atoms()[&right.local]
        };
        let local = draft.add_bond(local_a, local_b, order)?;
        self.commit_merge(target, other, draft, &mapping);
        Ok(self.register_bond(target, local, a, b))
    }

    /// Deletes one bond, returning it. Deleting a bridge splits its molecule at
    /// publication.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidBond`] when `id` is deleted or foreign.
    pub fn delete_bond(&mut self, id: EditBondId) -> Result<EditBond, TopologyEditError> {
        let previous = self.require_bond(id)?;
        let loc = self.bond_location(id)?;
        self.draft_mut(loc.group).delete_bond(loc.local)?;
        let source = self.source.as_deref();
        self.bonds.remove(&id, |index| source_bond(source, index));
        self.group_mut(loc.group).bonds.remove(&loc.local);
        self.properties.bonds_mut().clear_index(loc.slot);
        self.deleted_slots.1.push(loc.slot);
        self.mark_group_changed(loc.group);
        self.invalidate_component_classes([previous.a, previous.b]);
        self.invalidate_residue_for_atom(previous.a);
        self.invalidate_residue_for_atom(previous.b);
        self.changed();
        Ok(previous)
    }
    /// Deletes several bonds, returning them in handle order.
    ///
    /// Validates every handle first; duplicates are ignored.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidBond`] when any handle is deleted or foreign;
    /// nothing is deleted then.
    pub fn delete_bonds(
        &mut self,
        ids: impl IntoIterator<Item = EditBondId>,
    ) -> Result<Vec<(EditBondId, EditBond)>, TopologyEditError> {
        let ids = ids.into_iter().collect::<BTreeSet<_>>();
        for &id in &ids {
            self.bond_location(id)?;
        }
        ids.into_iter()
            .map(|id| self.delete_bond(id).map(|b| (id, b)))
            .collect()
    }
    /// Changes represented order, removing stereo assertions focused on this bond.
    /// Assigning the current order leaves the draft unchanged.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidBond`] when `id` is deleted or foreign, or
    /// [`TopologyEditError::Molecule`] when the molecular draft rejects the order.
    pub fn set_bond_order(
        &mut self,
        id: EditBondId,
        order: BondOrder,
    ) -> Result<(), TopologyEditError> {
        let previous = self.require_bond(id)?;
        if previous.order == order {
            return Ok(());
        }
        let loc = self.bond_location(id)?;
        self.draft_mut(loc.group).set_bond_order(loc.local, order)?;
        self.mark_group_changed(loc.group);
        self.invalidate_component_classes([previous.a]);
        self.invalidate_residue_for_atom(previous.a);
        self.invalidate_residue_for_atom(previous.b);
        self.changed();
        Ok(())
    }

    /// Rewires a bond while retaining its editing handle and static annotations.
    /// The transaction may merge groups and split the old component at publication.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidBond`] or [`TopologyEditError::InvalidAtom`]
    /// for deleted or foreign handles, or [`TopologyEditError::Molecule`] for a
    /// self bond or an existing bond between the new endpoints. The draft is
    /// unchanged on error.
    pub fn set_bond_endpoints(
        &mut self,
        id: EditBondId,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<(), TopologyEditError> {
        let old = self.require_bond(id)?;
        self.atom_location(a)?;
        self.atom_location(b)?;
        if (old.a == a && old.b == b) || (old.a == b && old.b == a) {
            return Ok(());
        }
        let mut staged = self.clone();
        staged.set_bond_endpoints_staged(id, a, b)?;
        *self = staged;
        Ok(())
    }

    // The caller owns rollback; group merging may precede a failed rewire.
    fn set_bond_endpoints_staged(
        &mut self,
        id: EditBondId,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<(), TopologyEditError> {
        let old = self.require_bond(id)?;
        self.atom_location(a)?;
        self.atom_location(b)?;
        if (old.a == a && old.b == b) || (old.a == b && old.b == a) {
            return Ok(());
        }
        let old_loc = self.bond_location(id)?;
        // A checked replacement on the combined draft preserves definition-level
        // bond annotations too. Combine groups before rewiring, without new bonds.
        for atom in [a, b] {
            let bond_group = self.bond_location(id)?.group;
            let atom_group = self.atom_location(atom)?.group;
            if bond_group != atom_group {
                let target = bond_group.min(atom_group);
                let other = bond_group.max(atom_group);
                let mut draft = self.molecule(target).to_editor();
                let map = draft.append_working(self.molecule(other))?;
                self.commit_merge(target, other, draft, &map);
            }
        }
        let loc = self.bond_location(id)?;
        let aa = self.atom_location(a)?.local;
        let bb = self.atom_location(b)?.local;
        self.draft_mut(loc.group)
            .set_bond_endpoints(loc.local, aa, bb)?;
        debug_assert_eq!(old_loc.slot, loc.slot);
        self.mark_group_changed(loc.group);
        self.invalidate_component_classes([old.a, old.b, a, b]);
        for atom in [old.a, old.b, a, b] {
            self.invalidate_residue_for_atom(atom);
        }
        self.changed();
        Ok(())
    }

    /// Replaces one bond's endpoints and order, returning the previous bond.
    ///
    /// Transactional: combines [`Self::set_bond_endpoints`] and
    /// [`Self::set_bond_order`].
    ///
    /// # Errors
    ///
    /// As for [`Self::set_bond_endpoints`] and [`Self::set_bond_order`]; the draft
    /// is unchanged on error.
    pub fn replace_bond(
        &mut self,
        id: EditBondId,
        replacement: EditBond,
    ) -> Result<EditBond, TopologyEditError> {
        let previous = self.require_bond(id)?;
        let mut staged = self.clone();
        staged.set_bond_endpoints_staged(id, replacement.a, replacement.b)?;
        staged.set_bond_order(id, replacement.order)?;
        *self = staged;
        Ok(previous)
    }

    /// Explicitly assigns the current connected component's class. Assigning its
    /// existing inferred class still records an override for subsequent metadata
    /// edits. A later chemical change invalidates this component's override.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when `atom` is deleted or foreign.
    pub fn set_molecule_class(
        &mut self,
        atom: EditAtomId,
        class: MoleculeClass,
    ) -> Result<(), TopologyEditError> {
        let loc = self.atom_location(atom)?;
        let local_component = self
            .molecule(loc.group)
            .connected_components()
            .into_iter()
            .find(|atoms| atoms.contains(&loc.local))
            .expect("live component");
        let component = local_component
            .into_iter()
            .map(|id| self.group_atom(loc.group, id))
            .collect::<BTreeSet<_>>();
        let (group_class, group_explicit) = self.group_class(loc.group);
        if component.len() == self.group_atom_count(loc.group)
            && (group_explicit
                || self
                    .molecule_classes
                    .keys()
                    .any(|id| component.contains(id)))
            && self.component_class(&component, group_class) == Some(class)
        {
            return Ok(());
        }
        // Class overrides are occurrence-local too: detach a reused definition.
        self.draft_mut(loc.group);
        self.molecule_classes
            .retain(|id, _| !component.contains(id));
        self.molecule_classes.insert(atom, class);
        self.revision += 1;
        Ok(())
    }

    /// Replaces stereo on one intact source occurrence. Elements use that
    /// occurrence's molecule-local IDs. Call on a fresh editor before structural
    /// edits to this occurrence; other occurrences may still be edited normally.
    /// Existing enhanced stereo groups are deliberately removed. Atomic on error.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidSourceInstance`] when the draft has no source,
    /// `instance` is not a source instance, or the occurrence was already
    /// edited; [`TopologyEditError::InvalidSourceAtom`] when its atoms were
    /// deleted; [`TopologyEditError::Molecule`] or
    /// [`TopologyEditError::Publication`] when the stereo elements are invalid.
    pub fn replace_source_instance_stereo(
        &mut self,
        instance: MoleculeInstanceId,
        elements: &[crate::stereo::StereoElement],
    ) -> Result<(), TopologyEditError> {
        let source = self
            .source
            .as_ref()
            .ok_or(TopologyEditError::InvalidSourceInstance(instance))?;
        let molecule = source
            .instance(instance)
            .ok_or(TopologyEditError::InvalidSourceInstance(instance))?;
        let atom = molecule.atoms().next().unwrap().id();
        let handle = self.require_atom_handle(atom)?;
        let group = self.atom_location(handle)?.group;
        if self.group_changed(group)
            || self.group_atom_count(group) != molecule.molecule().atom_count()
        {
            return Err(TopologyEditError::InvalidSourceInstance(instance));
        }
        let mut draft = self.molecule(group).to_editor();
        let groups = draft.stereo_groups().map(|(id, _)| id).collect::<Vec<_>>();
        for id in groups {
            draft.remove_stereo_group(id)?;
        }
        let ids = draft.stereo_element_ids().collect::<Vec<_>>();
        for id in ids {
            draft.remove_stereo_element(id)?;
        }
        for element in elements {
            draft.add_stereo_element(element.clone())?;
        }
        draft.clone().finish()?;
        self.group_mut(group).chemistry = GroupChemistry::Draft(Box::new(draft));
        self.mark_group_changed(group);
        self.changed();
        Ok(())
    }

    fn component_class(
        &self,
        atoms: &BTreeSet<EditAtomId>,
        fallback: Option<MoleculeClass>,
    ) -> Option<MoleculeClass> {
        self.molecule_classes
            .iter()
            .find_map(|(id, &class)| atoms.contains(id).then_some(class))
            .or(fallback)
    }

    pub(crate) fn atom_slot(&self, id: EditAtomId) -> Result<usize, TopologyEditError> {
        Ok(self.atom_location(id)?.slot)
    }
    pub(crate) fn bond_slot(&self, id: EditBondId) -> Result<usize, TopologyEditError> {
        Ok(self.bond_location(id)?.slot)
    }
    /// How many atom and bond deletions this draft has recorded, a mark for
    /// [`Self::deleted_slots_since`].
    pub(crate) fn deletion_mark(&self) -> (usize, usize) {
        (self.deleted_slots.0.len(), self.deleted_slots.1.len())
    }
    /// The slots of the atoms and bonds deleted since `mark`.
    pub(crate) fn deleted_slots_since(&self, mark: (usize, usize)) -> (&[usize], &[usize]) {
        (
            self.deleted_slots.0.get(mark.0..).unwrap_or_default(),
            self.deleted_slots.1.get(mark.1..).unwrap_or_default(),
        )
    }
    /// The slots of the live atoms in handle order, which is publication order.
    pub(crate) fn live_atom_slots(&self) -> Vec<usize> {
        self.atom_locations().map(|(_, loc)| loc.slot).collect()
    }
    pub(crate) fn atom_slot_count(&self) -> usize {
        self.properties.atoms().len()
    }
    pub(crate) fn bond_slot_count(&self) -> usize {
        self.properties.bonds().len()
    }
    fn atom_location(&self, id: EditAtomId) -> Result<Location<AtomId>, TopologyEditError> {
        let source = self.source.as_deref();
        self.atoms
            .get(&id, |index| source_atom(source, index))
            .ok_or(TopologyEditError::InvalidAtom(id))
    }
    fn bond_location(&self, id: EditBondId) -> Result<Location<BondId>, TopologyEditError> {
        let source = self.source.as_deref();
        self.bonds
            .get(&id, |index| source_bond(source, index))
            .ok_or(TopologyEditError::InvalidBond(id))
    }
    /// The live atom locations in handle order.
    fn atom_locations(&self) -> impl ExactSizeIterator<Item = (EditAtomId, Location<AtomId>)> + '_ {
        let source = self.source.as_deref();
        self.atoms.iter(move |index| source_atom(source, index))
    }
    fn source_ref(&self) -> &Topology {
        self.source.as_deref().expect("source groups have a source")
    }
    fn group_is_live(&self, group: usize) -> bool {
        !matches!(self.groups[group], GroupSlot::Merged)
    }
    /// The stored group, copying an untouched source occurrence's handles
    /// into storage first.
    fn group_mut(&mut self, group: usize) -> &mut Group {
        if matches!(self.groups[group], GroupSlot::Source) {
            let stored = self.source_group(group);
            self.groups[group] = GroupSlot::Live(Box::new(stored));
        }
        match &mut self.groups[group] {
            GroupSlot::Live(group) => group,
            _ => unreachable!("live groups are stored"),
        }
    }
    /// The stored form of an untouched source occurrence.
    fn source_group(&self, group: usize) -> Group {
        let (class, class_explicit) = self.group_class(group);
        Group {
            chemistry: GroupChemistry::Source(self.source_definition(group)),
            atoms: self.group_atoms(group).into_iter().collect(),
            bonds: self.group_bonds(group).into_iter().collect(),
            instance_slot: Some(group),
            changed: false,
            class,
            class_explicit,
        }
    }
    /// The source definition of an untouched source occurrence.
    fn source_definition(&self, group: usize) -> MoleculeDefinitionId {
        self.source_ref()
            .instance_record(MoleculeInstanceId::new(group as u32))
            .expect("source group")
            .definition()
    }
    /// The handle of a group's local atom.
    fn group_atom(&self, group: usize, local: AtomId) -> EditAtomId {
        match &self.groups[group] {
            GroupSlot::Source => {
                let instance = MoleculeInstanceId::new(group as u32);
                let index = self
                    .source_ref()
                    .atom_index(InstanceAtomId::new(instance, local))
                    .expect("source group atom");
                self.atoms.block().handle(index.index())
            }
            GroupSlot::Live(stored) => stored.atoms[&local],
            GroupSlot::Merged => unreachable!("merged groups have no atoms"),
        }
    }
    /// The handle of a group's local bond.
    fn group_bond(&self, group: usize, local: BondId) -> EditBondId {
        match &self.groups[group] {
            GroupSlot::Source => {
                let instance = MoleculeInstanceId::new(group as u32);
                let index = self
                    .source_ref()
                    .bond_index(InstanceBondId::new(instance, local))
                    .expect("source group bond");
                self.bonds.block().handle(index.index())
            }
            GroupSlot::Live(stored) => stored.bonds[&local],
            GroupSlot::Merged => unreachable!("merged groups have no bonds"),
        }
    }
    /// A live group's atoms with their handles, in local ID order.
    fn group_atoms(&self, group: usize) -> Vec<(AtomId, EditAtomId)> {
        match &self.groups[group] {
            GroupSlot::Live(stored) => stored.atoms.iter().map(|(&l, &h)| (l, h)).collect(),
            _ => self
                .molecule(group)
                .atom_ids()
                .map(|local| (local, self.group_atom(group, local)))
                .collect(),
        }
    }
    /// A live group's bonds with their handles, in local ID order.
    fn group_bonds(&self, group: usize) -> Vec<(BondId, EditBondId)> {
        match &self.groups[group] {
            GroupSlot::Live(stored) => stored.bonds.iter().map(|(&l, &h)| (l, h)).collect(),
            _ => self
                .molecule(group)
                .bond_ids()
                .map(|local| (local, self.group_bond(group, local)))
                .collect(),
        }
    }
    /// The number of live atoms of a live group.
    fn group_atom_count(&self, group: usize) -> usize {
        match &self.groups[group] {
            GroupSlot::Live(stored) => stored.atoms.len(),
            _ => self.molecule(group).atom_count(),
        }
    }
    /// Whether a chemical edit changed a live group.
    fn group_changed(&self, group: usize) -> bool {
        match &self.groups[group] {
            GroupSlot::Live(stored) => stored.changed,
            _ => false,
        }
    }
    /// A live group's retained class and whether it was assigned explicitly.
    fn group_class(&self, group: usize) -> (Option<MoleculeClass>, bool) {
        match &self.groups[group] {
            GroupSlot::Live(stored) => (stored.class, stored.class_explicit),
            _ => {
                let source = self.source_ref();
                let definition = self.source_definition(group);
                (
                    Some(source.definitions[definition.index()].class()),
                    source.molecule_class_overrides.contains_key(&definition),
                )
            }
        }
    }
    fn molecule(&self, group: usize) -> &Molecule {
        match &self.groups[group] {
            GroupSlot::Source => self
                .source_ref()
                .definition_molecule(MoleculeInstanceId::new(group as u32)),
            GroupSlot::Live(stored) => match &stored.chemistry {
                GroupChemistry::Source(id) => self.source_ref().definition(*id).unwrap().molecule(),
                GroupChemistry::Added(molecule, _) => molecule,
                GroupChemistry::Draft(editor) => editor.working(),
            },
            GroupSlot::Merged => unreachable!("merged groups have no molecule"),
        }
    }
    fn draft_mut(&mut self, group: usize) -> &mut MoleculeEditor {
        if !matches!(self.group_mut(group).chemistry, GroupChemistry::Draft(_)) {
            let draft = self.molecule(group).to_editor();
            self.group_mut(group).chemistry = GroupChemistry::Draft(Box::new(draft));
        }
        match &mut self.group_mut(group).chemistry {
            GroupChemistry::Draft(draft) => draft,
            _ => unreachable!(),
        }
    }
    fn mark_group_changed(&mut self, group: usize) {
        let group = self.group_mut(group);
        group.changed = true;
        group.class = None;
        group.class_explicit = false;
    }

    // Historical staging groups can contain several disconnected components.
    // Invalidate overrides only in the components touched by a successful edit.
    // After a split, callers supply seeds on each side; after an atom deletion,
    // they supply the deleted atom's surviving neighbors.
    fn invalidate_component_classes(&mut self, atoms: impl IntoIterator<Item = EditAtomId>) {
        if self.molecule_classes.is_empty() {
            return;
        }
        let mut visited = BTreeSet::new();
        let mut pending = atoms.into_iter().collect::<Vec<_>>();
        while let Some(atom) = pending.pop() {
            if !visited.insert(atom) {
                continue;
            }
            self.molecule_classes.remove(&atom);
            if self.molecule_classes.is_empty() {
                return;
            }
            pending.extend(
                self.neighbors(atom)
                    .expect("affected component contains live atoms")
                    .filter(|neighbor| !visited.contains(neighbor)),
            );
        }
    }
    fn changed(&mut self) {
        self.revision += 1;
        self.structural_revision += 1;
        self.properties.owner_mut().clear();
    }
    /// Stores an added group whose atoms and bonds take the given handles and
    /// slots, and returns its index.
    fn register_group(
        &mut self,
        chemistry: GroupChemistry,
        molecule: &Molecule,
        class: Option<MoleculeClass>,
        class_explicit: bool,
        identity: GroupIdentity,
    ) -> usize {
        let group = self.groups.len();
        let mut atoms = DenseMap::new();
        let mut bonds = DenseMap::new();
        for (local, (id, slot)) in molecule.atom_ids().zip(identity.atoms) {
            self.atoms.insert(id, Location { group, local, slot });
            atoms.insert(local, id);
        }
        for (local, (id, slot)) in molecule.bond_ids().zip(identity.bonds) {
            self.bonds.insert(id, Location { group, local, slot });
            bonds.insert(local, id);
        }
        self.groups.push(GroupSlot::Live(Box::new(Group {
            chemistry,
            atoms,
            bonds,
            instance_slot: None,
            changed: false,
            class,
            class_explicit,
        })));
        group
    }
    /// The handles of one group's atoms and bonds by local ID.
    fn group_handles(&self, group: usize) -> EditMolecule {
        EditMolecule {
            atoms: self.group_atoms(group).into_iter().collect(),
            bonds: self.group_bonds(group).into_iter().collect(),
        }
    }
    fn register_bond(
        &mut self,
        group: usize,
        local: BondId,
        a: EditAtomId,
        b: EditAtomId,
    ) -> EditBondId {
        let id = EditBondId::new();
        let slot = self.properties.bonds().len();
        self.properties.resize_bonds(slot + 1);
        self.bonds.insert(id, Location { group, local, slot });
        self.group_mut(group).bonds.insert(local, id);
        self.mark_group_changed(group);
        self.invalidate_component_classes([a]);
        self.invalidate_residue_for_atom(a);
        self.invalidate_residue_for_atom(b);
        self.changed();
        id
    }
    fn commit_merge(
        &mut self,
        target: usize,
        other: usize,
        draft: MoleculeEditor,
        map: &crate::molecule::MoleculeCorrespondence,
    ) {
        self.group_mut(other);
        let GroupSlot::Live(removed) =
            std::mem::replace(&mut self.groups[other], GroupSlot::Merged)
        else {
            unreachable!("merged groups are stored");
        };
        self.group_mut(target).chemistry = GroupChemistry::Draft(Box::new(draft));
        let source = self.source.as_deref();
        let GroupSlot::Live(group) = &mut self.groups[target] else {
            unreachable!("the merge target is stored");
        };
        for (old, handle) in removed.atoms {
            let local = map.atoms()[&old];
            group.atoms.insert(local, handle);
            let loc = self
                .atoms
                .get_mut(&handle, |index| source_atom(source, index))
                .unwrap();
            loc.group = target;
            loc.local = local;
        }
        for (old, handle) in removed.bonds {
            let local = map.bonds()[&old];
            group.bonds.insert(local, handle);
            let loc = self
                .bonds
                .get_mut(&handle, |index| source_bond(source, index))
                .unwrap();
            loc.group = target;
            loc.local = local;
        }
        self.mark_group_changed(target);
    }
}

/// Borrowed chemistry inspection translated to stable editing endpoints.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EditBond {
    a: EditAtomId,
    b: EditAtomId,
    /// The represented bond order.
    pub order: BondOrder,
}
impl EditBond {
    /// A bond between handles `a` and `b`. Nothing is validated until it is applied.
    pub fn new(a: EditAtomId, b: EditAtomId, order: BondOrder) -> Self {
        Self { a, b, order }
    }
    /// The first endpoint handle.
    pub fn a(self) -> EditAtomId {
        self.a
    }
    /// The second endpoint handle.
    pub fn b(self) -> EditAtomId {
        self.b
    }
    /// Both endpoint handles as `(a, b)`.
    pub fn endpoints(self) -> (EditAtomId, EditAtomId) {
        (self.a, self.b)
    }
}

/// Definition-local IDs translated to newly allocated occurrence editing handles.
#[derive(Debug, Clone, Default)]
pub struct EditMolecule {
    atoms: BTreeMap<AtomId, EditAtomId>,
    bonds: BTreeMap<BondId, EditBondId>,
}
impl EditMolecule {
    /// New atom handles keyed by the added molecule's local atom IDs.
    pub fn atoms(&self) -> &BTreeMap<AtomId, EditAtomId> {
        &self.atoms
    }
    /// New bond handles keyed by the added molecule's local bond IDs.
    pub fn bonds(&self) -> &BTreeMap<BondId, EditBondId> {
        &self.bonds
    }
}
