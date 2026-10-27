//! Transactional system edits, with molecular chemistry staged per affected occurrence.
mod append;
mod error;
mod hierarchy;
mod identity;
mod placement;
mod properties;
mod publication;
pub use publication::{DenseCorrespondence, EditCorrespondence};

use super::dense_map::{DenseKey, DenseMap, QualifiedMap};
use super::{
    AtomSiteId, AtomSiteMetadata, ChainId, HierarchyError, InstanceAtomId, InstanceBondId,
    MoleculeClass, MoleculeDefinitionId, MoleculeInstanceId, ResidueClass, ResidueId, Topology,
    TopologyBuildError, TopologyBuilder,
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
use hierarchy::EditHierarchy;
pub use hierarchy::{EditAtomSite, EditChain, EditResidue};
pub use identity::{EditAtomId, EditAtomSiteId, EditBondId, EditChainId, EditResidueId};
pub(crate) use placement::{AttachedInOccurrence, DraftHydrogenError};
pub use properties::{EditProperties, EditPropertiesMut, EditPropertyTable, EditPropertyTableMut};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Arc;

#[derive(Debug, Clone)]
enum GroupChemistry {
    Source(MoleculeDefinitionId),
    Added(Arc<Molecule>),
    Draft(Box<MoleculeEditor>),
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

    /// One instance of a dense source layout, with handles allocated in
    /// source dense order from the given first handles and slots offset by the
    /// given bases.
    fn dense(
        source: &Topology,
        instance: MoleculeInstanceId,
        handles: (u64, u64),
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
            (
                EditAtomId::from_position(handles.0 + dense as u64),
                bases.0 + dense,
            )
        };
        let bond = |local| {
            let dense = source
                .bond_index(InstanceBondId::new(instance, local))
                .expect("source instance bond has a dense index")
                .index();
            (
                EditBondId::from_position(handles.1 + dense as u64),
                bases.1 + dense,
            )
        };
        Self {
            atoms: molecule.atom_ids().map(atom).collect(),
            bonds: molecule.bond_ids().map(bond).collect(),
        }
    }

    /// Reserves one consecutive handle per dense source atom and bond, in
    /// dense order, and returns the first of each.
    fn source_handles(source: &Topology) -> (u64, u64) {
        (
            identity::fresh_block(source.atom_count()),
            identity::fresh_block(source.bond_count()),
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
    groups: Vec<Option<Group>>,
    atoms: DenseMap<EditAtomId, Location<AtomId>>,
    bonds: DenseMap<EditBondId, Location<BondId>>,
    source_atoms: QualifiedMap<InstanceAtomId, EditAtomId>,
    source_bonds: QualifiedMap<InstanceBondId, EditBondId>,
    molecule_classes: BTreeMap<EditAtomId, MoleculeClass>,
    hierarchy: EditHierarchy,
    properties: TopologyProperties,
    revision: u64,
    pub(crate) structural_revision: u64,
    // Storage slots of atoms deleted from this draft, so a slot's liveness is
    // known without walking every live atom.
    dead_atom_slots: BTreeSet<usize>,
    // The handles of source entities: consecutive blocks in source order.
    source_blocks: SourceBlocks,
}

/// One block of consecutive handles, one per source entity in source order.
#[derive(Debug, Clone, Copy, Default)]
struct HandleBlock {
    base: u64,
    len: usize,
}

impl HandleBlock {
    /// The source index of a handle in this block.
    fn source(self, handle: impl DenseKey) -> Option<u32> {
        let offset = handle.position().checked_sub(self.base)?;
        (offset < self.len as u64).then_some(offset as u32)
    }
}

#[derive(Debug, Clone, Copy, Default)]
struct SourceBlocks {
    atoms: HandleBlock,
    bonds: HandleBlock,
    chains: HandleBlock,
    residues: HandleBlock,
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
        self.hierarchy.chains.clear();
        self.hierarchy.residues.clear();
        self.hierarchy.sites.clear();
        self.hierarchy.atom_sites.clear();
        self.hierarchy.residue_atoms.clear();
        self.properties = TopologyProperties::default();
        self.dead_atom_slots.clear();
        self.changed();
    }

    /// Retains the supplied immutable source snapshot without cloning chemistry.
    pub fn from_topology(topology: impl Into<Arc<Topology>>) -> Self {
        let source = topology.into();
        let mut editor = Self {
            properties: source.properties().clone(),
            source: Some(Arc::clone(&source)),
            ..Self::default()
        };
        // Source slots are source dense indices, so an unedited draft
        // republishes the source order exactly. Source handles are one block
        // in dense order, so vectors index them.
        let handles = GroupIdentity::source_handles(&source);
        editor.atoms = DenseMap::block(handles.0, source.atom_count());
        editor.bonds = DenseMap::block(handles.1, source.bond_count());
        editor.source_blocks.atoms = HandleBlock {
            base: handles.0,
            len: source.atom_count(),
        };
        editor.source_blocks.bonds = HandleBlock {
            base: handles.1,
            len: source.bond_count(),
        };
        for (instance, value) in source.instance_records() {
            let molecule = source
                .definition(value.definition())
                .expect("published definition")
                .molecule();
            let identity = GroupIdentity::dense(&source, instance, handles, (0, 0));
            editor.register_group(
                GroupChemistry::Source(value.definition()),
                molecule,
                Some(instance),
                Some(source.definition(value.definition()).unwrap().class()),
                source
                    .molecule_class_overrides
                    .contains_key(&value.definition()),
                identity,
            );
        }
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
        self.atoms.keys().copied()
    }
    /// Live bond handles in handle order.
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = EditBondId> + '_ {
        self.bonds.keys().copied()
    }
    /// Live atoms with their handles, in publication order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (EditAtomId, &Atom)> + '_ {
        self.atoms.iter().map(|(&id, loc)| {
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
        self.source_atoms
            .get(&source)
            .copied()
            .filter(|id| self.atoms.contains_key(id))
    }
    /// The draft handle of a source bond, or `None` when the bond is not in the
    /// source or was deleted.
    pub fn bond_handle(&self, source: InstanceBondId) -> Option<EditBondId> {
        self.source_bonds
            .get(&source)
            .copied()
            .filter(|id| self.bonds.contains_key(id))
    }
    /// One draft atom, or `None` when `id` names no live atom of this draft.
    pub fn atom(&self, id: EditAtomId) -> Option<&Atom> {
        let loc = self.atoms.get(&id)?;
        self.molecule(loc.group).atom(loc.local)
    }
    /// One draft bond with handle endpoints, or `None` when `id` names no live
    /// bond of this draft.
    pub fn bond(&self, id: EditBondId) -> Option<EditBond> {
        let loc = self.bonds.get(&id)?;
        let group = self.groups[loc.group]
            .as_ref()
            .expect("bond locations name live groups");
        let bond = self.molecule(loc.group).bond(loc.local)?;
        Some(EditBond {
            a: group.atoms[&bond.a()],
            b: group.atoms[&bond.b()],
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
        let group = self.groups[loc.group].as_ref().unwrap();
        Ok(self
            .molecule(loc.group)
            .require_neighbors(loc.local)?
            .map(|id| group.atoms[&id]))
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
        let group = self.groups[loc.group].as_ref().unwrap();
        Ok(self
            .molecule(loc.group)
            .require_incident_bonds(loc.local)?
            .map(|(id, _)| group.bonds[&id]))
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
            .map(|id| self.groups[left.group].as_ref().unwrap().bonds[&id]))
    }
    /// The current connected components as atom-handle lists.
    ///
    /// Each becomes one molecule instance at publication.
    pub fn connected_components(&self) -> Vec<Vec<EditAtomId>> {
        self.groups
            .iter()
            .enumerate()
            .filter_map(|(i, group)| group.as_ref().map(|g| (i, g)))
            .flat_map(|(i, g)| {
                self.molecule(i)
                    .connected_components()
                    .into_iter()
                    .map(|atoms| atoms.into_iter().map(|a| g.atoms[&a]).collect())
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
        self.groups.push(Some(Group {
            chemistry: GroupChemistry::Draft(Box::new(draft)),
            atoms: DenseMap::from_iter([(local, id)]),
            bonds: DenseMap::new(),
            instance_slot: None,
            changed: true,
            class: None,
            class_explicit: false,
        }));
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
        let identity = GroupIdentity::appended(
            &owned,
            self.properties.atoms().len(),
            self.properties.bonds().len(),
        );
        let group = self.register_group(
            GroupChemistry::Added(Arc::clone(&owned)),
            &owned,
            None,
            None,
            false,
            identity,
        );
        let result = self.group_handles(group);
        self.properties.resize_atoms(
            self.atoms
                .values()
                .map(|l| l.slot + 1)
                .max()
                .unwrap_or(0)
                .max(self.properties.atoms().len()),
        );
        self.properties.resize_bonds(
            self.bonds
                .values()
                .map(|l| l.slot + 1)
                .max()
                .unwrap_or(0)
                .max(self.properties.bonds().len()),
        );
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
            let bond_loc = self.bonds.remove(&bond).unwrap();
            self.groups[loc.group]
                .as_mut()
                .unwrap()
                .bonds
                .remove(&bond_loc.local);
            self.properties.bonds_mut().clear_index(bond_loc.slot);
        }
        self.atoms.remove(&id);
        self.dead_atom_slots.insert(loc.slot);
        self.groups[loc.group]
            .as_mut()
            .unwrap()
            .atoms
            .remove(&loc.local);
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
        if self
            .source
            .as_ref()
            .is_none_or(|t| t.instance_record(source).is_none())
        {
            return Err(TopologyEditError::InvalidSourceInstance(source));
        }
        let ids = self
            .source_atoms
            .instance(source.index())
            .filter(|(_, handle)| self.atoms.contains_key(handle))
            .map(|(_, &h)| h)
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
        self.bonds.remove(&id);
        self.groups[loc.group]
            .as_mut()
            .unwrap()
            .bonds
            .remove(&loc.local);
        self.properties.bonds_mut().clear_index(loc.slot);
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
        let group = self.groups[loc.group].as_ref().unwrap();
        let component = local_component
            .into_iter()
            .map(|id| group.atoms[&id])
            .collect::<BTreeSet<_>>();
        if component.len() == group.atoms.len()
            && (group.class_explicit
                || self
                    .molecule_classes
                    .keys()
                    .any(|id| component.contains(id)))
            && self.component_class(&component, group.class) == Some(class)
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
        let stored = self.groups[group].as_ref().unwrap();
        if stored.changed || stored.atoms.len() != molecule.molecule().atom_count() {
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
        self.groups[group].as_mut().unwrap().chemistry = GroupChemistry::Draft(Box::new(draft));
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
    pub(crate) fn atom_slot_count(&self) -> usize {
        self.properties.atoms().len()
    }
    pub(crate) fn bond_slot_count(&self) -> usize {
        self.properties.bonds().len()
    }
    fn atom_location(&self, id: EditAtomId) -> Result<Location<AtomId>, TopologyEditError> {
        self.atoms
            .get(&id)
            .copied()
            .ok_or(TopologyEditError::InvalidAtom(id))
    }
    fn bond_location(&self, id: EditBondId) -> Result<Location<BondId>, TopologyEditError> {
        self.bonds
            .get(&id)
            .copied()
            .ok_or(TopologyEditError::InvalidBond(id))
    }
    fn molecule(&self, group: usize) -> &Molecule {
        match &self.groups[group].as_ref().unwrap().chemistry {
            GroupChemistry::Source(id) => self
                .source
                .as_ref()
                .unwrap()
                .definition(*id)
                .unwrap()
                .molecule(),
            GroupChemistry::Added(molecule) => molecule,
            GroupChemistry::Draft(editor) => editor.working(),
        }
    }
    fn draft_mut(&mut self, group: usize) -> &mut MoleculeEditor {
        if !matches!(
            self.groups[group].as_ref().unwrap().chemistry,
            GroupChemistry::Draft(_)
        ) {
            let draft = self.molecule(group).to_editor();
            self.groups[group].as_mut().unwrap().chemistry = GroupChemistry::Draft(Box::new(draft));
        }
        match &mut self.groups[group].as_mut().unwrap().chemistry {
            GroupChemistry::Draft(draft) => draft,
            _ => unreachable!(),
        }
    }
    fn mark_group_changed(&mut self, group: usize) {
        let group = self.groups[group].as_mut().unwrap();
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
    fn register_group(
        &mut self,
        chemistry: GroupChemistry,
        molecule: &Molecule,
        source: Option<MoleculeInstanceId>,
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
            if let Some(source) = source {
                self.source_atoms
                    .insert(InstanceAtomId::new(source, local), id);
            }
        }
        for (local, (id, slot)) in molecule.bond_ids().zip(identity.bonds) {
            self.bonds.insert(id, Location { group, local, slot });
            bonds.insert(local, id);
            if let Some(source) = source {
                self.source_bonds
                    .insert(InstanceBondId::new(source, local), id);
            }
        }
        self.groups.push(Some(Group {
            chemistry,
            atoms,
            bonds,
            instance_slot: source.map(MoleculeInstanceId::index),
            changed: false,
            class,
            class_explicit,
        }));
        group
    }
    /// The handles of one group's atoms and bonds by local ID.
    fn group_handles(&self, group: usize) -> EditMolecule {
        let group = self.groups[group].as_ref().expect("live group");
        EditMolecule {
            atoms: group
                .atoms
                .iter()
                .map(|(&local, &id)| (local, id))
                .collect(),
            bonds: group
                .bonds
                .iter()
                .map(|(&local, &id)| (local, id))
                .collect(),
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
        self.groups[group].as_mut().unwrap().bonds.insert(local, id);
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
        let removed = self.groups[other].take().unwrap();
        let group = self.groups[target].as_mut().unwrap();
        group.chemistry = GroupChemistry::Draft(Box::new(draft));
        group.changed = true;
        group.class = None;
        group.class_explicit = false;
        for (old, handle) in removed.atoms {
            let local = map.atoms()[&old];
            group.atoms.insert(local, handle);
            let loc = self.atoms.get_mut(&handle).unwrap();
            loc.group = target;
            loc.local = local;
        }
        for (old, handle) in removed.bonds {
            let local = map.bonds()[&old];
            group.bonds.insert(local, handle);
            let loc = self.bonds.get_mut(&handle).unwrap();
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
