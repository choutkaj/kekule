//! Transactional system edits, with molecular chemistry staged per affected occurrence.
mod append;
mod error;
mod hierarchy;
mod identity;
mod properties;
mod publication;

use super::{
    AtomSiteId, AtomSiteMetadata, ChainId, HierarchyError, InstanceAtomId, InstanceBondId,
    MoleculeClass, MoleculeDefinitionId, MoleculeInstanceId, ResidueClass, ResidueId, Topology,
    TopologyBuildError, TopologyBuilder,
};
use crate::core::{
    Atom, AtomId, BondId, BondOrder, Molecule, MoleculeEditor, MoleculeError,
    MoleculePublicationError,
};
use crate::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
pub(crate) use append::AppendMapping;
pub use error::*;
use hierarchy::EditHierarchy;
pub use hierarchy::{EditAtomSite, EditChain, EditResidue};
pub use identity::*;
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
    atoms: BTreeMap<AtomId, EditAtomId>,
    bonds: BTreeMap<BondId, EditBondId>,
    instance_slot: Option<usize>,
    changed: bool,
    class: Option<MoleculeClass>,
    class_explicit: bool,
}
#[derive(Debug, Clone, Copy)]
struct Location<Id> {
    group: usize,
    local: Id,
    slot: usize,
}

/// Detached coordinate-free structural editing state.
///
/// Chemical edits affect individual occurrences, even when definitions are reused.
/// Deleting bonds can split molecules; adding bonds can merge them. Publication
/// constructs valid connected definitions and one immutable topology snapshot.
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
/// use kekule::{smiles, topology::{Topology, TopologyEditor}};
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Topology::from_molecule(&molecule)?;
/// let source_bond = topology.bond_ids()[0];
/// let mut editor = TopologyEditor::from_topology(topology);
/// editor.delete_bond(editor.bond_handle(source_bond)?)?;
/// let edited = editor.finish()?;
/// assert_eq!(edited.instance_count(), 2);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct TopologyEditor {
    source: Option<Arc<Topology>>,
    groups: Vec<Option<Group>>,
    atoms: BTreeMap<EditAtomId, Location<AtomId>>,
    bonds: BTreeMap<EditBondId, Location<BondId>>,
    source_atoms: BTreeMap<InstanceAtomId, EditAtomId>,
    source_bonds: BTreeMap<InstanceBondId, EditBondId>,
    molecule_classes: BTreeMap<EditAtomId, MoleculeClass>,
    hierarchy: EditHierarchy,
    properties: Properties,
    revision: u64,
    pub(crate) structural_revision: u64,
}

impl Topology {
    /// Starts a detached draft from a shared topology, retaining its exact snapshot.
    ///
    /// This method borrows an `Arc<Topology>` so no-op publication can return
    /// that same allocation. For an owned `Topology`, use [`Self::into_editor`]
    /// to transfer ownership, or [`TopologyEditor::from_topology`] with either
    /// an owned value or an `Arc`. Clone the `Arc` when another owner must retain
    /// the original snapshot; molecular data is copied only when edited.
    pub fn edit(self: &Arc<Self>) -> TopologyEditor {
        TopologyEditor::from_topology(Arc::clone(self))
    }
    /// Transfers an owned topology into a structural draft. Successful consuming
    /// publication can move its unchanged molecular definitions into the result.
    pub fn into_editor(self) -> TopologyEditor {
        TopologyEditor::from_topology(self)
    }
}

impl TopologyEditor {
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
        self.properties = Properties::new();
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
        for (instance, value) in source.instances() {
            let molecule = source
                .definition(value.definition())
                .expect("published definition")
                .molecule();
            editor.register_group(
                GroupChemistry::Source(value.definition()),
                molecule,
                Some(instance),
                Some(source.definition(value.definition()).unwrap().class()),
                source
                    .molecule_class_overrides
                    .contains_key(&value.definition()),
            );
        }
        editor.import_hierarchy(&source);
        editor
    }

    pub fn source_topology(&self) -> Option<&Topology> {
        self.source.as_deref()
    }
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }
    pub fn atom_ids(&self) -> impl ExactSizeIterator<Item = EditAtomId> + '_ {
        self.atoms.keys().copied()
    }
    pub fn bond_ids(&self) -> impl ExactSizeIterator<Item = EditBondId> + '_ {
        self.bonds.keys().copied()
    }
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
    pub fn atom_handle(&self, source: InstanceAtomId) -> Result<EditAtomId, TopologyEditError> {
        self.source_atoms
            .get(&source)
            .copied()
            .filter(|id| self.atoms.contains_key(id))
            .ok_or(TopologyEditError::InvalidSourceAtom(source))
    }
    pub fn bond_handle(&self, source: InstanceBondId) -> Result<EditBondId, TopologyEditError> {
        self.source_bonds
            .get(&source)
            .copied()
            .filter(|id| self.bonds.contains_key(id))
            .ok_or(TopologyEditError::InvalidSourceBond(source))
    }
    pub fn atom(&self, id: EditAtomId) -> Result<&Atom, TopologyEditError> {
        let loc = self.atom_location(id)?;
        Ok(self.molecule(loc.group).atom(loc.local)?)
    }
    pub fn bond(&self, id: EditBondId) -> Result<EditBond, TopologyEditError> {
        let loc = self.bond_location(id)?;
        let group = self.groups[loc.group].as_ref().unwrap();
        let bond = self.molecule(loc.group).bond(loc.local)?;
        Ok(EditBond {
            a: group.atoms[&bond.a()],
            b: group.atoms[&bond.b()],
            order: bond.order,
        })
    }
    pub fn neighbors(
        &self,
        id: EditAtomId,
    ) -> Result<impl Iterator<Item = EditAtomId> + '_, TopologyEditError> {
        let loc = self.atom_location(id)?;
        let group = self.groups[loc.group].as_ref().unwrap();
        Ok(self
            .molecule(loc.group)
            .neighbors(loc.local)?
            .map(|id| group.atoms[&id]))
    }
    pub fn incident_bonds(
        &self,
        id: EditAtomId,
    ) -> Result<impl Iterator<Item = EditBondId> + '_, TopologyEditError> {
        let loc = self.atom_location(id)?;
        let group = self.groups[loc.group].as_ref().unwrap();
        Ok(self
            .molecule(loc.group)
            .incident_bonds(loc.local)?
            .map(|(id, _)| group.bonds[&id]))
    }
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
            .bond_between(left.local, right.local)?
            .map(|id| self.groups[left.group].as_ref().unwrap().bonds[&id]))
    }
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
    pub fn add_atom(&mut self, atom: Atom) -> Result<EditAtomId, TopologyEditError> {
        let mut draft = MoleculeEditor::new();
        let local = draft.add_atom(atom)?;
        let id = EditAtomId::new();
        let slot = self.properties.atoms().len();
        let group = self.groups.len();
        self.atoms.insert(id, Location { group, local, slot });
        self.groups.push(Some(Group {
            chemistry: GroupChemistry::Draft(Box::new(draft)),
            atoms: BTreeMap::from([(local, id)]),
            bonds: BTreeMap::new(),
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
    pub fn add_molecule(&mut self, molecule: &Molecule) -> Result<EditMolecule, TopologyEditError> {
        let owned = Arc::new(molecule.clone());
        let result = self.register_group(
            GroupChemistry::Added(Arc::clone(&owned)),
            &owned,
            None,
            None,
            false,
        );
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

    pub fn replace_atom(&mut self, id: EditAtomId, atom: Atom) -> Result<Atom, TopologyEditError> {
        let loc = self.atom_location(id)?;
        if self.atom(id)? == &atom {
            return Ok(atom);
        }
        let previous = self.draft_mut(loc.group).replace_atom(loc.local, atom)?;
        self.mark_group_changed(loc.group);
        self.invalidate_component_classes([id]);
        self.invalidate_residue_for_atom(id);
        self.changed();
        Ok(previous)
    }

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

    /// Validates every handle before deleting any atom; duplicates are ignored.
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
    pub fn delete_instance(
        &mut self,
        source: MoleculeInstanceId,
    ) -> Result<Vec<(EditAtomId, Atom)>, TopologyEditError> {
        if self
            .source
            .as_ref()
            .is_none_or(|t| t.instance(source).is_err())
        {
            return Err(TopologyEditError::InvalidSourceInstance(source));
        }
        let ids = self
            .source_atoms
            .iter()
            .filter(|(id, handle)| id.molecule() == source && self.atoms.contains_key(handle))
            .map(|(_, &h)| h)
            .collect::<Vec<_>>();
        self.delete_atoms(ids)
    }

    /// Adds an asserted bond. Different occurrences are merged transactionally.
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
        let mut draft = self.molecule(target).edit();
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

    pub fn delete_bond(&mut self, id: EditBondId) -> Result<EditBond, TopologyEditError> {
        let previous = self.bond(id)?;
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
    pub fn set_bond_order(
        &mut self,
        id: EditBondId,
        order: BondOrder,
    ) -> Result<(), TopologyEditError> {
        let previous = self.bond(id)?;
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
    pub fn set_bond_endpoints(
        &mut self,
        id: EditBondId,
        a: EditAtomId,
        b: EditAtomId,
    ) -> Result<(), TopologyEditError> {
        let old = self.bond(id)?;
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
        let old = self.bond(id)?;
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
                let mut draft = self.molecule(target).edit();
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

    pub fn replace_bond(
        &mut self,
        id: EditBondId,
        replacement: EditBond,
    ) -> Result<EditBond, TopologyEditError> {
        let previous = self.bond(id)?;
        let mut staged = self.clone();
        staged.set_bond_endpoints_staged(id, replacement.a, replacement.b)?;
        staged.set_bond_order(id, replacement.order)?;
        *self = staged;
        Ok(previous)
    }

    /// Explicitly assigns the current connected component's class. Assigning its
    /// existing inferred class still records an override for subsequent metadata
    /// edits. A later chemical change invalidates this component's override.
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
            let draft = self.molecule(group).edit();
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
        self.properties.clear_owner();
    }
    fn register_group(
        &mut self,
        chemistry: GroupChemistry,
        molecule: &Molecule,
        source: Option<MoleculeInstanceId>,
        class: Option<MoleculeClass>,
        class_explicit: bool,
    ) -> EditMolecule {
        let group = self.groups.len();
        let mut result = EditMolecule::default();
        let atom_start = if source.is_some() {
            self.atoms.len()
        } else {
            self.properties.atoms().len()
        };
        let bond_start = if source.is_some() {
            self.bonds.len()
        } else {
            self.properties.bonds().len()
        };
        for (index, local) in molecule.atom_ids().enumerate() {
            let id = EditAtomId::new();
            self.atoms.insert(
                id,
                Location {
                    group,
                    local,
                    slot: atom_start + index,
                },
            );
            result.atoms.insert(local, id);
            if let Some(source) = source {
                self.source_atoms
                    .insert(InstanceAtomId::new(source, local), id);
            }
        }
        for (index, local) in molecule.bond_ids().enumerate() {
            let id = EditBondId::new();
            self.bonds.insert(
                id,
                Location {
                    group,
                    local,
                    slot: bond_start + index,
                },
            );
            result.bonds.insert(local, id);
            if let Some(source) = source {
                self.source_bonds
                    .insert(InstanceBondId::new(source, local), id);
            }
        }
        self.groups.push(Some(Group {
            chemistry,
            atoms: result.atoms.clone(),
            bonds: result.bonds.clone(),
            instance_slot: source.map(MoleculeInstanceId::index),
            changed: false,
            class,
            class_explicit,
        }));
        result
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
        map: &crate::core::MoleculeAppendMapping,
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
    pub order: BondOrder,
}
impl EditBond {
    pub fn new(a: EditAtomId, b: EditAtomId, order: BondOrder) -> Self {
        Self { a, b, order }
    }
    pub fn a(self) -> EditAtomId {
        self.a
    }
    pub fn b(self) -> EditAtomId {
        self.b
    }
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
    pub fn atoms(&self) -> &BTreeMap<AtomId, EditAtomId> {
        &self.atoms
    }
    pub fn bonds(&self) -> &BTreeMap<BondId, EditBondId> {
        &self.bonds
    }
}
