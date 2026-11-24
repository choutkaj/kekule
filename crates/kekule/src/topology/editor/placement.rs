//! Placement edits inside the draft of the instance that holds the parent:
//! hydrogen materialization and collapse, attached atoms, and the atoms of an
//! attached fragment, so represented stereo follows the hydrogens and the work
//! scales with the edited atoms rather than the system.
use super::{EditAtomId, EditBondId, GroupGraph, Location, TopologyEditError, TopologyEditor};
use crate::algorithms::{
    assign_collapsed_hydrogen_counts, collapse_stereo_hydrogen_carrier,
    materialize_stereo_hydrogen_carriers, rdkit_inferred_hydrogen_count, rdkit_valence_violation,
    removable_hydrogen, stereo_hydrogen_bindings, stereo_hydrogen_is_collapsible,
    stereo_hydrogen_parents, AddedHydrogenOrigin, RetainedHydrogenReason,
};
use crate::molecule::{Atom, AtomId, BondOrder, Element, ImplicitHydrogens};
use crate::stereo::{StereoCarrier, StereoElementKind};
use std::collections::{BTreeMap, BTreeSet};

/// The explicit hydrogens
/// [`StructureEditor::remove_hydrogens`](crate::structure::StructureEditor::remove_hydrogens)
/// collapsed on one atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedHydrogens {
    atom: EditAtomId,
    removed: Vec<EditAtomId>,
    retained: Vec<(EditAtomId, RetainedHydrogenReason)>,
    hydrogens: ImplicitHydrogens,
    implicit_hydrogens: usize,
}

impl RemovedHydrogens {
    /// The atom whose hydrogens were collapsed.
    pub const fn atom(&self) -> EditAtomId {
        self.atom
    }
    /// The deleted hydrogens, which no longer name draft atoms.
    pub fn removed(&self) -> &[EditAtomId] {
        &self.removed
    }
    /// Explicit hydrogens kept, with why, in handle order.
    pub fn retained(&self) -> &[(EditAtomId, RetainedHydrogenReason)] {
        &self.retained
    }
    /// The atom's resulting declaration: inference is kept only when it
    /// reproduces the complete implicit count, which is otherwise fixed.
    pub const fn declaration(&self) -> ImplicitHydrogens {
        self.hydrogens
    }
    /// The atom's complete implicit hydrogen count afterwards.
    pub const fn implicit_hydrogens(&self) -> usize {
        self.implicit_hydrogens
    }
}

/// The outcome of [`TopologyEditor::collapse_hydrogens`]: one entry per
/// parent, and the storage slots of the deleted atoms and bonds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CollapsedHydrogenSet {
    pub(crate) parents: Vec<RemovedHydrogens>,
    pub(crate) dead_atom_slots: Vec<usize>,
    pub(crate) dead_bond_slots: Vec<usize>,
}

/// The atom [`TopologyEditor::attach_in_instance`] added.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AttachedInInstance {
    pub(crate) atom: EditAtomId,
    pub(crate) bond: EditBondId,
    /// The deleted hydrogen's atom slot and bond slots.
    pub(crate) dead: Option<(usize, Vec<usize>)>,
}

/// One atom of the bond a ring fuses onto, in a [`FusionPlan`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FusionEndPlan {
    pub(crate) atom: EditAtomId,
    /// The explicit hydrogen whose stereo place the atom's new ring bond
    /// takes and which is deleted; `None` takes an implicit hydrogen's.
    pub(crate) replaced: Option<EditAtomId>,
    /// Further explicit hydrogens the atom gives up.
    pub(crate) removed: Vec<EditAtomId>,
}

/// The edit [`TopologyEditor::fuse_in_instance`] applies.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FusionPlan {
    /// The atoms of the shared bond.
    pub(crate) ends: [FusionEndPlan; 2],
    /// The shared bond and its order afterwards.
    pub(crate) shared: EditBondId,
    pub(crate) shared_order: BondOrder,
    /// The new atoms.
    pub(crate) atoms: Vec<Atom>,
    /// The new bonds, between indices into `atoms`; `atoms.len() + i` names
    /// end `i`, which has exactly one new bond.
    pub(crate) bonds: Vec<(usize, usize, BondOrder)>,
}

/// What [`TopologyEditor::fuse_in_instance`] added and deleted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FusedInInstance {
    /// The new atoms, in plan order.
    pub(crate) atoms: Vec<EditAtomId>,
    /// The deleted hydrogens' atom slots and bond slots.
    pub(crate) dead_atom_slots: Vec<usize>,
    pub(crate) dead_bond_slots: Vec<usize>,
}

/// Why a hydrogen operation inside a topology draft failed.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum DraftHydrogenError {
    /// A handle is invalid or the molecular draft rejected an edit.
    Topology(TopologyEditError),
    /// The atom's represented stereo cannot follow its hydrogens.
    Stereo(EditAtomId),
    /// The atom's implicit count would overflow.
    Overflow(EditAtomId, usize),
}

impl From<TopologyEditError> for DraftHydrogenError {
    fn from(error: TopologyEditError) -> Self {
        Self::Topology(error)
    }
}

impl From<crate::molecule::MoleculeError> for DraftHydrogenError {
    fn from(error: crate::molecule::MoleculeError) -> Self {
        Self::Topology(error.into())
    }
}

/// A planned collapse of one parent's hydrogens: locals within its group.
struct ParentCollapse {
    parent: EditAtomId,
    local: AtomId,
    removable: Vec<AtomId>,
    retained: Vec<(EditAtomId, RetainedHydrogenReason)>,
}

/// The planned collapse of one group.
struct GroupCollapse {
    group: usize,
    entries: Vec<ParentCollapse>,
    totals: BTreeMap<AtomId, usize>,
    by_parent: BTreeMap<AtomId, Vec<AtomId>>,
    fixed: BTreeSet<AtomId>,
    bindings: BTreeMap<AtomId, Vec<crate::stereo::StereoElementId>>,
    handles: BTreeMap<AtomId, EditAtomId>,
}

impl TopologyEditor {
    /// The implicit hydrogens `atom` would materialize, with their origin: its
    /// fixed count, or the valence model's inference on the draft graph, which
    /// needs no installed perception.
    pub(crate) fn materializable_hydrogens(
        &self,
        atom: EditAtomId,
    ) -> Result<(usize, AddedHydrogenOrigin), TopologyEditError> {
        let loc = self.atom_location(atom)?;
        Ok(with_group_graph!(self.molecule(loc.group), |molecule| {
            let value = molecule.atom(loc.local).expect("live draft atom");
            match value.hydrogens {
                ImplicitHydrogens::Fixed(count) => (usize::from(count), AddedHydrogenOrigin::Fixed),
                ImplicitHydrogens::Inferred => (
                    usize::from(rdkit_inferred_hydrogen_count(molecule, loc.local, value)),
                    AddedHydrogenOrigin::Inferred,
                ),
            }
        }))
    }

    /// Checks that materializing `count` hydrogens on `atom` keeps its
    /// represented stereo complete: an implicit-hydrogen carrier needs exactly
    /// one materialized hydrogen, and an axis cannot name one.
    pub(crate) fn check_hydrogen_stereo(
        &self,
        atom: EditAtomId,
        count: usize,
    ) -> Result<(), DraftHydrogenError> {
        let loc = self.atom_location(atom)?;
        for (_, value) in self.molecule(loc.group).stereo_elements() {
            let implicit = match &value.kind {
                StereoElementKind::Tetrahedral(stereo) if stereo.center == loc.local => stereo
                    .carriers
                    .iter()
                    .filter(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen))
                    .count(),
                StereoElementKind::DoubleBond(stereo) => {
                    usize::from(
                        stereo.left == loc.local
                            && stereo.left_carrier == StereoCarrier::ImplicitHydrogen,
                    ) + usize::from(
                        stereo.right == loc.local
                            && stereo.right_carrier == StereoCarrier::ImplicitHydrogen,
                    )
                }
                StereoElementKind::Axis(stereo)
                    if count > 0 && stereo.carriers.contains(&StereoCarrier::ImplicitHydrogen) =>
                {
                    return Err(DraftHydrogenError::Stereo(atom));
                }
                _ => 0,
            };
            if implicit > 1 || (implicit == 1 && count > 0 && count != 1) {
                return Err(DraftHydrogenError::Stereo(atom));
            }
        }
        Ok(())
    }

    /// Adds `count` hydrogens bonded to `parent` inside its instance's
    /// draft, naming them in its implicit-hydrogen stereo carriers, and zeroes
    /// a fixed count. The caller has checked [`Self::check_hydrogen_stereo`].
    pub(crate) fn materialize_hydrogens(
        &mut self,
        parent: EditAtomId,
        count: usize,
        origin: AddedHydrogenOrigin,
    ) -> Result<Vec<(EditAtomId, EditBondId)>, TopologyEditError> {
        let loc = self.atom_location(parent)?;
        if count == 0 {
            return Ok(Vec::new());
        }
        let hydrogen = Element::from_atomic_number(1).expect("hydrogen is an element");
        let draft = self.draft_mut(loc.group);
        let mut locals = Vec::with_capacity(count);
        for _ in 0..count {
            let atom = draft.add_atom(Atom::new(hydrogen))?;
            // The parent's stereo lists an implicit hydrogen until the carrier
            // is renamed below, so insert without pruning.
            let bond = draft.add_bond_unpruned(loc.local, atom, BondOrder::Single)?;
            locals.push((atom, bond));
        }
        if origin == AddedHydrogenOrigin::Fixed {
            draft.require_atom_mut(loc.local)?.hydrogens = ImplicitHydrogens::Fixed(0);
        }
        let added = BTreeMap::from([(loc.local, locals.iter().map(|(atom, _)| *atom).collect())]);
        materialize_stereo_hydrogen_carriers(draft, &added)?;
        Ok(locals
            .into_iter()
            .map(|(local, bond)| self.register_bonded_atom(loc.group, local, bond, parent))
            .collect())
    }

    /// Bonds a new atom to `parent` inside its instance's draft, as
    /// [`Self::add_atom`] and [`Self::add_bond`] would without copying the
    /// instance. With `replaced`, an explicit hydrogen of the parent, the new
    /// atom takes that hydrogen's place in represented stereo and the
    /// hydrogen is deleted; otherwise a single bond takes the place of an
    /// implicit hydrogen there.
    pub(crate) fn attach_in_instance(
        &mut self,
        parent: EditAtomId,
        atom: Atom,
        order: BondOrder,
        replaced: Option<EditAtomId>,
    ) -> Result<AttachedInInstance, TopologyEditError> {
        let loc = self.atom_location(parent)?;
        let replaced_local = replaced
            .map(|hydrogen| self.atom_location(hydrogen).map(|loc| loc.local))
            .transpose()?;
        // The parent's count follows the new bond and the bond to the
        // hydrogen it replaces, as edits make it follow.
        let mut changes = super::ValenceChanges::default();
        changes.change(parent, None, Some(order));
        if let Some(hydrogen) = replaced {
            let bond = self
                .bond_between(parent, hydrogen)?
                .ok_or(TopologyEditError::InvalidAtom(hydrogen))?;
            let bond = self.require_bond(bond)?.order;
            changes.change(parent, Some(bond), None);
        }
        let counts = self.fixed_counts_after(&changes)?;
        let local = self.draft_mut(loc.group).add_atom(atom)?;
        let bond = bond_in_place_of_hydrogen(
            self.draft_mut(loc.group),
            loc.local,
            local,
            order,
            replaced_local,
        )?;
        let (atom, bond) = self.register_bonded_atom(loc.group, local, bond, parent);
        self.record_pruned(loc.group);
        let mut dead = None;
        if let Some(hydrogen) = replaced {
            let atom_slot = self.atom_location(hydrogen)?.slot;
            let bond_slots = self
                .incident_bonds(hydrogen)?
                .map(|bond| self.bond_location(bond).map(|loc| loc.slot))
                .collect::<Result<Vec<_>, _>>()?;
            self.remove_atom(hydrogen)?;
            dead = Some((atom_slot, bond_slots));
        }
        self.set_fixed_counts(counts)?;
        Ok(AttachedInInstance { atom, bond, dead })
    }

    /// Fuses new ring atoms onto the bond `plan.shared` inside its instance's
    /// draft: adds `plan.atoms`, bonds them as `plan.bonds` lists, gives the
    /// shared bond `plan.shared_order`, and deletes each end's replaced and
    /// removed explicit hydrogens. An end's one new bond takes the replaced
    /// hydrogen's place in represented stereo, or an implicit hydrogen's
    /// without one; stereo of the shared bond is dropped, as the ring fixes
    /// it, and listed in [`Self::pruned_stereo`]; fixed counts follow the
    /// whole edit, as edits make them follow. The caller has checked the
    /// plan, so nothing changes on error.
    pub(crate) fn fuse_in_instance(
        &mut self,
        plan: &FusionPlan,
    ) -> Result<FusedInInstance, TopologyEditError> {
        let shared = self.require_bond(plan.shared)?;
        let group = self.atom_location(plan.ends[0].atom)?.group;
        let count = plan.atoms.len();
        let end_of = |index: usize| index.checked_sub(count);
        let mut changes = super::ValenceChanges::default();
        for &(a, b, order) in &plan.bonds {
            for end in [a, b].into_iter().filter_map(end_of) {
                changes.change(plan.ends[end].atom, None, Some(order));
            }
        }
        if shared.order != plan.shared_order {
            changes.bond(
                [shared.a, shared.b],
                Some(shared.order),
                Some(plan.shared_order),
            );
        }
        for end in &plan.ends {
            for &hydrogen in end.replaced.iter().chain(&end.removed) {
                let bond = self
                    .bond_between(end.atom, hydrogen)?
                    .ok_or(TopologyEditError::InvalidAtom(hydrogen))?;
                changes.change(end.atom, Some(self.require_bond(bond)?.order), None);
            }
        }
        let counts = self.fixed_counts_after(&changes)?;
        self.draft_mut(group).check_room(count, plan.bonds.len())?;
        self.mark_group_changed(group);
        self.changed();
        let mut atoms = Vec::with_capacity(count);
        for atom in &plan.atoms {
            let local = self.draft_mut(group).add_atom(atom.clone())?;
            atoms.push(self.register_atom(group, local));
        }
        // The ring fixes the configuration of the double bond it shares, so
        // the bond's stereo is dropped and reported as pruned.
        let shared_local = self.bond_location(plan.shared)?.local;
        let fixed = self
            .draft_mut(group)
            .stereo_elements()
            .filter(|(_, element)| {
                matches!(&element.kind, StereoElementKind::DoubleBond(stereo) if stereo.bond == shared_local)
            })
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        for id in fixed {
            let element = self.draft_mut(group).remove_stereo_element(id)?;
            let entry = self.pruned_entry(group, &element);
            self.pruned_stereo.push(entry);
        }
        for &(a, b, order) in &plan.bonds {
            if end_of(a).is_none() && end_of(b).is_none() {
                // The fragment's counts already describe its own bonds.
                self.connect(atoms[a], atoms[b], order)?;
            }
        }
        // Each end's hydrogens go as soon as its ring bond is in, so no edit
        // sees the end over its valence and prunes its stereo.
        let mut dead_atom_slots = Vec::new();
        let mut dead_bond_slots = Vec::new();
        for (end, plan_end) in plan.ends.iter().enumerate() {
            let &(a, b, order) = plan
                .bonds
                .iter()
                .find(|(a, b, _)| end_of(*a) == Some(end) || end_of(*b) == Some(end))
                .expect("each end has one new bond");
            let new = atoms[a.min(b)];
            let parent_local = self.atom_location(plan_end.atom)?.local;
            let new_local = self.atom_location(new)?.local;
            let replaced = plan_end
                .replaced
                .map(|hydrogen| self.atom_location(hydrogen).map(|loc| loc.local))
                .transpose()?;
            let bond = bond_in_place_of_hydrogen(
                self.draft_mut(group),
                parent_local,
                new_local,
                order,
                replaced,
            )?;
            self.register_bond(group, bond, plan_end.atom, new);
            self.record_pruned(group);
            for &hydrogen in plan_end.replaced.iter().chain(&plan_end.removed) {
                dead_atom_slots.push(self.atom_location(hydrogen)?.slot);
                for bond in self.incident_bonds(hydrogen)?.collect::<Vec<_>>() {
                    dead_bond_slots.push(self.bond_location(bond)?.slot);
                }
                self.remove_atom(hydrogen)?;
            }
        }
        if shared.order != plan.shared_order {
            self.reorder_bond(plan.shared, plan.shared_order)?;
        }
        self.set_fixed_counts(counts)?;
        Ok(FusedInInstance {
            atoms,
            dead_atom_slots,
            dead_bond_slots,
        })
    }

    /// Whether `atom`'s represented valence breaks the valence model's
    /// strict rules.
    pub(crate) fn exceeds_valence(&self, atom: EditAtomId) -> Result<bool, TopologyEditError> {
        let loc = self.atom_location(atom)?;
        Ok(with_group_graph!(self.molecule(loc.group), |molecule| {
            let value = molecule.atom(loc.local).expect("live draft atom");
            rdkit_valence_violation(molecule, loc.local, value).is_some()
        }))
    }

    /// Adds an unbonded atom to the instance draft that holds `anchor`, as
    /// [`Self::add_atom`] would without starting an instance of its own; bond
    /// it with [`Self::add_bond`], which then merges nothing.
    pub(crate) fn add_atom_beside(
        &mut self,
        anchor: EditAtomId,
        atom: Atom,
    ) -> Result<EditAtomId, TopologyEditError> {
        let loc = self.atom_location(anchor)?;
        let local = self.draft_mut(loc.group).add_atom(atom)?;
        let id = self.register_atom(loc.group, local);
        self.mark_group_changed(loc.group);
        self.changed();
        Ok(id)
    }

    /// Whether the atom storage slot `slot` holds a live atom.
    pub(crate) fn is_live_atom_slot(&self, slot: usize) -> bool {
        slot < self.atom_slot_count() && !self.dead_atom_slots.contains(&slot)
    }

    fn register_atom(&mut self, group: usize, local: AtomId) -> EditAtomId {
        let id = EditAtomId::new();
        let slot = self.properties.atoms().len();
        self.properties.resize_atoms(slot + 1);
        self.atoms.insert(id, Location { group, local, slot });
        self.group_mut(group).atoms.insert(local, id);
        id
    }

    fn register_bonded_atom(
        &mut self,
        group: usize,
        local: AtomId,
        bond: crate::molecule::BondId,
        parent: EditAtomId,
    ) -> (EditAtomId, EditBondId) {
        let id = self.register_atom(group, local);
        (id, self.register_bond(group, bond, parent, id))
    }

    /// Collapses the ordinary explicit hydrogens bonded to each parent into its
    /// implicit count, as [`Molecule::remove_hydrogens`] does, and returns one
    /// entry per parent in input order, leaving out selected hydrogens that
    /// the collapse deletes. `annotated` tells, from a hydrogen's
    /// atom slot and its bond's slot, whether the caller holds annotations
    /// on either that the topology does not see; such hydrogens are kept.
    /// Nothing changes on error.
    ///
    /// [`Molecule::remove_hydrogens`]: crate::molecule::Molecule::remove_hydrogens
    pub(crate) fn collapse_hydrogens(
        &mut self,
        parents: &[EditAtomId],
        annotated: impl Fn(usize, usize) -> Option<RetainedHydrogenReason>,
    ) -> Result<CollapsedHydrogenSet, DraftHydrogenError> {
        let mut groups = BTreeMap::<usize, Vec<ParentCollapse>>::new();
        let mut seen = BTreeSet::new();
        for &parent in parents {
            let loc = self.atom_location(parent)?;
            if seen.insert(parent) {
                groups.entry(loc.group).or_default().push(ParentCollapse {
                    parent,
                    local: loc.local,
                    removable: Vec::new(),
                    retained: Vec::new(),
                });
            }
        }
        let mut plans = Vec::with_capacity(groups.len());
        for (group, mut entries) in groups {
            let properties = self.molecule(group).properties();
            let perception = self.molecule(group).perception();
            with_group_graph!(self.molecule(group), |molecule| {
                let mut candidates = BTreeSet::new();
                for entry in &mut entries {
                    for hydrogen in molecule.require_neighbors(entry.local)? {
                        let atom = molecule.atom(hydrogen).expect("neighbors are live");
                        if atom.element.atomic_number() != 1 {
                            continue;
                        }
                        let handle = self.group_atom(group, hydrogen);
                        let outcome = removable_hydrogen(molecule, properties, hydrogen, atom)
                            .expect("neighbors are live")
                            .and_then(|(_, bond)| {
                                let atom_slot = self.atom_slot(handle).expect("live atom");
                                let bond_slot = self
                                    .bond_slot(self.group_bond(group, bond))
                                    .expect("live bond");
                                let has_data =
                                    |table: &crate::properties::RawPropertyTable, slot| {
                                        table.row_has_data(slot).expect("live rows")
                                    };
                                if has_data(self.properties.atoms().raw(), atom_slot) {
                                    return Err(RetainedHydrogenReason::AtomProperties);
                                }
                                if has_data(self.properties.bonds().raw(), bond_slot) {
                                    return Err(RetainedHydrogenReason::BondProperties);
                                }
                                annotated(atom_slot, bond_slot).map_or(Ok(()), Err)
                            });
                        match outcome {
                            Ok(()) => {
                                candidates.insert(hydrogen);
                                entry.removable.push(hydrogen);
                            }
                            Err(reason) => entry.retained.push((handle, reason)),
                        }
                    }
                }
                let mut removable = BTreeSet::new();
                for entry in &mut entries {
                    let mut kept = Vec::new();
                    for &hydrogen in &entry.removable {
                        if stereo_hydrogen_is_collapsible(
                            molecule,
                            perception,
                            hydrogen,
                            entry.local,
                            &candidates,
                        ) {
                            removable.insert(hydrogen);
                            kept.push(hydrogen);
                        } else {
                            entry.retained.push((
                                self.group_atom(group, hydrogen),
                                RetainedHydrogenReason::UnsupportedStereoRole,
                            ));
                        }
                    }
                    entry.removable = kept;
                    entry.retained.sort_by_key(|(handle, _)| *handle);
                }
                let mut totals = BTreeMap::new();
                let mut by_parent = BTreeMap::new();
                for entry in entries.iter().filter(|entry| !entry.removable.is_empty()) {
                    let value = molecule.atom(entry.local).expect("live parent");
                    let implicit = match value.hydrogens {
                        ImplicitHydrogens::Fixed(count) => usize::from(count),
                        ImplicitHydrogens::Inferred => {
                            usize::from(rdkit_inferred_hydrogen_count(molecule, entry.local, value))
                        }
                    };
                    let total = implicit + entry.removable.len();
                    if u8::try_from(total).is_err() {
                        return Err(DraftHydrogenError::Overflow(entry.parent, total));
                    }
                    totals.insert(entry.local, total);
                    by_parent.insert(entry.local, entry.removable.clone());
                }
                plans.push(GroupCollapse {
                    group,
                    fixed: stereo_hydrogen_parents(molecule, &removable),
                    bindings: stereo_hydrogen_bindings(molecule, &removable),
                    handles: removable
                        .iter()
                        .map(|hydrogen| (*hydrogen, self.group_atom(group, *hydrogen)))
                        .collect(),
                    entries,
                    totals,
                    by_parent,
                });
            });
        }
        let mut results = BTreeMap::new();
        let mut dead_atoms = Vec::new();
        let mut dead_bonds = Vec::new();
        for plan in plans {
            let group = plan.group;
            for (&local, &handle) in &plan.handles {
                collapse_stereo_hydrogen_carrier(
                    self.draft_mut(group),
                    local,
                    plan.bindings.get(&local).map(Vec::as_slice).unwrap_or(&[]),
                )
                .expect("planned stereo elements are live");
                dead_atoms.push(self.atom_location(handle)?.slot);
                for bond in self.incident_bonds(handle)?.collect::<Vec<_>>() {
                    dead_bonds.push(self.bond_location(bond)?.slot);
                }
                self.remove_atom(handle)?;
            }
            let adjustments = assign_collapsed_hydrogen_counts(
                self.draft_mut(group),
                &plan.totals,
                &plan.by_parent,
                &plan.fixed,
            )
            .expect("planned counts fit and name live parents")
            .into_iter()
            .map(|adjustment| (adjustment.parent, adjustment))
            .collect::<BTreeMap<_, _>>();
            // A selected hydrogen collapsed into its selected parent is
            // reported under that parent only.
            let deleted = plan.handles.values().copied().collect::<BTreeSet<_>>();
            for entry in plan.entries {
                if deleted.contains(&entry.parent) {
                    continue;
                }
                let (hydrogens, implicit_hydrogens) = match adjustments.get(&entry.local) {
                    Some(adjustment) => (adjustment.hydrogens, adjustment.implicit_hydrogens),
                    None => with_group_graph!(self.molecule(group), |molecule| {
                        let value = molecule.atom(entry.local).expect("live parent");
                        let implicit = match value.hydrogens {
                            ImplicitHydrogens::Fixed(count) => usize::from(count),
                            ImplicitHydrogens::Inferred => usize::from(
                                rdkit_inferred_hydrogen_count(molecule, entry.local, value),
                            ),
                        };
                        (value.hydrogens, implicit)
                    }),
                };
                results.insert(
                    entry.parent,
                    RemovedHydrogens {
                        atom: entry.parent,
                        removed: entry
                            .removable
                            .iter()
                            .map(|local| plan.handles[local])
                            .collect(),
                        retained: entry.retained,
                        hydrogens,
                        implicit_hydrogens,
                    },
                );
            }
            if !plan.totals.is_empty() {
                self.mark_group_changed(group);
                self.changed();
            }
        }
        Ok(CollapsedHydrogenSet {
            parents: parents
                .iter()
                .filter_map(|parent| results.remove(parent))
                .collect(),
            dead_atom_slots: dead_atoms,
            dead_bond_slots: dead_bonds,
        })
    }

    /// The atoms placed in `residue`, or `None` when it names no live residue.
    pub(crate) fn residue_atom_ids(
        &self,
        residue: super::EditResidueId,
    ) -> Option<impl Iterator<Item = EditAtomId> + '_> {
        self.hierarchy
            .residues
            .contains_key(&residue)
            .then(|| self.residue_atom_list(residue).into_iter())
    }
}

/// Bonds `parent` to `atom` in `draft`. A single bond takes the place of
/// `replaced`, an explicit hydrogen of the parent about to be deleted, or
/// else of an implicit hydrogen in the parent's represented stereo; another
/// order prunes the stereo it invalidates.
fn bond_in_place_of_hydrogen(
    draft: &mut crate::molecule::MoleculeEditor,
    parent: AtomId,
    atom: AtomId,
    order: BondOrder,
    replaced: Option<AtomId>,
) -> Result<crate::molecule::BondId, crate::molecule::MoleculeError> {
    if order != BondOrder::Single {
        return draft.add_bond(parent, atom, order);
    }
    // The new atom takes over a hydrogen's stereo role below, so insert
    // without pruning.
    let bond = draft.add_bond_unpruned(parent, atom, order)?;
    let from = replaced.map_or(StereoCarrier::ImplicitHydrogen, StereoCarrier::Atom);
    let replacements = draft
        .stereo_elements()
        .filter_map(|(id, element)| {
            let mut element = element.clone();
            let mut changed = false;
            let mut swap = |owner: AtomId, carrier: &mut StereoCarrier| {
                if owner == parent && *carrier == from && !changed {
                    *carrier = StereoCarrier::Atom(atom);
                    changed = true;
                }
            };
            match &mut element.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    let center = stereo.center;
                    stereo.carriers.iter_mut().for_each(|c| swap(center, c));
                }
                StereoElementKind::DoubleBond(stereo) => {
                    swap(stereo.left, &mut stereo.left_carrier);
                    swap(stereo.right, &mut stereo.right_carrier);
                }
                StereoElementKind::Axis(_) => {}
            }
            changed.then_some((id, element))
        })
        .collect::<Vec<_>>();
    // The replaced hydrogen is still bonded until it is deleted, so the new
    // carrier is written without validation, as hydrogen collapse does;
    // publication validates the result.
    for (id, element) in replacements {
        draft.set_stereo_element_unchecked(id, element);
    }
    Ok(bond)
}
