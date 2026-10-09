//! Compound editor operations. None exposes the unfinished graph as a Molecule.
use std::collections::{BTreeMap, BTreeSet};

use super::*;

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

    /// Maps every live entity of `molecule` to itself.
    pub(crate) fn identity(molecule: &Molecule) -> Self {
        Self {
            atoms: molecule.atom_ids().map(|id| (id, id)).collect(),
            bonds: molecule.bond_ids().map(|id| (id, id)).collect(),
            stereo_elements: molecule.stereo_element_ids().map(|id| (id, id)).collect(),
            stereo_groups: molecule.stereo_groups().map(|(id, _)| (id, id)).collect(),
        }
    }

    pub fn atom(&self, source: AtomId) -> Option<AtomId> {
        self.atoms.get(&source).copied()
    }
    pub fn bond(&self, source: BondId) -> Option<BondId> {
        self.bonds.get(&source).copied()
    }
    pub fn stereo_element(&self, source: StereoElementId) -> Option<StereoElementId> {
        self.stereo_elements.get(&source).copied()
    }
    pub fn stereo_group(&self, source: StereoGroupId) -> Option<StereoGroupId> {
        self.stereo_groups.get(&source).copied()
    }
    pub fn atoms(&self) -> &BTreeMap<AtomId, AtomId> {
        &self.atoms
    }
    pub fn bonds(&self) -> &BTreeMap<BondId, BondId> {
        &self.bonds
    }
    pub fn stereo_elements(&self) -> &BTreeMap<StereoElementId, StereoElementId> {
        &self.stereo_elements
    }
    pub fn stereo_groups(&self) -> &BTreeMap<StereoGroupId, StereoGroupId> {
        &self.stereo_groups
    }
}

impl MoleculeEditor {
    pub fn is_empty(&self) -> bool {
        self.atom_count() == 0
    }

    /// Live connected components in stable atom order, including isolated atoms.
    pub fn connected_components(&self) -> Vec<Vec<AtomId>> {
        self.working.connected_components()
    }

    /// Empty editing state is not a connected molecule.
    pub fn is_connected(&self) -> bool {
        self.working.validate_connected().is_ok()
    }

    /// Replaces represented atom state while retaining its ID and properties.
    pub fn replace_atom(&mut self, id: AtomId, atom: Atom) -> Result<Atom> {
        if self.atom(id)? == &atom {
            return Ok(atom);
        }
        Ok(std::mem::replace(&mut *self.atom_mut(id)?, atom))
    }

    /// Replaces a bond's endpoints and order, retaining its ID and properties.
    /// Invalid endpoints, self-bonds, and duplicate bonds leave state unchanged.
    /// Rewiring removes stereo assertions referencing the bond or changed
    /// endpoints, whose chemical neighborhoods have changed.
    /// Changing only the order removes assertions focused on this bond.
    pub fn replace_bond(&mut self, id: BondId, replacement: Bond) -> Result<Bond> {
        let previous = self.bond(id)?.clone();
        let (a, b) = replacement.endpoints();
        self.atom(a)?;
        self.atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.bond_between(a, b)?.is_some_and(|other| other != id) {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        if previous == replacement {
            return Ok(previous);
        }
        if (previous.a() == a && previous.b() == b) || (previous.a() == b && previous.b() == a) {
            self.bond_mut(id)?.set_order(replacement.order);
            return Ok(previous);
        }
        for atom in [previous.a(), previous.b()] {
            self.working.graph.adjacency[atom.index()].retain(|bond| *bond != id);
        }
        self.working.graph.adjacency[a.index()].push(id);
        self.working.graph.adjacency[b.index()].push(id);
        self.working.graph.bonds[id.index()] = Some(replacement);
        let invalid = self
            .stereo_elements()
            .filter_map(|(element_id, element)| {
                (element.references_bond(id)
                    || [previous.a(), previous.b(), a, b]
                        .into_iter()
                        .any(|atom| element.references_atom(atom)))
                .then_some(element_id)
            })
            .collect::<Vec<_>>();
        for element in invalid {
            self.remove_stereo_element(element)?;
        }
        self.working.clear_perception();
        self.working.properties.owner_mut().clear();
        Ok(previous)
    }

    /// Changes represented order, removing stereo assertions focused on this bond.
    /// Assigning the current order leaves chemistry, stereo, and annotations intact.
    pub fn set_bond_order(&mut self, id: BondId, order: BondOrder) -> Result<()> {
        if self.bond(id)?.order == order {
            return Ok(());
        }
        self.bond_mut(id)?.set_order(order);
        Ok(())
    }

    pub fn set_bond_endpoints(&mut self, id: BondId, a: AtomId, b: AtomId) -> Result<()> {
        self.replace_bond(id, Bond::new(a, b, self.bond(id)?.order))?;
        Ok(())
    }

    /// Deletes a set of live atoms and their incident bonds. IDs are validated
    /// before mutation; duplicates are ignored. Surviving IDs are never renumbered.
    pub fn delete_atoms(
        &mut self,
        atoms: impl IntoIterator<Item = AtomId>,
    ) -> Result<Vec<(AtomId, Atom)>> {
        let atoms = atoms.into_iter().collect::<BTreeSet<_>>();
        for &id in &atoms {
            self.atom(id)?;
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
    ) -> Result<Vec<(BondId, Bond)>> {
        let bonds = bonds.into_iter().collect::<BTreeSet<_>>();
        for &id in &bonds {
            self.bond(id)?;
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
    ) -> Result<Vec<(AtomId, Atom)>> {
        let retained = atoms.into_iter().collect::<BTreeSet<_>>();
        for &id in &retained {
            self.atom(id)?;
        }
        let removed = self
            .atom_ids()
            .filter(|id| !retained.contains(id))
            .collect::<Vec<_>>();
        self.delete_atoms(removed)
    }

    /// Clears graph, properties, and perception, restarting the local ID space.
    pub fn clear(&mut self) {
        *self = Self::new();
    }

    /// Replaces a relation group without changing its ID. All membership checks
    /// precede mutation, and members already in another group are rejected.
    pub fn replace_stereo_group(
        &mut self,
        id: StereoGroupId,
        replacement: StereoGroup,
    ) -> Result<StereoGroup> {
        let previous = self.stereo_group(id)?.clone();
        if replacement.members.is_empty()
            || replacement.members.iter().collect::<BTreeSet<_>>().len()
                != replacement.members.len()
        {
            return Err(MoleculeError::InvalidStereoReference(
                "stereo group requires nonempty unique members",
            ));
        }
        for &member in &replacement.members {
            if self
                .stereo_element(member)?
                .group
                .is_some_and(|group| group != id)
            {
                return Err(MoleculeError::InvalidStereoReference(
                    "stereo element already belongs to another group",
                ));
            }
        }
        if previous == replacement {
            return Ok(previous);
        }
        for &member in &previous.members {
            self.working.graph.stereo_elements[member.index()]
                .as_mut()
                .expect("validated group member")
                .group = None;
        }
        for &member in &replacement.members {
            self.working.graph.stereo_elements[member.index()]
                .as_mut()
                .expect("validated replacement member")
                .group = Some(id);
        }
        self.working.graph.stereo_groups[id.index()] = Some(replacement);
        self.working.properties.owner_mut().clear();
        self.working.invalidate_stereo();
        Ok(previous)
    }

    /// Appends a published fragment transactionally and returns semantic ID maps.
    /// Copies live atom/bond properties and represented stereo, including groups.
    /// Source owner properties and perception are not transferred. Successful
    /// append clears target owner properties and perception as a structural edit.
    /// A temporary disconnected result is allowed; connect it before finishing.
    /// This compound operation stages a clone of the target for rollback on error.
    pub fn append_molecule(&mut self, source: &Molecule) -> Result<MoleculeCorrespondence> {
        self.append_working(source)
    }

    // System edits may combine disconnected private drafts before partitioning.
    pub(crate) fn append_working(&mut self, source: &Molecule) -> Result<MoleculeCorrespondence> {
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
        // Published sources are compact, so each source row maps to one draft row.
        let property_error = |error| MoleculeError::Property(Box::new(error));
        let atom_rows = map.atoms.values().map(|id| id.index()).collect::<Vec<_>>();
        staged
            .working
            .properties
            .atoms_mut()
            .copy_rows_from(source.properties().atoms().raw(), &atom_rows)
            .map_err(property_error)?;
        let bond_rows = map.bonds.values().map(|id| id.index()).collect::<Vec<_>>();
        staged
            .working
            .properties
            .bonds_mut()
            .copy_rows_from(source.properties().bonds().raw(), &bond_rows)
            .map_err(property_error)?;
        let carrier = |value| match value {
            StereoCarrier::Atom(id) => StereoCarrier::Atom(map.atoms[&id]),
            other => other,
        };
        for (id, element) in source.stereo_elements() {
            let mut kind = element.kind.clone();
            match &mut kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    stereo.center = map.atoms[&stereo.center];
                    for c in &mut stereo.carriers {
                        *c = carrier(*c);
                    }
                }
                StereoElementKind::DoubleBond(stereo) => {
                    stereo.bond = map.bonds[&stereo.bond];
                    stereo.left = map.atoms[&stereo.left];
                    stereo.right = map.atoms[&stereo.right];
                    stereo.left_carrier = carrier(stereo.left_carrier);
                    stereo.right_carrier = carrier(stereo.right_carrier);
                }
                StereoElementKind::Axis(stereo) => {
                    stereo.axis = map.bonds[&stereo.axis];
                    for c in &mut stereo.carriers {
                        *c = carrier(*c);
                    }
                }
            }
            map.stereo_elements
                .insert(id, staged.add_stereo_element(StereoElement::new(kind))?);
        }
        for (id, group) in source.stereo_groups() {
            let mut group = group.clone();
            group.members = group
                .members
                .iter()
                .map(|member| map.stereo_elements[member])
                .collect();
            map.stereo_groups
                .insert(id, staged.add_stereo_group(group)?);
        }
        *self = staged;
        Ok(map)
    }
}
