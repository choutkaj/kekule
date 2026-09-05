//! Compound editor operations. None exposes the unfinished graph as a Molecule.
use std::collections::{BTreeMap, BTreeSet};

use super::*;
use crate::properties::{PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue};

/// Correspondence from source fragment IDs to newly allocated editor IDs.
/// Deleted source slots have no entry. Existing target IDs remain unchanged.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MoleculeAppendMapping {
    atoms: BTreeMap<AtomId, AtomId>,
    bonds: BTreeMap<BondId, BondId>,
    stereo_elements: BTreeMap<StereoElementId, StereoElementId>,
    stereo_groups: BTreeMap<StereoGroupId, StereoGroupId>,
}

impl MoleculeAppendMapping {
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
        Ok(std::mem::replace(&mut *self.atom_mut(id)?, atom))
    }

    /// Replaces a bond's endpoints and order, retaining its ID and properties.
    /// Invalid endpoints, self-bonds, and duplicate bonds leave state unchanged.
    /// Rewiring removes stereo assertions referencing the bond or changed
    /// endpoints, whose chemical neighborhoods have changed.
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
        self.working.properties.clear_owner();
        Ok(previous)
    }

    pub fn set_bond_order(&mut self, id: BondId, order: BondOrder) -> Result<()> {
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

    /// Applies a batch of values for one atom property transactionally. Repeated
    /// IDs are applied in input order; `None` removes a value. Only the property
    /// table is staged, without cloning graph or perception state.
    pub fn set_atom_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (AtomId, Option<PropertyValue>)>,
    ) -> Result<()> {
        let mut table = self.atom_properties().clone();
        for (id, value) in values {
            self.atom(id)?;
            table
                .set_value(key.clone(), id.index(), value)
                .map_err(|e| MoleculeError::Property(Box::new(e)))?;
        }
        *self.working.properties.atoms_mut() = table;
        Ok(())
    }

    /// Bond counterpart of [`Self::set_atom_properties`].
    pub fn set_bond_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (BondId, Option<PropertyValue>)>,
    ) -> Result<()> {
        let mut table = self.bond_properties().clone();
        for (id, value) in values {
            self.bond(id)?;
            table
                .set_value(key.clone(), id.index(), value)
                .map_err(|e| MoleculeError::Property(Box::new(e)))?;
        }
        *self.working.properties.bonds_mut() = table;
        Ok(())
    }

    pub fn remove_atom_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.working.properties.atoms_mut().remove(key)
    }

    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.working.properties.bonds_mut().remove(key)
    }

    /// Replaces a complete column in current live [`Self::atom_ids`] order.
    /// Deleted slots are filled with missing values automatically. Counts,
    /// values, and units are checked before mutation; unrelated columns survive.
    pub fn set_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<()> {
        let slots = live_column_slots(self.working.graph.atoms.iter().map(Option::is_some));
        set_live_column(self.working.properties.atoms_mut(), key, column, &slots)
    }

    /// Bond counterpart of [`Self::set_atom_property_column`], in live bond order.
    pub fn set_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<()> {
        let slots = live_column_slots(self.working.graph.bonds.iter().map(Option::is_some));
        set_live_column(self.working.properties.bonds_mut(), key, column, &slots)
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
        self.working.properties.clear_owner();
        self.working.invalidate_stereo();
        Ok(previous)
    }

    /// Appends a published fragment transactionally and returns semantic ID maps.
    /// Copies live atom/bond properties and represented stereo, including groups.
    /// Source owner properties and perception are not transferred. Successful
    /// append clears target owner properties and perception as a structural edit.
    /// A temporary disconnected result is allowed; connect it before finishing.
    /// This compound operation stages a clone of the target for rollback on error.
    pub fn append_molecule(&mut self, source: &Molecule) -> Result<MoleculeAppendMapping> {
        let mut staged = self.clone();
        let mut map = MoleculeAppendMapping::default();
        for (id, atom) in source.atoms() {
            map.atoms.insert(id, staged.add_atom(atom.clone())?);
        }
        for (id, bond) in source.bonds() {
            map.bonds.insert(
                id,
                staged.add_bond(map.atoms[&bond.a()], map.atoms[&bond.b()], bond.order)?,
            );
        }
        for (old, new) in &map.atoms {
            for (key, _) in source.atom_properties().iter() {
                staged.set_atom_property(*new, key.clone(), source.atom_property(*old, key)?)?;
            }
        }
        for (old, new) in &map.bonds {
            for (key, _) in source.bond_properties().iter() {
                staged.set_bond_property(*new, key.clone(), source.bond_property(*old, key)?)?;
            }
        }
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

fn live_column_slots(live: impl Iterator<Item = bool>) -> Vec<Option<usize>> {
    let mut next = 0;
    live.map(|live| {
        live.then(|| {
            let index = next;
            next += 1;
            index
        })
    })
    .collect()
}

fn set_live_column(
    table: &mut PropertyTable,
    key: PropertyKey,
    column: PropertyColumn,
    slots: &[Option<usize>],
) -> Result<()> {
    let mut dense = PropertyTable::new(slots.iter().flatten().count());
    let map_error = |error: PropertyError| MoleculeError::Property(Box::new(error));
    dense.insert(key.clone(), column).map_err(map_error)?;
    let mut projected = dense.select_optional_indices(slots).map_err(map_error)?;
    if let Some(column) = projected.remove(&key) {
        table.insert(key, column).map_err(map_error)?;
    } else {
        table.remove(&key);
    }
    Ok(())
}
