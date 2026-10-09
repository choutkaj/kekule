use super::*;

impl TopologyEditor {
    /// Staged topology-level owner annotations. Entity annotations are
    /// addressed by edit handles because draft storage slots are private.
    pub fn owner_properties(&self) -> &OwnerProperties {
        self.properties.owner()
    }
    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        let previous = self.properties.owner_mut().insert(key.clone(), value)?;
        if previous.as_ref() != self.properties.owner().get(&key) {
            self.revision += 1;
        }
        Ok(previous)
    }
    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        let previous = self.properties.owner_mut().remove(key);
        if previous.is_some() {
            self.revision += 1;
        }
        previous
    }
    pub fn clear_properties(&mut self) {
        if !self.properties.owner().is_empty() {
            self.properties.owner_mut().clear();
            self.revision += 1;
        }
    }
    pub fn atom_property(
        &self,
        id: EditAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self
            .properties
            .atoms()
            .raw()
            .value(key, self.atom_slot(id)?)?)
    }
    pub fn bond_property(
        &self,
        id: EditBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self
            .properties
            .bonds()
            .raw()
            .value(key, self.bond_slot(id)?)?)
    }
    pub fn set_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.atom_slot(id)?;
        let previous = self.properties.atoms().raw().value(&key, slot)?;
        self.properties
            .atoms_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.atoms().raw().value(&key, slot)? {
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_bond_property(
        &mut self,
        id: EditBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.bond_slot(id)?;
        let previous = self.properties.bonds().raw().value(&key, slot)?;
        self.properties
            .bonds_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.bonds().raw().value(&key, slot)? {
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_atom_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditAtomId, Option<PropertyValue>)>,
    ) -> Result<(), TopologyEditError> {
        let mut staged = self.properties.atoms().raw().stage_column(&key);
        for (id, value) in values {
            staged.set_value(key.clone(), self.atom_slot(id)?, value)?;
        }
        if staged.get(&key) != self.properties.atoms().raw().get(&key) {
            self.properties.atoms_mut().commit_column(key, staged);
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_bond_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditBondId, Option<PropertyValue>)>,
    ) -> Result<(), TopologyEditError> {
        let mut staged = self.properties.bonds().raw().stage_column(&key);
        for (id, value) in values {
            staged.set_value(key.clone(), self.bond_slot(id)?, value)?;
        }
        if staged.get(&key) != self.properties.bonds().raw().get(&key) {
            self.properties.bonds_mut().commit_column(key, staged);
            self.revision += 1;
        }
        Ok(())
    }
    /// One static atom column in live atom-handle order.
    pub fn atom_property_column(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.properties.atoms().get(key)?;
        self.properties
            .atoms()
            .raw()
            .select_indices(&self.atoms.values().map(|l| l.slot).collect::<Vec<_>>())
            .expect("live slots are allocated rows")
            .remove(key)
    }
    /// One static bond column in live bond-handle order.
    pub fn bond_property_column(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.properties.bonds().get(key)?;
        self.properties
            .bonds()
            .raw()
            .select_indices(&self.bonds.values().map(|l| l.slot).collect::<Vec<_>>())
            .expect("live slots are allocated rows")
            .remove(key)
    }
    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        let previous = self.atom_property_column(&key);
        let slots = self.atoms.values().map(|l| l.slot).collect::<Vec<_>>();
        install_column(self.properties.atoms_mut(), &slots, key.clone(), column)?;
        if previous != self.atom_property_column(&key) {
            self.revision += 1;
        }
        Ok(previous)
    }
    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        let previous = self.bond_property_column(&key);
        let slots = self.bonds.values().map(|l| l.slot).collect::<Vec<_>>();
        install_column(self.properties.bonds_mut(), &slots, key.clone(), column)?;
        if previous != self.bond_property_column(&key) {
            self.revision += 1;
        }
        Ok(previous)
    }
    pub fn remove_atom_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.atom_property_column(key);
        self.properties.atoms_mut().remove(key);
        if previous.is_some() {
            self.revision += 1;
        }
        previous
    }
    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.bond_property_column(key);
        self.properties.bonds_mut().remove(key);
        if previous.is_some() {
            self.revision += 1;
        }
        previous
    }
    pub fn chain_property(
        &self,
        id: EditChainId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self
            .properties
            .chains()
            .raw()
            .value(key, self.chain(id)?.slot)?)
    }
    pub fn residue_property(
        &self,
        id: EditResidueId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self
            .properties
            .residues()
            .raw()
            .value(key, self.residue(id)?.slot)?)
    }
    pub fn atom_site_property(
        &self,
        id: EditAtomSiteId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self
            .properties
            .atom_sites()
            .raw()
            .value(key, self.atom_site(id)?.slot)?)
    }
    pub fn set_chain_property(
        &mut self,
        id: EditChainId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.chain(id)?.slot;
        let previous = self.properties.chains().raw().value(&key, slot)?;
        self.properties
            .chains_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.chains().raw().value(&key, slot)? {
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_residue_property(
        &mut self,
        id: EditResidueId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.residue(id)?.slot;
        let previous = self.properties.residues().raw().value(&key, slot)?;
        self.properties
            .residues_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.residues().raw().value(&key, slot)? {
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_atom_site_property(
        &mut self,
        id: EditAtomSiteId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.atom_site(id)?.slot;
        let previous = self.properties.atom_sites().raw().value(&key, slot)?;
        self.properties
            .atom_sites_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.atom_sites().raw().value(&key, slot)? {
            self.revision += 1;
        }
        Ok(())
    }
    /// Definition-scoped annotation for this occurrence, distinct from static
    /// system atom properties. Other occurrences of the definition are unchanged.
    pub fn set_definition_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let loc = self.atom_location(id)?;
        let mut staged = self.molecule(loc.group).edit();
        let row = loc.local;
        staged
            .properties_mut()
            .atoms_mut()
            .set_value(key.clone(), row, value)
            .map_err(|error| MoleculeError::Property(Box::new(error)))?;
        if staged.properties().atoms().value_ref(&key, row)
            != self
                .molecule(loc.group)
                .properties()
                .atoms()
                .value_ref(&key, row)
        {
            self.groups[loc.group].as_mut().unwrap().chemistry =
                GroupChemistry::Draft(Box::new(staged));
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_definition_bond_property(
        &mut self,
        id: EditBondId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let loc = self.bond_location(id)?;
        let mut staged = self.molecule(loc.group).edit();
        let row = loc.local;
        staged
            .properties_mut()
            .bonds_mut()
            .set_value(key.clone(), row, value)
            .map_err(|error| MoleculeError::Property(Box::new(error)))?;
        if staged.properties().bonds().value_ref(&key, row)
            != self
                .molecule(loc.group)
                .properties()
                .bonds()
                .value_ref(&key, row)
        {
            self.groups[loc.group].as_mut().unwrap().chemistry =
                GroupChemistry::Draft(Box::new(staged));
            self.revision += 1;
        }
        Ok(())
    }
}

fn install_column(
    table: &mut RawPropertyTable,
    live: &[usize],
    key: PropertyKey,
    column: PropertyColumn,
) -> Result<(), PropertyError> {
    let mut dense = RawPropertyTable::new(live.len());
    dense.insert(key.clone(), column)?;
    let mut slots = vec![None; table.len()];
    for (index, &slot) in live.iter().enumerate() {
        slots[slot] = Some(index);
    }
    let mut projected = dense.select_optional_indices(&slots)?;
    if let Some(column) = projected.remove(&key) {
        table.insert(key, column)?;
    } else {
        table.remove(&key);
    }
    Ok(())
}

#[cfg(test)]
mod allocation_tests {
    use super::*;
    use crate::properties::{PropertyColumn, RawPropertyTable};

    fn key(name: &str) -> PropertyKey {
        PropertyKey::new(name).unwrap()
    }

    fn int_column(table: &RawPropertyTable, key: &PropertyKey) -> *const Option<i64> {
        match table.get(key).unwrap() {
            PropertyColumn::Int(values) => values.as_ptr(),
            _ => panic!("expected an integer column"),
        }
    }

    // A batch stages and commits only its own column; unrelated columns keep
    // their allocations through successful and failed batches.
    #[test]
    fn property_batches_keep_unrelated_column_allocations() {
        let mut editor = TopologyEditor::from_topology(crate::smiles::to_topology("CCC").unwrap());
        let atoms = editor.atom_ids().collect::<Vec<_>>();
        let bonds = editor.bond_ids().collect::<Vec<_>>();
        let untouched = key("untouched");
        editor
            .set_atom_property(atoms[0], untouched.clone(), Some(PropertyValue::Int(1)))
            .unwrap();
        editor
            .set_bond_property(bonds[0], untouched.clone(), Some(PropertyValue::Int(2)))
            .unwrap();
        let atom_column = int_column(editor.properties.atoms().raw(), &untouched);
        let bond_column = int_column(editor.properties.bonds().raw(), &untouched);
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
            int_column(editor.properties.atoms().raw(), &untouched),
            atom_column
        );
        assert_eq!(
            int_column(editor.properties.bonds().raw(), &untouched),
            bond_column
        );
    }
}
