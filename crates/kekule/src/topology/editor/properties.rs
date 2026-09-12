use super::*;

impl TopologyEditor {
    /// Static properties in stable draft-slot order. Use semantic cell access or
    /// live-order column helpers rather than interpreting slots as live indices.
    pub fn properties(&self) -> &Properties {
        &self.properties
    }
    pub fn atom_properties(&self) -> &PropertyTable {
        self.properties.atoms()
    }
    pub fn bond_properties(&self) -> &PropertyTable {
        self.properties.bonds()
    }
    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        let previous = self.properties.insert(key.clone(), value)?;
        if previous.as_ref() != self.properties.get(&key) {
            self.revision += 1;
        }
        Ok(previous)
    }
    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        let previous = self.properties.remove(key);
        if previous.is_some() {
            self.revision += 1;
        }
        previous
    }
    pub fn clear_properties(&mut self) {
        if !self.properties.owner_is_empty() {
            self.properties.clear_owner();
            self.revision += 1;
        }
    }
    pub fn atom_property(
        &self,
        id: EditAtomId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self.properties.atoms().value(key, self.atom_slot(id)?)?)
    }
    pub fn bond_property(
        &self,
        id: EditBondId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self.properties.bonds().value(key, self.bond_slot(id)?)?)
    }
    pub fn set_atom_property(
        &mut self,
        id: EditAtomId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.atom_slot(id)?;
        let previous = self.properties.atoms().value(&key, slot)?;
        self.properties
            .atoms_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.atoms().value(&key, slot)? {
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
        let previous = self.properties.bonds().value(&key, slot)?;
        self.properties
            .bonds_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.bonds().value(&key, slot)? {
            self.revision += 1;
        }
        Ok(())
    }
    pub fn set_atom_properties(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (EditAtomId, Option<PropertyValue>)>,
    ) -> Result<(), TopologyEditError> {
        let mut staged = self.properties.atoms().stage_column(&key);
        for (id, value) in values {
            staged.set_value(key.clone(), self.atom_slot(id)?, value)?;
        }
        if staged.get(&key) != self.properties.atoms().get(&key) {
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
        let mut staged = self.properties.bonds().stage_column(&key);
        for (id, value) in values {
            staged.set_value(key.clone(), self.bond_slot(id)?, value)?;
        }
        if staged.get(&key) != self.properties.bonds().get(&key) {
            self.properties.bonds_mut().commit_column(key, staged);
            self.revision += 1;
        }
        Ok(())
    }
    pub fn atom_property_column(
        &self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        Ok(self
            .atom_properties()
            .select_indices(&self.atoms.values().map(|l| l.slot).collect::<Vec<_>>())?
            .remove(key))
    }
    pub fn bond_property_column(
        &self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        Ok(self
            .bond_properties()
            .select_indices(&self.bonds.values().map(|l| l.slot).collect::<Vec<_>>())?
            .remove(key))
    }
    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        let previous = self.atom_property_column(&key)?;
        let slots = self.atoms.values().map(|l| l.slot).collect::<Vec<_>>();
        install_column(self.properties.atoms_mut(), &slots, key.clone(), column)?;
        if previous != self.atom_property_column(&key)? {
            self.revision += 1;
        }
        Ok(previous)
    }
    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        let previous = self.bond_property_column(&key)?;
        let slots = self.bonds.values().map(|l| l.slot).collect::<Vec<_>>();
        install_column(self.properties.bonds_mut(), &slots, key.clone(), column)?;
        if previous != self.bond_property_column(&key)? {
            self.revision += 1;
        }
        Ok(previous)
    }
    pub fn remove_atom_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.atom_property_column(key).expect("live property slots");
        self.properties.atoms_mut().remove(key);
        if previous.is_some() {
            self.revision += 1;
        }
        previous
    }
    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.bond_property_column(key).expect("live property slots");
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
        Ok(self.properties.chains().value(key, self.chain(id)?.slot)?)
    }
    pub fn residue_property(
        &self,
        id: EditResidueId,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        Ok(self
            .properties
            .residues()
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
            .value(key, self.atom_site(id)?.slot)?)
    }
    pub fn set_chain_property(
        &mut self,
        id: EditChainId,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let slot = self.chain(id)?.slot;
        let previous = self.properties.chains().value(&key, slot)?;
        self.properties
            .chains_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.chains().value(&key, slot)? {
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
        let previous = self.properties.residues().value(&key, slot)?;
        self.properties
            .residues_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.residues().value(&key, slot)? {
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
        let previous = self.properties.atom_sites().value(&key, slot)?;
        self.properties
            .atom_sites_mut()
            .set_value(key.clone(), slot, value)?;
        if previous != self.properties.atom_sites().value(&key, slot)? {
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
        staged.set_atom_property(loc.local, key.clone(), value)?;
        if staged.atom_property(loc.local, &key)?
            != self.molecule(loc.group).atom_property(loc.local, &key)?
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
        staged.set_bond_property(loc.local, key.clone(), value)?;
        if staged.bond_property(loc.local, &key)?
            != self.molecule(loc.group).bond_property(loc.local, &key)?
        {
            self.groups[loc.group].as_mut().unwrap().chemistry =
                GroupChemistry::Draft(Box::new(staged));
            self.revision += 1;
        }
        Ok(())
    }
}

fn install_column(
    table: &mut PropertyTable,
    live: &[usize],
    key: PropertyKey,
    column: PropertyColumn,
) -> Result<(), PropertyError> {
    let mut dense = PropertyTable::new(live.len());
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
