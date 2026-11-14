use std::marker::PhantomData;

use super::*;
use crate::properties::{rows::RowKey, HandleSlots, LiveRows, PropertyTableMut, TableRow};

macro_rules! handle_rows {
    ($($handle:ident => $storage:ty, $invalid:ident;)*) => {$(
        impl RowKey for $handle {
            type Storage = $storage;

            fn slot(self, rows: LiveRows<'_>) -> Result<usize, TopologyEditError> {
                rows.handle_slot(self.position())
                    .ok_or(TopologyEditError::$invalid(self))
            }
        }
        impl TableRow for $handle {
            type Error = TopologyEditError;
        }
    )*};
}

handle_rows! {
    EditAtomId => crate::topology::TopologyAtomIndex, InvalidAtom;
    EditBondId => crate::topology::TopologyBondIndex, InvalidBond;
    EditChainId => ChainId, InvalidChain;
    EditResidueId => ResidueId, InvalidResidue;
    EditAtomSiteId => AtomSiteId, InvalidAtomSite;
}

/// Read access to one draft property table, addressed by editing handle
/// `H` (`EditAtomId`, `EditBondId`, `EditChainId`, `EditResidueId`, or
/// `EditAtomSiteId`).
///
/// Columns read in live-handle order, which is the dense order of the
/// published result. A deleted or foreign handle fails with the matching
/// [`TopologyEditError`].
pub struct EditPropertyTable<'a, H: TableRow<Error = TopologyEditError>> {
    table: &'a RawPropertyTable,
    handles: &'a dyn HandleSlots,
    handle: PhantomData<fn() -> H>,
}

impl<H: TableRow<Error = TopologyEditError>> fmt::Debug for EditPropertyTable<'_, H> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.table.fmt(formatter)
    }
}

impl<'a, H: TableRow<Error = TopologyEditError>> EditPropertyTable<'a, H> {
    pub(crate) fn new(table: &'a RawPropertyTable, handles: &'a dyn HandleSlots) -> Self {
        Self {
            table,
            handles,
            handle: PhantomData,
        }
    }

    /// Whether any column is present.
    pub fn has_data(&self) -> bool {
        self.table.has_data()
    }

    /// Column keys in deterministic order.
    pub fn keys(&self) -> impl ExactSizeIterator<Item = &'a PropertyKey> + 'a {
        self.table.iter().map(|(key, _)| key)
    }

    /// One column in live-handle order.
    pub fn get(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        crate::properties::live_column(self.table, self.handles, key)
    }

    /// Reads one cell of a live row.
    pub fn value(
        &self,
        key: &PropertyKey,
        row: H,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        let slot = row.slot(LiveRows::Handles(self.handles))?;
        Ok(self.table.value(key, slot)?)
    }
}

/// The reads of [`EditPropertyTable`] on a draft table's guard.
impl<H: TableRow<Error = TopologyEditError>> PropertyTableMut<'_, H> {
    fn view(&self) -> EditPropertyTable<'_, H> {
        let (table, live) = self.parts();
        let LiveRows::Handles(handles) = live else {
            unreachable!("draft tables address rows by editing handle")
        };
        EditPropertyTable::new(table, handles)
    }

    /// Whether any column is present.
    pub fn has_data(&self) -> bool {
        self.view().has_data()
    }

    /// Column keys in deterministic order.
    pub fn keys(&self) -> impl ExactSizeIterator<Item = &PropertyKey> {
        self.parts().0.iter().map(|(key, _)| key)
    }

    /// One column in live-handle order.
    pub fn get(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.view().get(key)
    }

    /// Reads one cell of a live row.
    pub fn value(
        &self,
        key: &PropertyKey,
        row: H,
    ) -> Result<Option<PropertyValue>, TopologyEditError> {
        self.view().value(key, row)
    }
}

/// Read access to a draft's static topology annotations, addressed by
/// editing handle.
#[derive(Debug, Clone, Copy)]
pub struct EditProperties<'a> {
    editor: &'a TopologyEditor,
}

impl<'a> EditProperties<'a> {
    /// Topology-level owner annotations.
    pub fn owner(&self) -> &'a OwnerProperties {
        self.editor.properties.owner()
    }

    /// Atom annotations.
    pub fn atoms(&self) -> EditPropertyTable<'a, EditAtomId> {
        EditPropertyTable::new(self.editor.properties.atoms().raw(), &self.editor.atoms)
    }

    /// Bond annotations.
    pub fn bonds(&self) -> EditPropertyTable<'a, EditBondId> {
        EditPropertyTable::new(self.editor.properties.bonds().raw(), &self.editor.bonds)
    }

    /// Chain annotations.
    pub fn chains(&self) -> EditPropertyTable<'a, EditChainId> {
        EditPropertyTable::new(
            self.editor.properties.chains().raw(),
            &self.editor.hierarchy.chains,
        )
    }

    /// Residue annotations.
    pub fn residues(&self) -> EditPropertyTable<'a, EditResidueId> {
        EditPropertyTable::new(
            self.editor.properties.residues().raw(),
            &self.editor.hierarchy.residues,
        )
    }

    /// Atom-site annotations.
    pub fn atom_sites(&self) -> EditPropertyTable<'a, EditAtomSiteId> {
        EditPropertyTable::new(
            self.editor.properties.atom_sites().raw(),
            &self.editor.hierarchy.sites,
        )
    }
}

/// Mutable access to a draft's static topology annotations, addressed by
/// editing handle; see [`PropertyTableMut`] for the table methods.
///
/// ```
/// use kekule::{properties::{PropertyKey, PropertyValue}, smiles};
/// let mut editor = smiles::to_topology("CO")?.to_editor();
/// let oxygen = editor.atom_ids().nth(1).unwrap();
/// let key = PropertyKey::new("label")?;
/// editor
///     .properties_mut()
///     .atoms_mut()
///     .set_value(key.clone(), oxygen, Some(PropertyValue::String("O1".into())))?;
/// let topology = editor.finish()?;
/// assert_eq!(
///     topology.atoms().nth(1).unwrap().property(&key),
///     Some(PropertyValue::String("O1".into()))
/// );
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct EditPropertiesMut<'a> {
    editor: &'a mut TopologyEditor,
}

impl EditPropertiesMut<'_> {
    /// Topology-level owner annotations.
    pub fn owner(&self) -> &OwnerProperties {
        self.editor.properties.owner()
    }

    /// Topology-level owner annotations. Publication compares them with the
    /// source's, so an unchanged draft still publishes its source.
    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        self.editor.properties.owner_mut()
    }

    /// Atom annotations.
    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, EditAtomId> {
        let editor = &mut *self.editor;
        PropertyTableMut::new(
            editor.properties.atom_table_mut(),
            LiveRows::Handles(&editor.atoms),
            Some(&mut editor.revision),
        )
    }

    /// Bond annotations.
    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, EditBondId> {
        let editor = &mut *self.editor;
        PropertyTableMut::new(
            editor.properties.bond_table_mut(),
            LiveRows::Handles(&editor.bonds),
            Some(&mut editor.revision),
        )
    }

    /// Chain annotations.
    pub fn chains_mut(&mut self) -> PropertyTableMut<'_, EditChainId> {
        let editor = &mut *self.editor;
        PropertyTableMut::new(
            editor.properties.chain_table_mut(),
            LiveRows::Handles(&editor.hierarchy.chains),
            Some(&mut editor.revision),
        )
    }

    /// Residue annotations.
    pub fn residues_mut(&mut self) -> PropertyTableMut<'_, EditResidueId> {
        let editor = &mut *self.editor;
        PropertyTableMut::new(
            editor.properties.residue_table_mut(),
            LiveRows::Handles(&editor.hierarchy.residues),
            Some(&mut editor.revision),
        )
    }

    /// Atom-site annotations.
    pub fn atom_sites_mut(&mut self) -> PropertyTableMut<'_, EditAtomSiteId> {
        let editor = &mut *self.editor;
        PropertyTableMut::new(
            editor.properties.atom_site_table_mut(),
            LiveRows::Handles(&editor.hierarchy.sites),
            Some(&mut editor.revision),
        )
    }

    /// Sets a definition-scoped atom annotation for this instance, distinct
    /// from static topology atom properties. Other instances of the
    /// definition are unchanged.
    pub fn set_definition_atom_value(
        &mut self,
        key: PropertyKey,
        atom: EditAtomId,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let editor = &mut *self.editor;
        let loc = editor.atom_location(atom)?;
        let mut staged = editor.molecule(loc.group).to_editor();
        let row = loc.local;
        staged
            .properties_mut()
            .atoms_mut()
            .set_value(key.clone(), row, value)
            .map_err(MoleculeError::Property)?;
        if staged.properties().atoms().value_ref(&key, row)
            != editor
                .molecule(loc.group)
                .properties()
                .atoms()
                .value_ref(&key, row)
        {
            editor.group_mut(loc.group).chemistry = GroupChemistry::Draft(Box::new(staged));
            editor.revision += 1;
        }
        Ok(())
    }

    /// Sets a definition-scoped bond annotation for this instance; see
    /// [`Self::set_definition_atom_value`].
    pub fn set_definition_bond_value(
        &mut self,
        key: PropertyKey,
        bond: EditBondId,
        value: Option<PropertyValue>,
    ) -> Result<(), TopologyEditError> {
        let editor = &mut *self.editor;
        let loc = editor.bond_location(bond)?;
        let mut staged = editor.molecule(loc.group).to_editor();
        let row = loc.local;
        staged
            .properties_mut()
            .bonds_mut()
            .set_value(key.clone(), row, value)
            .map_err(MoleculeError::Property)?;
        if staged.properties().bonds().value_ref(&key, row)
            != editor
                .molecule(loc.group)
                .properties()
                .bonds()
                .value_ref(&key, row)
        {
            editor.group_mut(loc.group).chemistry = GroupChemistry::Draft(Box::new(staged));
            editor.revision += 1;
        }
        Ok(())
    }
}

impl TopologyEditor {
    /// Staged static topology annotations, addressed by editing handle.
    pub fn properties(&self) -> EditProperties<'_> {
        EditProperties { editor: self }
    }

    /// Mutable staged static topology annotations, addressed by editing
    /// handle.
    pub fn properties_mut(&mut self) -> EditPropertiesMut<'_> {
        EditPropertiesMut { editor: self }
    }

    /// The live rows of the draft's atoms, for structure tables kept beside
    /// this draft.
    pub(crate) fn atom_rows(&self) -> &dyn HandleSlots {
        &self.atoms
    }

    /// The live rows of the draft's bonds.
    pub(crate) fn bond_rows(&self) -> &dyn HandleSlots {
        &self.bonds
    }
}

#[cfg(test)]
mod tests;
