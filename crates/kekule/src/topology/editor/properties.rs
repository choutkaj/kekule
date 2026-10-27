use std::marker::PhantomData;

use super::*;

/// The live rows of one draft property domain: which handles are live and
/// which storage slot each one addresses.
#[derive(Debug, Clone, Copy)]
pub(crate) enum EditRows<'a> {
    Atoms(&'a DenseMap<EditAtomId, Location<AtomId>>),
    Bonds(&'a DenseMap<EditBondId, Location<BondId>>),
    Chains(&'a DenseMap<EditChainId, EditChain>),
    Residues(&'a DenseMap<EditResidueId, EditResidue>),
    AtomSites(&'a DenseMap<EditAtomSiteId, EditAtomSite>),
}

impl EditRows<'_> {
    /// Live slots in handle order, which is publication order.
    fn live_slots(self) -> Vec<usize> {
        match self {
            Self::Atoms(rows) => rows.values().map(|row| row.slot).collect(),
            Self::Bonds(rows) => rows.values().map(|row| row.slot).collect(),
            Self::Chains(rows) => rows.values().map(|row| row.slot).collect(),
            Self::Residues(rows) => rows.values().map(|row| row.slot).collect(),
            Self::AtomSites(rows) => rows.values().map(|row| row.slot).collect(),
        }
    }
}

/// A draft handle that addresses one row of an editor property table.
trait EditRow: Copy {
    fn slot(self, rows: EditRows<'_>) -> Result<usize, TopologyEditError>;
}

macro_rules! edit_rows {
    ($($handle:ty => $variant:ident, $invalid:ident;)*) => {$(
        impl EditRow for $handle {
            fn slot(self, rows: EditRows<'_>) -> Result<usize, TopologyEditError> {
                match rows {
                    EditRows::$variant(rows) => rows
                        .get(&self)
                        .map(|row| row.slot)
                        .ok_or(TopologyEditError::$invalid(self)),
                    _ => unreachable!("a handle reads only its own domain's rows"),
                }
            }
        }
    )*};
}

edit_rows! {
    EditAtomId => Atoms, InvalidAtom;
    EditBondId => Bonds, InvalidBond;
    EditChainId => Chains, InvalidChain;
    EditResidueId => Residues, InvalidResidue;
    EditAtomSiteId => AtomSites, InvalidAtomSite;
}

/// One property column of the live rows, in handle order.
fn live_column(
    table: &RawPropertyTable,
    rows: EditRows<'_>,
    key: &PropertyKey,
) -> Option<PropertyColumn> {
    table.get(key)?;
    table
        .select_indices(&rows.live_slots())
        .expect("live slots are allocated rows")
        .remove(key)
}

/// Read access to one draft property table, addressed by editing handle
/// `H` (`EditAtomId`, `EditBondId`, `EditChainId`, `EditResidueId`, or
/// `EditAtomSiteId`).
///
/// Columns read in live-handle order, which is the dense order of the
/// published result.
pub struct EditPropertyTable<'a, H> {
    table: &'a RawPropertyTable,
    rows: EditRows<'a>,
    handle: PhantomData<fn() -> H>,
}

impl<H> fmt::Debug for EditPropertyTable<'_, H> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.table.fmt(formatter)
    }
}

/// Mutable access to one draft property table, addressed by editing handle
/// `H`, with the method set of [`PropertyTableMut`](crate::properties::PropertyTableMut).
///
/// Every mutation is transactional. A deleted or foreign handle fails with
/// the matching [`TopologyEditError`] (`InvalidAtom`, `InvalidBond`, ...).
/// Columns are read and inserted in live-handle order.
pub struct EditPropertyTableMut<'a, H> {
    table: &'a mut RawPropertyTable,
    rows: EditRows<'a>,
    revision: Option<&'a mut u64>,
    handle: PhantomData<fn() -> H>,
}

impl<H> fmt::Debug for EditPropertyTableMut<'_, H> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.table.fmt(formatter)
    }
}

impl<'a, H> EditPropertyTable<'a, H> {
    pub(crate) const fn new(table: &'a RawPropertyTable, rows: EditRows<'a>) -> Self {
        Self {
            table,
            rows,
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
        live_column(self.table, self.rows, key)
    }
}

impl<'a, H> EditPropertyTableMut<'a, H> {
    pub(crate) fn new(
        table: &'a mut RawPropertyTable,
        rows: EditRows<'a>,
        revision: Option<&'a mut u64>,
    ) -> Self {
        Self {
            table,
            rows,
            revision,
            handle: PhantomData,
        }
    }

    fn changed(&mut self) {
        if let Some(revision) = self.revision.as_deref_mut() {
            *revision += 1;
        }
    }

    /// Whether any column is present.
    pub fn has_data(&self) -> bool {
        self.table.has_data()
    }

    /// Column keys in deterministic order.
    pub fn keys(&self) -> impl ExactSizeIterator<Item = &PropertyKey> {
        self.table.iter().map(|(key, _)| key)
    }

    /// One column in live-handle order.
    pub fn get(&self, key: &PropertyKey) -> Option<PropertyColumn> {
        live_column(self.table, self.rows, key)
    }

    /// Inserts or replaces a complete column given in live-handle order,
    /// returning the previous one.
    pub fn insert(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, TopologyEditError> {
        let previous = self.get(&key);
        install_column(self.table, &self.rows.live_slots(), key.clone(), column)?;
        if previous != self.get(&key) {
            self.changed();
        }
        Ok(previous)
    }

    /// Removes a column, returning it in live-handle order.
    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = self.get(key);
        self.table.remove(key);
        if previous.is_some() {
            self.changed();
        }
        previous
    }

    /// Removes every column.
    pub fn clear(&mut self) {
        if self.table.has_data() {
            self.table.clear_columns();
            self.changed();
        }
    }
}

macro_rules! edit_table_rows {
    ($($handle:ty),*) => {$(
        impl EditPropertyTable<'_, $handle> {
            /// Reads one cell of a live row.
            pub fn value(
                &self,
                key: &PropertyKey,
                row: $handle,
            ) -> Result<Option<PropertyValue>, TopologyEditError> {
                Ok(self.table.value(key, row.slot(self.rows)?)?)
            }
        }

        impl EditPropertyTableMut<'_, $handle> {
            /// Reads one cell of a live row.
            pub fn value(
                &self,
                key: &PropertyKey,
                row: $handle,
            ) -> Result<Option<PropertyValue>, TopologyEditError> {
                Ok(self.table.value(key, row.slot(self.rows)?)?)
            }

            /// Sets or clears one cell, removing a column that becomes
            /// all-missing.
            pub fn set_value(
                &mut self,
                key: PropertyKey,
                row: $handle,
                value: Option<PropertyValue>,
            ) -> Result<(), TopologyEditError> {
                let slot = row.slot(self.rows)?;
                let previous = self.table.value(&key, slot)?;
                self.table.set_value(key.clone(), slot, value)?;
                if previous != self.table.value(&key, slot)? {
                    self.changed();
                }
                Ok(())
            }

/// Clears one cell; the same as [`Self::set_value`] with `None`.
            pub fn clear_value(
                &mut self,
                key: PropertyKey,
                row: $handle,
            ) -> Result<(), TopologyEditError> {
                self.set_value(key, row, None)
            }

            /// Applies several cell updates to one column transactionally;
            /// repeated handles use the last value. Other columns are never
            /// copied.
            pub fn set_values(
                &mut self,
                key: PropertyKey,
                values: impl IntoIterator<Item = ($handle, Option<PropertyValue>)>,
            ) -> Result<(), TopologyEditError> {
                let mut staged = self.table.stage_column(&key);
                for (row, value) in values {
                    staged.set_value(key.clone(), row.slot(self.rows)?, value)?;
                }
                // An equal column may still differ in signed zeros or unit
                // symbols; a table without publication tracking keeps them.
                let changed = staged.get(&key) != self.table.get(&key);
                if changed || self.revision.is_none() {
                    self.table.commit_column(key, staged);
                }
                if changed {
                    self.changed();
                }
                Ok(())
            }
        }
    )*};
}

edit_table_rows!(
    EditAtomId,
    EditBondId,
    EditChainId,
    EditResidueId,
    EditAtomSiteId
);

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

    /// Atom annotations, addressed by editing handle.
    pub fn atoms(&self) -> EditPropertyTable<'a, EditAtomId> {
        EditPropertyTable::new(
            self.editor.properties.atoms().raw(),
            EditRows::Atoms(&self.editor.atoms),
        )
    }

    /// Bond annotations, addressed by editing handle.
    pub fn bonds(&self) -> EditPropertyTable<'a, EditBondId> {
        EditPropertyTable::new(
            self.editor.properties.bonds().raw(),
            EditRows::Bonds(&self.editor.bonds),
        )
    }

    /// Chain annotations, addressed by editing handle.
    pub fn chains(&self) -> EditPropertyTable<'a, EditChainId> {
        EditPropertyTable::new(
            self.editor.properties.chains().raw(),
            EditRows::Chains(&self.editor.hierarchy.chains),
        )
    }

    /// Residue annotations, addressed by editing handle.
    pub fn residues(&self) -> EditPropertyTable<'a, EditResidueId> {
        EditPropertyTable::new(
            self.editor.properties.residues().raw(),
            EditRows::Residues(&self.editor.hierarchy.residues),
        )
    }

    /// Atom-site annotations, addressed by editing handle.
    pub fn atom_sites(&self) -> EditPropertyTable<'a, EditAtomSiteId> {
        EditPropertyTable::new(
            self.editor.properties.atom_sites().raw(),
            EditRows::AtomSites(&self.editor.hierarchy.sites),
        )
    }
}

/// Mutable access to a draft's static topology annotations, addressed by
/// editing handle; see [`EditPropertyTableMut`] for the table methods.
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

    /// Atom annotations, addressed by editing handle.
    pub fn atoms_mut(&mut self) -> EditPropertyTableMut<'_, EditAtomId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.properties.atoms_mut(),
            EditRows::Atoms(&editor.atoms),
            Some(&mut editor.revision),
        )
    }

    /// Bond annotations, addressed by editing handle.
    pub fn bonds_mut(&mut self) -> EditPropertyTableMut<'_, EditBondId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.properties.bonds_mut(),
            EditRows::Bonds(&editor.bonds),
            Some(&mut editor.revision),
        )
    }

    /// Chain annotations, addressed by editing handle.
    pub fn chains_mut(&mut self) -> EditPropertyTableMut<'_, EditChainId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.properties.chains_mut(),
            EditRows::Chains(&editor.hierarchy.chains),
            Some(&mut editor.revision),
        )
    }

    /// Residue annotations, addressed by editing handle.
    pub fn residues_mut(&mut self) -> EditPropertyTableMut<'_, EditResidueId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.properties.residues_mut(),
            EditRows::Residues(&editor.hierarchy.residues),
            Some(&mut editor.revision),
        )
    }

    /// Atom-site annotations, addressed by editing handle.
    pub fn atom_sites_mut(&mut self) -> EditPropertyTableMut<'_, EditAtomSiteId> {
        let editor = &mut *self.editor;
        EditPropertyTableMut::new(
            editor.properties.atom_sites_mut(),
            EditRows::AtomSites(&editor.hierarchy.sites),
            Some(&mut editor.revision),
        )
    }

    /// Sets a definition-scoped atom annotation for this occurrence, distinct
    /// from static topology atom properties. Other occurrences of the
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
            .map_err(|error| MoleculeError::Property(Box::new(error)))?;
        if staged.properties().atoms().value_ref(&key, row)
            != editor
                .molecule(loc.group)
                .properties()
                .atoms()
                .value_ref(&key, row)
        {
            editor.groups[loc.group].as_mut().unwrap().chemistry =
                GroupChemistry::Draft(Box::new(staged));
            editor.revision += 1;
        }
        Ok(())
    }

    /// Sets a definition-scoped bond annotation for this occurrence; see
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
            .map_err(|error| MoleculeError::Property(Box::new(error)))?;
        if staged.properties().bonds().value_ref(&key, row)
            != editor
                .molecule(loc.group)
                .properties()
                .bonds()
                .value_ref(&key, row)
        {
            editor.groups[loc.group].as_mut().unwrap().chemistry =
                GroupChemistry::Draft(Box::new(staged));
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
    pub(crate) fn atom_rows(&self) -> EditRows<'_> {
        EditRows::Atoms(&self.atoms)
    }

    /// The live rows of the draft's bonds.
    pub(crate) fn bond_rows(&self) -> EditRows<'_> {
        EditRows::Bonds(&self.bonds)
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
        let mut properties = editor.properties_mut();
        properties
            .atoms_mut()
            .set_value(untouched.clone(), atoms[0], Some(PropertyValue::Int(1)))
            .unwrap();
        properties
            .bonds_mut()
            .set_value(untouched.clone(), bonds[0], Some(PropertyValue::Int(2)))
            .unwrap();
        let atom_column = int_column(editor.properties.atoms().raw(), &untouched);
        let bond_column = int_column(editor.properties.bonds().raw(), &untouched);
        let mut properties = editor.properties_mut();
        properties
            .atoms_mut()
            .set_values(
                key("edited"),
                [
                    (atoms[0], Some(PropertyValue::Int(3))),
                    (atoms[1], Some(PropertyValue::Int(4))),
                ],
            )
            .unwrap();
        assert!(properties
            .atoms_mut()
            .set_values(
                key("edited"),
                [(atoms[1], Some(PropertyValue::String("bad type".into())))],
            )
            .is_err());
        properties
            .bonds_mut()
            .set_values(key("edited"), [(bonds[1], Some(PropertyValue::Int(5)))])
            .unwrap();
        properties
            .bonds_mut()
            .set_values(key("edited"), [(bonds[1], None)])
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
