//! Unified scalar and columnar annotations for canonical Kekule objects.
//!
//! Property scope follows the object that owns it, and each scope has its own
//! type: definition-invariant annotations are [`MoleculeProperties`] of a
//! [`crate::core::Molecule`], system annotations are [`TopologyProperties`] of a
//! [`crate::topology::Topology`], and coordinate-dependent annotations are
//! [`RealizationProperties`] of a [`crate::structure::Conformation`]. Each scope
//! holds [`OwnerProperties`] plus one [`PropertyTable`] per entity domain,
//! addressed by that domain's identifier ([`PropertyRow`]).
//!
//! Owners hand out read access through `properties()` and length-preserving
//! mutable access through `properties_mut()`. Row counts are fixed by the
//! owner and can never be changed through a property guard.
//!
//! Properties are annotations, not represented graph chemistry. Changing a
//! generic property does not change molecular identity or trigger chemical
//! perception.

use std::collections::BTreeMap;
use std::fmt;
use std::marker::PhantomData;
use std::str::FromStr;

use crate::core::{AtomId, BondId};
use crate::topology::{
    AtomSiteId, ChainId, MoleculeInstanceId, ResidueId, TopologyAtomIndex, TopologyBondIndex,
};
use crate::units::{Unit, UnitError};

pub const MAX_PROPERTY_KEY_LEN: usize = 128;

/// A validated, deterministic property identifier.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PropertyKey(String);

impl PropertyKey {
    pub fn new(value: impl Into<String>) -> Result<Self, PropertyError> {
        let value = value.into();
        if !valid_property_key(&value) {
            return Err(PropertyError::InvalidKey(value));
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl AsRef<str> for PropertyKey {
    fn as_ref(&self) -> &str {
        self.as_str()
    }
}

impl fmt::Display for PropertyKey {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl TryFrom<&str> for PropertyKey {
    type Error = PropertyError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl TryFrom<String> for PropertyKey {
    type Error = PropertyError;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        Self::new(value)
    }
}

impl FromStr for PropertyKey {
    type Err = PropertyError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::new(value)
    }
}

fn valid_property_key(value: &str) -> bool {
    if value.is_empty() || value.len() > MAX_PROPERTY_KEY_LEN {
        return false;
    }
    let mut bytes = value.bytes();
    let Some(first) = bytes.next() else {
        return false;
    };
    (first.is_ascii_alphabetic() || first == b'_')
        && bytes.all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.'))
}

/// One generic scalar annotation.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyValue {
    Bool(bool),
    Int(i64),
    Real { value: f64, unit: Unit },
    String(String),
}

impl PropertyValue {
    /// Borrows this scalar without copying its string storage.
    pub fn as_ref(&self) -> PropertyValueRef<'_> {
        match self {
            Self::Bool(value) => PropertyValueRef::Bool(*value),
            Self::Int(value) => PropertyValueRef::Int(*value),
            Self::Real { value, unit } => PropertyValueRef::Real {
                value: *value,
                unit: *unit,
            },
            Self::String(value) => PropertyValueRef::String(value),
        }
    }

    pub fn real(value: f64, unit: Unit) -> Result<Self, PropertyError> {
        let value = Self::Real { value, unit };
        value.validate()?;
        Ok(value)
    }

    fn validate(&self) -> Result<(), PropertyError> {
        if let Self::Real { value, .. } = self {
            if !value.is_finite() {
                return Err(PropertyError::NonFiniteValue { index: None });
            }
        }
        Ok(())
    }
}

/// A scalar property view that borrows string storage from its owner or column.
///
/// Numerical values and units are copied; strings remain borrowed. Use
/// [`Self::to_value`] when an independently owned value is needed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PropertyValueRef<'a> {
    Bool(bool),
    Int(i64),
    Real { value: f64, unit: Unit },
    String(&'a str),
}

impl PropertyValueRef<'_> {
    /// Materializes an owned scalar, copying string storage when present.
    pub fn to_value(self) -> PropertyValue {
        match self {
            Self::Bool(value) => PropertyValue::Bool(value),
            Self::Int(value) => PropertyValue::Int(value),
            Self::Real { value, unit } => PropertyValue::Real { value, unit },
            Self::String(value) => PropertyValue::String(value.to_owned()),
        }
    }
}

/// One homogeneous optional-valued property column.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyColumn {
    Bool(Vec<Option<bool>>),
    Int(Vec<Option<i64>>),
    Real {
        unit: Unit,
        values: Vec<Option<f64>>,
    },
    String(Vec<Option<String>>),
}

impl PropertyColumn {
    /// Checked expansion from live entity order into a private editor slot space.
    pub(crate) fn into_editor_slots(
        self,
        live_slots: &[usize],
        slot_count: usize,
    ) -> Result<Self, PropertyError> {
        if self.len() != live_slots.len() {
            return Err(PropertyError::LengthMismatch {
                expected: live_slots.len(),
                actual: self.len(),
            });
        }
        self.validate()?;
        let mut indices = vec![None; slot_count];
        for (index, &slot) in live_slots.iter().enumerate() {
            let destination = indices.get_mut(slot).ok_or(PropertyError::InvalidIndex {
                len: slot_count,
                index: slot,
            })?;
            *destination = Some(index);
        }
        Ok(self.select_optional_indices(&indices))
    }
    pub fn len(&self) -> usize {
        match self {
            Self::Bool(values) => values.len(),
            Self::Int(values) => values.len(),
            Self::Real { values, .. } => values.len(),
            Self::String(values) => values.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn is_all_missing(&self) -> bool {
        match self {
            Self::Bool(values) => values.iter().all(Option::is_none),
            Self::Int(values) => values.iter().all(Option::is_none),
            Self::Real { values, .. } => values.iter().all(Option::is_none),
            Self::String(values) => values.iter().all(Option::is_none),
        }
    }

    /// Reads a scalar, copying string storage. See [`Self::value_ref`] to borrow.
    pub fn value(&self, index: usize) -> Result<Option<PropertyValue>, PropertyError> {
        self.value_ref(index)
            .map(|value| value.map(PropertyValueRef::to_value))
    }

    /// Borrows one cell without allocating, preserving missing-value semantics.
    pub fn value_ref(&self, index: usize) -> Result<Option<PropertyValueRef<'_>>, PropertyError> {
        if index >= self.len() {
            return Err(PropertyError::InvalidIndex {
                len: self.len(),
                index,
            });
        }
        Ok(match self {
            Self::Bool(values) => values[index].map(PropertyValueRef::Bool),
            Self::Int(values) => values[index].map(PropertyValueRef::Int),
            Self::Real { unit, values } => {
                values[index].map(|value| PropertyValueRef::Real { value, unit: *unit })
            }
            Self::String(values) => values[index].as_deref().map(PropertyValueRef::String),
        })
    }

    fn validate(&self) -> Result<(), PropertyError> {
        if let Self::Real { values, .. } = self {
            if let Some(index) = values
                .iter()
                .position(|value| value.is_some_and(|value| !value.is_finite()))
            {
                return Err(PropertyError::NonFiniteValue { index: Some(index) });
            }
        }
        Ok(())
    }

    fn converted_to(self, unit: Unit) -> Result<Self, PropertyError> {
        let Self::Real {
            unit: source_unit,
            values,
        } = self
        else {
            return Ok(self);
        };
        let factor = source_unit.conversion_factor_to(unit)?;
        let values = values
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                let converted = value.map(|value| value * factor);
                if converted.is_some_and(|value| !value.is_finite()) {
                    return Err(PropertyError::NonFiniteValue { index: Some(index) });
                }
                Ok(converted)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::Real { unit, values })
    }

    fn select_indices(&self, indices: &[usize]) -> Self {
        match self {
            Self::Bool(values) => Self::Bool(indices.iter().map(|index| values[*index]).collect()),
            Self::Int(values) => Self::Int(indices.iter().map(|index| values[*index]).collect()),
            Self::Real { unit, values } => Self::Real {
                unit: *unit,
                values: indices.iter().map(|index| values[*index]).collect(),
            },
            Self::String(values) => {
                Self::String(indices.iter().map(|index| values[*index].clone()).collect())
            }
        }
    }

    fn select_optional_indices(&self, indices: &[Option<usize>]) -> Self {
        match self {
            Self::Bool(values) => Self::Bool(
                indices
                    .iter()
                    .map(|index| index.and_then(|index| values[index]))
                    .collect(),
            ),
            Self::Int(values) => Self::Int(
                indices
                    .iter()
                    .map(|index| index.and_then(|index| values[index]))
                    .collect(),
            ),
            Self::Real { unit, values } => Self::Real {
                unit: *unit,
                values: indices
                    .iter()
                    .map(|index| index.and_then(|index| values[index]))
                    .collect(),
            },
            Self::String(values) => Self::String(
                indices
                    .iter()
                    .map(|index| index.and_then(|index| values[index].clone()))
                    .collect(),
            ),
        }
    }

    fn populated_len(&self) -> usize {
        match self {
            Self::Bool(values) => values.iter().filter(|value| value.is_some()).count(),
            Self::Int(values) => values.iter().filter(|value| value.is_some()).count(),
            Self::Real { values, .. } => values.iter().filter(|value| value.is_some()).count(),
            Self::String(values) => values.iter().filter(|value| value.is_some()).count(),
        }
    }

    pub(crate) fn is_populated_at(&self, index: usize) -> bool {
        match self {
            Self::Bool(values) => values[index].is_some(),
            Self::Int(values) => values[index].is_some(),
            Self::Real { values, .. } => values[index].is_some(),
            Self::String(values) => values[index].is_some(),
        }
    }

    fn resize_missing(&mut self, len: usize) {
        match self {
            Self::Bool(values) => values.resize(len, None),
            Self::Int(values) => values.resize(len, None),
            Self::Real { values, .. } => values.resize(len, None),
            Self::String(values) => values.resize(len, None),
        }
    }

    fn clear(&mut self, index: usize) {
        match self {
            Self::Bool(values) => values[index] = None,
            Self::Int(values) => values[index] = None,
            Self::Real { values, .. } => values[index] = None,
            Self::String(values) => values[index] = None,
        }
    }
}

mod sealed {
    pub trait Sealed {}
}

/// Typed row address of one property table.
///
/// Each table is indexed by the identifier of the domain it annotates, so a
/// molecule atom table cannot be read with a topology index by accident.
pub trait PropertyRow: Copy + sealed::Sealed {
    /// Zero-based row in the owning table.
    fn row(self) -> usize;
}

macro_rules! property_rows {
    ($($row:ty),* $(,)?) => {$(
        impl sealed::Sealed for $row {}
        impl PropertyRow for $row {
            fn row(self) -> usize {
                self.index()
            }
        }
    )*};
}

property_rows!(
    crate::core::AtomId,
    crate::core::BondId,
    crate::topology::MoleculeInstanceId,
    crate::topology::TopologyAtomIndex,
    crate::topology::TopologyBondIndex,
    crate::topology::ChainId,
    crate::topology::ResidueId,
    crate::topology::AtomSiteId,
);

/// Columnar properties for one homogeneous entity domain, addressed by `R`.
///
/// The row count is fixed by the owning object. Read through this table and
/// mutate through the owner's [`PropertyTableMut`] guard.
pub struct PropertyTable<R> {
    raw: RawPropertyTable,
    row: PhantomData<fn() -> R>,
}

impl<R> fmt::Debug for PropertyTable<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.raw.fmt(formatter)
    }
}

impl<R> Clone for PropertyTable<R> {
    fn clone(&self) -> Self {
        Self::from_raw(self.raw.clone())
    }
}

impl<R> PartialEq for PropertyTable<R> {
    fn eq(&self, other: &Self) -> bool {
        self.raw == other.raw
    }
}

impl<R> Default for PropertyTable<R> {
    fn default() -> Self {
        Self::from_raw(RawPropertyTable::default())
    }
}

impl<R> PropertyTable<R> {
    pub(crate) fn new(len: usize) -> Self {
        Self::from_raw(RawPropertyTable::new(len))
    }

    pub(crate) const fn from_raw(raw: RawPropertyTable) -> Self {
        Self {
            raw,
            row: PhantomData,
        }
    }

    pub(crate) const fn raw(&self) -> &RawPropertyTable {
        &self.raw
    }

    /// Number of rows fixed by the owner.
    pub fn len(&self) -> usize {
        self.raw.len()
    }

    pub fn is_empty(&self) -> bool {
        self.raw.is_empty()
    }

    /// Whether any column is present. All-missing columns are never stored.
    pub fn has_data(&self) -> bool {
        self.raw.has_data()
    }

    pub fn get(&self, key: &PropertyKey) -> Option<&PropertyColumn> {
        self.raw.get(key)
    }

    /// Columns in deterministic key order.
    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyColumn)> + DoubleEndedIterator {
        self.raw.iter()
    }

    pub fn keys(&self) -> impl ExactSizeIterator<Item = &PropertyKey> + DoubleEndedIterator {
        self.raw.columns.keys()
    }
}

impl<R: PropertyRow> PropertyTable<R> {
    pub fn row_has_data(&self, row: R) -> Result<bool, PropertyError> {
        self.raw.row_has_data(row.row())
    }

    /// Reads a cell, copying string storage. See [`Self::value_ref`] to borrow.
    pub fn value(&self, key: &PropertyKey, row: R) -> Result<Option<PropertyValue>, PropertyError> {
        self.raw.value(key, row.row())
    }

    /// Borrows one cell without allocating.
    ///
    /// Missing keys and missing cells return `None`; a row outside the table
    /// is an error even when the key is absent.
    pub fn value_ref(
        &self,
        key: &PropertyKey,
        row: R,
    ) -> Result<Option<PropertyValueRef<'_>>, PropertyError> {
        self.raw.value_ref(key, row.row())
    }
}

/// Mutable columns of an owner-sized table. The row count cannot change.
///
/// Reads are available through `Deref`. Every mutation is transactional.
///
/// ```compile_fail
/// use kekule::{properties::PropertyTable, topology::TopologyBuilder};
/// let mut builder = TopologyBuilder::new();
/// let mut properties = builder.properties_mut();
/// *properties.atoms_mut() = PropertyTable::default();
/// ```
pub struct PropertyTableMut<'a, R> {
    table: &'a mut PropertyTable<R>,
    live: LiveRows<'a>,
}

/// Which rows of a table may receive values; editor drafts keep the rows of
/// removed atoms and bonds allocated but writable only as missing.
#[derive(Debug, Clone, Copy)]
pub(crate) enum LiveRows<'a> {
    All,
    Atoms(&'a [Option<crate::core::Atom>]),
    Bonds(&'a [Option<crate::core::Bond>]),
}

impl LiveRows<'_> {
    fn check(self, row: usize) -> Result<(), PropertyError> {
        let live = match self {
            Self::All => true,
            Self::Atoms(slots) => slots.get(row).is_none_or(Option::is_some),
            Self::Bonds(slots) => slots.get(row).is_none_or(Option::is_some),
        };
        if live {
            Ok(())
        } else {
            Err(PropertyError::RemovedRow { index: row })
        }
    }

    fn check_column(self, column: &PropertyColumn) -> Result<(), PropertyError> {
        if matches!(self, Self::All) {
            return Ok(());
        }
        (0..column.len())
            .filter(|row| column.is_populated_at(*row))
            .try_for_each(|row| self.check(row))
    }
}

impl<R> fmt::Debug for PropertyTableMut<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.table.fmt(formatter)
    }
}

impl<R> std::ops::Deref for PropertyTableMut<'_, R> {
    type Target = PropertyTable<R>;

    fn deref(&self) -> &Self::Target {
        self.table
    }
}

impl<'a, R> PropertyTableMut<'a, R> {
    pub(crate) fn new(table: &'a mut PropertyTable<R>) -> Self {
        Self {
            table,
            live: LiveRows::All,
        }
    }

    pub(crate) fn with_live_rows(table: &'a mut PropertyTable<R>, live: LiveRows<'a>) -> Self {
        Self { table, live }
    }

    /// Inserts or replaces a complete column in row order.
    ///
    /// Replacing a real column preserves its storage unit and converts
    /// compatible input values. An all-missing input removes the key.
    pub fn insert(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        if column.len() == self.table.len() {
            self.live.check_column(&column)?;
        }
        self.table.raw.insert(key, column)
    }

    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.table.raw.remove(key)
    }

    /// Removes every column, keeping the row count.
    pub fn clear(&mut self) {
        self.table.raw.clear_columns();
    }
}

impl<R: PropertyRow> PropertyTableMut<'_, R> {
    /// Sets or clears one cell, removing a column that becomes all-missing.
    pub fn set_value(
        &mut self,
        key: PropertyKey,
        row: R,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        self.check_row(row.row())?;
        self.table.raw.set_value(key, row.row(), value)
    }

    pub fn clear_value(&mut self, key: PropertyKey, row: R) -> Result<(), PropertyError> {
        self.check_row(row.row())?;
        self.table.raw.set_value(key, row.row(), None)
    }

    fn check_row(&self, row: usize) -> Result<(), PropertyError> {
        if row < self.table.len() {
            self.live.check(row)?;
        }
        Ok(())
    }

    /// Applies several cell updates to one column transactionally: either
    /// every update succeeds or the column is unchanged. Other columns are
    /// never copied.
    pub fn set_values(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (R, Option<PropertyValue>)>,
    ) -> Result<(), PropertyError> {
        let mut staged = self.table.raw.stage_column(&key);
        for (row, value) in values {
            self.check_row(row.row())?;
            staged.set_value(key.clone(), row.row(), value)?;
        }
        self.table.raw.commit_column(key, staged);
        Ok(())
    }
}

/// Untyped row storage shared by every typed [`PropertyTable`].
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct RawPropertyTable {
    len: usize,
    columns: BTreeMap<PropertyKey, PropertyColumn>,
    populated: BTreeMap<PropertyKey, usize>,
}

impl RawPropertyTable {
    pub(crate) fn new(len: usize) -> Self {
        Self {
            len,
            columns: BTreeMap::new(),
            populated: BTreeMap::new(),
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.len
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub(crate) fn has_data(&self) -> bool {
        !self.columns.is_empty()
    }

    pub(crate) fn row_has_data(&self, index: usize) -> Result<bool, PropertyError> {
        self.validate_index(index)?;
        Ok(self
            .columns
            .values()
            .any(|column| column.is_populated_at(index)))
    }

    pub(crate) fn get(&self, key: &PropertyKey) -> Option<&PropertyColumn> {
        self.columns.get(key)
    }

    pub(crate) fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyColumn)> + DoubleEndedIterator {
        self.columns.iter()
    }

    /// Inserts or replaces a complete column transactionally.
    ///
    /// Replacing a real column preserves its existing storage unit and converts
    /// compatible input values into it. An all-missing input removes the key.
    pub(crate) fn insert(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        self.insert_validated(key, column)
    }

    fn insert_validated(
        &mut self,
        key: PropertyKey,
        mut column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, PropertyError> {
        if column.len() != self.len {
            return Err(PropertyError::LengthMismatch {
                expected: self.len,
                actual: column.len(),
            });
        }
        column.validate()?;
        if column.is_all_missing() {
            self.populated.remove(&key);
            return Ok(self.columns.remove(&key));
        }
        if let Some(PropertyColumn::Real { unit, .. }) = self.columns.get(&key) {
            column = match column {
                PropertyColumn::Real { .. } => column.converted_to(*unit)?,
                _ => return Err(PropertyError::TypeMismatch { key }),
            };
        } else if let Some(existing) = self.columns.get(&key) {
            if std::mem::discriminant(existing) != std::mem::discriminant(&column) {
                return Err(PropertyError::TypeMismatch { key });
            }
        }
        let populated = column.populated_len();
        self.populated.insert(key.clone(), populated);
        Ok(self.columns.insert(key, column))
    }

    pub(crate) fn remove(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.populated.remove(key);
        self.columns.remove(key)
    }

    /// Stages one column for a transactional sequence of ordinary cell updates.
    /// Other columns are not copied. Commit only after every update succeeds.
    pub(crate) fn stage_column(&self, key: &PropertyKey) -> Self {
        let mut staged = Self::new(self.len);
        if let Some(column) = self.columns.get(key) {
            staged.columns.insert(key.clone(), column.clone());
            staged.populated.insert(key.clone(), self.populated[key]);
        }
        staged
    }

    /// Installs the final staged column, preserving its exact stored values.
    /// A sequence may clear a column and recreate it with a new type or unit.
    pub(crate) fn commit_column(&mut self, key: PropertyKey, mut staged: Self) {
        debug_assert_eq!(self.len, staged.len);
        self.remove(&key);
        if let Some(column) = staged.columns.remove(&key) {
            let populated = staged.populated.remove(&key).expect("staged column count");
            self.populated.insert(key.clone(), populated);
            self.columns.insert(key, column);
        }
    }

    pub(crate) fn value(
        &self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        self.value_ref(key, index)
            .map(|value| value.map(PropertyValueRef::to_value))
    }

    /// Borrows one property cell without allocating.
    ///
    /// Missing keys and missing cells return `None`; an invalid row is an error
    /// even when the key is absent.
    pub(crate) fn value_ref(
        &self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValueRef<'_>>, PropertyError> {
        self.validate_index(index)?;
        self.columns
            .get(key)
            .map(|column| column.value_ref(index))
            .transpose()
            .map(Option::flatten)
    }

    /// Sets or clears one cell transactionally, removing an all-missing column.
    pub(crate) fn set_value(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        self.set_value_validated(key, index, value)
    }

    fn set_value_validated(
        &mut self,
        key: PropertyKey,
        index: usize,
        value: Option<PropertyValue>,
    ) -> Result<(), PropertyError> {
        self.validate_index(index)?;
        if let Some(value) = &value {
            value.validate()?;
        }
        let Some(column) = self.columns.get(&key) else {
            let Some(value) = value else {
                return Ok(());
            };
            let mut column = match &value {
                PropertyValue::Bool(_) => PropertyColumn::Bool(vec![None; self.len]),
                PropertyValue::Int(_) => PropertyColumn::Int(vec![None; self.len]),
                PropertyValue::Real { unit, .. } => PropertyColumn::Real {
                    unit: *unit,
                    values: vec![None; self.len],
                },
                PropertyValue::String(_) => PropertyColumn::String(vec![None; self.len]),
            };
            assign_value(&mut column, index, Some(value));
            self.populated.insert(key.clone(), 1);
            self.columns.insert(key, column);
            return Ok(());
        };

        let value = validate_value_for_column(column, &key, index, value)?;
        let was_populated = column.is_populated_at(index);
        let will_be_populated = value.is_some();
        let column = self
            .columns
            .get_mut(&key)
            .expect("validated property column must remain present");
        assign_value(column, index, value);

        let populated = self
            .populated
            .get_mut(&key)
            .expect("published property column must have a populated count");
        match (was_populated, will_be_populated) {
            (false, true) => *populated += 1,
            (true, false) => *populated -= 1,
            _ => {}
        }
        if *populated == 0 {
            self.populated.remove(&key);
            self.columns.remove(&key);
        }
        Ok(())
    }

    pub(crate) fn select_indices(&self, indices: &[usize]) -> Result<Self, PropertyError> {
        for index in indices {
            self.validate_index(*index)?;
        }
        let columns: BTreeMap<_, _> = self
            .columns
            .iter()
            .filter_map(|(key, column)| {
                let selected = column.select_indices(indices);
                (!selected.is_all_missing()).then(|| (key.clone(), selected))
            })
            .collect();
        let populated = columns
            .iter()
            .map(|(key, column)| (key.clone(), column.populated_len()))
            .collect();
        Ok(Self {
            len: indices.len(),
            columns,
            populated,
        })
    }

    pub(crate) fn select_optional_indices(
        &self,
        indices: &[Option<usize>],
    ) -> Result<Self, PropertyError> {
        for index in indices.iter().flatten() {
            self.validate_index(*index)?;
        }
        let columns: BTreeMap<_, _> = self
            .columns
            .iter()
            .filter_map(|(key, column)| {
                let selected = column.select_optional_indices(indices);
                (!selected.is_all_missing()).then(|| (key.clone(), selected))
            })
            .collect();
        let populated = columns
            .iter()
            .map(|(key, column)| (key.clone(), column.populated_len()))
            .collect();
        Ok(Self {
            len: indices.len(),
            columns,
            populated,
        })
    }

    pub(crate) fn clear_columns(&mut self) {
        self.columns.clear();
        self.populated.clear();
    }

    pub(crate) fn resize_missing(&mut self, len: usize) {
        let previous_len = self.len;
        self.len = len;
        if len >= previous_len {
            for column in self.columns.values_mut() {
                column.resize_missing(len);
            }
            return;
        }
        let populated = &mut self.populated;
        self.columns.retain(|key, column| {
            column.resize_missing(len);
            let count = column.populated_len();
            if count == 0 {
                populated.remove(key);
                false
            } else {
                populated.insert(key.clone(), count);
                true
            }
        });
    }

    /// Copies source rows into explicitly allocated destination slots. The caller
    /// stages the containing operation, so errors may leave this private draft changed.
    pub(crate) fn copy_rows_from(
        &mut self,
        source: &Self,
        rows: &[usize],
    ) -> Result<(), PropertyError> {
        if rows.len() != source.len() {
            return Err(PropertyError::LengthMismatch {
                expected: source.len(),
                actual: rows.len(),
            });
        }
        for &row in rows {
            self.validate_index(row)?;
        }
        for (key, column) in source.iter() {
            for (index, &row) in rows.iter().enumerate() {
                if let Some(value) = column.value(index)? {
                    self.set_value(key.clone(), row, Some(value))?;
                }
            }
        }
        Ok(())
    }

    pub(crate) fn clear_index(&mut self, index: usize) {
        if index >= self.len {
            return;
        }
        let populated = &mut self.populated;
        self.columns.retain(|key, column| {
            let was_populated = column.is_populated_at(index);
            column.clear(index);
            let count = populated
                .get_mut(key)
                .expect("published property column must have a populated count");
            if was_populated {
                *count -= 1;
            }
            if *count == 0 {
                populated.remove(key);
                false
            } else {
                true
            }
        });
    }

    fn validate_index(&self, index: usize) -> Result<(), PropertyError> {
        if index >= self.len {
            return Err(PropertyError::InvalidIndex {
                len: self.len,
                index,
            });
        }
        Ok(())
    }
}

/// Owner-level scalar annotations, keyed deterministically.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct OwnerProperties {
    values: BTreeMap<PropertyKey, PropertyValue>,
}

impl OwnerProperties {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn get(&self, key: &PropertyKey) -> Option<&PropertyValue> {
        self.values.get(key)
    }

    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyValue)> + DoubleEndedIterator {
        self.values.iter()
    }

    pub fn keys(&self) -> impl ExactSizeIterator<Item = &PropertyKey> + DoubleEndedIterator {
        self.values.keys()
    }

    /// Inserts a validated value, returning the one it replaces.
    pub fn insert(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        value.validate()?;
        Ok(self.values.insert(key, value))
    }

    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.values.remove(key)
    }

    pub fn clear(&mut self) {
        self.values.clear();
    }
}

/// Annotations scoped to one molecule definition: owner values plus one row
/// per atom and bond slot, addressed by [`AtomId`] and [`BondId`].
///
/// Published molecules have dense IDs, so every row is live. In a
/// [`crate::core::MoleculeEditor`] draft, rows of removed atoms and bonds stay
/// allocated and missing.
#[derive(Debug, Clone, PartialEq)]
pub struct MoleculeProperties {
    owner: OwnerProperties,
    atoms: PropertyTable<AtomId>,
    bonds: PropertyTable<BondId>,
}

impl MoleculeProperties {
    pub(crate) fn new(atoms: usize, bonds: usize) -> Self {
        Self {
            owner: OwnerProperties::new(),
            atoms: PropertyTable::new(atoms),
            bonds: PropertyTable::new(bonds),
        }
    }

    pub fn owner(&self) -> &OwnerProperties {
        &self.owner
    }

    pub fn atoms(&self) -> &PropertyTable<AtomId> {
        &self.atoms
    }

    pub fn bonds(&self) -> &PropertyTable<BondId> {
        &self.bonds
    }

    pub fn is_empty(&self) -> bool {
        self.owner.is_empty() && !self.atoms.has_data() && !self.bonds.has_data()
    }

    pub(crate) fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.owner
    }

    pub(crate) fn atoms_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.atoms.raw
    }

    pub(crate) fn bonds_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.bonds.raw
    }
}

/// Length-preserving mutable access to [`MoleculeProperties`].
#[derive(Debug)]
pub struct MoleculePropertiesMut<'a> {
    properties: &'a mut MoleculeProperties,
    atoms_live: LiveRows<'a>,
    bonds_live: LiveRows<'a>,
}

impl std::ops::Deref for MoleculePropertiesMut<'_> {
    type Target = MoleculeProperties;

    fn deref(&self) -> &Self::Target {
        self.properties
    }
}

impl<'a> MoleculePropertiesMut<'a> {
    /// Rows of removed draft atoms and bonds reject values.
    pub(crate) fn new(
        properties: &'a mut MoleculeProperties,
        atoms: &'a [Option<crate::core::Atom>],
        bonds: &'a [Option<crate::core::Bond>],
    ) -> Self {
        Self {
            properties,
            atoms_live: LiveRows::Atoms(atoms),
            bonds_live: LiveRows::Bonds(bonds),
        }
    }

    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties.owner
    }

    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, AtomId> {
        PropertyTableMut::with_live_rows(&mut self.properties.atoms, self.atoms_live)
    }

    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, BondId> {
        PropertyTableMut::with_live_rows(&mut self.properties.bonds, self.bonds_live)
    }

    /// Removes every owner value and column, keeping row counts.
    pub fn clear(&mut self) {
        self.properties.owner.clear();
        self.properties.atoms.raw.clear_columns();
        self.properties.bonds.raw.clear_columns();
    }
}

/// Annotations scoped to one topology: owner values plus dense rows for
/// molecule instances, atoms, bonds, chains, residues, and atom sites.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TopologyProperties {
    owner: OwnerProperties,
    molecule_instances: PropertyTable<MoleculeInstanceId>,
    atoms: PropertyTable<TopologyAtomIndex>,
    bonds: PropertyTable<TopologyBondIndex>,
    chains: PropertyTable<ChainId>,
    residues: PropertyTable<ResidueId>,
    atom_sites: PropertyTable<AtomSiteId>,
}

impl TopologyProperties {
    pub fn owner(&self) -> &OwnerProperties {
        &self.owner
    }
    pub fn molecule_instances(&self) -> &PropertyTable<MoleculeInstanceId> {
        &self.molecule_instances
    }
    pub fn atoms(&self) -> &PropertyTable<TopologyAtomIndex> {
        &self.atoms
    }
    pub fn bonds(&self) -> &PropertyTable<TopologyBondIndex> {
        &self.bonds
    }
    pub fn chains(&self) -> &PropertyTable<ChainId> {
        &self.chains
    }
    pub fn residues(&self) -> &PropertyTable<ResidueId> {
        &self.residues
    }
    pub fn atom_sites(&self) -> &PropertyTable<AtomSiteId> {
        &self.atom_sites
    }

    pub fn is_empty(&self) -> bool {
        self.owner.is_empty() && self.tables().iter().all(|table| !table.has_data())
    }

    pub(crate) fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.owner
    }
    pub(crate) fn molecule_instances_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.molecule_instances.raw
    }
    pub(crate) fn atoms_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.atoms.raw
    }
    pub(crate) fn bonds_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.bonds.raw
    }
    pub(crate) fn chains_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.chains.raw
    }
    pub(crate) fn residues_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.residues.raw
    }
    pub(crate) fn atom_sites_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.atom_sites.raw
    }

    pub(crate) fn resize_atoms(&mut self, len: usize) {
        self.atoms.raw.resize_missing(len);
    }

    pub(crate) fn resize_bonds(&mut self, len: usize) {
        self.bonds.raw.resize_missing(len);
    }

    pub(crate) fn resize_domains(&mut self, dimensions: [usize; 6]) {
        for (table, len) in self.tables_mut().into_iter().zip(dimensions) {
            // Extension appends missing values. A shorter hierarchy replacement
            // must not truncate populated annotations before publication checks.
            if !table.has_data() || len >= table.len() {
                table.resize_missing(len);
            }
        }
    }

    pub(crate) fn validate_dimensions(&self, dimensions: [usize; 6]) -> Result<(), PropertyError> {
        for (table, expected) in self.tables().into_iter().zip(dimensions) {
            if table.len() != expected {
                return Err(PropertyError::LengthMismatch {
                    expected,
                    actual: table.len(),
                });
            }
        }
        Ok(())
    }

    /// Projects entity rows through one operation's correspondence; owner
    /// values are dropped because the projected system is a new owner.
    pub(crate) fn project(
        &self,
        molecule_instances: &[Option<usize>],
        atoms: &[usize],
        bonds: &[usize],
        chains: &[usize],
        residues: &[usize],
        atom_sites: &[usize],
    ) -> Result<Self, PropertyError> {
        Ok(Self {
            owner: OwnerProperties::new(),
            molecule_instances: PropertyTable::from_raw(
                self.molecule_instances
                    .raw
                    .select_optional_indices(molecule_instances)?,
            ),
            atoms: PropertyTable::from_raw(self.atoms.raw.select_indices(atoms)?),
            bonds: PropertyTable::from_raw(self.bonds.raw.select_indices(bonds)?),
            chains: PropertyTable::from_raw(self.chains.raw.select_indices(chains)?),
            residues: PropertyTable::from_raw(self.residues.raw.select_indices(residues)?),
            atom_sites: PropertyTable::from_raw(self.atom_sites.raw.select_indices(atom_sites)?),
        })
    }

    fn tables(&self) -> [&RawPropertyTable; 6] {
        [
            &self.molecule_instances.raw,
            &self.atoms.raw,
            &self.bonds.raw,
            &self.chains.raw,
            &self.residues.raw,
            &self.atom_sites.raw,
        ]
    }

    fn tables_mut(&mut self) -> [&mut RawPropertyTable; 6] {
        [
            &mut self.molecule_instances.raw,
            &mut self.atoms.raw,
            &mut self.bonds.raw,
            &mut self.chains.raw,
            &mut self.residues.raw,
            &mut self.atom_sites.raw,
        ]
    }
}

/// Length-preserving mutable access to staged [`TopologyProperties`].
#[derive(Debug)]
pub struct TopologyPropertiesMut<'a> {
    properties: &'a mut TopologyProperties,
}

impl std::ops::Deref for TopologyPropertiesMut<'_> {
    type Target = TopologyProperties;

    fn deref(&self) -> &Self::Target {
        self.properties
    }
}

impl<'a> TopologyPropertiesMut<'a> {
    pub(crate) fn new(properties: &'a mut TopologyProperties) -> Self {
        Self { properties }
    }
    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties.owner
    }
    pub fn molecule_instances_mut(&mut self) -> PropertyTableMut<'_, MoleculeInstanceId> {
        PropertyTableMut::new(&mut self.properties.molecule_instances)
    }
    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, TopologyAtomIndex> {
        PropertyTableMut::new(&mut self.properties.atoms)
    }
    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, TopologyBondIndex> {
        PropertyTableMut::new(&mut self.properties.bonds)
    }
    pub fn chains_mut(&mut self) -> PropertyTableMut<'_, ChainId> {
        PropertyTableMut::new(&mut self.properties.chains)
    }
    pub fn residues_mut(&mut self) -> PropertyTableMut<'_, ResidueId> {
        PropertyTableMut::new(&mut self.properties.residues)
    }
    pub fn atom_sites_mut(&mut self) -> PropertyTableMut<'_, AtomSiteId> {
        PropertyTableMut::new(&mut self.properties.atom_sites)
    }

    /// Removes every owner value and column, keeping row counts.
    pub fn clear(&mut self) {
        self.properties.owner.clear();
        for table in self.properties.tables_mut() {
            table.clear_columns();
        }
    }
}

/// Coordinate-dependent annotations of one realization: owner values plus one
/// row per topology atom and bond, in dense topology order.
///
/// A detached realization payload has atom rows only; bond rows are allocated
/// when the payload is bound to a topology.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct RealizationProperties {
    owner: OwnerProperties,
    atoms: PropertyTable<TopologyAtomIndex>,
    bonds: PropertyTable<TopologyBondIndex>,
}

impl RealizationProperties {
    pub(crate) fn new(atoms: usize) -> Self {
        Self {
            owner: OwnerProperties::new(),
            atoms: PropertyTable::new(atoms),
            bonds: PropertyTable::new(0),
        }
    }

    pub fn owner(&self) -> &OwnerProperties {
        &self.owner
    }

    pub fn atoms(&self) -> &PropertyTable<TopologyAtomIndex> {
        &self.atoms
    }

    pub fn bonds(&self) -> &PropertyTable<TopologyBondIndex> {
        &self.bonds
    }

    pub fn is_empty(&self) -> bool {
        self.owner.is_empty() && !self.atoms.has_data() && !self.bonds.has_data()
    }

    pub(crate) fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.owner
    }

    /// Whether bond rows can be bound to a topology with `bond_count` bonds:
    /// either they already match or none were allocated yet.
    pub(crate) fn validate_bond_rows(&self, bond_count: usize) -> Result<(), PropertyError> {
        if self.bonds.len() == bond_count || self.bonds.is_empty() {
            return Ok(());
        }
        Err(PropertyError::LengthMismatch {
            expected: bond_count,
            actual: self.bonds.len(),
        })
    }

    /// Allocates bond rows after [`Self::validate_bond_rows`] succeeded.
    pub(crate) fn bind_bond_rows(&mut self, bond_count: usize) {
        debug_assert!(self.validate_bond_rows(bond_count).is_ok());
        self.bonds.raw.resize_missing(bond_count);
    }

    pub(crate) fn atoms_raw_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.atoms.raw
    }

    pub(crate) fn bonds_raw_mut(&mut self) -> &mut RawPropertyTable {
        &mut self.bonds.raw
    }

    /// Grows staged rows with missing values.
    pub(crate) fn resize(&mut self, atoms: usize, bonds: usize) {
        self.atoms.raw.resize_missing(atoms);
        self.bonds.raw.resize_missing(bonds);
    }

    /// Reorders atom rows; `previous[new]` is the old row.
    pub(crate) fn reorder_atoms(&mut self, previous: &[usize]) -> Result<(), PropertyError> {
        self.atoms.raw = self.atoms.raw.select_indices(previous)?;
        Ok(())
    }

    /// Projects atom and bond rows; owner values are dropped because the
    /// projected realization is a new owner.
    pub(crate) fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, PropertyError> {
        let bonds = if self.bonds.is_empty() {
            RawPropertyTable::new(0)
        } else {
            self.bonds.raw.select_indices(bonds)?
        };
        Ok(Self {
            owner: OwnerProperties::new(),
            atoms: PropertyTable::from_raw(self.atoms.raw.select_indices(atoms)?),
            bonds: PropertyTable::from_raw(bonds),
        })
    }
}

/// Length-preserving mutable access to [`RealizationProperties`].
#[derive(Debug)]
pub struct RealizationPropertiesMut<'a> {
    properties: &'a mut RealizationProperties,
}

impl std::ops::Deref for RealizationPropertiesMut<'_> {
    type Target = RealizationProperties;

    fn deref(&self) -> &Self::Target {
        self.properties
    }
}

impl<'a> RealizationPropertiesMut<'a> {
    pub(crate) fn new(properties: &'a mut RealizationProperties) -> Self {
        Self { properties }
    }

    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties.owner
    }

    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, TopologyAtomIndex> {
        PropertyTableMut::new(&mut self.properties.atoms)
    }

    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, TopologyBondIndex> {
        PropertyTableMut::new(&mut self.properties.bonds)
    }

    /// Removes every owner value and column, keeping row counts.
    pub fn clear(&mut self) {
        self.properties.owner.clear();
        self.properties.atoms.raw.clear_columns();
        self.properties.bonds.raw.clear_columns();
    }
}

fn validate_value_for_column(
    column: &PropertyColumn,
    key: &PropertyKey,
    index: usize,
    value: Option<PropertyValue>,
) -> Result<Option<PropertyValue>, PropertyError> {
    let Some(value) = value else {
        return Ok(None);
    };
    match (column, value) {
        (PropertyColumn::Bool(_), value @ PropertyValue::Bool(_))
        | (PropertyColumn::Int(_), value @ PropertyValue::Int(_))
        | (PropertyColumn::String(_), value @ PropertyValue::String(_)) => Ok(Some(value)),
        (
            PropertyColumn::Real { unit, .. },
            PropertyValue::Real {
                value,
                unit: source_unit,
            },
        ) => {
            let value = value * source_unit.conversion_factor_to(*unit)?;
            if !value.is_finite() {
                return Err(PropertyError::NonFiniteValue { index: Some(index) });
            }
            Ok(Some(PropertyValue::Real { value, unit: *unit }))
        }
        _ => Err(PropertyError::TypeMismatch { key: key.clone() }),
    }
}

fn assign_value(column: &mut PropertyColumn, index: usize, value: Option<PropertyValue>) {
    match (column, value) {
        (PropertyColumn::Bool(values), None) => values[index] = None,
        (PropertyColumn::Bool(values), Some(PropertyValue::Bool(value))) => {
            values[index] = Some(value)
        }
        (PropertyColumn::Int(values), None) => values[index] = None,
        (PropertyColumn::Int(values), Some(PropertyValue::Int(value))) => {
            values[index] = Some(value)
        }
        (PropertyColumn::String(values), None) => values[index] = None,
        (PropertyColumn::String(values), Some(PropertyValue::String(value))) => {
            values[index] = Some(value)
        }
        (PropertyColumn::Real { values, .. }, None) => values[index] = None,
        (PropertyColumn::Real { values, .. }, Some(PropertyValue::Real { value, unit: _ })) => {
            values[index] = Some(value)
        }
        _ => unreachable!("property cell assignment must be prevalidated"),
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PropertyError {
    InvalidKey(String),
    LengthMismatch {
        expected: usize,
        actual: usize,
    },
    InvalidIndex {
        len: usize,
        index: usize,
    },
    /// The row belongs to an atom or bond removed from an editing draft.
    RemovedRow {
        index: usize,
    },
    TypeMismatch {
        key: PropertyKey,
    },
    NonFiniteValue {
        index: Option<usize>,
    },
    Unit(UnitError),
}

impl fmt::Display for PropertyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey(key) => write!(
                formatter,
                "invalid property key {key:?}; use a 1-{MAX_PROPERTY_KEY_LEN} character ASCII identifier"
            ),
            Self::LengthMismatch { expected, actual } => {
                write!(formatter, "property column requires {expected} values, received {actual}")
            }
            Self::InvalidIndex { len, index } => {
                write!(formatter, "property index {index} is outside table length {len}")
            }
            Self::RemovedRow { index } => {
                write!(formatter, "property row {index} belongs to a removed draft entity")
            }
            Self::TypeMismatch { key } => write!(formatter, "property {key:?} has a different type"),
            Self::NonFiniteValue { index: Some(index) } => {
                write!(formatter, "real property value at index {index} must be finite")
            }
            Self::NonFiniteValue { index: None } => {
                formatter.write_str("real property value must be finite")
            }
            Self::Unit(error) => write!(formatter, "invalid property unit: {error}"),
        }
    }
}

impl std::error::Error for PropertyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for PropertyError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::AtomId;
    use crate::units::{ANGSTROM, DIMENSIONLESS, KELVIN, NANOMETER};

    fn key(value: &str) -> PropertyKey {
        PropertyKey::new(value).unwrap()
    }

    #[test]
    fn detached_realization_bond_rows_bind_once_and_then_must_match() {
        let score = key("score");
        let mut properties = RealizationProperties::new(3);
        properties
            .owner_mut()
            .insert(score.clone(), PropertyValue::Int(9))
            .unwrap();
        assert_eq!(properties.atoms().len(), 3);
        assert!(properties.bonds().is_empty());
        assert!(properties.validate_bond_rows(2).is_ok());
        properties.bind_bond_rows(2);
        RealizationPropertiesMut::new(&mut properties)
            .bonds_mut()
            .set_value(
                score.clone(),
                TopologyBondIndex::new(1),
                Some(PropertyValue::Int(3)),
            )
            .unwrap();
        assert_eq!(
            properties.validate_bond_rows(1),
            Err(PropertyError::LengthMismatch {
                expected: 1,
                actual: 2
            })
        );
        let projected = properties.project(&[2, 0], &[1]).unwrap();
        assert!(projected.owner().is_empty());
        assert_eq!(
            projected
                .bonds()
                .value(&score, TopologyBondIndex::new(0))
                .unwrap(),
            Some(PropertyValue::Int(3))
        );
        assert_eq!(properties.owner().get(&score), Some(&PropertyValue::Int(9)));
    }

    #[test]
    fn typed_tables_address_rows_and_batch_updates_are_transactional() {
        let tag = key("tag");
        let mut properties = MoleculeProperties::new(3, 0);
        let mut guard = MoleculePropertiesMut::new(&mut properties, &[], &[]);
        guard
            .atoms_mut()
            .set_values(
                tag.clone(),
                [
                    (AtomId::new(0), Some(PropertyValue::Int(1))),
                    (AtomId::new(2), Some(PropertyValue::Int(3))),
                ],
            )
            .unwrap();
        let before = guard.atoms().clone();
        assert!(matches!(
            guard.atoms_mut().set_values(
                tag.clone(),
                [
                    (AtomId::new(1), Some(PropertyValue::Int(2))),
                    (AtomId::new(3), Some(PropertyValue::Int(4))),
                ],
            ),
            Err(PropertyError::InvalidIndex { len: 3, index: 3 })
        ));
        assert_eq!(guard.atoms(), &before);
        assert_eq!(
            properties.atoms().value(&tag, AtomId::new(2)).unwrap(),
            Some(PropertyValue::Int(3))
        );
        assert!(!properties.atoms().row_has_data(AtomId::new(1)).unwrap());
        assert_eq!(properties.atoms().keys().collect::<Vec<_>>(), [&tag]);
    }

    #[test]
    fn property_keys_use_one_conservative_validation_path() {
        for valid in ["x", "_private", "partial_charge", "force-field.v1", "a-2"] {
            assert_eq!(PropertyKey::new(valid).unwrap().as_str(), valid);
            assert_eq!(valid.parse::<PropertyKey>().unwrap().as_str(), valid);
        }
        for invalid in ["", "2bad", "has space", "slash/name", "unicode_µ"] {
            assert!(PropertyKey::new(invalid).is_err(), "accepted {invalid:?}");
        }
        assert!(PropertyKey::new("x".repeat(MAX_PROPERTY_KEY_LEN + 1)).is_err());
    }

    #[test]
    fn scalar_values_cover_the_complete_domain_and_reject_non_finite_reals() {
        let mut properties = OwnerProperties::new();
        properties
            .insert(key("bool"), PropertyValue::Bool(true))
            .unwrap();
        properties
            .insert(key("int"), PropertyValue::Int(-7))
            .unwrap();
        properties
            .insert(
                key("real"),
                PropertyValue::Real {
                    value: 2.5,
                    unit: DIMENSIONLESS,
                },
            )
            .unwrap();
        properties
            .insert(key("string"), PropertyValue::String("value".into()))
            .unwrap();
        assert_eq!(properties.iter().count(), 4);
        assert_eq!(
            properties
                .iter()
                .map(|(key, _)| key.as_str())
                .collect::<Vec<_>>(),
            ["bool", "int", "real", "string"]
        );
        assert!(matches!(
            properties.insert(
                key("nan"),
                PropertyValue::Real {
                    value: f64::NAN,
                    unit: DIMENSIONLESS,
                },
            ),
            Err(PropertyError::NonFiniteValue { index: None })
        ));
    }

    #[test]
    fn every_column_type_supports_missing_values_and_deterministic_iteration() {
        let mut table = RawPropertyTable::new(3);
        table
            .insert(
                key("z_string"),
                PropertyColumn::String(vec![None, Some("x".into()), None]),
            )
            .unwrap();
        table
            .insert(
                key("a_bool"),
                PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
            )
            .unwrap();
        table
            .insert(key("m_int"), PropertyColumn::Int(vec![None, Some(3), None]))
            .unwrap();
        table
            .insert(
                key("r_real"),
                PropertyColumn::Real {
                    unit: DIMENSIONLESS,
                    values: vec![Some(1.0), None, Some(2.0)],
                },
            )
            .unwrap();
        assert_eq!(
            table
                .iter()
                .map(|(key, _)| key.as_str())
                .collect::<Vec<_>>(),
            ["a_bool", "m_int", "r_real", "z_string"]
        );
        assert_eq!(table.value(&key("m_int"), 0).unwrap(), None);
        assert_eq!(
            table.value(&key("m_int"), 1).unwrap(),
            Some(PropertyValue::Int(3))
        );
    }

    #[test]
    fn table_updates_are_transactional_unit_aware_and_normalize_all_missing() {
        let mut table = RawPropertyTable::new(2);
        let length_error = table.insert(key("bad_len"), PropertyColumn::Bool(vec![Some(true)]));
        assert!(matches!(
            length_error,
            Err(PropertyError::LengthMismatch { .. })
        ));
        assert!(!table.has_data());

        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: NANOMETER,
                    values: vec![Some(1.0), None],
                },
            )
            .unwrap();
        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: ANGSTROM,
                    values: vec![Some(5.0), Some(10.0)],
                },
            )
            .unwrap();
        let Some(PropertyValue::Real { value, unit }) = table.value(&key("distance"), 0).unwrap()
        else {
            panic!("distance should be a real value");
        };
        assert_eq!(unit, NANOMETER);
        assert!((value - 0.5).abs() < 1.0e-12);

        let before = table.clone();
        assert!(matches!(
            table.insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: KELVIN,
                    values: vec![Some(1.0), Some(2.0)],
                },
            ),
            Err(PropertyError::Unit(_))
        ));
        assert_eq!(table, before);
        assert!(matches!(
            table.insert(
                key("finite"),
                PropertyColumn::Real {
                    unit: DIMENSIONLESS,
                    values: vec![Some(f64::INFINITY), None],
                },
            ),
            Err(PropertyError::NonFiniteValue { index: Some(0) })
        ));

        table.set_value(key("distance"), 0, None).unwrap();
        table.set_value(key("distance"), 1, None).unwrap();
        assert!(table.get(&key("distance")).is_none());
        table
            .insert(key("missing"), PropertyColumn::Int(vec![None, None]))
            .unwrap();
        assert!(table.get(&key("missing")).is_none());
    }

    #[test]
    fn failed_single_cell_updates_leave_the_existing_column_unchanged() {
        let mut table = RawPropertyTable::new(3);
        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: NANOMETER,
                    values: vec![Some(1.0), Some(2.0), None],
                },
            )
            .unwrap();

        for invalid in [
            PropertyValue::Int(4),
            PropertyValue::Real {
                value: 3.0,
                unit: KELVIN,
            },
            PropertyValue::Real {
                value: f64::INFINITY,
                unit: NANOMETER,
            },
        ] {
            let before = table.clone();
            assert!(table.set_value(key("distance"), 1, Some(invalid)).is_err());
            assert_eq!(table, before);
        }
    }

    #[test]
    fn compatible_single_cell_updates_convert_without_changing_column_unit() {
        let mut table = RawPropertyTable::new(2);
        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: NANOMETER,
                    values: vec![Some(1.0), None],
                },
            )
            .unwrap();
        table
            .set_value(
                key("distance"),
                1,
                Some(PropertyValue::Real {
                    value: 5.0,
                    unit: ANGSTROM,
                }),
            )
            .unwrap();

        let Some(PropertyColumn::Real { unit, values }) = table.get(&key("distance")) else {
            panic!("distance should remain a real column");
        };
        assert_eq!(*unit, NANOMETER);
        assert_eq!(values[0], Some(1.0));
        assert!((values[1].unwrap() - 0.5).abs() < 1.0e-12);
    }

    #[test]
    fn successful_single_cell_update_keeps_the_column_allocation() {
        let mut table = RawPropertyTable::new(4);
        table
            .insert(
                key("labels"),
                PropertyColumn::String(vec![Some("a".into()), None, None, None]),
            )
            .unwrap();
        let before = match table.get(&key("labels")).unwrap() {
            PropertyColumn::String(values) => values.as_ptr(),
            _ => unreachable!(),
        };

        table
            .set_value(key("labels"), 3, Some(PropertyValue::String("d".into())))
            .unwrap();

        let after = match table.get(&key("labels")).unwrap() {
            PropertyColumn::String(values) => values.as_ptr(),
            _ => unreachable!(),
        };
        assert_eq!(after, before);
    }

    #[test]
    fn clearing_the_final_populated_cell_removes_the_column() {
        let mut table = RawPropertyTable::new(4);
        table
            .set_value(key("flag"), 2, Some(PropertyValue::Bool(true)))
            .unwrap();
        table.set_value(key("flag"), 2, None).unwrap();
        assert_eq!(table.get(&key("flag")), None);
        assert!(!table.has_data());
    }

    #[test]
    fn checked_projection_preserves_columns_and_missing_cells() {
        let mut table = RawPropertyTable::new(3);
        table
            .insert(
                key("flag"),
                PropertyColumn::Bool(vec![Some(true), None, Some(false)]),
            )
            .unwrap();
        let selected = table.select_indices(&[2, 1]).unwrap();
        assert_eq!(selected.len(), 2);
        assert_eq!(
            selected.value(&key("flag"), 0).unwrap(),
            Some(PropertyValue::Bool(false))
        );
        assert_eq!(selected.value(&key("flag"), 1).unwrap(), None);
        assert!(matches!(
            table.select_indices(&[3]),
            Err(PropertyError::InvalidIndex { .. })
        ));
    }

    #[test]
    fn optional_projection_is_column_wise_and_preserves_real_units() {
        let mut table = RawPropertyTable::new(3);
        table
            .insert(
                key("distance"),
                PropertyColumn::Real {
                    unit: ANGSTROM,
                    values: vec![Some(1.0), None, Some(3.0)],
                },
            )
            .unwrap();
        table
            .insert(
                key("label"),
                PropertyColumn::String(vec![Some("a".into()), None, Some("c".into())]),
            )
            .unwrap();

        let selected = table
            .select_optional_indices(&[Some(2), None, Some(0), Some(1)])
            .unwrap();
        assert_eq!(
            selected.get(&key("distance")),
            Some(&PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(3.0), None, Some(1.0), None],
            })
        );
        assert_eq!(
            selected.get(&key("label")),
            Some(&PropertyColumn::String(vec![
                Some("c".into()),
                None,
                Some("a".into()),
                None,
            ]))
        );
        assert!(matches!(
            table.select_optional_indices(&[Some(3)]),
            Err(PropertyError::InvalidIndex { .. })
        ));
    }

    #[test]
    fn generic_tables_do_not_reserve_realization_semantic_names() {
        let mut table = RawPropertyTable::new(1);
        table
            .set_value(key("occupancy"), 0, Some(PropertyValue::Int(7)))
            .unwrap();
        table
            .insert(
                key("b_factor"),
                PropertyColumn::String(vec![Some("generic".into())]),
            )
            .unwrap();
        assert_eq!(
            table.value(&key("occupancy"), 0).unwrap(),
            Some(PropertyValue::Int(7))
        );
        assert_eq!(
            table.value(&key("b_factor"), 0).unwrap(),
            Some(PropertyValue::String("generic".into()))
        );
    }
}
