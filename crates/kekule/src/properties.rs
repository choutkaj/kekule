//! Unified scalar and columnar annotations for canonical Kekule objects.
//!
//! Property scope follows the object that owns it, and each scope has its own
//! type: definition-invariant annotations are [`MoleculeProperties`] of a
//! [`crate::molecule::Molecule`], system annotations are [`TopologyProperties`] of a
//! [`crate::topology::Topology`], and coordinate-dependent annotations are
//! [`GeometryProperties`] of a [`crate::structure::Structure`], geometry,
//! or trajectory frame. Each scope holds [`OwnerProperties`] plus one
//! [`PropertyTable`] per entity domain, addressed by that domain's identifier
//! ([`PropertyRow`]).
//!
//! Owners hand out read access through `properties()` and length-preserving
//! mutable access through `properties_mut()`. Row counts are fixed by the
//! owner and can never be changed through a property guard. Every table,
//! published or draft, is written through one guard, [`PropertyTableMut`],
//! addressed by a [`TableRow`]: an identifier, or the editing handle of a
//! topology or structure draft.
//!
//! Properties are annotations, not represented graph chemistry. Changing a
//! generic property does not change molecular identity or trigger chemical
//! perception.
//!
//! Every column of a table has one type, and a real column one storage unit;
//! [`PropertyValue::Real`] values written to it are converted to that unit.
//! Columns whose cells are all missing are removed rather than stored.
//!
//! # Example
//!
//! ```
//! use kekule::molecule::AtomId;
//! use kekule::properties::{PropertyKey, PropertyValue};
//! use kekule::units::{ELEMENTARY_CHARGE};
//!
//! let mut molecule = kekule::smiles::to_molecules("CO")?.remove(0);
//! let charge = PropertyKey::new("partial_charge")?;
//! let oxygen = AtomId::new(1);
//! molecule.properties_mut().atoms_mut().set_value(
//!     charge.clone(),
//!     oxygen,
//!     Some(PropertyValue::real(-0.5, ELEMENTARY_CHARGE)?),
//! )?;
//!
//! let atoms = molecule.properties().atoms();
//! assert_eq!(atoms.len(), 2);
//! assert_eq!(
//!     atoms.value(&charge, oxygen)?,
//!     Some(PropertyValue::Real { value: -0.5, unit: ELEMENTARY_CHARGE }),
//! );
//! assert_eq!(atoms.value(&charge, AtomId::new(0))?, None);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use std::collections::BTreeMap;
use std::fmt;
use std::marker::PhantomData;
use std::str::FromStr;

use crate::molecule::{AtomId, BondId};
use crate::topology::{
    AtomSiteId, ChainId, MoleculeInstanceId, ResidueId, TopologyAtomIndex, TopologyBondIndex,
};
use crate::units::{Unit, UnitError};

/// Maximum length of a [`PropertyKey`], in bytes.
pub const MAX_PROPERTY_KEY_LEN: usize = 128;

/// A validated, deterministic property identifier.
///
/// A key is 1 to [`MAX_PROPERTY_KEY_LEN`] ASCII bytes: a letter or `_`,
/// followed by letters, digits, `_`, `-` or `.`. Keys order and compare
/// byte-wise.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PropertyKey(String);

impl PropertyKey {
    /// Validates `value` as a key.
    ///
    /// # Errors
    ///
    /// [`PropertyError::InvalidKey`] when `value` is not a valid key.
    pub fn new(value: impl Into<String>) -> Result<Self, PropertyError> {
        let value = value.into();
        if !valid_property_key(&value) {
            return Err(PropertyError::InvalidKey(value));
        }
        Ok(Self(value))
    }

    /// The key text.
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
///
/// Stored real values are finite; [`Self::real`] validates one, and storing
/// a value validates it.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PropertyValue {
    /// A boolean.
    Bool(bool),
    /// A signed integer.
    Int(i64),
    /// A real number with its physical unit.
    Real {
        /// The magnitude, in `unit`.
        value: f64,
        /// The unit of `value`; [`DIMENSIONLESS`](crate::units::DIMENSIONLESS)
        /// for a pure number.
        unit: Unit,
    },
    /// A UTF-8 string.
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

    /// A [`Self::Real`] value.
    ///
    /// # Errors
    ///
    /// [`PropertyError::NonFiniteValue`] with `index: None` when `value` is NaN
    /// or infinite.
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
#[non_exhaustive]
pub enum PropertyValueRef<'a> {
    /// A boolean.
    Bool(bool),
    /// A signed integer.
    Int(i64),
    /// A real number with its physical unit.
    Real {
        /// The magnitude, in `unit`.
        value: f64,
        /// The unit of `value`.
        unit: Unit,
    },
    /// A borrowed UTF-8 string.
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
///
/// Element `i` is the cell of row `i`; `None` is a missing cell. Storing a
/// column requires every real value to be finite.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PropertyColumn {
    /// Boolean cells.
    Bool(Vec<Option<bool>>),
    /// Signed integer cells.
    Int(Vec<Option<i64>>),
    /// Real cells sharing one unit.
    Real {
        /// The unit of every value.
        unit: Unit,
        /// The magnitudes, in `unit`.
        values: Vec<Option<f64>>,
    },
    /// UTF-8 string cells.
    String(Vec<Option<String>>),
}

impl PropertyColumn {
    /// The number of cells, missing ones included.
    pub fn len(&self) -> usize {
        match self {
            Self::Bool(values) => values.len(),
            Self::Int(values) => values.len(),
            Self::Real { values, .. } => values.len(),
            Self::String(values) => values.len(),
        }
    }

    /// Whether the column has no cells.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether every cell is missing; true for an empty column.
    pub fn is_all_missing(&self) -> bool {
        match self {
            Self::Bool(values) => values.iter().all(Option::is_none),
            Self::Int(values) => values.iter().all(Option::is_none),
            Self::Real { values, .. } => values.iter().all(Option::is_none),
            Self::String(values) => values.iter().all(Option::is_none),
        }
    }

    /// Reads a scalar, copying string storage. See [`Self::value_ref`] to borrow.
    ///
    /// # Errors
    ///
    /// As [`Self::value_ref`].
    pub fn value(&self, index: usize) -> Result<Option<PropertyValue>, PropertyError> {
        self.value_ref(index)
            .map(|value| value.map(PropertyValueRef::to_value))
    }

    /// Borrows one cell without allocating; a missing cell is `None`.
    ///
    /// # Errors
    ///
    /// [`PropertyError::InvalidIndex`] when `index` is outside the column.
    pub fn value_ref(&self, index: usize) -> Result<Option<PropertyValueRef<'_>>, PropertyError> {
        if index >= self.len() {
            return Err(PropertyError::InvalidIndex(IndexOutOfRange::new(
                Countable::Row,
                index,
                self.len(),
            )));
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

    /// Checks that `value` can be stored in this column under `key`, as
    /// [`PropertyTableMut::set_value`] would store it.
    pub(crate) fn check_value(
        &self,
        key: &PropertyKey,
        value: &PropertyValue,
    ) -> Result<(), PropertyError> {
        validate_value_for_column(self, key, 0, Some(value.clone())).map(drop)
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

    fn remove(&mut self, index: usize) {
        match self {
            Self::Bool(values) => {
                values.remove(index);
            }
            Self::Int(values) => {
                values.remove(index);
            }
            Self::Real { values, .. } => {
                values.remove(index);
            }
            Self::String(values) => {
                values.remove(index);
            }
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

/// Row addressing shared by [`PropertyTableMut`] and the draft tables beside
/// it. The items are nominally public so the sealed [`TableRow`] can name them.
pub(crate) mod rows {
    use super::{PropertyError, TableRow};

    /// How a row key reaches its storage row.
    pub trait RowKey: Copy {
        /// The row type of the table that stores these rows.
        type Storage;

        /// The storage row of `self` among `rows`.
        fn slot(self, rows: LiveRows<'_>) -> Result<usize, <Self as TableRow>::Error>
        where
            Self: TableRow;
    }

    /// The row-liveness policy of a table: every row of a published owner,
    /// the live slots of a molecule draft, whose removed rows accept missing
    /// values alone, or the live handles of a topology or structure draft.
    #[derive(Clone, Copy)]
    pub enum LiveRows<'a> {
        All,
        Slots(&'a dyn RowLiveness),
        Handles(&'a dyn HandleSlots),
    }

    /// Which slots of a draft identifier space are live.
    pub trait RowLiveness {
        fn is_live(&self, row: usize) -> bool;
    }

    /// The storage slots of a draft's editing handles.
    pub trait HandleSlots {
        /// The slot of the live handle at `position`, or `None`.
        fn slot(&self, position: u64) -> Option<usize>;
        /// Live slots in handle order. Publication keeps it, except that an
        /// added atom joins its residue's atoms.
        fn live_slots(&self) -> Vec<usize>;
    }

    impl std::fmt::Debug for LiveRows<'_> {
        fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            formatter.write_str(match self {
                Self::All => "All",
                Self::Slots(_) => "Slots",
                Self::Handles(_) => "Handles",
            })
        }
    }

    impl LiveRows<'_> {
        pub(crate) fn handle_slot(self, position: u64) -> Option<usize> {
            match self {
                Self::Handles(handles) => handles.slot(position),
                Self::All | Self::Slots(_) => None,
            }
        }

        pub(super) fn check(self, row: usize) -> Result<(), PropertyError> {
            match self {
                Self::Slots(slots) if !slots.is_live(row) => {
                    Err(PropertyError::RemovedRow { index: row })
                }
                _ => Ok(()),
            }
        }

        pub(super) fn check_column(
            self,
            column: &super::PropertyColumn,
        ) -> Result<(), PropertyError> {
            if !matches!(self, Self::Slots(_)) {
                return Ok(());
            }
            (0..column.len())
                .filter(|row| column.is_populated_at(*row))
                .try_for_each(|row| self.check(row))
        }
    }
}
pub(crate) use rows::{HandleSlots, LiveRows, RowLiveness};

/// The row key of a [`PropertyTableMut`]: a [`PropertyRow`] identifier, or an
/// editing handle of a topology or structure draft (`EditAtomId`, ...).
pub trait TableRow: rows::RowKey {
    /// What a write or a read of one row fails with: [`PropertyError`] for
    /// identifiers, the editor's error for editing handles.
    type Error: std::error::Error + From<PropertyError>;
}

/// Typed row address of one property table.
///
/// Each table is indexed by the identifier of the domain it annotates, so a
/// molecule atom table cannot be read with a topology index by accident.
pub trait PropertyRow: TableRow<Error = PropertyError> + rows::RowKey<Storage = Self> {
    /// Zero-based row in the owning table.
    fn row(self) -> usize;
}

macro_rules! property_rows {
    ($($row:ty),* $(,)?) => {$(
        impl rows::RowKey for $row {
            type Storage = Self;

            fn slot(self, _: LiveRows<'_>) -> Result<usize, PropertyError> {
                Ok(self.index())
            }
        }
        impl TableRow for $row {
            type Error = PropertyError;
        }
        impl PropertyRow for $row {
            fn row(self) -> usize {
                self.index()
            }
        }
    )*};
}

property_rows!(
    crate::molecule::AtomId,
    crate::molecule::BondId,
    crate::topology::MoleculeInstanceId,
    crate::topology::TopologyAtomIndex,
    crate::topology::TopologyBondIndex,
    crate::topology::ChainId,
    crate::topology::ResidueId,
    crate::topology::AtomSiteId,
    crate::structure::GeometryIndex,
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

    /// Grows or shrinks the table to `len` rows; new rows are missing.
    pub(crate) fn resize_missing(&mut self, len: usize) {
        self.raw.resize_missing(len);
    }

    /// Removes row `index`, shifting later rows down. A column left without
    /// values is dropped.
    pub(crate) fn remove_row(&mut self, index: usize) {
        self.raw.remove_row(index);
    }

    /// The rows `indices`, each checked, in the given order; rows may repeat.
    pub(crate) fn select_rows(&self, indices: &[usize]) -> Result<Self, PropertyError> {
        self.raw.select_indices(indices).map(Self::from_raw)
    }

    /// Number of rows fixed by the owner.
    pub fn len(&self) -> usize {
        self.raw.len()
    }

    /// Whether the table has no rows.
    pub fn is_empty(&self) -> bool {
        self.raw.is_empty()
    }

    /// Whether any column is present. All-missing columns are never stored.
    pub fn has_data(&self) -> bool {
        self.raw.has_data()
    }

    /// The column stored under `key`, or `None` when no cell of it is
    /// populated.
    pub fn get(&self, key: &PropertyKey) -> Option<&PropertyColumn> {
        self.raw.get(key)
    }

    /// Columns in deterministic key order.
    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyColumn)> + DoubleEndedIterator {
        self.raw.iter()
    }

    /// Column keys in deterministic order.
    pub fn keys(&self) -> impl ExactSizeIterator<Item = &PropertyKey> + DoubleEndedIterator {
        self.raw.columns.keys()
    }
}

impl<R: PropertyRow> PropertyTable<R> {
    /// Whether any column has a value in `row`.
    ///
    /// # Errors
    ///
    /// [`PropertyError::InvalidIndex`] when `row` is outside the table.
    pub fn row_has_data(&self, row: R) -> Result<bool, PropertyError> {
        self.raw.row_has_data(row.row())
    }

    /// Reads a cell, copying string storage. See [`Self::value_ref`] to borrow.
    ///
    /// # Errors
    ///
    /// As [`Self::value_ref`].
    pub fn value(&self, key: &PropertyKey, row: R) -> Result<Option<PropertyValue>, PropertyError> {
        self.raw.value(key, row.row())
    }

    /// Borrows one cell without allocating.
    ///
    /// Missing keys and missing cells return `None`; a row outside the table
    /// is an error even when the key is absent.
    ///
    /// # Errors
    ///
    /// [`PropertyError::InvalidIndex`] when `row` is outside the table.
    pub fn value_ref(
        &self,
        key: &PropertyKey,
        row: R,
    ) -> Result<Option<PropertyValueRef<'_>>, PropertyError> {
        self.raw.value_ref(key, row.row())
    }
}

/// Mutable columns of an owner-sized table, addressed by row key `R`. The row
/// count cannot change.
///
/// With a [`PropertyRow`] identifier, columns are in row order and reads are
/// available through `Deref`; rows of atoms and bonds removed from a molecule
/// draft accept missing values alone and fail with
/// [`PropertyError::RemovedRow`]. With an editing handle of a topology or
/// structure draft, columns are read and inserted in live-handle order, which
/// publication keeps except that an added atom joins its residue's atoms, and
/// a deleted or foreign handle fails with the editor's matching error
/// (`InvalidAtom`, `InvalidBond`, ...).
///
/// Every mutation is transactional: on error the table is unchanged.
///
/// ```compile_fail
/// use kekule::{properties::PropertyTable, topology::TopologyBuilder};
/// let mut builder = TopologyBuilder::new();
/// let mut properties = builder.properties_mut();
/// *properties.atoms_mut() = PropertyTable::default();
/// ```
pub struct PropertyTableMut<'a, R: TableRow> {
    table: &'a mut PropertyTable<R::Storage>,
    live: LiveRows<'a>,
    // A draft's edit counter, advanced only when a write changes the table, so
    // an unchanged draft still publishes its source.
    revision: Option<&'a mut u64>,
    row: PhantomData<fn() -> R>,
}

impl<R: TableRow> fmt::Debug for PropertyTableMut<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.table.fmt(formatter)
    }
}

impl<R: PropertyRow> std::ops::Deref for PropertyTableMut<'_, R> {
    type Target = PropertyTable<R>;

    fn deref(&self) -> &Self::Target {
        self.table
    }
}

impl<'a, R: TableRow> PropertyTableMut<'a, R> {
    pub(crate) fn new(
        table: &'a mut PropertyTable<R::Storage>,
        live: LiveRows<'a>,
        revision: Option<&'a mut u64>,
    ) -> Self {
        Self {
            table,
            live,
            revision,
            row: PhantomData,
        }
    }

    /// The raw table and its row policy, for the draft reads beside this guard.
    pub(crate) fn parts(&self) -> (&RawPropertyTable, LiveRows<'_>) {
        (&self.table.raw, self.live)
    }

    fn changed(&mut self) {
        if let Some(revision) = self.revision.as_deref_mut() {
            *revision += 1;
        }
    }

    fn slot(&self, row: R) -> Result<usize, R::Error> {
        let slot = row.slot(self.live)?;
        if slot < self.table.len() {
            self.live.check(slot)?;
        }
        Ok(slot)
    }

    /// Inserts or replaces a complete column, returning the replaced one.
    ///
    /// Replacing a real column preserves its storage unit and converts
    /// compatible input values. An all-missing input removes the key.
    ///
    /// # Errors
    ///
    /// [`PropertyError::LengthMismatch`] when the column length differs from
    /// the number of rows (live rows, for editing handles),
    /// [`PropertyError::NonFiniteValue`] for a NaN or
    /// infinite real value, [`PropertyError::TypeMismatch`] when an existing
    /// column under `key` has another type, [`PropertyError::Unit`] when real
    /// values cannot be converted to the existing column's unit, and
    /// [`PropertyError::RemovedRow`] when a value lands in a removed row.
    pub fn insert(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, R::Error> {
        let LiveRows::Handles(handles) = self.live else {
            if column.len() == self.table.len() {
                self.live.check_column(&column)?;
            }
            let previous = self.table.raw.insert(key.clone(), column)?;
            if self.revision.is_some() && previous.as_ref() != self.table.raw.get(&key) {
                self.changed();
            }
            return Ok(previous);
        };
        let previous = live_column(&self.table.raw, handles, &key);
        install_column(
            &mut self.table.raw,
            &handles.live_slots(),
            key.clone(),
            column,
        )?;
        if previous != live_column(&self.table.raw, handles, &key) {
            self.changed();
        }
        Ok(previous)
    }

    /// Removes and returns the column under `key`, if present.
    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        let previous = match self.live {
            LiveRows::Handles(handles) => {
                let previous = live_column(&self.table.raw, handles, key);
                self.table.raw.remove(key);
                previous
            }
            LiveRows::All | LiveRows::Slots(_) => self.table.raw.remove(key),
        };
        if previous.is_some() {
            self.changed();
        }
        previous
    }

    /// Removes every column, keeping the row count.
    pub fn clear(&mut self) {
        if self.table.raw.has_data() {
            self.table.raw.clear_columns();
            self.changed();
        }
    }

    /// Sets or clears one cell, removing a column that becomes all-missing.
    ///
    /// Setting a cell of an absent key creates a column of the value's type
    /// (and unit, for a real value); a real value written to an existing real
    /// column is converted to its unit.
    ///
    /// # Errors
    ///
    /// [`PropertyError::InvalidIndex`] when `row` is outside the table,
    /// [`PropertyError::RemovedRow`] when it is a removed row,
    /// [`PropertyError::NonFiniteValue`] for a NaN or infinite real value,
    /// before or after conversion, [`PropertyError::TypeMismatch`] when the
    /// existing column has another type, and [`PropertyError::Unit`] when the
    /// value's unit is incompatible with the column's.
    pub fn set_value(
        &mut self,
        key: PropertyKey,
        row: R,
        value: Option<PropertyValue>,
    ) -> Result<(), R::Error> {
        let slot = self.slot(row)?;
        let tracked = self.revision.is_some();
        let previous = tracked
            .then(|| self.table.raw.value(&key, slot))
            .transpose()?;
        self.table.raw.set_value(key.clone(), slot, value)?;
        if tracked && previous != Some(self.table.raw.value(&key, slot)?) {
            self.changed();
        }
        Ok(())
    }

    /// Clears one cell; the same as [`Self::set_value`] with `None`. Clearing
    /// a cell of an absent key is a no-op.
    ///
    /// # Errors
    ///
    /// [`PropertyError::InvalidIndex`] when `row` is outside the table, and
    /// [`PropertyError::RemovedRow`] when it is a removed row.
    pub fn clear_value(&mut self, key: PropertyKey, row: R) -> Result<(), R::Error> {
        self.set_value(key, row, None)
    }

    /// Applies several cell updates to one column transactionally: either
    /// every update succeeds or the column is unchanged. Other columns are
    /// never copied. Updates apply in order, so a later update of a row wins.
    ///
    /// # Errors
    ///
    /// The first error [`Self::set_value`] would return for an update.
    pub fn set_values(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (R, Option<PropertyValue>)>,
    ) -> Result<(), R::Error> {
        let mut staged = self.table.raw.stage_column(&key);
        for (row, value) in values {
            staged.set_value(key.clone(), self.slot(row)?, value)?;
        }
        // An equal column may still differ in signed zeros or unit symbols;
        // only a draft that tracks publication keeps its stored column then.
        let changed = staged.get(&key) != self.table.raw.get(&key);
        if changed || self.revision.is_none() {
            self.table.raw.commit_column(key, staged);
        }
        if changed {
            self.changed();
        }
        Ok(())
    }
}

/// One column of a draft's live rows, in handle order.
pub(crate) fn live_column(
    table: &RawPropertyTable,
    handles: &dyn HandleSlots,
    key: &PropertyKey,
) -> Option<PropertyColumn> {
    table.get(key)?;
    table
        .select_indices(&handles.live_slots())
        .expect("live slots are allocated rows")
        .remove(key)
}

/// Stores a column given in live-handle order at the live slots `live`.
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
            return Err(PropertyError::LengthMismatch(LengthMismatch::new(
                Countable::Value,
                self.len,
                column.len(),
            )));
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

    /// Removes row `index` (checked by the caller), shifting later rows down
    /// and dropping a column left without values.
    pub(crate) fn remove_row(&mut self, index: usize) {
        debug_assert!(index < self.len);
        self.len -= 1;
        let populated = &mut self.populated;
        self.columns.retain(|key, column| {
            let was_populated = column.is_populated_at(index);
            column.remove(index);
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
            return Err(PropertyError::LengthMismatch(LengthMismatch::new(
                Countable::Value,
                source.len(),
                rows.len(),
            )));
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

    /// Whether [`Self::copy_rows_from`] takes every value of `source`: each
    /// source column whose key this table has matches its type, and its real
    /// values convert to its unit and stay finite.
    pub(crate) fn accepts_rows_of(&self, source: &Self) -> bool {
        source.columns.iter().all(|(key, column)| {
            let Some(own) = self.columns.get(key) else {
                return true;
            };
            match (own, column) {
                (PropertyColumn::Bool(_), PropertyColumn::Bool(_))
                | (PropertyColumn::Int(_), PropertyColumn::Int(_))
                | (PropertyColumn::String(_), PropertyColumn::String(_)) => true,
                (
                    PropertyColumn::Real { unit, .. },
                    PropertyColumn::Real { unit: from, values },
                ) => from.conversion_factor_to(*unit).is_ok_and(|factor| {
                    values
                        .iter()
                        .flatten()
                        .all(|value| (value * factor).is_finite())
                }),
                _ => false,
            }
        })
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
            return Err(PropertyError::InvalidIndex(IndexOutOfRange::new(
                Countable::Row,
                index,
                self.len,
            )));
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
    /// An empty set of owner values.
    pub fn new() -> Self {
        Self::default()
    }

    /// The number of stored values.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether no value is stored.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// The value stored under `key`, if any.
    pub fn get(&self, key: &PropertyKey) -> Option<&PropertyValue> {
        self.values.get(key)
    }

    /// Values in deterministic key order.
    pub fn iter(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyKey, &PropertyValue)> + DoubleEndedIterator {
        self.values.iter()
    }

    /// Keys in deterministic order.
    pub fn keys(&self) -> impl ExactSizeIterator<Item = &PropertyKey> + DoubleEndedIterator {
        self.values.keys()
    }

    /// Inserts a validated value, returning the one it replaces. Unlike a
    /// table column, an owner value may change type or unit on replacement.
    ///
    /// # Errors
    ///
    /// [`PropertyError::NonFiniteValue`] with `index: None` for a NaN or
    /// infinite real value; nothing is stored.
    pub fn insert(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        value.validate()?;
        Ok(self.values.insert(key, value))
    }

    /// Removes and returns the value under `key`, if present.
    pub fn remove(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.values.remove(key)
    }

    /// Removes every value.
    pub fn clear(&mut self) {
        self.values.clear();
    }
}

/// Annotations scoped to one molecule definition: owner values plus one row
/// per atom and bond slot, addressed by [`AtomId`] and [`BondId`].
///
/// Published molecules have dense IDs, so every row is live. In a
/// [`crate::molecule::MoleculeEditor`] draft, rows of removed atoms and bonds stay
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

    /// Owner-level values.
    pub fn owner(&self) -> &OwnerProperties {
        &self.owner
    }

    /// One row per atom slot, addressed by [`AtomId`].
    pub fn atoms(&self) -> &PropertyTable<AtomId> {
        &self.atoms
    }

    /// One row per bond slot, addressed by [`BondId`].
    pub fn bonds(&self) -> &PropertyTable<BondId> {
        &self.bonds
    }

    /// Whether no owner value and no column is stored.
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
    pub(crate) fn new(
        properties: &'a mut MoleculeProperties,
        atoms_live: LiveRows<'a>,
        bonds_live: LiveRows<'a>,
    ) -> Self {
        Self {
            properties,
            atoms_live,
            bonds_live,
        }
    }

    /// Mutable owner-level values.
    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties.owner
    }

    /// Mutable atom columns; rows of removed draft atoms accept only missing
    /// values.
    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, AtomId> {
        PropertyTableMut::new(&mut self.properties.atoms, self.atoms_live, None)
    }

    /// Mutable bond columns; rows of removed draft bonds accept only missing
    /// values.
    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, BondId> {
        PropertyTableMut::new(&mut self.properties.bonds, self.bonds_live, None)
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
    /// Owner-level values.
    pub fn owner(&self) -> &OwnerProperties {
        &self.owner
    }
    /// One row per molecule instance, addressed by [`MoleculeInstanceId`].
    pub fn molecule_instances(&self) -> &PropertyTable<MoleculeInstanceId> {
        &self.molecule_instances
    }
    /// One row per topology atom, addressed by its dense [`TopologyAtomIndex`].
    pub fn atoms(&self) -> &PropertyTable<TopologyAtomIndex> {
        &self.atoms
    }
    /// One row per topology bond, addressed by its dense [`TopologyBondIndex`].
    pub fn bonds(&self) -> &PropertyTable<TopologyBondIndex> {
        &self.bonds
    }
    /// One row per hierarchy chain, addressed by [`ChainId`].
    pub fn chains(&self) -> &PropertyTable<ChainId> {
        &self.chains
    }
    /// One row per hierarchy residue, addressed by [`ResidueId`].
    pub fn residues(&self) -> &PropertyTable<ResidueId> {
        &self.residues
    }
    /// One row per hierarchy atom site, addressed by [`AtomSiteId`].
    pub fn atom_sites(&self) -> &PropertyTable<AtomSiteId> {
        &self.atom_sites
    }

    /// Whether no owner value and no column is stored.
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

    // Typed tables, for the guards of an editing draft.
    pub(crate) fn atom_table_mut(&mut self) -> &mut PropertyTable<TopologyAtomIndex> {
        &mut self.atoms
    }
    pub(crate) fn bond_table_mut(&mut self) -> &mut PropertyTable<TopologyBondIndex> {
        &mut self.bonds
    }
    pub(crate) fn chain_table_mut(&mut self) -> &mut PropertyTable<ChainId> {
        &mut self.chains
    }
    pub(crate) fn residue_table_mut(&mut self) -> &mut PropertyTable<ResidueId> {
        &mut self.residues
    }
    pub(crate) fn atom_site_table_mut(&mut self) -> &mut PropertyTable<AtomSiteId> {
        &mut self.atom_sites
    }

    pub(crate) fn resize_atoms(&mut self, len: usize) {
        self.atoms.raw.resize_missing(len);
    }

    pub(crate) fn resize_bonds(&mut self, len: usize) {
        self.bonds.raw.resize_missing(len);
    }

    /// Whether every table takes the rows of `source`'s; see
    /// [`RawPropertyTable::accepts_rows_of`].
    pub(crate) fn accepts_rows_of(&self, source: &Self) -> bool {
        self.tables()
            .into_iter()
            .zip(source.tables())
            .all(|(own, from)| own.accepts_rows_of(from))
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
                return Err(PropertyError::LengthMismatch(LengthMismatch::new(
                    Countable::Value,
                    expected,
                    table.len(),
                )));
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

/// Length-preserving mutable access to staged or published [`TopologyProperties`].
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
    /// Mutable owner-level values.
    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties.owner
    }
    /// Mutable molecule-instance columns.
    pub fn molecule_instances_mut(&mut self) -> PropertyTableMut<'_, MoleculeInstanceId> {
        PropertyTableMut::new(&mut self.properties.molecule_instances, LiveRows::All, None)
    }
    /// Mutable topology-atom columns.
    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, TopologyAtomIndex> {
        PropertyTableMut::new(&mut self.properties.atoms, LiveRows::All, None)
    }
    /// Mutable topology-bond columns.
    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, TopologyBondIndex> {
        PropertyTableMut::new(&mut self.properties.bonds, LiveRows::All, None)
    }
    /// Mutable chain columns.
    pub fn chains_mut(&mut self) -> PropertyTableMut<'_, ChainId> {
        PropertyTableMut::new(&mut self.properties.chains, LiveRows::All, None)
    }
    /// Mutable residue columns.
    pub fn residues_mut(&mut self) -> PropertyTableMut<'_, ResidueId> {
        PropertyTableMut::new(&mut self.properties.residues, LiveRows::All, None)
    }
    /// Mutable atom-site columns.
    pub fn atom_sites_mut(&mut self) -> PropertyTableMut<'_, AtomSiteId> {
        PropertyTableMut::new(&mut self.properties.atom_sites, LiveRows::All, None)
    }

    /// Removes every owner value and column, keeping row counts.
    pub fn clear(&mut self) {
        self.properties.owner.clear();
        for table in self.properties.tables_mut() {
            table.clear_columns();
        }
    }
}

/// Coordinate-dependent annotations of one structure: owner values plus one
/// row per topology atom and bond, in dense topology order.
///
/// A detached structure payload has atom rows only; bond rows are allocated
/// when the payload is bound to a topology.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct GeometryProperties {
    owner: OwnerProperties,
    atoms: PropertyTable<TopologyAtomIndex>,
    bonds: PropertyTable<TopologyBondIndex>,
}

impl GeometryProperties {
    pub(crate) fn new(atoms: usize) -> Self {
        Self {
            owner: OwnerProperties::new(),
            atoms: PropertyTable::new(atoms),
            bonds: PropertyTable::new(0),
        }
    }

    /// Owner-level values.
    pub fn owner(&self) -> &OwnerProperties {
        &self.owner
    }

    /// One row per topology atom, addressed by its dense [`TopologyAtomIndex`].
    pub fn atoms(&self) -> &PropertyTable<TopologyAtomIndex> {
        &self.atoms
    }

    /// One row per topology bond, addressed by its dense
    /// [`TopologyBondIndex`]; empty until a detached payload is bound to a
    /// topology.
    pub fn bonds(&self) -> &PropertyTable<TopologyBondIndex> {
        &self.bonds
    }

    /// Whether no owner value and no column is stored.
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
        Err(PropertyError::LengthMismatch(LengthMismatch::new(
            Countable::Value,
            bond_count,
            self.bonds.len(),
        )))
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

    /// Whether the atom and bond tables take the rows of `source`'s, which
    /// has a row for each of `atoms` and `bonds`; see
    /// [`RawPropertyTable::accepts_rows_of`].
    pub(crate) fn accepts_rows_of(&self, source: &Self, atoms: usize, bonds: usize) -> bool {
        source.atoms.len() == atoms
            && source.bonds.len() == bonds
            && self.atoms.raw.accepts_rows_of(&source.atoms.raw)
            && self.bonds.raw.accepts_rows_of(&source.bonds.raw)
    }

    // Typed tables, for the guards of an editing draft.
    pub(crate) fn atom_table_mut(&mut self) -> &mut PropertyTable<TopologyAtomIndex> {
        &mut self.atoms
    }
    pub(crate) fn bond_table_mut(&mut self) -> &mut PropertyTable<TopologyBondIndex> {
        &mut self.bonds
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
    /// projected structure is a new owner.
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

/// Length-preserving mutable access to [`GeometryProperties`].
#[derive(Debug)]
pub struct GeometryPropertiesMut<'a> {
    properties: &'a mut GeometryProperties,
}

impl std::ops::Deref for GeometryPropertiesMut<'_> {
    type Target = GeometryProperties;

    fn deref(&self) -> &Self::Target {
        self.properties
    }
}

impl<'a> GeometryPropertiesMut<'a> {
    pub(crate) fn new(properties: &'a mut GeometryProperties) -> Self {
        Self { properties }
    }

    /// Mutable owner-level values.
    pub fn owner_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties.owner
    }

    /// Mutable topology-atom columns.
    pub fn atoms_mut(&mut self) -> PropertyTableMut<'_, TopologyAtomIndex> {
        PropertyTableMut::new(&mut self.properties.atoms, LiveRows::All, None)
    }

    /// Mutable topology-bond columns.
    pub fn bonds_mut(&mut self) -> PropertyTableMut<'_, TopologyBondIndex> {
        PropertyTableMut::new(&mut self.properties.bonds, LiveRows::All, None)
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

/// A table of per-item properties: the items each row describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PropertyDomain {
    /// Molecule instances of a topology.
    MoleculeInstances,
    /// Atoms of a topology.
    TopologyAtoms,
    /// Bonds of a topology.
    TopologyBonds,
    /// Hierarchy chains.
    Chains,
    /// Hierarchy residues.
    Residues,
    /// Hierarchy atom sites.
    AtomSites,
    /// Atoms of a geometry.
    GeometryAtoms,
    /// Bonds of a geometry.
    GeometryBonds,
}

impl fmt::Display for PropertyDomain {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::MoleculeInstances => "molecule instance",
            Self::TopologyAtoms => "topology atom",
            Self::TopologyBonds => "topology bond",
            Self::Chains => "chain",
            Self::Residues => "residue",
            Self::AtomSites => "atom site",
            Self::GeometryAtoms => "geometry atom",
            Self::GeometryBonds => "geometry bond",
        })
    }
}

/// A property key, column or value is invalid for its table.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PropertyError {
    /// The text is not a 1 to `MAX_PROPERTY_KEY_LEN` character ASCII identifier.
    InvalidKey(String),
    /// A column or row set has a different length than the table.
    LengthMismatch(LengthMismatch),
    /// A row index is outside the table.
    InvalidIndex(IndexOutOfRange),
    /// The row belongs to an atom or bond removed from an editing draft.
    RemovedRow {
        /// The removed row's index.
        index: usize,
    },
    /// A column or value has a different type than the existing column.
    TypeMismatch {
        /// The key of the existing column.
        key: PropertyKey,
    },
    /// A real value is NaN or infinite.
    NonFiniteValue {
        /// The row of the value, or `None` for a single value.
        index: Option<usize>,
    },
    /// A real column could not be converted to the existing column's unit; see
    /// the source.
    Unit(UnitError),
}

impl fmt::Display for PropertyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidKey(key) => write!(
                formatter,
                "The property key {key:?} is invalid. Use an ASCII identifier of 1 to {MAX_PROPERTY_KEY_LEN} characters."
            ),
            Self::LengthMismatch(error) => error.fmt(formatter),
            Self::InvalidIndex(error) => error.fmt(formatter),
            Self::RemovedRow { index } => write!(
                formatter,
                "Property row {index} belongs to a removed atom or bond."
            ),
            Self::TypeMismatch { key } => write!(
                formatter,
                "The property \"{key}\" already has another type."
            ),
            Self::NonFiniteValue { index: Some(index) } => {
                write!(formatter, "The property value in row {index} is not finite.")
            }
            Self::NonFiniteValue { index: None } => {
                formatter.write_str("The property value is not finite.")
            }
            Self::Unit(_) => formatter.write_str("Cannot convert the property to the column's unit."),
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
mod tests;
