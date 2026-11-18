//! What a frame table holds: [`ColumnType`], [`ColumnSpec`], and
//! [`FrameTableHeader`].

use std::sync::Arc;

use crate::properties::PropertyKey;
use crate::structure::sampling::is_energy;
use crate::structure::Sampling;
use crate::units::{named_unit, Unit, CANONICAL_TIME_UNIT};

use super::error::{FrameTableError, Reason};

/// The version of the frame-table format this crate writes and reads.
pub const FRAME_TABLE_VERSION: u32 = 1;

/// The key of a frame table's step column, always its first.
pub(super) const STEP: &str = "step";

/// The type of a frame-table column's cells.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum ColumnType {
    /// `true` or `false`.
    Bool,
    /// A signed 64-bit integer.
    Int,
    /// A finite real number in this unit, which must be one of kekule's
    /// named unit constants, such as kJ/mol.
    Real(Unit),
    /// UTF-8 text.
    Text,
}

/// One column of a frame table: its key and the type of its cells.
#[derive(Debug, Clone, PartialEq)]
pub struct ColumnSpec {
    key: PropertyKey,
    kind: ColumnType,
}

impl ColumnSpec {
    /// The column `key` of `kind` cells.
    pub const fn new(key: PropertyKey, kind: ColumnType) -> Self {
        Self { key, kind }
    }

    /// The column's key.
    pub const fn key(&self) -> &PropertyKey {
        &self.key
    }

    /// The type of the column's cells.
    pub const fn kind(&self) -> ColumnType {
        self.kind
    }
}

/// The header of a frame table: its columns after the step, the column of
/// each row's time, and the sampling record of the run.
///
/// Every table starts with the column `step`, the MD step of each row, which
/// the header does not list. Steps strictly increase.
///
/// ```
/// use kekule::properties::PropertyKey;
/// use kekule::structure::{Sampling, ThermodynamicState};
/// use kekule::traj::frames::{ColumnSpec, ColumnType, FrameTableHeader};
/// use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE, PICOSECOND};
///
/// let key = |name: &str| PropertyKey::new(name).unwrap();
/// let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN))?;
/// let header = FrameTableHeader::new(vec![
///     ColumnSpec::new(key("time"), ColumnType::Real(PICOSECOND)),
///     ColumnSpec::new(key("potential_energy"), ColumnType::Real(KILOJOULE_PER_MOLE)),
///     ColumnSpec::new(key("bias"), ColumnType::Real(KILOJOULE_PER_MOLE)),
/// ])?
/// .with_time(key("time"))?
/// .with_sampling(Sampling::new(state).with_bias([key("bias")])?)?;
/// assert_eq!(header.columns().len(), 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct FrameTableHeader {
    columns: Arc<[ColumnSpec]>,
    time: Option<PropertyKey>,
    sampling: Option<Sampling>,
}

impl FrameTableHeader {
    /// A header of `columns`, after the step column, without a time column
    /// or a sampling record.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Invalid`](super::FrameTableErrorKind::Invalid)
    /// for a column named `step`, two columns of one key, and a real column
    /// whose unit is not one of kekule's named unit constants.
    pub fn new(columns: Vec<ColumnSpec>) -> Result<Self, FrameTableError> {
        for (index, column) in columns.iter().enumerate() {
            if column.key.as_str() == STEP {
                return Err(FrameTableError::invalid(Reason::ReservedStep));
            }
            if columns[..index].iter().any(|other| other.key == column.key) {
                return Err(FrameTableError::invalid(Reason::RepeatedColumn(
                    column.key.clone(),
                )));
            }
            if let ColumnType::Real(unit) = column.kind {
                if unit_symbol(unit).is_none() {
                    return Err(FrameTableError::invalid(Reason::UnnamedUnit(
                        column.key.clone(),
                    )));
                }
            }
        }
        Ok(Self {
            columns: columns.into(),
            time: None,
            sampling: None,
        })
    }

    /// Returns this header whose `column` holds each row's time.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Invalid`](super::FrameTableErrorKind::Invalid)
    /// when the header has no such column or it is not a real column of
    /// times.
    pub fn with_time(mut self, column: PropertyKey) -> Result<Self, FrameTableError> {
        match self.kind(&column) {
            None => return Err(FrameTableError::invalid(Reason::NoColumn(column))),
            Some(ColumnType::Real(unit)) if unit.is_compatible(CANONICAL_TIME_UNIT) => {}
            Some(_) => return Err(FrameTableError::invalid(Reason::NoTimes(column))),
        }
        self.time = Some(column);
        Ok(self)
    }

    /// Returns this header with the run's sampling record.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Invalid`](super::FrameTableErrorKind::Invalid)
    /// when a column the record names is not a real column of energies of
    /// the header.
    pub fn with_sampling(mut self, sampling: Sampling) -> Result<Self, FrameTableError> {
        for key in sampling.columns() {
            match self.kind(key) {
                None => return Err(FrameTableError::invalid(Reason::NoColumn(key.clone()))),
                Some(ColumnType::Real(unit)) if is_energy(unit) => {}
                Some(_) => return Err(FrameTableError::invalid(Reason::NoEnergies(key.clone()))),
            }
        }
        self.sampling = Some(sampling);
        Ok(self)
    }

    /// The columns after the step column, in file order.
    pub fn columns(&self) -> &[ColumnSpec] {
        &self.columns
    }

    /// The column of each row's time, when named.
    pub const fn time(&self) -> Option<&PropertyKey> {
        self.time.as_ref()
    }

    /// The run's sampling record, when the table has one.
    pub const fn sampling(&self) -> Option<&Sampling> {
        self.sampling.as_ref()
    }

    /// The position of `key` among the columns.
    pub(super) fn position(&self, key: &PropertyKey) -> Option<usize> {
        self.columns.iter().position(|column| &column.key == key)
    }

    /// The shared columns, which rows of this header hold.
    pub(super) fn shared_columns(&self) -> &Arc<[ColumnSpec]> {
        &self.columns
    }

    fn kind(&self, key: &PropertyKey) -> Option<ColumnType> {
        self.position(key).map(|index| self.columns[index].kind)
    }
}

/// The symbol of `unit` when it is a named unit constant that the symbol
/// restores.
pub(super) fn unit_symbol(unit: Unit) -> Option<&'static str> {
    let symbol = unit.symbol()?;
    (named_unit(symbol) == Some(unit)).then_some(symbol)
}
