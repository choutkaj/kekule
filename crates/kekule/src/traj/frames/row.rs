//! One row of a frame table being written: [`FrameRow`].

use std::sync::Arc;

use crate::properties::{PropertyKey, PropertyValue};
use crate::units::Quantity;

use super::error::{FrameTableError, Reason};
use super::header::{ColumnSpec, ColumnType};

/// One row of a frame table being written: its step and a cell per column,
/// each missing until set.
///
/// [`FrameTableWriter::row`](super::FrameTableWriter::row) gives an empty
/// row of its header. A real value is stored in its column's unit, and a
/// value that is not finite, before or after conversion, is stored as
/// missing, so a run never stops on a diverging energy.
#[derive(Debug, Clone, PartialEq)]
pub struct FrameRow {
    columns: Arc<[ColumnSpec]>,
    step: Option<u64>,
    cells: Vec<Option<PropertyValue>>,
}

impl FrameRow {
    pub(super) fn new(columns: Arc<[ColumnSpec]>) -> Self {
        let cells = vec![None; columns.len()];
        Self {
            columns,
            step: None,
            cells,
        }
    }

    /// Sets the row's MD step.
    pub fn set_step(&mut self, step: u64) {
        self.step = Some(step);
    }

    /// The row's MD step, once set.
    pub const fn step(&self) -> Option<u64> {
        self.step
    }

    /// Sets the cell of the column `key`.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Invalid`](super::FrameTableErrorKind::Invalid)
    /// for a key the header has no column of, a value of another type than
    /// the column's, and a real value in a unit of another dimension. The
    /// cell is then unchanged.
    pub fn set(&mut self, key: &PropertyKey, value: PropertyValue) -> Result<(), FrameTableError> {
        let index = self
            .columns
            .iter()
            .position(|column| column.key() == key)
            .ok_or_else(|| FrameTableError::invalid(Reason::NoColumn(key.clone())))?;
        let cell = match (self.columns[index].kind(), value) {
            (ColumnType::Bool, value @ PropertyValue::Bool(_))
            | (ColumnType::Int, value @ PropertyValue::Int(_))
            | (ColumnType::Text, value @ PropertyValue::String(_)) => Some(value),
            (ColumnType::Real(column), PropertyValue::Real { value, unit }) => {
                if !unit.is_compatible(column) {
                    return Err(FrameTableError::invalid(Reason::CellUnit(key.clone())));
                }
                Quantity::new(value, unit)
                    .value_in(column)
                    .ok()
                    .filter(|value| value.is_finite())
                    .map(|value| PropertyValue::Real {
                        value,
                        unit: column,
                    })
            }
            _ => return Err(FrameTableError::invalid(Reason::CellType(key.clone()))),
        };
        self.cells[index] = cell;
        Ok(())
    }

    /// The cell of the column `key`, when the header has the column and the
    /// cell is set.
    pub fn value(&self, key: &PropertyKey) -> Option<&PropertyValue> {
        let index = self.columns.iter().position(|column| column.key() == key)?;
        self.cells[index].as_ref()
    }

    /// Clears the step and every cell, for the next row.
    pub fn clear(&mut self) {
        self.step = None;
        self.cells.iter_mut().for_each(|cell| *cell = None);
    }

    pub(super) fn columns(&self) -> &Arc<[ColumnSpec]> {
        &self.columns
    }

    pub(super) fn cells(&self) -> &[Option<PropertyValue>] {
        &self.cells
    }
}
