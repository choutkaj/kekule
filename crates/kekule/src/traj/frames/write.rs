//! Writing a frame table: [`FrameTableWriter`] and [`create_frame_table`].

use std::fmt;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use super::error::{FrameTableError, Reason};
use super::format::{header_text, row_text};
use super::header::FrameTableHeader;
use super::row::FrameRow;

/// Writes the rows of a frame table, one complete line per row, in step
/// order.
///
/// The header is written when the writer is made. Each row is written and
/// flushed whole, so a reader following the file sees complete rows, and a
/// process that dies leaves at most one unfinished line, which readers
/// ignore and [`append_frame_table`](super::append_frame_table) removes.
/// For the rows to survive the failure of a node, sync the file
/// (`writer.get_ref().get_ref().sync_data()` for a file writer) before a
/// checkpoint counts them. After an I/O error the file may end in an
/// unfinished line; continue it with
/// [`append_frame_table`](super::append_frame_table), not with this writer.
///
/// ```
/// use kekule::properties::{PropertyKey, PropertyValue};
/// use kekule::traj::frames::{
///     read_frame_table_from, ColumnSpec, ColumnType, FrameTableHeader, FrameTableWriter,
/// };
/// use kekule::units::KILOJOULE_PER_MOLE;
///
/// let energy = PropertyKey::new("potential_energy")?;
/// let header = FrameTableHeader::new(vec![ColumnSpec::new(
///     energy.clone(),
///     ColumnType::Real(KILOJOULE_PER_MOLE),
/// )])?;
/// let mut writer = FrameTableWriter::new(Vec::new(), header)?;
/// let mut row = writer.row();
/// for (step, value) in [(0, -30521.25), (500, f64::NAN)] {
///     row.clear();
///     row.set_step(step);
///     row.set(&energy, PropertyValue::Real { value, unit: KILOJOULE_PER_MOLE })?;
///     writer.write_row(&row)?;
/// }
/// let table = read_frame_table_from(writer.into_inner()?.as_slice())?;
/// assert_eq!(table.steps(), [0, 500]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct FrameTableWriter<W: Write> {
    out: W,
    header: FrameTableHeader,
    last_step: Option<u64>,
    line: String,
}

impl<W: Write> fmt::Debug for FrameTableWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FrameTableWriter")
            .field("header", &self.header)
            .field("last_step", &self.last_step)
            .finish_non_exhaustive()
    }
}

impl<W: Write> FrameTableWriter<W> {
    /// Writes the header of a new table to `out` and returns a writer of its
    /// rows.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Io`](super::FrameTableErrorKind::Io) when
    /// writing fails.
    pub fn new(mut out: W, header: FrameTableHeader) -> Result<Self, FrameTableError> {
        out.write_all(header_text(&header).as_bytes())
            .and_then(|()| out.flush())
            .map_err(FrameTableError::write)?;
        Ok(Self::resumed(out, header, None))
    }

    /// A writer of rows after `last_step` in a table whose header `out`
    /// already holds.
    pub(super) fn resumed(out: W, header: FrameTableHeader, last_step: Option<u64>) -> Self {
        Self {
            out,
            header,
            last_step,
            line: String::new(),
        }
    }

    /// An empty row of this table: no step, every cell missing.
    pub fn row(&self) -> FrameRow {
        FrameRow::new(Arc::clone(self.header.shared_columns()))
    }

    /// Writes `row` as the next line and flushes it.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Invalid`](super::FrameTableErrorKind::Invalid)
    /// for a row of another header, without a step, or whose step does not
    /// exceed the last one, before anything is written, and
    /// [`FrameTableErrorKind::Io`](super::FrameTableErrorKind::Io) when
    /// writing fails.
    pub fn write_row(&mut self, row: &FrameRow) -> Result<(), FrameTableError> {
        let columns = self.header.shared_columns();
        if !Arc::ptr_eq(row.columns(), columns) && row.columns() != columns {
            return Err(FrameTableError::invalid(Reason::ForeignRow));
        }
        let step = row
            .step()
            .ok_or_else(|| FrameTableError::invalid(Reason::NoStep))?;
        if let Some(last) = self.last_step.filter(|&last| step <= last) {
            return Err(FrameTableError::invalid(Reason::StepOrder { step, last }));
        }
        self.line.clear();
        row_text(&mut self.line, step, columns, row.cells());
        self.out
            .write_all(self.line.as_bytes())
            .and_then(|()| self.out.flush())
            .map_err(FrameTableError::write)?;
        self.last_step = Some(step);
        Ok(())
    }

    /// Flushes the output.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Io`](super::FrameTableErrorKind::Io) when
    /// flushing fails.
    pub fn flush(&mut self) -> Result<(), FrameTableError> {
        self.out.flush().map_err(FrameTableError::write)
    }

    /// The output, for syncing a file.
    pub const fn get_ref(&self) -> &W {
        &self.out
    }

    /// Flushes and returns the output.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Io`](super::FrameTableErrorKind::Io) when
    /// flushing fails.
    pub fn into_inner(mut self) -> Result<W, FrameTableError> {
        self.flush()?;
        Ok(self.out)
    }

    /// The table's header.
    pub const fn header(&self) -> &FrameTableHeader {
        &self.header
    }

    /// The step of the last row the table holds, when it holds one.
    pub const fn last_step(&self) -> Option<u64> {
        self.last_step
    }
}

/// Creates the frame table at `path`, replacing an existing file, writes its
/// header, and returns a writer of its rows. Continue a table with
/// [`append_frame_table`](super::append_frame_table) instead.
///
/// # Errors
///
/// [`FrameTableErrorKind::Io`](super::FrameTableErrorKind::Io) when the file
/// cannot be created or written.
pub fn create_frame_table(
    path: impl AsRef<Path>,
    header: FrameTableHeader,
) -> Result<FrameTableWriter<BufWriter<File>>, FrameTableError> {
    let file = File::create(path).map_err(FrameTableError::write)?;
    FrameTableWriter::new(BufWriter::new(file), header)
}
