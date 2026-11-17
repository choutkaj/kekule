//! Frame tables: a run's per-step records, such as energies, temperature,
//! and bias, in a text file beside its trajectory, keyed by MD step.
//!
//! A frame table is a tab-separated text file with one line per recorded
//! step; `docs/frame-table.md` specifies it. Its header lists the columns
//! after the step ([`ColumnSpec`], [`ColumnType`]), the column of each row's
//! time, and the run's [`Sampling`](crate::structure::Sampling) record
//! ([`FrameTableHeader`]). Real numbers keep every bit, a missing cell is
//! explicit, and a value that is not finite is written as missing.
//!
//! - Writing: [`create_frame_table`] or [`FrameTableWriter::new`] write the
//!   header, and [`FrameTableWriter::write_row`] writes and flushes one
//!   [`FrameRow`] per line, so a crash leaves at most one unfinished line.
//! - Continuing: [`append_frame_table`] keeps the rows a checkpoint covers
//!   ([`KeepRows`]), removes the rest, and returns a writer after them, so a
//!   continued run's file equals one written straight through
//!   ([`FrameTableAppendReport`]).
//! - Reading: [`read_frame_table`] and [`read_frame_table_from`] read the
//!   complete rows into a [`FrameTable`] and leave out an unfinished line;
//!   [`FrameTableReader`] follows a growing file and reads only new rows.
//! - Joining: [`StructureSet::join_frames`](crate::structure::StructureSet::join_frames)
//!   puts the real and integer columns into a set's columns by step and
//!   installs the sampling record ([`JoinReport`]). The frames of a
//!   trajectory are usually a subset of the table's rows.
//!
//! Every failure is a [`FrameTableError`], classified by
//! [`FrameTableErrorKind`].
//!
//! ```
//! use kekule::properties::{PropertyKey, PropertyValue};
//! use kekule::traj::frames::{
//!     append_frame_table, create_frame_table, read_frame_table, ColumnSpec, ColumnType,
//!     FrameTableHeader, KeepRows,
//! };
//! use kekule::units::KILOJOULE_PER_MOLE;
//!
//! let path = std::env::temp_dir().join("kekule-frames-doc.tsv");
//! let energy = PropertyKey::new("potential_energy")?;
//! let header = FrameTableHeader::new(vec![ColumnSpec::new(
//!     energy.clone(),
//!     ColumnType::Real(KILOJOULE_PER_MOLE),
//! )])?;
//! let mut writer = create_frame_table(&path, header.clone())?;
//! let mut row = writer.row();
//! for step in [0, 500, 1000] {
//!     row.set_step(step);
//!     row.set(&energy, PropertyValue::Real { value: -1.0, unit: KILOJOULE_PER_MOLE })?;
//!     writer.write_row(&row)?;
//! }
//! drop(writer);
//!
//! // A run restarted from its checkpoint at step 500 continues the table.
//! let (writer, report) = append_frame_table(&path, &header, KeepRows::ThroughStep(500))?;
//! assert_eq!(report.last_step(), Some(500));
//! assert_eq!(report.dropped_rows(), 1);
//! drop(writer);
//! assert_eq!(read_frame_table(&path)?.steps(), [0, 500]);
//! # std::fs::remove_file(&path)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod append;
mod error;
mod format;
mod header;
mod join;
mod read;
mod row;
mod write;

pub use append::{append_frame_table, FrameTableAppendReport, KeepRows};
pub use error::{FrameTableError, FrameTableErrorKind};
pub use header::{ColumnSpec, ColumnType, FrameTableHeader, FRAME_TABLE_VERSION};
pub use join::JoinReport;
pub use read::{read_frame_table, read_frame_table_from, FrameTable, FrameTableReader};
pub use row::FrameRow;
pub use write::{create_frame_table, FrameTableWriter};
