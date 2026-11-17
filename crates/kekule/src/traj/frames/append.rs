//! Continuing a frame table after a restart: [`append_frame_table`].

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use super::error::{FrameTableError, FrameTableErrorKind, Reason};
use super::format::header_text;
use super::header::FrameTableHeader;
use super::read::parse;
use super::write::FrameTableWriter;

/// Which rows of a frame table [`append_frame_table`] keeps; it removes every
/// row after them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum KeepRows {
    /// The rows whose step is at most this one: those a run had written when
    /// it checkpointed at this step, which a run restarted from the
    /// checkpoint does not write again.
    ThroughStep(u64),
    /// Every complete row.
    All,
}

/// What [`append_frame_table`] kept of a frame table, and what it removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrameTableAppendReport {
    kept_rows: usize,
    last_step: Option<u64>,
    dropped_rows: usize,
    dropped_bytes: u64,
}

impl FrameTableAppendReport {
    /// The rows the file still holds, which the writer appends after.
    pub const fn kept_rows(&self) -> usize {
        self.kept_rows
    }

    /// The step of the last kept row, when a row was kept; a run compares
    /// it with its checkpoint to see that no row it counts is missing.
    pub const fn last_step(&self) -> Option<u64> {
        self.last_step
    }

    /// The complete rows removed after the kept ones.
    pub const fn dropped_rows(&self) -> usize {
        self.dropped_rows
    }

    /// The bytes removed: those of the dropped rows and of the pending tail,
    /// or of the whole file when its header was not complete.
    pub const fn dropped_bytes(&self) -> u64 {
        self.dropped_bytes
    }
}

/// Opens the frame table at `path` to continue it, as a run restarted from a
/// checkpoint does: keeps the rows `keep` selects, removes the rest and the
/// pending tail, and returns a writer that appends after the kept rows, with
/// a report of what was removed.
///
/// A run checkpointed at step N and restarted from there passes
/// [`KeepRows::ThroughStep`]`(N)` and writes the rows after the checkpoint
/// as it computes them; the file then holds, byte for byte, what one
/// uninterrupted run would have written. A file whose header is not complete
/// yet is started again with `header`. Otherwise its header must equal
/// `header`, and every complete row must be valid; the scan reads the whole
/// file.
///
/// The file is cut after the kept rows before the writer is returned.
/// Appending survives a process that dies at any point: before the cut, the
/// file is as it was; after it, it holds the kept rows and maybe some
/// appended ones, possibly with an unfinished line, which another call takes
/// as this one took the first file. A [`FrameTableReader`](super::FrameTableReader)
/// following the file whose rows were removed gets
/// [`FrameTableErrorKind::Changed`] from its next refresh. No other process
/// may write the file meanwhile.
///
/// # Errors
///
/// [`FrameTableErrorKind::Io`] when the file cannot be read, cut, or
/// written; [`FrameTableErrorKind::Mismatch`] when its header differs from
/// `header`; and the errors of [`read_frame_table`](super::read_frame_table)
/// for a file that is not a valid frame table. The file is unchanged then.
pub fn append_frame_table(
    path: impl AsRef<Path>,
    header: &FrameTableHeader,
    keep: KeepRows,
) -> Result<(FrameTableWriter<BufWriter<File>>, FrameTableAppendReport), FrameTableError> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)
        .map_err(FrameTableError::read)?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .map_err(FrameTableError::read)?;
    let parsed = match parse(&bytes) {
        Ok(parsed) => parsed,
        Err(error) if error.kind() == FrameTableErrorKind::Pending => {
            cut(&mut file, 0)?;
            file.write_all(header_text(header).as_bytes())
                .and_then(|()| file.flush())
                .map_err(FrameTableError::write)?;
            let report = FrameTableAppendReport {
                kept_rows: 0,
                last_step: None,
                dropped_rows: 0,
                dropped_bytes: bytes.len() as u64,
            };
            let writer = FrameTableWriter::resumed(BufWriter::new(file), header.clone(), None);
            return Ok((writer, report));
        }
        Err(error) => return Err(error),
    };
    if parsed.table.header() != header {
        return Err(FrameTableError::new(
            FrameTableErrorKind::Mismatch,
            Reason::Mismatch,
        ));
    }
    let steps = parsed.table.steps();
    let kept = match keep {
        KeepRows::ThroughStep(step) => steps.partition_point(|&row| row <= step),
        KeepRows::All => steps.len(),
    };
    let end = kept
        .checked_sub(1)
        .map_or(parsed.header_end, |last| parsed.row_ends[last]);
    cut(&mut file, end as u64)?;
    let last_step = kept.checked_sub(1).map(|last| steps[last]);
    let report = FrameTableAppendReport {
        kept_rows: kept,
        last_step,
        dropped_rows: steps.len() - kept,
        dropped_bytes: (bytes.len() - end) as u64,
    };
    let writer = FrameTableWriter::resumed(BufWriter::new(file), header.clone(), last_step);
    Ok((writer, report))
}

/// Cuts `file` to `length` bytes and moves to its end.
fn cut(file: &mut File, length: u64) -> Result<(), FrameTableError> {
    file.set_len(length)
        .and_then(|()| file.seek(SeekFrom::Start(length)).map(drop))
        .map_err(FrameTableError::write)
}
