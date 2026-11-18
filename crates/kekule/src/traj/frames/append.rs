//! Continuing a frame table after a restart: [`append_frame_table`].

use std::fs::{File, OpenOptions};
use std::io::{Read, Seek, SeekFrom, Write};
use std::path::Path;

use super::error::{FrameTableError, FrameTableErrorKind, Reason};
use super::format::{header_text, parse_header, ParsedHeader, Rows};
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
    /// checkpoint does not write again. The rows are kept in order up to the
    /// first row of a later step or the first line that is not a valid row,
    /// such as one a crash damaged; that line and everything after it are
    /// removed.
    ThroughStep(u64),
    /// Every complete row. Every complete line must be a valid row.
    All,
}

/// What [`append_frame_table`] kept of a frame table, and what it removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrameTableAppendReport {
    kept_rows: usize,
    last_step: Option<u64>,
    dropped_lines: usize,
    dropped_bytes: u64,
}

impl FrameTableAppendReport {
    /// The rows the file still holds, which the writer appends after.
    pub const fn kept_rows(&self) -> usize {
        self.kept_rows
    }

    /// The step of the last kept row, when a row was kept. A run compares it
    /// with its checkpoint: an earlier step means rows the checkpoint counts
    /// are missing, as when a damaged line ended the kept rows.
    pub const fn last_step(&self) -> Option<u64> {
        self.last_step
    }

    /// The complete lines removed after the kept rows.
    pub const fn dropped_lines(&self) -> usize {
        self.dropped_lines
    }

    /// The bytes removed: those of the dropped lines and of the pending
    /// tail, or of the whole file when its header was not complete.
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
/// `header`: the same records with the same values, in any order, whatever
/// its `##` comments.
///
/// The file is cut after the kept rows, and the cut is synced to the disk,
/// before the writer is returned.
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
/// for a file whose header is not valid, or under [`KeepRows::All`] for a
/// complete line that is not a valid row. The file is unchanged then.
pub fn append_frame_table(
    path: impl AsRef<Path>,
    header: &FrameTableHeader,
    keep: KeepRows,
) -> Result<(FrameTableWriter<File>, FrameTableAppendReport), FrameTableError> {
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)
        .map_err(FrameTableError::read)?;
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .map_err(FrameTableError::read)?;
    let found = match parse_header(&bytes) {
        Ok(ParsedHeader { header, .. }) => header,
        Err(error) if error.kind() == FrameTableErrorKind::Pending => {
            cut(&mut file, 0)?;
            file.write_all(header_text(header).as_bytes())
                .and_then(|()| file.sync_data())
                .map_err(FrameTableError::write)?;
            let report = FrameTableAppendReport {
                kept_rows: 0,
                last_step: None,
                dropped_lines: 0,
                dropped_bytes: bytes.len() as u64,
            };
            let writer = FrameTableWriter::resumed(file, header.clone(), None);
            return Ok((writer, report));
        }
        Err(error) => return Err(error),
    };
    if &found != header {
        return Err(FrameTableError::new(
            FrameTableErrorKind::Mismatch,
            Reason::Mismatch,
        ));
    }
    let Kept {
        rows: kept,
        last_step,
        end,
        lines,
    } = match keep {
        KeepRows::ThroughStep(step) => kept_through(&bytes, header, step)?,
        KeepRows::All => {
            let parsed = parse(&bytes)?;
            Kept {
                rows: parsed.table.len(),
                last_step: parsed.table.steps().last().copied(),
                end: parsed.end,
                lines: 0,
            }
        }
    };
    cut(&mut file, end as u64)?;
    let report = FrameTableAppendReport {
        kept_rows: kept,
        last_step,
        dropped_lines: lines,
        dropped_bytes: (bytes.len() - end) as u64,
    };
    let writer = FrameTableWriter::resumed(file, header.clone(), last_step);
    Ok((writer, report))
}

/// The rows [`append_frame_table`] keeps and where they end.
struct Kept {
    rows: usize,
    last_step: Option<u64>,
    /// The end of the kept rows.
    end: usize,
    /// The complete lines after them.
    lines: usize,
}

/// The rows of `bytes`, in order, up to the first row of a step after `step`
/// or the first complete line that is not a valid row.
fn kept_through(
    bytes: &[u8],
    header: &FrameTableHeader,
    step: u64,
) -> Result<Kept, FrameTableError> {
    let ParsedHeader {
        end: header_end,
        lines: header_lines,
        ..
    } = parse_header(bytes)?;
    let mut rows = Rows::new(header);
    let mut end = header_end;
    let mut complete = bytes[header_end..]
        .split_inclusive(|&byte| byte == b'\n')
        .filter(|line| line.ends_with(b"\n"));
    let mut lines = 0;
    for line in complete.by_ref() {
        let first = line.split(|&byte| byte == b'\t').next().unwrap_or_default();
        let later = std::str::from_utf8(first)
            .ok()
            .and_then(|text| text.parse::<u64>().ok())
            .is_some_and(|row| row > step);
        if later
            || rows
                .extend(header, line, header_lines + rows.steps.len() + 1)
                .is_err()
        {
            lines = 1;
            break;
        }
        end += line.len();
    }
    Ok(Kept {
        rows: rows.steps.len(),
        last_step: rows.steps.last().copied(),
        end,
        lines: lines + complete.count(),
    })
}

/// Cuts `file` to `length` bytes, syncs the cut, and moves to its end.
fn cut(file: &mut File, length: u64) -> Result<(), FrameTableError> {
    file.set_len(length)
        .and_then(|()| file.sync_data())
        .and_then(|()| file.seek(SeekFrom::Start(length)).map(drop))
        .map_err(FrameTableError::write)
}
