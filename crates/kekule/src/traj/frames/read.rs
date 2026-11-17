//! Reading a frame table: [`FrameTable`], [`read_frame_table`],
//! [`read_frame_table_from`], and the growing [`FrameTableReader`].

use std::fs::File;
use std::io::{BufRead, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::properties::{PropertyColumn, PropertyKey};

use super::error::{FrameTableError, FrameTableErrorKind, Reason};
use super::format::{parse_header, ParsedHeader, Rows};
use super::header::FrameTableHeader;

/// A frame table read from a file: the header, and the complete rows as a
/// step and a column of cells per header column.
///
/// A final line without a line feed, and NUL bytes after the last one, are
/// the pending tail: a row its writer has not finished, or blocks lost at
/// power loss. The table leaves it out ([`Self::pending_bytes`]).
///
/// ```
/// use kekule::traj::frames::read_frame_table_from;
///
/// let text = "#kekule-frames\t1\n#column\tstep\tint\n#column\tbias\treal\tkJ/mol\nstep\tbias\n\
///             0\t1.5\n500\t\n1000\t2";
/// let table = read_frame_table_from(text.as_bytes())?;
/// assert_eq!(table.steps(), [0, 500]);
/// assert_eq!(table.pending_bytes(), 6);
/// assert_eq!(table.row_of_step(500), Some(1));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct FrameTable {
    header: FrameTableHeader,
    rows: Rows,
    pending_bytes: u64,
}

impl FrameTable {
    /// The table's header.
    pub const fn header(&self) -> &FrameTableHeader {
        &self.header
    }

    /// The number of complete rows.
    pub fn len(&self) -> usize {
        self.rows.steps.len()
    }

    /// Whether the table has no complete row.
    pub fn is_empty(&self) -> bool {
        self.rows.steps.is_empty()
    }

    /// The step of every row, strictly increasing.
    pub fn steps(&self) -> &[u64] {
        &self.rows.steps
    }

    /// The cells of the column `key`, one per row, in the column's type:
    /// text columns are string columns, and real cells are in the column's
    /// unit. `None` for a key the header does not list.
    pub fn column(&self, key: &PropertyKey) -> Option<&PropertyColumn> {
        self.header
            .position(key)
            .map(|index| &self.rows.columns[index])
    }

    /// The row of `step`, when the table has one.
    pub fn row_of_step(&self, step: u64) -> Option<usize> {
        self.rows.steps.binary_search(&step).ok()
    }

    /// The bytes after the last complete row, which the table leaves out.
    pub const fn pending_bytes(&self) -> u64 {
        self.pending_bytes
    }
}

/// A frame table parsed from bytes, with where its parts end.
pub(super) struct Parsed {
    pub(super) table: FrameTable,
    /// The bytes of the header, through the line naming the columns.
    pub(super) header_end: usize,
    /// The bytes of the complete rows, through the last line feed.
    pub(super) end: usize,
    /// The start of the last complete row, when there is one.
    pub(super) last_row: Option<usize>,
    /// The number of complete lines.
    pub(super) lines: usize,
    /// The end of each complete row.
    pub(super) row_ends: Vec<usize>,
}

/// Parses a whole frame table.
pub(super) fn parse(bytes: &[u8]) -> Result<Parsed, FrameTableError> {
    let ParsedHeader { header, end, lines } = parse_header(bytes)?;
    let complete = bytes[end..]
        .iter()
        .rposition(|&byte| byte == b'\n')
        .map_or(end, |last| end + last + 1);
    let mut rows = Rows::new(&header);
    rows.extend(&header, &bytes[end..complete], lines + 1)?;
    let mut row_ends = Vec::with_capacity(rows.steps.len());
    for (offset, &byte) in bytes[end..complete].iter().enumerate() {
        if byte == b'\n' {
            row_ends.push(end + offset + 1);
        }
    }
    let last_row = match row_ends.len() {
        0 => None,
        1 => Some(end),
        count => Some(row_ends[count - 2]),
    };
    Ok(Parsed {
        lines: lines + row_ends.len(),
        table: FrameTable {
            header,
            rows,
            pending_bytes: (bytes.len() - complete) as u64,
        },
        header_end: end,
        end: complete,
        last_row,
        row_ends,
    })
}

/// Reads the frame table at `path`: its header and complete rows.
///
/// # Errors
///
/// [`FrameTableErrorKind::Io`] when the file cannot be read,
/// [`FrameTableErrorKind::Pending`] when its header is not complete,
/// [`FrameTableErrorKind::Unsupported`] for another version or an unknown
/// header record, and [`FrameTableErrorKind::Corrupt`] for a complete line
/// that is not valid.
pub fn read_frame_table(path: impl AsRef<Path>) -> Result<FrameTable, FrameTableError> {
    let bytes = std::fs::read(path).map_err(FrameTableError::read)?;
    parse(&bytes).map(|parsed| parsed.table)
}

/// Reads a frame table from `input`, as [`read_frame_table`] reads a file.
///
/// # Errors
///
/// As [`read_frame_table`].
pub fn read_frame_table_from(mut input: impl BufRead) -> Result<FrameTable, FrameTableError> {
    let mut bytes = Vec::new();
    input
        .read_to_end(&mut bytes)
        .map_err(FrameTableError::read)?;
    parse(&bytes).map(|parsed| parsed.table)
}

/// A frame table that its writer is still growing, read again as rows are
/// appended.
///
/// [`Self::refresh`] parses only the rows appended since the last read. When
/// bytes already read change, as when a run is cut at a checkpoint and
/// continued, it fails with [`FrameTableErrorKind::Changed`]; open the file
/// again.
#[derive(Debug, Clone)]
pub struct FrameTableReader {
    path: PathBuf,
    table: FrameTable,
    header: Vec<u8>,
    /// The start of the last complete row, and its bytes.
    last_row: Option<(u64, Vec<u8>)>,
    /// The end of the complete rows.
    end: u64,
    lines: usize,
}

impl FrameTableReader {
    /// Opens the frame table at `path` and reads its complete rows.
    ///
    /// # Errors
    ///
    /// As [`read_frame_table`]; retry later on
    /// [`FrameTableErrorKind::Pending`].
    pub fn open(path: impl AsRef<Path>) -> Result<Self, FrameTableError> {
        let path = path.as_ref().to_path_buf();
        let bytes = std::fs::read(&path).map_err(FrameTableError::read)?;
        let parsed = parse(&bytes)?;
        Ok(Self {
            path,
            header: bytes[..parsed.header_end].to_vec(),
            last_row: parsed
                .last_row
                .map(|start| (start as u64, bytes[start..parsed.end].to_vec())),
            end: parsed.end as u64,
            lines: parsed.lines,
            table: parsed.table,
        })
    }

    /// The rows read so far.
    pub const fn table(&self) -> &FrameTable {
        &self.table
    }

    /// Reads the rows appended since the last read and returns how many
    /// there are. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`FrameTableErrorKind::Changed`] when the header or the last row read
    /// changed or the file is shorter, [`FrameTableErrorKind::Corrupt`] for a
    /// new complete line that is not valid, and
    /// [`FrameTableErrorKind::Io`] when the file cannot be read.
    pub fn refresh(&mut self) -> Result<usize, FrameTableError> {
        let changed = || FrameTableError::new(FrameTableErrorKind::Changed, Reason::Changed);
        let mut file = File::open(&self.path).map_err(FrameTableError::read)?;
        let mut header = vec![0; self.header.len()];
        if file.read_exact(&mut header).is_err() || header != self.header {
            return Err(changed());
        }
        if let Some((start, row)) = &self.last_row {
            let mut bytes = vec![0; row.len()];
            file.seek(SeekFrom::Start(*start))
                .map_err(FrameTableError::read)?;
            if file.read_exact(&mut bytes).is_err() || &bytes != row {
                return Err(changed());
            }
        }
        file.seek(SeekFrom::Start(self.end))
            .map_err(FrameTableError::read)?;
        let mut appended = Vec::new();
        file.read_to_end(&mut appended)
            .map_err(FrameTableError::read)?;
        let complete = appended
            .iter()
            .rposition(|&byte| byte == b'\n')
            .map_or(0, |last| last + 1);
        let added =
            self.table
                .rows
                .extend(&self.table.header, &appended[..complete], self.lines + 1)?;
        if added > 0 {
            let last = appended[..complete - 1]
                .iter()
                .rposition(|&byte| byte == b'\n')
                .map_or(0, |last| last + 1);
            self.last_row = Some((self.end + last as u64, appended[last..complete].to_vec()));
        }
        self.end += complete as u64;
        self.lines += added;
        self.table.pending_bytes = (appended.len() - complete) as u64;
        Ok(added)
    }
}
