//! The public Molfile entry points that `crate::molfile` re-exports.

use crate::diagnostics::FormatError;
use crate::io::{MolfileDocument, MolfileParseOptions, MolfileSource, MolfileWriteOptions};

/// Parses one V2000 or V3000 Molfile into its typed connection table
/// without assigning canonical meaning. `options` bounds the input and the
/// parsed table.
///
/// # Errors
///
/// Returns a [`FormatError`] of kind
/// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
/// when the input exceeds `options`, and `Syntax`, `InvalidValue`,
/// `InvalidReference`, `Inconsistent` or `Unsupported` when the connection
/// table is malformed or uses content Kekule does not read.
pub fn parse_str(
    input: &str,
    options: MolfileParseOptions,
) -> Result<MolfileDocument, FormatError> {
    crate::io::parse_molfile_document(input, options)
}

/// Reads and parses one Molfile from `reader`.
///
/// # Errors
///
/// As [`parse_str`]. The configured byte limit is enforced while reading;
/// a read failure or non-UTF-8 input fails with
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn parse_reader(
    reader: impl std::io::Read,
    options: MolfileParseOptions,
) -> Result<MolfileDocument, FormatError> {
    crate::io::parse_molfile_reader(reader, options)
}

/// Writes one CTAB from a molecule (zero coordinates) or a structure, with
/// the [`WriteReport`](crate::diagnostics::WriteReport) of what it left
/// out under [`UnrepresentablePolicy::Omit`](crate::diagnostics::UnrepresentablePolicy::Omit):
/// a structure's cell, a geometry's state, occupancies, B factors, and
/// properties ([`MolfileWriteOptions::cell`],
/// [`MolfileWriteOptions::geometry_metadata`],
/// [`MolfileWriteOptions::properties`]).
///
/// [`MolfileWriteVersion::Auto`](crate::molfile::MolfileWriteVersion::Auto) uses V2000 when it can represent the
/// chemistry and V3000 otherwise. V2000 rounds coordinates to four decimal
/// places in angstroms; V3000 preserves round-trip decimal text. Stereo is
/// projected against the coordinates emitted by the chosen version, so
/// specified stereo requires a structure.
///
/// # Errors
///
/// Returns a [`FormatError`] of kind
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// when the selected version cannot represent the chemistry, or
/// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
/// when the source is incomplete.
pub fn write<'a>(
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<(String, crate::diagnostics::WriteReport), FormatError> {
    crate::io::write_molfile(source, options)
}

/// Writes the CTAB [`write()`] would return and returns its report. It
/// renders the whole CTAB first, since automatic version selection
/// learns whether V2000 can hold it only by rendering it, so nothing
/// reaches `writer` on failure.
///
/// # Errors
///
/// The errors of [`write()`]; a failing `writer` gives
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn write_to<'a>(
    writer: &mut impl std::io::Write,
    source: impl Into<MolfileSource<'a>>,
    options: MolfileWriteOptions,
) -> Result<crate::diagnostics::WriteReport, FormatError> {
    crate::io::write_molfile_to(writer, source, options)
}
