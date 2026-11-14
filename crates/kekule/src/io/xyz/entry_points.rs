//! The public XYZ entry points that `crate::xyz` re-exports.

use crate::diagnostics::FormatError;
use crate::io::{XyzDocument, XyzParseOptions};

/// Parses XYZ frames without assigning molecular meaning. `options` sets
/// the resource limits.
///
/// Blank lines may follow the last frame; any other content after them
/// is a syntax error. Empty input has no frames.
///
/// # Errors
///
/// Returns a [`FormatError`] of format
/// [`Format::Xyz`](crate::diagnostics::Format::Xyz) and kind
/// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
/// when the input exceeds `options`,
/// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
/// for a malformed count line or a truncated frame, and
/// [`DiagnosticKind::InvalidValue`](crate::diagnostics::DiagnosticKind::InvalidValue)
/// for a frame of zero atoms or a coordinate that is not a finite number.
pub fn parse_str(input: &str, options: XyzParseOptions) -> Result<XyzDocument, FormatError> {
    crate::io::parse_xyz_str(input, options)
}

/// Reads and parses XYZ frames from `reader`.
///
/// # Errors
///
/// As [`parse_str`]. The configured byte limit is enforced while reading,
/// before the input is held in memory; a read failure or non-UTF-8 input
/// fails with [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn parse_reader(
    reader: impl std::io::Read,
    options: XyzParseOptions,
) -> Result<XyzDocument, FormatError> {
    crate::io::parse_xyz_reader(reader, options)
}
