//! The public XYZ entry points that `crate::xyz` re-exports.

use crate::diagnostics::{FormatError, WriteReport};
use crate::io::{XyzDocument, XyzParseOptions, XyzSource, XyzWriteOptions};

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

/// Writes a structure as one frame, or a structure set as one frame per
/// geometry, and reports what the text leaves out.
///
/// XYZ carries each atom's element and position only. The report holds the
/// topology content every frame leaves out: bonds; atoms whose
/// implicit-hydrogen count is nonzero, or inferred and not yet perceived,
/// since XYZ readers fix every count at zero; atoms with a nonzero formal
/// charge, an isotope, radical electrons, or an atom-map number; atoms whose
/// atom site has an atom name, lost with its residue and chain, unlike the
/// unnamed `UNL` sites of an interpreted XYZ file; and molecule instances
/// with represented stereo. It also holds the state omitted under the policies of `options`,
/// and the occupancies and B factors, which are always left out. The whole
/// source is checked before the first byte is written.
///
/// ```
/// use kekule::xyz;
///
/// let text = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n";
/// let structure = xyz::parse_str(text, Default::default())?
///     .interpret(Default::default())?
///     .into_structure();
/// let options = xyz::XyzWriteOptions::default()
///     .with_decimal_places(3)?
///     .with_comment("water")?;
/// let (written, report) = xyz::write(&structure, options)?;
/// assert_eq!(written, text);
/// assert_eq!(report.omissions().len(), 1); // the two inferred bonds
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// A write [`FormatError`] of format
/// [`Format::Xyz`](crate::diagnostics::Format::Xyz) and kind
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// for a topology without atoms, a coordinate that is not finite in the
/// output unit, a structure set's own properties, and state `options` do not
/// omit, the last with
/// [`FormatErrorDetail::UnrepresentableState`](crate::diagnostics::FormatErrorDetail::UnrepresentableState);
/// and of kind
/// [`DiagnosticKind::Inconsistent`](crate::diagnostics::DiagnosticKind::Inconsistent)
/// for a structure set without geometries. An error about one geometry of a
/// set names it as its 1-based [`record`](FormatError::record).
pub fn write<'a>(
    source: impl Into<XyzSource<'a>>,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), FormatError> {
    crate::io::xyz::write::write_xyz(source, options)
}

/// Streams the text [`write()`] would return and returns its report.
///
/// # Errors
///
/// The errors of [`write()`]; a failing `writer` gives
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io).
pub fn write_to<'a>(
    writer: &mut impl std::io::Write,
    source: impl Into<XyzSource<'a>>,
    options: XyzWriteOptions,
) -> Result<WriteReport, FormatError> {
    crate::io::xyz::write::write_xyz_to(writer, source, options)
}
