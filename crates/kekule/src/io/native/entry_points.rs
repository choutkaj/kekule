//! The public native-format entry points that `crate::native` re-exports.

use crate::diagnostics::FormatError;
use crate::io::{NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions};

/// Encodes one object into a new byte vector.
///
/// # Errors
///
/// As [`write_to`], except that a vector never fails as a sink.
pub fn write<'a>(
    source: impl Into<NativeSource<'a>>,
    options: NativeWriteOptions,
) -> Result<Vec<u8>, FormatError> {
    crate::io::write_native(source, options)
}

/// Encodes one object into `writer`, staging output through a fixed
/// 64 KiB buffer rather than holding the whole encoding.
///
/// # Errors
///
/// A [`FormatError`] of kind
/// [`DiagnosticKind::Io`](crate::diagnostics::DiagnosticKind::Io) when
/// `writer` fails, or
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// when a list or string is too long for its 32-bit length field.
pub fn write_to<'a>(
    writer: &mut impl std::io::Write,
    source: impl Into<NativeSource<'a>>,
    options: NativeWriteOptions,
) -> Result<(), FormatError> {
    crate::io::write_native_to(writer, source, options)
}

/// Decodes one object from `reader`, reading exactly its bytes, so
/// further objects may follow in the same stream.
///
/// # Errors
///
/// A [`FormatError`] of format
/// [`Format::Native`](crate::diagnostics::Format::Native) with a
/// [`byte_offset`](FormatError::byte_offset):
/// [`DiagnosticKind::Syntax`](crate::diagnostics::DiagnosticKind::Syntax)
/// for truncated input or wrong magic bytes,
/// [`DiagnosticKind::Unsupported`](crate::diagnostics::DiagnosticKind::Unsupported)
/// with [`FormatErrorDetail::NativeVersion`](crate::diagnostics::FormatErrorDetail::NativeVersion)
/// for another format version,
/// [`DiagnosticKind::Limit`](crate::diagnostics::DiagnosticKind::Limit)
/// when `options` is exceeded, and `InvalidValue`, `InvalidReference`,
/// `Inconsistent`, or `Chemistry` when stored values do not reconstruct
/// a valid owner or the checksum does not match.
pub fn read(
    mut reader: impl std::io::Read,
    options: NativeReadOptions,
) -> Result<NativeObject, FormatError> {
    crate::io::read_native(&mut reader, options)
}
