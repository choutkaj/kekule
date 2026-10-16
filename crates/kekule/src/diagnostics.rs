//! Shared diagnostics of format parsing and interpretation.
//!
//! Every format reports failures through the same vocabulary: a structured
//! [`DiagnosticKind`] that callers can match on, and a [`SourceLocation`] in the
//! original input. Messages remain human-readable detail, not an API contract.
//! Writers share [`UnrepresentablePolicy`] for state a format cannot represent.

use std::fmt;

pub use crate::io::UnrepresentablePolicy;

/// A position in the original input text.
///
/// Lines are 1-based and count every line terminator of the original input
/// (`\n`, `\r\n`, or `\r`). Columns, when known, are 1-based byte offsets within
/// the line.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SourceLocation {
    line: usize,
    column: Option<usize>,
}

impl SourceLocation {
    /// A location identified by its line alone.
    pub const fn at_line(line: usize) -> Self {
        Self { line, column: None }
    }

    /// A location identified by line and byte column.
    pub const fn at_column(line: usize, column: usize) -> Self {
        Self {
            line,
            column: Some(column),
        }
    }

    pub const fn line(&self) -> usize {
        self.line
    }

    pub const fn column(&self) -> Option<usize> {
        self.column
    }
}

impl fmt::Display for SourceLocation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.column {
            Some(column) => write!(formatter, "line {}, column {column}", self.line),
            None => write!(formatter, "line {}", self.line),
        }
    }
}

/// What kind of problem a parse or interpretation diagnostic reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum DiagnosticKind {
    /// Reading the input failed, or it is not valid UTF-8.
    Io,
    /// The input exceeds a configured resource limit.
    Limit,
    /// The text does not follow the format grammar.
    Syntax,
    /// A well-formed field holds a value outside its domain, such as an unknown
    /// element symbol or a malformed number.
    InvalidValue,
    /// A field refers to something that does not exist, such as a bond
    /// endpoint outside the atom block.
    InvalidReference,
    /// Valid source content that Kekule does not support.
    Unsupported,
    /// Fields or records contradict each other, or a requested selection does
    /// not match the source.
    Inconsistent,
    /// The represented chemistry cannot be constructed or published.
    Chemistry,
}

impl fmt::Display for DiagnosticKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Io => "input error",
            Self::Limit => "resource limit",
            Self::Syntax => "syntax error",
            Self::InvalidValue => "invalid value",
            Self::InvalidReference => "invalid reference",
            Self::Unsupported => "unsupported content",
            Self::Inconsistent => "inconsistent content",
            Self::Chemistry => "chemistry error",
        })
    }
}

/// How interpretation treats source content it does not support.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum UnsupportedContent {
    /// Fail with a [`DiagnosticKind::Unsupported`] error.
    #[default]
    Reject,
    /// Omit the content and record a typed report entry naming it.
    Report,
}

/// Normalizes line terminators of `input` to `\n`, copying it once.
///
/// Every terminator (`\r\n`, `\r`, `\n`) counts as one line, so line numbers in
/// the normalized text equal line numbers in the original input.
pub(crate) fn normalize_line_endings(input: &str) -> String {
    if !input.contains('\r') {
        return input.to_owned();
    }
    let mut normalized = String::with_capacity(input.len());
    let mut chars = input.chars().peekable();
    while let Some(character) = chars.next() {
        if character == '\r' {
            if chars.peek() == Some(&'\n') {
                chars.next();
            }
            normalized.push('\n');
        } else {
            normalized.push(character);
        }
    }
    normalized
}

/// Reads at most `max_bytes` of UTF-8 text from `reader`.
pub(crate) fn read_bounded(
    reader: impl std::io::Read,
    max_bytes: usize,
) -> Result<String, (DiagnosticKind, String)> {
    use std::io::Read;

    let limit = u64::try_from(max_bytes)
        .unwrap_or(u64::MAX)
        .saturating_add(1);
    let mut bytes = Vec::new();
    reader
        .take(limit)
        .read_to_end(&mut bytes)
        .map_err(|error| (DiagnosticKind::Io, format!("cannot read input: {error}")))?;
    if bytes.len() > max_bytes {
        return Err((
            DiagnosticKind::Limit,
            "input exceeds configured byte limit".to_owned(),
        ));
    }
    String::from_utf8(bytes).map_err(|error| {
        (
            DiagnosticKind::Io,
            format!("input is not valid UTF-8: {error}"),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn line_endings_normalize_without_changing_line_numbers() {
        for (input, expected) in [
            ("a\nb", "a\nb"),
            ("a\r\nb\r\n", "a\nb\n"),
            ("a\rb\r\rc", "a\nb\n\nc"),
            ("\r\n\r", "\n\n"),
        ] {
            assert_eq!(normalize_line_endings(input), expected, "{input:?}");
        }
    }

    #[test]
    fn bounded_reads_report_limits_and_encoding_as_typed_kinds() {
        assert_eq!(read_bounded("abc".as_bytes(), 3).unwrap(), "abc");
        assert_eq!(
            read_bounded("abcd".as_bytes(), 3).unwrap_err().0,
            DiagnosticKind::Limit
        );
        assert_eq!(
            read_bounded([0xff_u8, 0xfe].as_slice(), 10).unwrap_err().0,
            DiagnosticKind::Io
        );
    }

    #[test]
    fn locations_display_line_and_optional_column() {
        assert_eq!(SourceLocation::at_line(4).to_string(), "line 4");
        assert_eq!(
            SourceLocation::at_column(4, 7).to_string(),
            "line 4, column 7"
        );
    }
}
