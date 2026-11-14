use super::{DiagnosticKind, Format, FormatError, LimitExceeded, Resource};

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

/// Reads at most `max_bytes` of UTF-8 text from `reader`, failing with a
/// parse-stage [`FormatError`] of `format` without a location.
pub(crate) fn read_bounded(
    reader: impl std::io::Read,
    max_bytes: usize,
    format: Format,
) -> Result<String, FormatError> {
    use std::io::Read;

    let limit = u64::try_from(max_bytes)
        .unwrap_or(u64::MAX)
        .saturating_add(1);
    let mut bytes = Vec::new();
    reader
        .take(limit)
        .read_to_end(&mut bytes)
        .map_err(|error| FormatError::parse(format, DiagnosticKind::Io, "").with_source(error))?;
    if bytes.len() > max_bytes {
        let limit = LimitExceeded::new(Resource::InputBytes, max_bytes, bytes.len());
        return Err(FormatError::parse(
            format,
            DiagnosticKind::Limit,
            format!("The input exceeds its {}.", limit.describe()),
        )
        .exceeding(limit));
    }
    String::from_utf8(bytes).map_err(|error| {
        FormatError::parse(format, DiagnosticKind::Io, "The input is not valid UTF-8.")
            .with_source(error)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagnostics::FormatStage;

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
        assert_eq!(
            read_bounded("abc".as_bytes(), 3, Format::Smiles).unwrap(),
            "abc"
        );
        let limit = read_bounded("abcd".as_bytes(), 3, Format::Smiles).unwrap_err();
        assert_eq!(
            (limit.format(), limit.stage(), limit.kind()),
            (Format::Smiles, FormatStage::Parse, DiagnosticKind::Limit)
        );
        let encoding = read_bounded([0xff_u8, 0xfe].as_slice(), 10, Format::Pdb).unwrap_err();
        assert_eq!(encoding.kind(), DiagnosticKind::Io);
        assert!(std::error::Error::source(&encoding)
            .and_then(|source| source.downcast_ref::<std::string::FromUtf8Error>())
            .is_some());
    }
}
