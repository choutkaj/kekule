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
mod tests;
