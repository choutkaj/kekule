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
