//! Focused regressions for the CIF source-syntax boundary.
use crate::mmcif::{self, MmcifParseOptions};

#[test]
fn mmcif_semicolon_text_preserves_opening_content_and_whitespace() {
    // CIF 1.1 syntax paragraphs 17–20: the closing delimiter's preceding
    // newline is excluded; other text (including initial blank lines) remains.
    for (field, expected) in [
        (";first\n;", "first"),
        (";first\nsecond\n;", "first\nsecond"),
        (";\n;", ""),
        (";\n\n  body  \n\n;", "\n\n  body  \n"),
        ("; foo\n  bar\n;", " foo\n  bar"),
        (";?\n;", "?"),
        (";.\n;", "."),
        (
            ";# literal\n ; not a delimiter\n;",
            "# literal\n ; not a delimiter",
        ),
    ] {
        let input = format!("data_x\n_x.text\n{field}\n_x.after yes\n");
        for newline in ["\n", "\r\n", "\r"] {
            let document = mmcif::parse_str(&input.replace('\n', newline)).unwrap();
            let value = document.blocks()[0].item("_x.text").unwrap();
            assert_eq!(value.text(), expected);
            assert_eq!(value.line(), 3);
            assert!(
                !value.is_missing(),
                "delimited placeholders are literal text"
            );
            assert_eq!(document.blocks()[0].item("_x.after").unwrap().text(), "yes");
        }
    }
}

#[test]
fn mmcif_semicolon_closing_line_preserves_following_tokens() {
    let doc = mmcif::parse_str("data_x\n_x.text\n;value\n; _x.after yes\n").unwrap();
    assert_eq!(doc.blocks()[0].item("_x.text").unwrap().text(), "value");
    assert_eq!(doc.blocks()[0].item("_x.after").unwrap().line(), 4);
    let error = mmcif::parse_str("data_x\n_x.text\n;value\n;not_delimited\n").unwrap_err();
    assert_eq!(error.line(), 4);
    assert!(error.message().contains("closing semicolon"));
    let error = mmcif::parse_str("data_x\n_x.text\n;value\nno end").unwrap_err();
    assert_eq!(error.line(), 3);
    assert!(error.message().contains("unterminated semicolon"));
}

#[test]
fn mmcif_semicolon_size_limit_counts_opening_line_and_newlines() {
    let options = MmcifParseOptions {
        max_token_bytes: 8,
        ..Default::default()
    };
    for field in [";12345678\n;", ";1234\n678\n;"] {
        let doc =
            mmcif::parse_str_with_options(&format!("data_x\n_x.v\n{field}\n"), options).unwrap();
        assert_eq!(doc.blocks()[0].item("_x.v").unwrap().text().len(), 8);
    }
    for field in [";123456789\n;", ";1234\n6789\n;", ";12345678\n\n;"] {
        let error = mmcif::parse_str_with_options(&format!("data_x\n_x.v\n{field}\n"), options)
            .unwrap_err();
        assert_eq!(error.line(), 3);
        assert!(error.message().contains("token limit"));
    }
}
