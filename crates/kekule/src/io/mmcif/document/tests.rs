use super::*;

#[test]
fn default_limits_admit_the_default_atom_site_row_limit() {
    let options = MmcifParseOptions::default();
    assert!(
        options.max_tokens >= options.max_atom_site_rows * (TYPICAL_ATOM_SITE_COLUMNS + 1),
        "token budget must cover a full default-sized atom-site table"
    );
    assert!(u32::try_from(options.max_input_bytes).is_ok());
}

#[test]
fn tokenizer_searches_each_line_end_once() {
    let repeats = 1_000;
    let source = format!(
        "data_x\nloop_\n_x.a\n_x.b\n_x.c\n{}# trailing comment\n;text\n; 2 3\n",
        "1 'a b' \"c\" ".repeat(repeats)
    );
    let mut tokenizer = Tokenizer::new(&source, MmcifParseOptions::default());
    let mut tokens = Vec::new();
    while let Some(span) = tokenizer.next_token().unwrap() {
        tokens.push((span.text(&source), span.line, span.bare));
    }
    let mut expected = vec![
        ("data_x", 1, true),
        ("loop_", 2, true),
        ("_x.a", 3, true),
        ("_x.b", 4, true),
        ("_x.c", 5, true),
    ];
    for _ in 0..repeats {
        expected.extend([("1", 6, true), ("a b", 6, false), ("c", 6, false)]);
    }
    expected.extend([("text", 7, false), ("2", 8, true), ("3", 8, true)]);
    assert_eq!(tokens, expected);
    // Eight lines plus the empty remainder after the final newline.
    assert_eq!(tokenizer.line_end_searches, 9);
}
