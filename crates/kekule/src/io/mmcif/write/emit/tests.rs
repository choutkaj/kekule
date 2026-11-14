use super::cif_value;

#[test]
fn cif_values_quote_exactly_the_text_cif_cannot_leave_bare() {
    for (value, expected) in [
        ("CA", "CA"),
        ("1-2", "1-2"),
        ("a#b", "'a#b'"),
        ("two words", "'two words'"),
        ("C1'", "\"C1'\""),
        ("say \"hi\"", "'say \"hi\"'"),
        (".", "'.'"),
        ("?", "'?'"),
        ("_tag", "'_tag'"),
        ("#comment", "'#comment'"),
        (";text", "';text'"),
        ("$frame", "'$frame'"),
        ("[1]", "'[1]'"),
        ("]", "']'"),
        ("loop_", "'loop_'"),
        ("DATA_block", "'DATA_block'"),
        ("save_frame", "'save_frame'"),
    ] {
        assert_eq!(cif_value(value, "_test").unwrap(), expected, "{value}");
    }
}
