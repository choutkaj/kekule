use super::*;
use crate::diagnostics::{LimitExceeded, Resource};

#[test]
fn counted_nouns_agree_with_their_count() {
    assert_eq!(counted(1, "atom", "atoms").to_string(), "1 atom");
    assert_eq!(counted(0, "atom", "atoms").to_string(), "0 atoms");
    assert_eq!(
        counted(1_024, "entry", "entries").to_string(),
        "1,024 entries"
    );
    assert_eq!(
        LimitExceeded::first_beyond(Resource::Rows, 0)
            .describe()
            .to_string(),
        "row limit: 1 row of 0"
    );
    assert_eq!(
        LimitExceeded::first_beyond(Resource::MatrixEntries, 0)
            .describe()
            .to_string(),
        "matrix entry limit: 1 entry of 0"
    );
}

#[test]
fn lists_join_with_and() {
    assert_eq!(listed::<&str>(&[]), "");
    assert_eq!(listed(&["a"]), "a");
    assert_eq!(listed(&["a", "b"]), "a and b");
    assert_eq!(listed(&["a", "b", "c"]), "a, b and c");
}

#[test]
fn counts_group_thousands() {
    for (count, text) in [
        (0, "0"),
        (999, "999"),
        (1_000, "1,000"),
        (12_345, "12,345"),
        (123_456, "123,456"),
        (1_234_567, "1,234,567"),
    ] {
        assert_eq!(Count(count).to_string(), text);
    }
}
