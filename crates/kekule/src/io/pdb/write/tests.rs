use super::*;

#[test]
fn hybrid36_encodes_what_the_reader_decodes() {
    for (value, width, text) in [
        (99_999, 5, "99999"),
        (100_000, 5, "A0000"),
        (100_001, 5, "A0001"),
        (87_440_031, 5, "zzzzz"),
        (9_999, 4, "9999"),
        (10_000, 4, "A000"),
        (-1, 4, "-1"),
        (-999, 4, "-999"),
        (2_436_111, 4, "zzzz"),
    ] {
        assert_eq!(hybrid36(value, width).as_deref(), Some(text), "{value}");
    }
    assert_eq!(hybrid36(87_440_032, 5), None);
    assert_eq!(hybrid36(-1000, 4), None);
}

#[test]
fn names_and_charges_follow_the_column_conventions() {
    assert_eq!(aligned_name("CA", "C"), " CA ");
    assert_eq!(aligned_name("FE", "FE"), "FE  ");
    assert_eq!(aligned_name("HG21", "H"), "HG21");
    assert_eq!(charge_field(-2), "2-");
    assert_eq!(charge_field(1), "1+");
    assert_eq!(charge_field(0), "0");
}
