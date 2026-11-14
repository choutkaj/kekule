use super::*;

#[test]
fn hybrid36_continues_the_decimal_range() {
    for (field, width, value) in [
        ("99999", 5, 99_999),
        ("A0000", 5, 100_000),
        ("A0001", 5, 100_001),
        ("ZZZZZ", 5, 100_000 + 26 * 36i64.pow(4) - 1),
        ("a0000", 5, 100_000 + 26 * 36i64.pow(4)),
        ("zzzzz", 5, 87_440_031),
        ("9999", 4, 9_999),
        ("A000", 4, 10_000),
        ("  -1", 4, -1),
    ] {
        assert_eq!(hybrid36(field, width), Ok(Some(value)), "{field}");
    }
    assert_eq!(hybrid36("    ", 4), Ok(None));
    for invalid in ["A0a0", " A00", "1A00", "!000"] {
        assert_eq!(hybrid36(invalid, 4), Err(()), "{invalid}");
    }
}
