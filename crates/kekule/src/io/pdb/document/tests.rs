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

#[test]
fn vmd_hexadecimal_numbers_need_the_sequence_to_prove_them() {
    let mut serials = Numbering::default();
    let mut read = |field: &str| serials.decode(field, 5);
    assert_eq!(read("99999"), Ok((Some(99_999), false)));
    // Hexadecimal 186a0 is 100,000, one past the decimal range.
    assert_eq!(read("186a0"), Ok((Some(100_000), true)));
    // All digits, yet hexadecimal while the sequence continues.
    assert_eq!(read("18700"), Ok((Some(100_096), false)));
    // A number below the last ends hexadecimal numbering.
    assert_eq!(read("00001"), Ok((Some(1), false)));
    assert_eq!(read("00002"), Ok((Some(2), false)));

    // Without a decimal predecessor, lowercase letters are hybrid-36.
    let mut serials = Numbering::default();
    assert_eq!(
        serials.decode("a0000", 5),
        Ok((Some(100_000 + 26 * 36i64.pow(4)), false))
    );
    assert_eq!(vmd_hexadecimal("0ffff", 5), None);
}
