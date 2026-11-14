use super::Crc32;

#[test]
fn matches_the_standard_check_value() {
    // The CRC-32/ISO-HDLC check value of the nine ASCII digits.
    let mut crc = Crc32::new();
    crc.update(b"123456789");
    assert_eq!(crc.finish(), 0xCBF4_3926);
    let mut split = Crc32::new();
    split.update(b"1234");
    split.update(b"56789");
    assert_eq!(split.finish(), 0xCBF4_3926);
    assert_eq!(Crc32::new().finish(), 0);
}
