use super::{Crc32, POLYNOMIAL};

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

/// The CRC of `bytes` one bit at a time.
fn bitwise(bytes: &[u8]) -> u32 {
    let mut crc = !0u32;
    for &byte in bytes {
        crc ^= u32::from(byte);
        for _ in 0..8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ POLYNOMIAL
            } else {
                crc >> 1
            };
        }
    }
    !crc
}

#[test]
fn matches_the_bitwise_crc_at_every_length_and_split() {
    let bytes = (0..300u32)
        .map(|index| (index.wrapping_mul(2_654_435_761) >> 13) as u8)
        .collect::<Vec<_>>();
    for len in 0..bytes.len() {
        let mut crc = Crc32::new();
        crc.update(&bytes[..len]);
        assert_eq!(crc.finish(), bitwise(&bytes[..len]), "length {len}");
    }
    for split in 0..bytes.len() {
        let mut crc = Crc32::new();
        crc.update(&bytes[..split]);
        crc.update(&bytes[split..]);
        assert_eq!(crc.finish(), bitwise(&bytes), "split at {split}");
    }
}
