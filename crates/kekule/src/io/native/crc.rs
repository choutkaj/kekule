//! CRC-32 (IEEE 802.3, reflected polynomial `0xEDB88320`), as used by zlib
//! and PNG, computed eight bytes at a time.

const POLYNOMIAL: u32 = 0xEDB8_8320;

const TABLES: [[u32; 256]; 8] = tables();

const fn tables() -> [[u32; 256]; 8] {
    let mut tables = [[0u32; 256]; 8];
    let mut byte = 0;
    while byte < 256 {
        let mut crc = byte as u32;
        let mut bit = 0;
        while bit < 8 {
            crc = if crc & 1 == 1 {
                (crc >> 1) ^ POLYNOMIAL
            } else {
                crc >> 1
            };
            bit += 1;
        }
        tables[0][byte] = crc;
        byte += 1;
    }
    let mut table = 1;
    while table < 8 {
        let mut byte = 0;
        while byte < 256 {
            let previous = tables[table - 1][byte];
            tables[table][byte] = (previous >> 8) ^ tables[0][(previous & 0xFF) as usize];
            byte += 1;
        }
        table += 1;
    }
    tables
}

/// A running CRC-32.
#[derive(Debug, Clone, Copy)]
pub(super) struct Crc32(u32);

impl Crc32 {
    pub(super) const fn new() -> Self {
        Self(!0)
    }

    pub(super) fn update(&mut self, mut bytes: &[u8]) {
        let mut crc = self.0;
        while let Some((chunk, rest)) = bytes.split_first_chunk::<8>() {
            let low = u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]) ^ crc;
            crc = TABLES[7][(low & 0xFF) as usize]
                ^ TABLES[6][((low >> 8) & 0xFF) as usize]
                ^ TABLES[5][((low >> 16) & 0xFF) as usize]
                ^ TABLES[4][(low >> 24) as usize]
                ^ TABLES[3][chunk[4] as usize]
                ^ TABLES[2][chunk[5] as usize]
                ^ TABLES[1][chunk[6] as usize]
                ^ TABLES[0][chunk[7] as usize];
            bytes = rest;
        }
        for &byte in bytes {
            crc = (crc >> 8) ^ TABLES[0][((crc ^ u32::from(byte)) & 0xFF) as usize];
        }
        self.0 = crc;
    }

    pub(super) const fn finish(self) -> u32 {
        !self.0
    }
}

#[cfg(test)]
mod tests {
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
}
