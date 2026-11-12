//! Fingerprints of recorded results: values the retired ensemble,
//! trajectory, and series paths gave, which the sample set paths must keep
//! bit for bit.

use std::fmt::Debug;

/// FNV-1a 64 over the `Debug` text of `value`. `Debug` writes every `f64`
/// as the shortest text that reads back to the same bits.
pub(crate) fn fingerprint(value: &impl Debug) -> u64 {
    format!("{value:?}")
        .bytes()
        .fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
            (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
        })
}

/// Checks a sequence of results against their recorded fingerprints, in
/// order, and that every recorded fingerprint was checked.
pub(crate) struct Recorded {
    expected: &'static [u64],
    next: usize,
}

impl Recorded {
    pub(crate) const fn new(expected: &'static [u64]) -> Self {
        Self { expected, next: 0 }
    }

    pub(crate) fn check(&mut self, value: &impl Debug) {
        assert_eq!(
            Some(&fingerprint(value)),
            self.expected.get(self.next),
            "recorded result {}",
            self.next
        );
        self.next += 1;
    }
}

impl Drop for Recorded {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            assert_eq!(self.next, self.expected.len(), "unchecked recorded results");
        }
    }
}
