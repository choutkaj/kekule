//! Distinct rows of string columns: the decoder keeps each distinct row of
//! atom-site metadata or residue labels once, and each node its row's key.

use std::collections::HashMap;
use std::rc::Rc;

use super::super::encode::ABSENT_STRING;

/// A string column read whole: its text and each entry's length, absent
/// entries as [`ABSENT_STRING`]. Every entry lies on character boundaries
/// of the text.
pub(super) struct StringColumn {
    pub(super) text: String,
    pub(super) lengths: Vec<u32>,
}

impl StringColumn {
    /// Each entry in order.
    pub(super) fn entries(&self) -> impl Iterator<Item = Option<&str>> {
        let mut start = 0;
        self.lengths.iter().map(move |&length| {
            if length == ABSENT_STRING {
                return None;
            }
            let end = start + length as usize;
            // The decoder checked every entry's bounds when it read them.
            let entry = &self.text[start..end];
            start = end;
            Some(entry)
        })
    }
}

/// A value stored beside a row's strings and compared with them.
pub(super) trait Extra: Copy {
    /// Appends this value's bytes, of a length fixed by its own leading
    /// bytes.
    fn pack(self, into: &mut Vec<u8>);
}

impl Extra for () {
    fn pack(self, _: &mut Vec<u8>) {}
}

impl Extra for Option<i32> {
    fn pack(self, into: &mut Vec<u8>) {
        match self {
            None => into.push(0),
            Some(value) => {
                into.push(1);
                into.extend_from_slice(&value.to_le_bytes());
            }
        }
    }
}

/// The rows of `columns`, of `count` entries each, with `extras` beside
/// them: each row's key, and the distinct rows in order of first
/// appearance, each made once by `make`.
pub(super) fn distinct_rows<'a, const N: usize, E: Extra, T>(
    columns: &'a [StringColumn; N],
    extras: impl IntoIterator<Item = E>,
    count: usize,
    mut make: impl FnMut([Option<&'a str>; N], E) -> T,
) -> (Vec<u32>, Vec<T>) {
    let mut entries = columns.each_ref().map(StringColumn::entries);
    let mut extras = extras.into_iter();
    let mut keys = Vec::with_capacity(count);
    let mut rows = Vec::new();
    // Rows are compared as one run of bytes: each entry's length code and
    // bytes, then the extra value. The codes make the runs of distinct rows
    // differ. The default hasher is keyed per process, so input cannot
    // choose colliding rows.
    let mut known = HashMap::<Rc<[u8]>, u32>::new();
    let mut runs = Vec::<Rc<[u8]>>::new();
    // Most nodes repeat a row met shortly before, as the sites of one atom
    // name in one chain do. A small table of recent rows, found by a cheap
    // fingerprint and confirmed by comparing runs, answers them before the
    // map does; any other row goes to the map.
    let bits = count
        .next_power_of_two()
        .trailing_zeros()
        .clamp(MIN_RECENT_BITS, MAX_RECENT_BITS);
    let mut recent = vec![u32::MAX; 1 << bits];
    let mut packed = Vec::new();
    for _ in 0..count {
        let texts = entries.each_mut().map(|column| column.next().flatten());
        let Some(extra) = extras.next() else { break };
        packed.clear();
        for text in texts {
            let code = text.map_or(ABSENT_STRING, |text| text.len() as u32);
            packed.extend_from_slice(&code.to_le_bytes());
            packed.extend_from_slice(text.unwrap_or_default().as_bytes());
        }
        extra.pack(&mut packed);
        let slot = &mut recent[fingerprint(&packed, bits)];
        let key = match runs.get(*slot as usize) {
            Some(run) if **run == *packed => *slot,
            _ => {
                let key = match known.get(packed.as_slice()) {
                    Some(&key) => key,
                    None => {
                        // A key per row: at most `count` rows, counted in
                        // `u32`.
                        let key = rows.len() as u32;
                        rows.push(make(texts, extra));
                        let run = Rc::<[u8]>::from(packed.as_slice());
                        runs.push(Rc::clone(&run));
                        known.insert(run, key);
                        key
                    }
                };
                *slot = key;
                key
            }
        };
        keys.push(key);
    }
    (keys, rows)
}

/// The table of recent rows has a slot per node, a power of two between
/// `2^MIN_RECENT_BITS` and `2^MAX_RECENT_BITS`.
const MIN_RECENT_BITS: u32 = 4;
const MAX_RECENT_BITS: u32 = 12;

/// The slot of the run `bytes` in a table of recent rows of `2^bits` slots:
/// a cheap mix that input may steer, which costs only table misses.
fn fingerprint(bytes: &[u8], bits: u32) -> usize {
    const MIX: u64 = 0x9E37_79B9_7F4A_7C15;
    let (words, rest) = bytes.as_chunks::<8>();
    let mut hash = bytes.len() as u64;
    for word in words {
        hash = (hash ^ u64::from_le_bytes(*word)).wrapping_mul(MIX);
    }
    for &byte in rest {
        hash = (hash ^ u64::from(byte)).wrapping_mul(MIX);
    }
    (hash >> (u64::BITS - bits)) as usize
}

#[cfg(test)]
mod tests;
