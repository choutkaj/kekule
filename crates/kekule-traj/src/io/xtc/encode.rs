//! XTC coordinate compression as GROMACS `xdr3dfcoord` writes it.
//!
//! The writer validates every frame before encoding, so each axis spans less
//! than [`XTC_MAX_SCALED`] and no integer arithmetic here can overflow.

use super::{XtcMagic, XTC_FIRST_SMALL_INDEX, XTC_MAGIC_INTS, XTC_MAX_SCALED};

const XTC_LAST_SMALL_INDEX: usize = XTC_MAGIC_INTS.len() - 1;
/// Axis sizes above this are written separately instead of packed together.
const XTC_MAX_PACKED_SIZE: u32 = 0x00ff_ffff;
/// Longest run of small coordinate deltas, in atoms.
const XTC_MAX_RUN: usize = 8;

/// Most-significant-bit-first writer of the compressed bit stream.
struct BitWriter<'a> {
    bytes: &'a mut Vec<u8>,
    pending: u64,
    pending_bits: u32,
}

impl<'a> BitWriter<'a> {
    fn new(bytes: &'a mut Vec<u8>) -> Self {
        Self {
            bytes,
            pending: 0,
            pending_bits: 0,
        }
    }

    /// Appends the low `bits` bits of `value`.
    fn write(&mut self, value: u32, bits: u32) {
        debug_assert!(bits <= 32 && u64::from(value) >> bits == 0);
        self.pending = (self.pending << bits) | u64::from(value);
        self.pending_bits += bits;
        while self.pending_bits >= 8 {
            self.pending_bits -= 8;
            self.bytes.push((self.pending >> self.pending_bits) as u8);
        }
        self.pending &= (1 << self.pending_bits) - 1;
    }

    /// Appends `value` in `bits` bits as GROMACS `sendints` does: whole bytes
    /// least significant first, then the remaining high bits.
    fn write_packed(&mut self, mut value: u128, bits: u32) {
        debug_assert!(bits >= 128 || value >> bits == 0);
        let mut remaining = bits;
        while remaining >= 8 {
            self.write((value & 0xff) as u32, 8);
            value >>= 8;
            remaining -= 8;
        }
        self.write(value as u32, remaining);
    }

    /// Pads the last partial byte with zero bits.
    fn finish(self) {
        if self.pending_bits > 0 {
            self.bytes
                .push((self.pending << (8 - self.pending_bits)) as u8);
        }
    }
}

const fn bit_length(value: u128) -> u32 {
    u128::BITS - value.leading_zeros()
}

/// Packs three values below their `sizes` into one integer, first value most
/// significant.
fn pack(values: [u32; 3], sizes: [u32; 3]) -> u128 {
    (u128::from(values[0]) * u128::from(sizes[1]) + u128::from(values[1])) * u128::from(sizes[2])
        + u128::from(values[2])
}

fn within(a: [i32; 3], b: [i32; 3], limit: i32) -> bool {
    (0..3).all(|axis| (a[axis] - b[axis]).abs() < limit)
}

/// How a whole coordinate is written relative to the frame's minimum.
enum FullCoordinate {
    /// All three axes packed into one integer of this many bits.
    Packed { sizes: [u32; 3], bits: u32 },
    /// Each axis in its own bit count, when an axis is too large to pack.
    Separate { bits: [u32; 3] },
}

impl FullCoordinate {
    fn write(&self, stream: &mut BitWriter<'_>, relative: [u32; 3]) {
        match self {
            Self::Packed { sizes, bits } => stream.write_packed(pack(relative, *sizes), *bits),
            Self::Separate { bits } => {
                for (value, bits) in relative.into_iter().zip(bits) {
                    stream.write(value, *bits);
                }
            }
        }
    }
}

/// Appends the compressed coordinate block of one frame: bounds, initial
/// small-delta index, byte count, bit stream, and XDR padding.
///
/// `coordinates` must hold at least two scaled positions whose per-axis span
/// is below [`XTC_MAX_SCALED`]. Like GROMACS, encoding reorders close pairs in
/// `coordinates` to improve compression of water molecules. `stream` is
/// reusable scratch for the bit stream. Fails only when the compressed byte
/// count does not fit the magic's count field.
pub(super) fn encode_coordinates(
    coordinates: &mut [[i32; 3]],
    magic: XtcMagic,
    stream: &mut Vec<u8>,
    output: &mut Vec<u8>,
) -> Result<(), &'static str> {
    debug_assert!(coordinates.len() >= 2);
    let mut minimum = [i32::MAX; 3];
    let mut maximum = [i32::MIN; 3];
    for coordinate in coordinates.iter() {
        for axis in 0..3 {
            minimum[axis] = minimum[axis].min(coordinate[axis]);
            maximum[axis] = maximum[axis].max(coordinate[axis]);
        }
    }
    let sizes: [u32; 3] = std::array::from_fn(|axis| {
        let span = i64::from(maximum[axis]) - i64::from(minimum[axis]);
        debug_assert!(span < XTC_MAX_SCALED);
        (span + 1) as u32
    });
    let full = if sizes.iter().any(|&size| size > XTC_MAX_PACKED_SIZE) {
        FullCoordinate::Separate {
            bits: sizes.map(|size| bit_length(u128::from(size))),
        }
    } else {
        FullCoordinate::Packed {
            sizes,
            bits: bit_length(sizes.iter().map(|&size| u128::from(size)).product()),
        }
    };
    // The smallest distance between consecutive atoms selects the initial
    // delta size. GROMACS sums in `int`; the sum is exact here.
    let minimum_difference = coordinates
        .windows(2)
        .map(|pair| {
            (0..3)
                .map(|axis| (i64::from(pair[0][axis]) - i64::from(pair[1][axis])).abs())
                .sum::<i64>()
        })
        .min()
        .unwrap_or(0);
    let mut small_index = (XTC_FIRST_SMALL_INDEX..XTC_LAST_SMALL_INDEX)
        .find(|&index| i64::from(XTC_MAGIC_INTS[index]) >= minimum_difference)
        .unwrap_or(XTC_LAST_SMALL_INDEX);

    for value in minimum.iter().chain(&maximum) {
        output.extend_from_slice(&value.to_be_bytes());
    }
    output.extend_from_slice(&(small_index as u32).to_be_bytes());

    stream.clear();
    let mut bits = BitWriter::new(stream);
    let max_index = XTC_LAST_SMALL_INDEX.min(small_index + 8);
    let min_index = max_index - 8;
    let mut smaller = XTC_MAGIC_INTS[XTC_FIRST_SMALL_INDEX.max(small_index - 1)] / 2;
    let mut small = XTC_MAGIC_INTS[small_index] / 2;
    let larger = XTC_MAGIC_INTS[max_index] / 2;
    let mut previous_run = None;
    let mut previous = [0; 3];
    let mut deltas = [[0_u32; 3]; XTC_MAX_RUN];
    let mut index = 0;
    while index < coordinates.len() {
        let mut size_change: i32 = if small_index < max_index
            && index >= 1
            && within(coordinates[index], previous, larger)
        {
            1
        } else if small_index > min_index {
            -1
        } else {
            0
        };
        let mut run_continues = index + 1 < coordinates.len()
            && within(coordinates[index], coordinates[index + 1], small);
        if run_continues {
            coordinates.swap(index, index + 1);
        }
        let current = coordinates[index];
        full.write(
            &mut bits,
            std::array::from_fn(|axis| (current[axis] - minimum[axis]) as u32),
        );
        previous = current;
        index += 1;
        if !run_continues && size_change == -1 {
            size_change = 0;
        }

        let mut run = 0;
        while run_continues && run < XTC_MAX_RUN {
            let current = coordinates[index];
            let squared_distance: i64 = (0..3)
                .map(|axis| {
                    let delta = i64::from(current[axis] - previous[axis]);
                    delta * delta
                })
                .sum();
            if size_change == -1 && squared_distance >= i64::from(smaller) * i64::from(smaller) {
                size_change = 0;
            }
            deltas[run] =
                std::array::from_fn(|axis| (current[axis] - previous[axis] + small) as u32);
            run += 1;
            previous = current;
            index += 1;
            run_continues =
                index < coordinates.len() && within(coordinates[index], previous, small);
        }

        // GROMACS counts a run in integers, three per atom.
        let run_integers = 3 * run as i32;
        if previous_run != Some(run_integers) || size_change != 0 {
            previous_run = Some(run_integers);
            bits.write(1, 1);
            bits.write((run_integers + size_change + 1) as u32, 5);
        } else {
            bits.write(0, 1);
        }
        let small_size = XTC_MAGIC_INTS[small_index] as u32;
        for delta in &deltas[..run] {
            bits.write_packed(pack(*delta, [small_size; 3]), small_index as u32);
        }
        if size_change < 0 {
            small_index -= 1;
            small = smaller;
            smaller = if small_index > XTC_FIRST_SMALL_INDEX {
                XTC_MAGIC_INTS[small_index - 1] / 2
            } else {
                0
            };
        } else if size_change > 0 {
            small_index += 1;
            smaller = small;
            small = XTC_MAGIC_INTS[small_index] / 2;
        }
    }
    bits.finish();

    match magic {
        XtcMagic::Xtc1995 => output.extend_from_slice(
            &i32::try_from(stream.len())
                .map_err(|_| "XTC 1995 compressed byte count exceeds a signed 32-bit integer")?
                .to_be_bytes(),
        ),
        XtcMagic::Xtc2023 => output.extend_from_slice(&(stream.len() as u64).to_be_bytes()),
    }
    output.extend_from_slice(stream);
    output.resize(output.len() + (4 - stream.len() % 4) % 4, 0);
    Ok(())
}
