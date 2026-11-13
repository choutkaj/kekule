//! Rejection of hostile structure set records, kinds 9 and 10: every stored
//! member, assignment, representative, weight, order value, and section tag
//! is checked against the set, and no truncation or corruption panics.

use super::native_reports::{find, parts, read, report, resealed};
use crate::diagnostics::{DiagnosticKind, Format};
use crate::math::Point3;
use crate::native::{self, NativeSource, NativeWriteOptions, StructureSetDescription};
use crate::structure::{
    Distribution, Geometry, Members, Order, Partition, Positions, State, StructureSet, TimeAxis,
};
use crate::units::{Quantity, DEGREE, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND};

/// Steps that mark where the order is stored.
const STEPS: [u64; 3] = [
    0x1111_2222_3333_4444,
    0x1111_2222_3333_4445,
    0x1111_2222_3333_4446,
];
const WEIGHTS: [f64; 2] = [0.123_456_789, 9.876_543_21];
const ENERGY: f64 = -12.345_678_9;
const TEMPERATURE: f64 = 312.345_678;
const X: [f64; 3] = [0.111_111, 0.222_222, 0.333_333];

/// Three one-atom geometries at [`X`] nm in `order`, without a distribution.
fn geometries(order: Order) -> StructureSet {
    let topology = crate::smiles::to_topology("O").unwrap();
    let geometries = X.map(|x| {
        Geometry::new(Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap())
    });
    let mut set = StructureSet::from_geometries(topology, geometries).unwrap();
    set.set_order(Some(order)).unwrap();
    set
}

/// [`geometries`] in a step order of [`STEPS`], with a distribution over
/// geometries 0 and 2 weighted [`WEIGHTS`], in two states: state 0 of geometry 2
/// at [`ENERGY`], state 1 of geometry 0, at [`TEMPERATURE`].
fn set() -> StructureSet {
    let mut set = geometries(Order::Step(STEPS.to_vec()));
    let distribution = Distribution::new(
        Members::Subset(vec![0, 2]),
        Partition::new(vec![1, 0], vec![1, 0]).unwrap(),
    )
    .unwrap()
    .with_weights(WEIGHTS.to_vec())
    .unwrap()
    .with_states(vec![
        State::new()
            .with_energy(Quantity::new(ENERGY, KILOJOULE_PER_MOLE))
            .unwrap(),
        State::new(),
    ])
    .unwrap()
    .with_temperature(Quantity::new(TEMPERATURE, KELVIN))
    .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    set
}

/// Where [`set`]'s encoding stores each part of its distribution.
struct Layout {
    listed: usize,
    members: usize,
    weights: usize,
    assignments: usize,
    representative: usize,
    geometries: usize,
}

impl Layout {
    fn of(bytes: &[u8]) -> Self {
        // The steps follow the order kind; the distribution follows them.
        let steps = find(bytes, &STEPS[0].to_le_bytes());
        assert_eq!(bytes[steps - 1], 2, "a step order");
        let start = steps + 24;
        assert_eq!(bytes[start..start + 2], [1, 1], "a subset distribution");
        let layout = Self {
            listed: start + 2,
            members: start + 6,
            weights: start + 15,
            assignments: start + 31,
            representative: start + 43,
            geometries: find(bytes, &TEMPERATURE.to_le_bytes()) + 8,
        };
        assert_eq!(find(bytes, &WEIGHTS[0].to_le_bytes()), layout.weights);
        assert_eq!(
            find(bytes, &ENERGY.to_le_bytes()),
            layout.representative + 5
        );
        assert_eq!(bytes[layout.geometries], 0, "geometries in the record");
        layout
    }
}

/// `bytes` with `patch` at `offset`, resealed.
fn patched(bytes: &[u8], offset: usize, patch: &[u8]) -> Vec<u8> {
    let mut corrupted = bytes.to_vec();
    corrupted[offset..offset + patch.len()].copy_from_slice(patch);
    resealed(corrupted)
}

/// `bytes` cut at `offset`, with the trailer, resealed.
fn cut(bytes: &[u8], offset: usize) -> Vec<u8> {
    let mut cut = bytes[..offset].to_vec();
    cut.extend_from_slice(&bytes[bytes.len() - 8..]);
    resealed(cut)
}

/// Asserts that `bytes` fail to read as an invalid value whose message
/// contains `message`.
fn assert_invalid(bytes: &[u8], message: &str, case: &str) {
    let error = read(bytes).unwrap_err();
    assert_eq!(error.format(), Format::Native, "{case}: {error}");
    assert_eq!(
        error.kind(),
        DiagnosticKind::InvalidValue,
        "{case}: {error}"
    );
    assert!(error.to_string().contains(message), "{case}: {error}");
}

fn u32s(values: &[u32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect()
}

#[test]
fn stored_distributions_are_checked_against_their_set() {
    let set = set();
    let bytes = native::write(&set, NativeWriteOptions::default()).unwrap();
    let decoded = read(&bytes).unwrap().into_structure_set().unwrap();
    assert_eq!(decoded.order(), set.order());
    assert_eq!(decoded.distribution(), set.distribution());
    let at = Layout::of(&bytes);
    let stored = "The stored structure is not valid.";
    let cases = [
        // Members unsorted, duplicated, and out of range.
        (at.members, u32s(&[2, 0]), stored, "unsorted members"),
        (at.members, u32s(&[0, 0]), stored, "duplicate members"),
        (at.members, u32s(&[0, 7]), stored, "a member out of range"),
        // More members than geometries.
        (
            at.listed,
            u32s(&[4]),
            "The distribution lists 4 members of 3 geometries.",
            "too many members",
        ),
        // An assignment to a state that does not exist.
        (at.assignments, u32s(&[2]), stored, "an assignment"),
        // State 0 represented by member 0, which is in state 1.
        (
            at.representative,
            u32s(&[0]),
            stored,
            "a representative outside its state",
        ),
        (
            at.representative,
            u32s(&[9]),
            stored,
            "a representative out of range",
        ),
        // Weights that are not finite and positive.
        (at.weights, f64::NAN.to_le_bytes().to_vec(), stored, "NaN"),
        (
            at.weights,
            f64::INFINITY.to_le_bytes().to_vec(),
            stored,
            "infinity",
        ),
        (at.weights, 0.0f64.to_le_bytes().to_vec(), stored, "zero"),
        (
            at.weights,
            (-1.0f64).to_le_bytes().to_vec(),
            stored,
            "negative",
        ),
        // An unknown geometries tag.
        (
            at.geometries,
            vec![7],
            "The geometries section tag 7 is unknown.",
            "a geometries tag",
        ),
    ];
    for (offset, patch, message, case) in cases {
        assert_invalid(&patched(&bytes, offset, &patch), message, case);
    }

    // Weights and assignments cut short.
    for (offset, case) in [
        (at.weights + 4, "weights"),
        (at.weights + 8, "weights"),
        (at.assignments + 2, "assignments"),
        (at.assignments + 4, "assignments"),
    ] {
        let error = read(&cut(&bytes, offset)).unwrap_err();
        assert_eq!(error.format(), Format::Native, "{case}: {error}");
        let error = read(&bytes[..offset]).unwrap_err();
        assert_eq!(error.format(), Format::Native, "{case}: {error}");
    }
}

#[test]
fn stored_order_values_and_positions_are_finite() {
    let times = [0.5, 1.5, 2.5];
    let timed = geometries(Order::Time(
        TimeAxis::new(Some(Quantity::new(times.to_vec(), PICOSECOND)), None).unwrap(),
    ));
    let scan = geometries(Order::Coordinate {
        unit: DEGREE,
        values: times.to_vec(),
    });
    for (set, case) in [(timed, "time"), (scan, "coordinate")] {
        let bytes = native::write(&set, NativeWriteOptions::default()).unwrap();
        let value = find(&bytes, &times[1].to_le_bytes());
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_invalid(
                &patched(&bytes, value, &bad.to_le_bytes()),
                "The stored structure is not valid.",
                case,
            );
        }
        let position = find(&bytes, &X[1].to_le_bytes());
        assert_invalid(
            &patched(&bytes, position, &f64::NAN.to_le_bytes()),
            "positions are not valid",
            case,
        );
    }
}

#[test]
fn interpreted_structure_sets_store_their_geometries_and_one_report_each() {
    let (_, mut set) = parts();
    set.set_order(Some(Order::Step(STEPS[..2].to_vec())))
        .unwrap();
    assert_eq!(super::debug_fingerprint(&set), 0xb5c1_674f_fb04_d6c5);
    let reports = vec![report(false), report(false)];
    let bytes = native::write(
        NativeSource::interpreted_structure_set(&set, &reports),
        NativeWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(super::debug_fingerprint(&bytes), 0x7e76_8e05_f793_9faf);
    assert_eq!(bytes[12], 10);
    let plain = native::write(&set, NativeWriteOptions::default()).unwrap();
    assert_eq!(bytes[13..plain.len() - 8], plain[13..plain.len() - 8]);

    // The geometries tag follows a distribution of every geometry, weighted, in
    // one state without an energy, and without a temperature.
    let steps = find(&bytes, &STEPS[0].to_le_bytes());
    let start = steps + 16;
    assert_eq!(bytes[start..start + 3], [1, 0, 1]);
    assert_eq!(bytes[start + 3..start + 11], 1.0f64.to_le_bytes());
    let geometries = start + 41;
    assert_eq!(bytes[geometries - 1..=geometries], [0, 0]);
    assert_invalid(
        &patched(&bytes, geometries, &[1]),
        "An interpreted structure set stores its geometries elsewhere.",
        "external geometries",
    );

    // A report fewer or more than the geometries.
    let reported = &bytes[plain.len() - 8..bytes.len() - 8];
    let one = &reported[..reported.len() / 2];
    assert_eq!(one, &reported[one.len()..]);
    let short = cut(&bytes, bytes.len() - 8 - one.len());
    assert_eq!(read(&short).unwrap_err().format(), Format::Native);
    let mut extra = bytes[..bytes.len() - 8].to_vec();
    extra.extend_from_slice(one);
    extra.extend_from_slice(&bytes[bytes.len() - 8..]);
    assert_eq!(read(&resealed(extra)).unwrap_err().format(), Format::Native);
}

/// Every prefix of `bytes` fails, and no single-byte corruption panics,
/// with or without a valid checksum.
fn sweep(bytes: &[u8]) {
    for len in 0..bytes.len() {
        let error = read(&bytes[..len]).unwrap_err();
        assert_eq!(error.format(), Format::Native);
    }
    for index in 0..bytes.len() {
        for flip in [0x01, 0x80, 0xFF] {
            let mut corrupted = bytes.to_vec();
            corrupted[index] ^= flip;
            assert!(
                read(&corrupted).is_err(),
                "corrupting byte {index} with {flip:#04x} went unnoticed"
            );
            let _ = read(&resealed(corrupted));
        }
    }
}

#[test]
fn every_truncation_and_corruption_of_a_structure_set_fails_cleanly() {
    let set = set();
    sweep(&native::write(&set, NativeWriteOptions::default()).unwrap());
    sweep(
        &native::write(
            &StructureSetDescription::of(&set),
            NativeWriteOptions::default(),
        )
        .unwrap(),
    );
}
