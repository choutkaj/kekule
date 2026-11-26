//! Rejection of hostile structure set, weights, and clustering records,
//! kinds 9, 10, 13, and 14: every stored log weight, cluster,
//! representative, order value, and section tag is checked, and no
//! truncation or corruption panics.

use super::reports::{find, parts, read, report, resealed};
use crate::clustering::Clustering;
use crate::diagnostics::{DiagnosticKind, Format};
use crate::math::Point3;
use crate::native::{self, NativeSource, NativeWriteOptions, StructureSetDescription};
use crate::properties::{PropertyColumn, PropertyKey};
use crate::statistics::{WeightSource, Weights};
use crate::structure::{
    Geometry, Order, Positions, Sampling, StructureSet, ThermodynamicState, TimeAxis,
};
use crate::units::{Quantity, DEGREE, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND};

/// Steps that mark where the order is stored.
const STEPS: [u64; 3] = [
    0x1111_2222_3333_4444,
    0x1111_2222_3333_4445,
    0x1111_2222_3333_4446,
];
/// Log weights of either sign, as any finite log weight may be.
const LOG_WEIGHTS: [f64; 2] = [-2.091_864_061, 2.290_163_202];
const TEMPERATURE: f64 = 312.345_678;
const X: [f64; 3] = [0.111_111, 0.222_222, 0.333_333];
/// Where a record's content starts: after the magic, the version, the kind,
/// the perception flag, and two reserved bytes.
const CONTENT: usize = 16;

/// Three one-atom geometries at [`X`] nm in `order`.
fn geometries(order: Order) -> StructureSet {
    let topology = crate::smiles::to_topology("O").unwrap();
    let geometries = X.map(|x| {
        Geometry::new(Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap())
    });
    let mut set = StructureSet::from_geometries(topology, geometries).unwrap();
    set.set_order(Some(order)).unwrap();
    set
}

/// Weights of geometries 0 and 2 of three, of [`LOG_WEIGHTS`], at
/// [`TEMPERATURE`].
fn weights() -> Weights {
    Weights::given(
        vec![Some(LOG_WEIGHTS[0]), None, Some(LOG_WEIGHTS[1])],
        Some(Quantity::new(TEMPERATURE, KELVIN)),
    )
    .unwrap()
}

/// Four geometries in two clusters: cluster 0 of geometries 0 and 3,
/// represented by 3, and cluster 1 of geometry 2; geometry 1 is in none.
fn clustering() -> Clustering {
    Clustering::new(vec![Some(0), None, Some(1), Some(0)], vec![3, 2]).unwrap()
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
fn stored_weights_are_finite_and_their_source_known() {
    let weights = weights();
    let bytes = native::write(&weights, NativeWriteOptions::default()).unwrap();
    assert_eq!(bytes[12], 13, "a weights record");
    assert_eq!(read(&bytes).unwrap().into_weights().unwrap(), weights);

    // The geometry count, the source tag, the temperature, and an optional
    // log weight per geometry.
    assert_eq!(bytes[CONTENT..CONTENT + 4], 3u32.to_le_bytes());
    let source = CONTENT + 4;
    assert_eq!(bytes[source], 1, "given weights");
    assert_eq!(bytes[source + 1], 1, "a temperature");
    assert_eq!(find(&bytes, &TEMPERATURE.to_le_bytes()), source + 2);
    let first = source + 11;
    assert_eq!(bytes[first - 1], 1, "geometry 0 is a sample");
    assert_eq!(find(&bytes, &LOG_WEIGHTS[0].to_le_bytes()), first);
    assert_eq!(bytes[first + 8], 0, "geometry 1 is no sample");

    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_invalid(
            &patched(&bytes, first, &bad.to_le_bytes()),
            "The log weight of geometry 0 is not finite.",
            "a log weight",
        );
    }
    assert_invalid(
        &patched(&bytes, source + 2, &(-1.0f64).to_le_bytes()),
        "The stored weights are not valid.",
        "a negative temperature",
    );
    assert_invalid(
        &patched(&bytes, source, &[4]),
        "The weight source tag 4 is unknown.",
        "a source tag",
    );
    assert_invalid(
        &patched(&bytes, first + 8, &[2]),
        "The optional value tag is 2, not 0 or 1.",
        "a sample tag",
    );
    // Zero and negative log weights are weights of one and below one.
    for log_weight in [0.0f64, -1.0] {
        let decoded = read(&patched(&bytes, first, &log_weight.to_le_bytes()))
            .unwrap()
            .into_weights()
            .unwrap();
        assert_eq!(decoded.log_weight(0), Some(log_weight));
    }
    // Log weights cut short.
    for offset in [first + 4, first + 8, first + 10] {
        let error = read(&cut(&bytes, offset)).unwrap_err();
        assert_eq!(error.format(), Format::Native);
        let error = read(&bytes[..offset]).unwrap_err();
        assert_eq!(error.format(), Format::Native);
    }
}

#[test]
fn every_weight_source_reads_back() {
    let mut set = geometries(Order::Step(vec![0, 1, 2]));
    let energy = PropertyKey::new("energy").unwrap();
    set.columns_mut()
        .insert(
            energy.clone(),
            PropertyColumn::Real {
                unit: KILOJOULE_PER_MOLE,
                values: vec![Some(1.5), Some(0.0), None],
            },
        )
        .unwrap();
    let state = ThermodynamicState::isothermal_isobaric(
        Quantity::new(TEMPERATURE, KELVIN),
        Quantity::new(1.0, crate::units::BAR),
    )
    .unwrap();
    set.set_sampling(Some(Sampling::new(state))).unwrap();
    let boltzmann = Weights::boltzmann(&set, &energy, Quantity::new(TEMPERATURE, KELVIN)).unwrap();
    let bias = Weights::from_bias(&set).unwrap().into_weights();
    assert_eq!(bias.source(), &WeightSource::Bias { state });
    let untimed = Weights::given(vec![Some(0.5), Some(-0.5)], None).unwrap();
    for weights in [Weights::equal(3), untimed, boltzmann, bias] {
        let bytes = native::write(&weights, NativeWriteOptions::default()).unwrap();
        assert_eq!(read(&bytes).unwrap().into_weights().unwrap(), weights);
    }
}

#[test]
fn stored_clusterings_are_checked() {
    let clustering = clustering();
    let bytes = native::write(&clustering, NativeWriteOptions::default()).unwrap();
    assert_eq!(bytes[12], 14, "a clustering record");
    assert_eq!(read(&bytes).unwrap().into_clustering().unwrap(), clustering);

    // The geometry count, the cluster count, an optional cluster per
    // geometry, and a representative per cluster.
    assert_eq!(bytes[CONTENT..CONTENT + 4], 4u32.to_le_bytes());
    assert_eq!(bytes[CONTENT + 4..CONTENT + 8], 2u32.to_le_bytes());
    let first = CONTENT + 8;
    assert_eq!(bytes[first], 1, "geometry 0 is in a cluster");
    assert_eq!(bytes[first + 5], 0, "geometry 1 is in none");
    let representatives = first + 16;
    assert_eq!(bytes[representatives..representatives + 8], u32s(&[3, 2]));
    let stored = "The stored clustering is not valid.";
    for (offset, patch, case) in [
        (first + 1, u32s(&[2]), "a cluster out of range"),
        (
            representatives,
            u32s(&[1]),
            "a representative outside its cluster",
        ),
        (representatives, u32s(&[9]), "a representative out of range"),
    ] {
        assert_invalid(&patched(&bytes, offset, &patch), stored, case);
    }
    assert_invalid(
        &patched(&bytes, first + 5, &[2]),
        "The cluster tag is 2, not 0 or 1.",
        "a cluster tag",
    );
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
    assert_eq!(set.len(), 2);
    assert_eq!(set.order(), Some(&Order::Step(STEPS[..2].to_vec())));
    let reports = vec![report(false), report(false)];
    let bytes = native::write(
        NativeSource::interpreted_structure_set(&set, &reports),
        NativeWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(bytes[8..12], 10u32.to_le_bytes(), "format version 10");
    assert_eq!(bytes[12], 10);
    let plain = native::write(&set, NativeWriteOptions::default()).unwrap();
    assert_eq!(plain[12], 9);
    assert_eq!(bytes[13..plain.len() - 8], plain[13..plain.len() - 8]);

    // An empty column table follows the two steps, the absent sampling
    // record's tag and the geometries tag follow it, and the first
    // geometry's positions follow that.
    let steps = find(&bytes, &STEPS[0].to_le_bytes());
    assert_eq!(bytes[steps + 8..steps + 16], STEPS[1].to_le_bytes());
    assert_eq!(bytes[steps + 16..steps + 20], 2u32.to_le_bytes(), "rows");
    assert_eq!(bytes[steps + 20..steps + 24], 0u32.to_le_bytes(), "columns");
    assert_eq!(bytes[steps + 24], 0, "no sampling record");
    let geometries = steps + 25;
    assert_eq!(bytes[geometries], 0);
    let first = set.get(0).unwrap().positions().values().value()[0];
    let stored = [first.x, first.y, first.z]
        .iter()
        .flat_map(|value| value.to_le_bytes())
        .collect::<Vec<_>>();
    assert_eq!(bytes[geometries + 1..geometries + 25], stored);
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
    let set = geometries(Order::Step(STEPS.to_vec()));
    sweep(&native::write(&set, NativeWriteOptions::default()).unwrap());
    sweep(
        &native::write(
            &StructureSetDescription::of(&set),
            NativeWriteOptions::default(),
        )
        .unwrap(),
    );
    sweep(&native::write(&weights(), NativeWriteOptions::default()).unwrap());
    sweep(&native::write(&clustering(), NativeWriteOptions::default()).unwrap());
}

/// [`geometries`] in a step order of [`STEPS`], whose geometries 0 and 1
/// hold the owner value `replica` as `first` and `second`, past the set's
/// checks.
fn with_owner_values(
    first: crate::properties::PropertyValue,
    second: crate::properties::PropertyValue,
) -> StructureSet {
    let mut set = geometries(Order::Step(STEPS.to_vec()));
    let key = crate::properties::PropertyKey::new("replica").unwrap();
    let stored = set.storage_mut().geometries_unchecked();
    for (geometry, value) in stored.iter_mut().zip([first, second]) {
        geometry
            .properties_mut()
            .owner_mut()
            .insert(key.clone(), value)
            .unwrap();
    }
    set
}

#[test]
fn a_set_geometry_with_values_of_its_own_fails_to_read() {
    use crate::properties::PropertyValue;
    let set = with_owner_values(PropertyValue::Int(4), PropertyValue::Int(5));
    let bytes = native::write(&set, NativeWriteOptions::default()).unwrap();
    assert_invalid(
        &bytes,
        "A geometry of the structure set has values of its own. The set's columns hold them.",
        "values of a geometry",
    );
}
