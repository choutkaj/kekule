//! XYZ writing of sample sets in a step or coordinate order, one frame per
//! sample.

use std::sync::Arc;

use kekule::diagnostics::{StructureField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::geometry::{PeriodicCell, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Order, Place, Sample, SampleSet};
use kekule::traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule::traj::{TrajectoryField, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryWriter};
use kekule::units::{Quantity, ANGSTROM, KILOJOULE_PER_MOLE, NANOMETER};

use crate::support::{self, Frame, Recorded};

const WATERS: &str = "3
first
O 0.0 0.0 0.0
H 0.957 0.0 0.0
H -0.24 0.927 0.0
3
second
O 0.0 0.0 0.1
H 0.957 0.0 0.1
H -0.24 0.927 0.1
3
third
O 0.0 0.0 0.2
H 0.957 0.0 0.2
H -0.24 0.927 0.2
";

/// The positions of the frames of [`WATERS`] in a set at `order`, each
/// sample at the place `place` gives its index.
fn ordered(order: Order, place: impl Fn(usize) -> Place) -> SampleSet {
    let read = kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_sample_set(Default::default())
        .unwrap()
        .into_sample_set();
    let mut set = SampleSet::new(read.shared_topology());
    set.set_order(Some(order)).unwrap();
    for (index, sample) in read.iter().enumerate() {
        let sample = Sample::new(sample.positions().clone());
        set.push(sample, Some(place(index))).unwrap();
    }
    set
}

/// The frames of [`WATERS`] at steps 0, 5, and 10, as an optimization
/// history.
fn history() -> SampleSet {
    ordered(Order::Step(Vec::new()), |index| {
        Place::Step(5 * index as u64)
    })
}

/// The frames of [`WATERS`] as a scan of a distance in angstroms.
fn scan() -> SampleSet {
    let order = Order::Coordinate {
        unit: ANGSTROM,
        values: Vec::new(),
    };
    ordered(order, |index| {
        Place::Coordinate(Quantity::new(0.9 + 0.1 * index as f64, ANGSTROM))
    })
}

/// [`history`] with every sample in [`cell`].
fn boxed() -> SampleSet {
    in_cells(history())
}

/// `set` with every sample in [`cell`].
fn in_cells(mut set: SampleSet) -> SampleSet {
    for index in 0..set.len() {
        set.get_mut(index).unwrap().set_cell(Some(cell()));
    }
    set
}

/// What any XYZ file of [`WATERS`] leaves out: the inferred bonds.
const BONDS: WriteOmission = WriteOmission::Bonds { bonds: 2 };

const TEXT: &str = "3\nitem\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n\
                    3\nitem\nO 0.000 0.000 0.100\nH 0.957 0.000 0.100\nH -0.240 0.927 0.100\n\
                    3\nitem\nO 0.000 0.000 0.200\nH 0.957 0.000 0.200\nH -0.240 0.927 0.200\n";

fn options() -> XyzWriteOptions {
    XyzWriteOptions::default()
        .with_decimal_places(3)
        .with_comment("item")
}

fn omit_places() -> XyzWriteOptions {
    options().with_sample_metadata(UnrepresentablePolicy::Omit)
}

/// Writes `set` and returns the text and report.
fn write(
    set: &SampleSet,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), TrajectoryIoError> {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "s.xyz")?;
    writer.write_sample_set(set)?;
    let report = writer.report().clone();
    Ok((String::from_utf8(writer.finish()?).unwrap(), report))
}

/// Asserts that writing `set` with `options` fails before any byte for
/// `field` at `frame`, and returns the error.
fn assert_rejects(
    set: &SampleSet,
    options: XyzWriteOptions,
    field: TrajectoryField,
    frame: usize,
) -> TrajectoryIoError {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "s.xyz").unwrap();
    let error = writer.write_sample_set(set).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::UnsupportedField);
    assert_eq!(error.field(), Some(field));
    assert_eq!(error.frame(), Some(frame));
    assert!(writer.writer().is_empty());
    assert_eq!(writer.report().omissions(), [BONDS]);
    error
}

fn state_omission(field: StructureField, structures: usize) -> WriteOmission {
    WriteOmission::StructureState { field, structures }
}

#[test]
fn a_step_or_coordinate_order_writes_one_frame_per_sample_once_its_places_are_omitted() {
    let history = history();
    let error = assert_rejects(&history, options(), TrajectoryField::OrderStep, 0);
    assert_eq!(error.detail(), Some("XYZ cannot store the order step."));
    let (text, report) = write(&history, omit_places()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::OrderStep, 3)]
    );

    let scan = scan();
    let error = assert_rejects(&scan, options(), TrajectoryField::OrderCoordinate, 0);
    assert_eq!(
        error.detail(),
        Some("XYZ cannot store the order coordinate.")
    );
    let (text, report) = write(&scan, omit_places()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::OrderCoordinate, 3)]
    );

    // The file reads back as the samples' positions in set order.
    let read = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret_sample_set(Default::default())
        .unwrap()
        .into_sample_set();
    assert_eq!(read.len(), history.len());
    for (sample, read) in history.iter().zip(&read) {
        let written = sample.positions().values();
        let read = read.positions().values();
        for (a, b) in written.value().iter().zip(read.value().iter()) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }
}

#[test]
fn sample_energies_and_set_values_are_rejected_or_omitted_with_a_report() {
    let energy = |value: f64| PropertyValue::real(value, KILOJOULE_PER_MOLE).unwrap();
    let key = || PropertyKey::new("energy").unwrap();
    let mut history = history();
    for index in 1..3 {
        history
            .get_mut(index)
            .unwrap()
            .properties_mut()
            .owner_mut()
            .insert(key(), energy(-1.0 - index as f64))
            .unwrap();
    }
    assert_rejects(&history, omit_places(), TrajectoryField::Properties, 1);
    let options = omit_places().with_frame_properties(UnrepresentablePolicy::Omit);
    let (text, report) = write(&history, options.clone()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(StructureField::OrderStep, 3),
            state_omission(StructureField::Properties, 2),
        ]
    );

    // Values of the whole set have nowhere to go.
    history
        .properties_mut()
        .insert(
            PropertyKey::new("method").unwrap(),
            PropertyValue::String("GFN2-xTB".into()),
        )
        .unwrap();
    assert_rejects(&history, options, TrajectoryField::CollectionProperties, 0);
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn cells_are_rejected_or_omitted_and_samples_follow_written_frames() {
    let history = history();
    let topology = history.shared_topology();
    let boxed = in_cells(history.clone());
    assert_rejects(&boxed, omit_places(), TrajectoryField::PeriodicCell, 0);
    let omitted = omit_places().with_cell(UnrepresentablePolicy::Omit);
    let (text, report) = write(&boxed, omitted.clone()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(StructureField::OrderStep, 3),
            WriteOmission::Cell { structures: 3 },
        ]
    );

    // After a frame, errors name the frame a sample would have become.
    let frame = Frame::new(history.get(0).unwrap().positions().clone()).bound(&topology);
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options(), "s.xyz").unwrap();
    writer.write_frame(frame.view(&topology)).unwrap();
    let error = writer.write_sample_set(&history).unwrap_err();
    assert_eq!(error.field(), Some(TrajectoryField::OrderStep));
    assert_eq!(error.frame(), Some(1));
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), omitted, "s.xyz").unwrap();
    writer.write_frame(frame.view(&topology)).unwrap();
    writer.write_sample_set(&boxed).unwrap();
    let written = String::from_utf8(writer.finish().unwrap()).unwrap();
    assert_eq!(written.lines().count(), 4 * 5);
}

/// Each set, and what it writes under each pair of place and other
/// policies, is what the series it replaces converted to and wrote, as
/// recorded before series became sample sets.
#[test]
fn ordered_sets_hold_and_write_what_the_series_they_replace_did() {
    let outcome = |result: Result<(String, WriteReport), TrajectoryIoError>| {
        result.map_err(|error| (error.kind(), error.field(), error.frame()))
    };
    let mut valued = history();
    valued
        .get_mut(1)
        .unwrap()
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("energy").unwrap(),
            PropertyValue::real(-1.0, KILOJOULE_PER_MOLE).unwrap(),
        )
        .unwrap();
    let fail = UnrepresentablePolicy::Fail;
    let omit = UnrepresentablePolicy::Omit;
    // Places, then sample properties and cells. A set checks every sample's
    // properties before any place, so sample properties are checked only
    // with places omitted.
    let cases = [
        (history(), vec![(fail, fail), (omit, fail)]),
        (scan(), vec![(fail, fail), (omit, fail)]),
        (valued, vec![(omit, fail), (omit, omit)]),
        (boxed(), vec![(fail, fail), (omit, fail), (omit, omit)]),
    ];
    let mut recorded = Recorded::new(&[
        0xa343_d6fe_8905_ac25, // history
        0x5d60_039d_d2ae_01d3, // history, places fail, rest fail
        0x557b_ee40_7701_4a6f, // history, places omit, rest fail
        0xc7a2_15ef_ec11_d93f, // scan
        0x7cb3_6b28_1b51_187f, // scan, places fail, rest fail
        0x779f_1eda_44dc_5c93, // scan, places omit, rest fail
        0x19f5_8772_5396_7c7b, // valued
        0x5001_e09a_cc79_1f23, // valued, places omit, rest fail
        0xd56e_4571_1877_e511, // valued, places omit, rest omit
        0xe59a_d794_16ac_794e, // boxed
        0x13e3_9b29_c744_e766, // boxed, places fail, rest fail
        0x13e3_9b29_c744_e766, // boxed, places omit, rest fail
        0xfc16_a448_989e_9f4e, // boxed, places omit, rest omit
    ]);
    for (set, policies) in cases {
        recorded.check(&support::contents(&set));
        for (places, rest) in policies {
            recorded.check(&outcome(write(
                &set,
                options()
                    .with_sample_metadata(places)
                    .with_frame_properties(rest)
                    .with_cell(rest),
            )));
        }
    }
}

#[test]
fn foreign_and_empty_ordered_sets_are_rejected() {
    let history = history();
    let other = kekule::smiles::to_topology("O").unwrap();
    let mut writer = XyzWriter::new(Vec::new(), other, omit_places(), "s.xyz").unwrap();
    let error = writer.write_sample_set(&history).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::TopologyMismatch);

    let mut empty = SampleSet::new(history.shared_topology());
    empty.set_order(Some(Order::Step(Vec::new()))).unwrap();
    let mut writer = XyzWriter::new(
        Vec::new(),
        history.shared_topology(),
        omit_places(),
        "s.xyz",
    )
    .unwrap();
    let error = writer.write_sample_set(&empty).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
    assert!(writer.writer().is_empty());
}
