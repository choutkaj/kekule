//! XYZ writing of series, one frame per item.

use std::sync::Arc;

use kekule::diagnostics::{StructureField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::geometry::{PeriodicCell, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Series, SeriesAxis, SeriesItem, TrajectoryFrame};
use kekule::units::{Quantity, ANGSTROM, KILOJOULE_PER_MOLE, NANOMETER};
use kekule_traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule_traj::{TrajectoryField, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryWriter};

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

/// The frames of [`WATERS`] at steps 0, 5, and 10, as an optimization
/// history.
fn history() -> Series {
    let ensemble = kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap()
        .into_ensemble();
    let items = ensemble.iter().enumerate().map(|(index, sample)| {
        let positions = sample.as_model_view().positions().clone();
        SeriesItem::at_step(positions, 5 * index as u64)
    });
    Series::from_items(ensemble.shared_topology(), SeriesAxis::Step, items).unwrap()
}

/// The frames of [`WATERS`] as a scan of a distance in angstroms.
fn scan() -> Series {
    let history = history();
    let items = history.iter().enumerate().map(|(index, item)| {
        let place = Quantity::new(0.9 + 0.1 * index as f64, ANGSTROM);
        SeriesItem::at_coordinate(item.positions().clone(), place).unwrap()
    });
    Series::from_items(
        history.shared_topology(),
        SeriesAxis::Coordinate(ANGSTROM),
        items,
    )
    .unwrap()
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
    options().with_series_metadata(UnrepresentablePolicy::Omit)
}

/// Writes `series` and returns the text and report.
fn write(
    series: &Series,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), TrajectoryIoError> {
    let mut writer = XyzWriter::new(Vec::new(), series.shared_topology(), options, "s.xyz")?;
    writer.write_series(series)?;
    let report = writer.report().clone();
    Ok((String::from_utf8(writer.finish()?).unwrap(), report))
}

/// Asserts that writing `series` with `options` fails before any byte for
/// `field` at `frame`.
fn assert_rejects(series: &Series, options: XyzWriteOptions, field: TrajectoryField, frame: usize) {
    let mut writer =
        XyzWriter::new(Vec::new(), series.shared_topology(), options, "s.xyz").unwrap();
    let error = writer.write_series(series).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::UnsupportedField);
    assert_eq!(error.field(), Some(field));
    assert_eq!(error.frame(), Some(frame));
    assert!(writer.writer().is_empty());
    assert_eq!(writer.report().omissions(), [BONDS]);
}

fn state_omission(field: StructureField, structures: usize) -> WriteOmission {
    WriteOmission::StructureState { field, structures }
}

#[test]
fn a_series_writes_one_frame_per_item_once_its_places_are_omitted() {
    let history = history();
    assert_rejects(&history, options(), TrajectoryField::ItemStep, 0);
    let (text, report) = write(&history, omit_places()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::ItemStep, 3)]
    );

    let scan = scan();
    assert_rejects(&scan, options(), TrajectoryField::ItemCoordinate, 0);
    let (text, report) = write(&scan, omit_places()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::ItemCoordinate, 3)]
    );

    // The file reads back as the items' positions in series order.
    let read = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap()
        .into_ensemble();
    assert_eq!(read.len(), history.len());
    for (item, read) in history.iter().zip(&read) {
        let written = item.positions().values();
        let read = read.as_model_view().positions().values();
        for (a, b) in written.value().iter().zip(read.value().iter()) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }
}

#[test]
fn item_energies_and_series_values_are_rejected_or_omitted_with_a_report() {
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
            state_omission(StructureField::ItemStep, 3),
            state_omission(StructureField::Properties, 2),
        ]
    );

    // Values of the whole series have nowhere to go.
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
fn cells_are_rejected_or_omitted_and_items_follow_written_frames() {
    let history = history();
    let topology = history.shared_topology();
    let items = history
        .iter()
        .map(|item| item.payload().clone().with_cell(cell()));
    let boxed = Series::from_items(Arc::clone(&topology), SeriesAxis::Step, items).unwrap();
    assert_rejects(&boxed, omit_places(), TrajectoryField::PeriodicCell, 0);
    let omitted = omit_places().with_cell(UnrepresentablePolicy::Omit);
    let (text, report) = write(&boxed, omitted.clone()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(StructureField::ItemStep, 3),
            WriteOmission::Cell { structures: 3 },
        ]
    );

    // After a frame, errors name the frame an item would have become.
    let frame = TrajectoryFrame::new(history.get(0).unwrap().positions().clone());
    let trajectory =
        kekule::structure::Trajectory::from_items(Arc::clone(&topology), [frame]).unwrap();
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options(), "s.xyz").unwrap();
    writer.write_frame(trajectory.get(0).unwrap()).unwrap();
    let error = writer.write_series(&history).unwrap_err();
    assert_eq!(error.field(), Some(TrajectoryField::ItemStep));
    assert_eq!(error.frame(), Some(1));
    let mut writer = XyzWriter::new(Vec::new(), topology, omitted, "s.xyz").unwrap();
    writer.write_frame(trajectory.get(0).unwrap()).unwrap();
    writer.write_series(&boxed).unwrap();
    let written = String::from_utf8(writer.finish().unwrap()).unwrap();
    assert_eq!(written.lines().count(), 4 * 5);
}

#[test]
fn foreign_and_empty_series_are_rejected() {
    let history = history();
    let other = kekule::smiles::to_topology("O").unwrap();
    let mut writer = XyzWriter::new(Vec::new(), other, omit_places(), "s.xyz").unwrap();
    let error = writer.write_series(&history).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::TopologyMismatch);

    let empty = Series::new(history.shared_topology(), SeriesAxis::Step);
    let mut writer = XyzWriter::new(
        Vec::new(),
        history.shared_topology(),
        omit_places(),
        "s.xyz",
    )
    .unwrap();
    let error = writer.write_series(&empty).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
    assert!(writer.writer().is_empty());
}
