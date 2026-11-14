//! XYZ writing of structure sets in a step or coordinate order, one frame per
//! geometry.

use std::sync::Arc;

use kekule::diagnostics::{GeometryField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::math::{PeriodicCell, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Geometry, Order, Place, StructureSet};
use kekule::traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule::traj::{TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryWriter};
use kekule::units::{Quantity, ANGSTROM, KILOJOULE_PER_MOLE, NANOMETER};

use crate::support::Frame;

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
/// geometry at the place `place` gives its index.
fn ordered(order: Order, place: impl Fn(usize) -> Place) -> StructureSet {
    let read = kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    let mut set = StructureSet::new(read.shared_topology());
    set.set_order(Some(order)).unwrap();
    for (index, geometry) in read.iter().enumerate() {
        let geometry = Geometry::new(geometry.positions().clone());
        set.push(geometry, Some(place(index))).unwrap();
    }
    set
}

/// The frames of [`WATERS`] at steps 0, 5, and 10, as an optimization
/// history.
fn history() -> StructureSet {
    ordered(Order::Step(Vec::new()), |index| {
        Place::Step(5 * index as u64)
    })
}

/// The frames of [`WATERS`] as a scan of a distance in angstroms.
fn scan() -> StructureSet {
    let order = Order::Coordinate {
        unit: ANGSTROM,
        values: Vec::new(),
    };
    ordered(order, |index| {
        Place::Coordinate(Quantity::new(0.9 + 0.1 * index as f64, ANGSTROM))
    })
}

/// [`history`] with every geometry in [`cell`].
fn boxed() -> StructureSet {
    in_cells(history())
}

/// `set` with every geometry in [`cell`].
fn in_cells(mut set: StructureSet) -> StructureSet {
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
    options().with_geometry_metadata(UnrepresentablePolicy::Omit)
}

/// Writes `set` and returns the text and report.
fn write(
    set: &StructureSet,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), TrajectoryIoError> {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "s.xyz")?;
    writer.write_structure_set(set)?;
    let report = writer.report().clone();
    Ok((String::from_utf8(writer.finish()?).unwrap(), report))
}

/// Asserts that writing `set` with `options` fails before any byte for
/// `field` at `frame`, and returns the error.
fn assert_rejects(
    set: &StructureSet,
    options: XyzWriteOptions,
    field: GeometryField,
    frame: usize,
) -> TrajectoryIoError {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "s.xyz").unwrap();
    let error = writer.write_structure_set(set).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::UnsupportedField);
    assert_eq!(error.field(), Some(field));
    assert_eq!(error.frame(), Some(frame));
    assert!(writer.writer().is_empty());
    assert_eq!(writer.report().omissions(), [BONDS]);
    error
}

fn state_omission(field: GeometryField, geometries: usize) -> WriteOmission {
    WriteOmission::GeometryState { field, geometries }
}

#[test]
fn a_step_or_coordinate_order_writes_one_frame_per_geometry_once_its_places_are_omitted() {
    let history = history();
    let error = assert_rejects(&history, options(), GeometryField::OrderStep, 0);
    assert_eq!(error.detail(), Some("XYZ cannot store the order step."));
    let (text, report) = write(&history, omit_places()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::OrderStep, 3)]
    );

    let scan = scan();
    let error = assert_rejects(&scan, options(), GeometryField::OrderCoordinate, 0);
    assert_eq!(
        error.detail(),
        Some("XYZ cannot store the order coordinate.")
    );
    let (text, report) = write(&scan, omit_places()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::OrderCoordinate, 3)]
    );

    // The file reads back as the geometries' positions in set order.
    let read = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    assert_eq!(read.len(), history.len());
    for (geometry, read) in history.iter().zip(&read) {
        let written = geometry.positions().values();
        let read = read.positions().values();
        for (a, b) in written.value().iter().zip(read.value().iter()) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }
}

#[test]
fn geometry_energies_and_set_values_are_rejected_or_omitted_with_a_report() {
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
    assert_rejects(&history, omit_places(), GeometryField::Properties, 1);
    let options = omit_places().with_frame_properties(UnrepresentablePolicy::Omit);
    let (text, report) = write(&history, options.clone()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(GeometryField::OrderStep, 3),
            state_omission(GeometryField::Properties, 2),
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
    assert_rejects(&history, options, GeometryField::SetProperties, 0);
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn cells_are_rejected_or_omitted_and_geometries_follow_written_frames() {
    let history = history();
    let topology = history.shared_topology();
    let boxed = in_cells(history.clone());
    assert_rejects(&boxed, omit_places(), GeometryField::Cell, 0);
    let omitted = omit_places().with_cell(UnrepresentablePolicy::Omit);
    let (text, report) = write(&boxed, omitted.clone()).unwrap();
    assert_eq!(text, TEXT);
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(GeometryField::OrderStep, 3),
            WriteOmission::Cell { structures: 3 },
        ]
    );

    // After a frame, errors name the frame a geometry would have become.
    let frame = Frame::new(history.get(0).unwrap().positions().clone()).bound(&topology);
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options(), "s.xyz").unwrap();
    writer.write_frame(frame.view(&topology)).unwrap();
    let error = writer.write_structure_set(&history).unwrap_err();
    assert_eq!(error.field(), Some(GeometryField::OrderStep));
    assert_eq!(error.frame(), Some(1));
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), omitted, "s.xyz").unwrap();
    writer.write_frame(frame.view(&topology)).unwrap();
    writer.write_structure_set(&boxed).unwrap();
    let written = String::from_utf8(writer.finish().unwrap()).unwrap();
    assert_eq!(written.lines().count(), 4 * 5);
}

/// What writing gives: the text and omissions, or the error's kind, field,
/// and frame.
type Outcome = Result<
    (String, Vec<WriteOmission>),
    (TrajectoryIoErrorKind, Option<GeometryField>, Option<usize>),
>;

fn refused(field: GeometryField, frame: usize) -> Outcome {
    Err((
        TrajectoryIoErrorKind::UnsupportedField,
        Some(field),
        Some(frame),
    ))
}

fn written(omissions: &[WriteOmission]) -> Outcome {
    Ok((TEXT.to_owned(), omissions.to_vec()))
}

/// Each set holds its order, and writes under each pair of place and other
/// policies as XYZ, which stores only elements and positions, allows.
#[test]
fn ordered_sets_write_what_xyz_can_store_under_each_policy() {
    let outcome = |result: Result<(String, WriteReport), TrajectoryIoError>| -> Outcome {
        result
            .map(|(text, report)| (text, report.omissions().to_vec()))
            .map_err(|error| (error.kind(), error.field(), error.frame()))
    };
    let history = history();
    assert_eq!(history.order(), Some(&Order::Step(vec![0, 5, 10])));
    let scan = scan();
    assert_eq!(
        scan.order(),
        Some(&Order::Coordinate {
            unit: ANGSTROM,
            values: (0..3).map(|index| 0.9 + 0.1 * index as f64).collect(),
        })
    );
    let mut valued = history.clone();
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
    let boxed = boxed();
    assert!(boxed
        .iter()
        .all(|geometry| geometry.cell() == Some(&cell())));
    let fail = UnrepresentablePolicy::Fail;
    let omit = UnrepresentablePolicy::Omit;
    let steps = state_omission(GeometryField::OrderStep, 3);
    // Places, then geometry properties and cells. A set checks every geometry's
    // properties before any place, so geometry properties are checked only
    // with places omitted. A cell is checked before any place.
    let cases = [
        (
            history,
            vec![
                (fail, fail, refused(GeometryField::OrderStep, 0)),
                (omit, fail, written(&[BONDS, steps])),
            ],
        ),
        (
            scan,
            vec![
                (fail, fail, refused(GeometryField::OrderCoordinate, 0)),
                (
                    omit,
                    fail,
                    written(&[BONDS, state_omission(GeometryField::OrderCoordinate, 3)]),
                ),
            ],
        ),
        (
            valued,
            vec![
                (omit, fail, refused(GeometryField::Properties, 1)),
                (
                    omit,
                    omit,
                    written(&[BONDS, steps, state_omission(GeometryField::Properties, 1)]),
                ),
            ],
        ),
        (
            boxed,
            vec![
                (fail, fail, refused(GeometryField::Cell, 0)),
                (omit, fail, refused(GeometryField::Cell, 0)),
                (
                    omit,
                    omit,
                    written(&[BONDS, steps, WriteOmission::Cell { structures: 3 }]),
                ),
            ],
        ),
    ];
    for (set, policies) in cases {
        for (places, rest, expected) in policies {
            let options = options()
                .with_geometry_metadata(places)
                .with_frame_properties(rest)
                .with_cell(rest);
            assert_eq!(
                outcome(write(&set, options)),
                expected,
                "places {places:?}, rest {rest:?}"
            );
        }
    }
}

#[test]
fn foreign_and_empty_ordered_sets_are_rejected() {
    let history = history();
    let other = kekule::smiles::to_topology("O").unwrap();
    let mut writer = XyzWriter::new(Vec::new(), other, omit_places(), "s.xyz").unwrap();
    let error = writer.write_structure_set(&history).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::TopologyMismatch);

    let mut empty = StructureSet::new(history.shared_topology());
    empty.set_order(Some(Order::Step(Vec::new()))).unwrap();
    let mut writer = XyzWriter::new(
        Vec::new(),
        history.shared_topology(),
        omit_places(),
        "s.xyz",
    )
    .unwrap();
    let error = writer.write_structure_set(&empty).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
    assert!(writer.writer().is_empty());
}
