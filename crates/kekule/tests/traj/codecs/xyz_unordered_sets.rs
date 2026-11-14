//! XYZ writing of unordered structure sets, one frame per geometry.

use std::sync::Arc;

use kekule::diagnostics::{GeometryField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::math::{PeriodicCell, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Order, StructureSet, TimeAxis, Velocities};
use kekule::traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule::traj::{TrajectoryField, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryWriter};
use kekule::units::{Quantity, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND};

use crate::support::{self, Frame};

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

/// The three unordered geometries `kekule::xyz` interprets from [`WATERS`].
fn waters() -> StructureSet {
    kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set()
}

/// What any XYZ file of [`waters`] leaves out: the inferred bonds.
const BONDS: WriteOmission = WriteOmission::Bonds { bonds: 2 };

fn options() -> XyzWriteOptions {
    XyzWriteOptions::default()
        .with_decimal_places(3)
        .with_comment("sample")
}

fn omit() -> XyzWriteOptions {
    options().with_geometry_metadata(UnrepresentablePolicy::Omit)
}

/// Writes `set` and returns the text and report.
fn write(
    set: &StructureSet,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), TrajectoryIoError> {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "e.xyz")?;
    writer.write_structure_set(set)?;
    let report = writer.report().clone();
    Ok((String::from_utf8(writer.finish()?).unwrap(), report))
}

/// Asserts that writing `set` with `options` fails before any byte for
/// `field` at `frame`, and returns the error.
fn assert_rejects(
    set: &StructureSet,
    options: XyzWriteOptions,
    field: TrajectoryField,
    frame: usize,
) -> TrajectoryIoError {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "e.xyz").unwrap();
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

/// The text of [`waters`] under [`options`].
fn plain() -> String {
    write(&waters(), options()).unwrap().0
}

#[test]
fn a_set_writes_one_frame_per_geometry_and_reads_back() {
    let set = waters();
    let (text, report) = write(&set, options()).unwrap();
    assert_eq!(report.omissions(), [BONDS]);
    assert_eq!(
        text,
        "3\nsample\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n\
         3\nsample\nO 0.000 0.000 0.100\nH 0.957 0.000 0.100\nH -0.240 0.927 0.100\n\
         3\nsample\nO 0.000 0.000 0.200\nH 0.957 0.000 0.200\nH -0.240 0.927 0.200\n"
    );

    // Reading the file back gives the geometries without an order.
    let read = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    assert_eq!(read.len(), set.len());
    assert_eq!(read.order(), None);
    for (written, read) in set.iter().zip(&read) {
        let written = written.positions().values();
        let read = read.positions().values();
        for (a, b) in written.value().iter().zip(read.value().iter()) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn cells_are_rejected_or_omitted_with_a_report_as_for_frames() {
    let mut set = waters();
    set.get_mut(1).unwrap().set_cell(Some(cell()));
    set.get_mut(2).unwrap().set_cell(Some(cell()));
    assert_rejects(&set, options(), TrajectoryField::PeriodicCell, 1);
    // A cell is not geometry metadata.
    assert_rejects(&set, omit(), TrajectoryField::PeriodicCell, 1);
    let (text, report) = write(&set, options().with_cell(UnrepresentablePolicy::Omit)).unwrap();
    assert_eq!(text, plain());
    assert_eq!(
        report.omissions(),
        [BONDS, WriteOmission::Cell { structures: 2 }]
    );

    // Frames follow the same policy, and the report counts every frame and
    // geometry written.
    let topology = set.shared_topology();
    let frame = Frame::new(set.get(0).unwrap().positions().clone())
        .with_cell(cell())
        .bound(&topology);
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options(), "f.xyz").unwrap();
    let error = writer.write_frame(frame.view(&topology)).unwrap_err();
    assert_eq!(error.field(), Some(TrajectoryField::PeriodicCell));
    let mut writer = XyzWriter::new(
        Vec::new(),
        Arc::clone(&topology),
        options().with_cell(UnrepresentablePolicy::Omit),
        "f.xyz",
    )
    .unwrap();
    writer.write_frame(frame.view(&topology)).unwrap();
    assert_eq!(
        writer.report().omissions(),
        [BONDS, WriteOmission::Cell { structures: 1 }]
    );
    writer.write_structure_set(&set).unwrap();
    assert_eq!(
        writer.report().omissions(),
        [BONDS, WriteOmission::Cell { structures: 3 }]
    );
    assert_eq!(
        writer
            .finish()
            .unwrap()
            .iter()
            .filter(|&&b| b == b'\n')
            .count(),
        20
    );
}

#[test]
fn occupancies_and_b_factors_are_always_omitted_with_a_report() {
    let mut set = waters();
    let mut geometry = set.get_mut(0).unwrap();
    geometry
        .set_occupancies(Some(vec![Some(1.0), Some(0.5), None]))
        .unwrap();
    geometry
        .set_b_factors(Some(Quantity::new(
            vec![Some(10.0), None, None],
            kekule::units::SQUARE_ANGSTROM,
        )))
        .unwrap();
    set.get_mut(1)
        .unwrap()
        .set_occupancies(Some(vec![Some(1.0); 3]))
        .unwrap();
    let (text, report) = write(&set, options()).unwrap();
    assert_eq!(text, plain());
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            WriteOmission::Occupancies { structures: 2 },
            WriteOmission::BFactors { structures: 1 },
        ]
    );
}

#[test]
fn properties_are_rejected_at_every_level() {
    let key = || PropertyKey::new("label").unwrap();
    let value = || PropertyValue::Bool(true);

    let mut set = waters();
    set.properties_mut().insert(key(), value()).unwrap();
    assert_rejects(&set, options(), TrajectoryField::CollectionProperties, 0);

    let mut set = waters();
    set.get_mut(1)
        .unwrap()
        .properties_mut()
        .owner_mut()
        .insert(key(), value())
        .unwrap();
    assert_rejects(&set, options(), TrajectoryField::Properties, 1);
}

#[test]
fn foreign_and_empty_sets_are_rejected() {
    let set = waters();
    let other = kekule::smiles::to_topology("O").unwrap();
    let mut writer = XyzWriter::new(Vec::new(), other, options(), "e.xyz").unwrap();
    let error = writer.write_structure_set(&set).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::TopologyMismatch);

    let empty = StructureSet::new(set.shared_topology());
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options(), "e.xyz").unwrap();
    let error = writer.write_structure_set(&empty).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
    assert!(writer.writer().is_empty());
}

/// `set` with zero velocities on its geometry `index`.
fn with_velocities(mut set: StructureSet, index: usize) -> StructureSet {
    set.get_mut(index)
        .unwrap()
        .set_velocities(Some(
            Velocities::from_vec(Quantity::new(
                vec![Vector3::new(0.0, 0.0, 0.0); 3],
                NANOMETER_PER_PICOSECOND,
            ))
            .unwrap(),
        ))
        .unwrap();
    set
}

#[test]
fn geometry_state_rules_match_the_pdb_writer() {
    // The XYZ writer's geometry-state rules are kekule's PDB and mmCIF ones.
    // PDB omits the fixed implicit-hydrogen counts XYZ interpretation gives.
    let set = with_velocities(waters(), 1);
    let pdb = kekule::pdb::write(
        &set,
        kekule::pdb::PdbWriteOptions::default()
            .with_geometry_metadata(UnrepresentablePolicy::Omit)
            .with_implicit_hydrogens(UnrepresentablePolicy::Omit),
    )
    .unwrap()
    .1;
    let (_, xyz) = write(&set, omit()).unwrap();
    let state = |report: &WriteReport| {
        report
            .omissions()
            .iter()
            .copied()
            .filter(|omission| matches!(omission, WriteOmission::GeometryState { .. }))
            .collect::<Vec<_>>()
    };
    assert_eq!(state(&xyz), state(&pdb));
    assert_eq!(state(&xyz), [state_omission(GeometryField::Velocities, 1)]);
}

/// The error's kind, field, and frame.
fn failure(
    error: TrajectoryIoError,
) -> (
    TrajectoryIoErrorKind,
    Option<TrajectoryField>,
    Option<usize>,
) {
    (error.kind(), error.field(), error.frame())
}

#[test]
fn a_policy_changes_only_whether_what_xyz_cannot_store_is_dropped() {
    let mut boxed = waters();
    boxed.get_mut(1).unwrap().set_cell(Some(cell()));
    let mut occupied = waters();
    occupied
        .get_mut(0)
        .unwrap()
        .set_occupancies(Some(vec![Some(0.5), None, None]))
        .unwrap();
    let mut labelled = waters();
    labelled
        .get_mut(1)
        .unwrap()
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("label").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();
    let plain = plain();
    for policy in [UnrepresentablePolicy::Fail, UnrepresentablePolicy::Omit] {
        let options = || options().with_geometry_metadata(policy).with_cell(policy);

        // Positions and elements are all XYZ stores, under either policy.
        let (text, report) = write(&waters(), options()).unwrap();
        assert_eq!(text, plain);
        assert_eq!(report.omissions(), [BONDS]);

        // Occupancies are always dropped and reported, never refused.
        let (text, report) = write(&occupied, options()).unwrap();
        assert_eq!(text, plain);
        assert_eq!(
            report.omissions(),
            [BONDS, WriteOmission::Occupancies { structures: 1 }]
        );

        // Geometry properties are refused under either policy, at their
        // geometry.
        assert_eq!(
            failure(write(&labelled, options()).unwrap_err()),
            (
                TrajectoryIoErrorKind::UnsupportedField,
                Some(TrajectoryField::Properties),
                Some(1)
            )
        );

        // A cell is refused or dropped by the cell policy.
        let written = write(&boxed, options());
        match policy {
            UnrepresentablePolicy::Fail => assert_eq!(
                failure(written.unwrap_err()),
                (
                    TrajectoryIoErrorKind::UnsupportedField,
                    Some(TrajectoryField::PeriodicCell),
                    Some(1)
                )
            ),
            _ => {
                let (text, report) = written.unwrap();
                assert_eq!(text, plain);
                assert_eq!(
                    report.omissions(),
                    [BONDS, WriteOmission::Cell { structures: 1 }]
                );
            }
        }
    }
}

#[test]
fn geometry_state_xyz_cannot_hold_is_rejected_or_omitted_with_a_report() {
    let plain = plain();

    // A time order.
    let mut timed = waters();
    timed
        .set_order(Some(Order::Time(
            TimeAxis::new(Some(Quantity::new(vec![0.0, 1.0, 2.0], PICOSECOND)), None).unwrap(),
        )))
        .unwrap();
    assert_rejects(&timed, options(), TrajectoryField::Time, 0);
    let (text, report) = write(&timed, omit()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Time, 3)]
    );

    // A time order without times or steps reads back unordered.
    let mut bare = waters();
    bare.set_order(Some(Order::Time(TimeAxis::new(None, None).unwrap())))
        .unwrap();
    assert_rejects(&bare, options(), TrajectoryField::Order, 0);
    let (text, report) = write(&bare, omit()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Order, 3)]
    );

    // Velocities on one geometry.
    let moving = with_velocities(waters(), 1);
    assert_rejects(&moving, options(), TrajectoryField::Velocities, 1);
    let (text, report) = write(&moving, omit()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Velocities, 1)]
    );
}

#[test]
fn path_writers_report_the_cells_their_frames_omit() {
    let set = waters();
    let frames = set
        .iter()
        .map(|geometry| Frame::new(geometry.positions().clone()).with_cell(cell()));
    let frames = support::time_ordered(set.shared_topology(), frames);
    let path = std::env::temp_dir().join(format!("kekule-traj-{}-cells.xyz", std::process::id()));
    let saved = kekule::traj::io::write_structure_set(
        &path,
        &frames,
        kekule::traj::io::TrajectoryWriteOptions::default()
            .with_xyz_options(options().with_cell(UnrepresentablePolicy::Omit)),
    );
    let _ = std::fs::remove_file(&path);
    assert_eq!(
        saved.unwrap().omissions(),
        [BONDS, WriteOmission::Cell { structures: 3 }]
    );
}
