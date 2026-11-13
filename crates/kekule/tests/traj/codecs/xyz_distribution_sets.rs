//! XYZ writing of structure sets with a distribution, one frame per geometry.

use std::sync::Arc;

use kekule::diagnostics::{GeometryField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::math::{PeriodicCell, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Members, Order, Partition, StructureSet, TimeAxis, Velocities,
};
use kekule::traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule::traj::{TrajectoryField, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryWriter};
use kekule::units::{
    Quantity, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND,
};

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

/// The three geometries `kekule::xyz` interprets from [`WATERS`] in one state
/// of equally weighted members.
fn waters() -> StructureSet {
    let mut set = kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    set.set_distribution(Some(
        Distribution::new(Members::All, Partition::one_state(3))
            .unwrap()
            .with_weights(vec![1.0; 3])
            .unwrap(),
    ))
    .unwrap();
    set
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
fn an_equally_weighted_set_writes_one_frame_per_geometry_and_reads_back() {
    let set = waters();
    let (text, report) = write(&set, options()).unwrap();
    assert_eq!(report.omissions(), [BONDS]);
    assert_eq!(
        text,
        "3\nsample\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n\
         3\nsample\nO 0.000 0.000 0.100\nH 0.957 0.000 0.100\nH -0.240 0.927 0.100\n\
         3\nsample\nO 0.000 0.000 0.200\nH 0.957 0.000 0.200\nH -0.240 0.927 0.200\n"
    );

    // Reading the file back gives the geometries without an order or a
    // distribution: one state of equally weighted geometries.
    let read = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    assert_eq!(read.len(), set.len());
    assert_eq!(read.distribution(), None);
    assert_eq!(read.order(), None);
    for (written, read) in set.iter().zip(&read) {
        let written = written.positions().values();
        let read = read.positions().values();
        for (a, b) in written.value().iter().zip(read.value().iter()) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }
}

#[test]
fn unequal_weights_are_rejected_or_omitted_with_a_report() {
    let mut set = waters();
    set.distribution_mut().unwrap().set_weight(2, 3.0).unwrap();
    let error = assert_rejects(&set, options(), TrajectoryField::Weight, 2);
    assert_eq!(error.detail(), Some("XYZ cannot store the weight."));
    let (text, report) = write(&set, omit()).unwrap();
    assert_eq!(text, plain());
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Weight, 3)]
    );
}

#[test]
fn several_states_are_rejected_or_omitted_with_a_report() {
    let mut set = waters();
    set.distribution_mut()
        .unwrap()
        .regroup(Partition::new(vec![0, 1, 1], vec![0, 1]).unwrap())
        .unwrap();
    let error = assert_rejects(&set, options(), TrajectoryField::State, 1);
    assert_eq!(error.detail(), Some("XYZ cannot store the state."));
    let (_, report) = write(&set, omit()).unwrap();
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::State, 3)]
    );
}

#[test]
fn energies_are_rejected_or_omitted_with_a_report() {
    let mut set = waters();
    let mut distribution = set.distribution_mut().unwrap();
    distribution
        .regroup(Partition::new(vec![0, 1, 1], vec![0, 1]).unwrap())
        .unwrap();
    // Only the second state, of the last two geometries, has an energy.
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-1.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    let error = assert_rejects(&set, options(), TrajectoryField::Energy, 1);
    assert_eq!(error.detail(), Some("XYZ cannot store the state energy."));
    let (_, report) = write(&set, omit()).unwrap();
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(GeometryField::Energy, 2),
            state_omission(GeometryField::State, 3),
        ]
    );
}

#[test]
fn a_temperature_is_rejected_or_omitted_with_a_report() {
    let mut set = waters();
    set.distribution_mut()
        .unwrap()
        .set_temperature(Some(Quantity::new(298.15, KELVIN)))
        .unwrap();
    let error = assert_rejects(&set, options(), TrajectoryField::Temperature, 0);
    assert_eq!(error.detail(), Some("XYZ cannot store the temperature."));
    let (_, report) = write(&set, omit()).unwrap();
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Temperature, 3)]
    );
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
fn every_kind_of_geometry_state_is_reported_together() {
    let mut set = waters();
    let mut distribution = set.distribution_mut().unwrap();
    distribution
        .regroup(Partition::state_per_sample(3))
        .unwrap();
    distribution.set_weight(0, 2.0).unwrap();
    distribution
        .state_mut(0)
        .unwrap()
        .set_energy(Some(Quantity::new(0.0, KILOJOULE_PER_MOLE)))
        .unwrap();
    distribution
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    // The first geometry's energy comes before its weight.
    assert_rejects(&set, options(), TrajectoryField::Energy, 0);
    let (_, report) = write(&set, omit()).unwrap();
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(GeometryField::Weight, 3),
            state_omission(GeometryField::Energy, 1),
            state_omission(GeometryField::State, 3),
            state_omission(GeometryField::Temperature, 3),
        ]
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
    let mut distribution = set.distribution_mut().unwrap();
    distribution
        .regroup(Partition::new(vec![0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key(), value())
        .unwrap();
    // Named by the state's representative.
    let error = assert_rejects(&set, options(), TrajectoryField::StateProperties, 2);
    assert_eq!(
        error.detail(),
        Some("XYZ cannot store the state properties.")
    );

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

#[test]
fn geometry_state_rules_match_the_pdb_writer() {
    // The XYZ writer's geometry-state rules are kekule's PDB and mmCIF ones.
    // PDB omits the fixed implicit-hydrogen counts XYZ interpretation gives.
    let mut set = waters();
    let mut distribution = set.distribution_mut().unwrap();
    distribution.set_weight(1, 0.5).unwrap();
    distribution
        .set_temperature(Some(Quantity::new(310.0, KELVIN)))
        .unwrap();
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
    assert_eq!(
        state(&xyz),
        [
            state_omission(GeometryField::Weight, 3),
            state_omission(GeometryField::Temperature, 3),
        ]
    );
}

/// The text and report, or the error's kind, field, and frame.
type Outcome = Result<
    (String, WriteReport),
    (
        TrajectoryIoErrorKind,
        Option<TrajectoryField>,
        Option<usize>,
    ),
>;

fn outcome(result: Result<(String, WriteReport), TrajectoryIoError>) -> Outcome {
    result.map_err(|error| (error.kind(), error.field(), error.frame()))
}

/// Each set, and what it writes under each policy, is what the ensemble it
/// replaces converted to and wrote, as recorded before ensembles became
/// structure sets.
#[test]
fn sets_with_a_distribution_hold_and_write_what_the_ensembles_they_replace_did() {
    let mut weighted = waters();
    weighted
        .distribution_mut()
        .unwrap()
        .set_weight(2, 3.0)
        .unwrap();
    let mut clustered = waters();
    let mut distribution = clustered.distribution_mut().unwrap();
    distribution
        .regroup(Partition::new(vec![0, 1, 1], vec![0, 1]).unwrap())
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-1.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    distribution
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
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
    let mut recorded = Recorded::new(&[
        0xf14d_cdd2_f26b_1088, // waters
        0x4939_b251_4e51_e9cb, // waters, fail
        0x4939_b251_4e51_e9cb, // waters, omit
        0x7e93_965f_759b_09ee, // weighted
        0x6c2d_94d0_1919_47c1, // weighted, fail
        0x9d0c_c4f1_8a4c_1fcd, // weighted, omit
        0x0ba9_867b_4244_e045, // clustered
        0x3ebe_6cee_6493_523a, // clustered, fail
        0x2c01_a7fc_aaa6_ee15, // clustered, omit
        0xf3a1_d233_9761_76d3, // boxed
        0x0ebd_2b21_c872_bfdf, // boxed, fail
        0x03c3_ae9e_91dc_8c0c, // boxed, omit
        0x77cc_a4b8_28e4_f383, // occupied
        0x447b_f888_5d4d_1835, // occupied, fail
        0x447b_f888_5d4d_1835, // occupied, omit
        0xd016_a79f_cb50_844c, // labelled
        0x5001_e09a_cc79_1f23, // labelled, fail
        0x5001_e09a_cc79_1f23, // labelled, omit
    ]);
    for set in [waters(), weighted, clustered, boxed, occupied, labelled] {
        recorded.check(&support::contents(&set));
        for policy in [UnrepresentablePolicy::Fail, UnrepresentablePolicy::Omit] {
            recorded.check(&outcome(write(
                &set,
                options().with_geometry_metadata(policy).with_cell(policy),
            )));
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

    // A lone state represented by its second member.
    let mut late = waters();
    late.set_distribution(Some(
        Distribution::new(
            Members::All,
            Partition::new(vec![0, 0, 0], vec![1]).unwrap(),
        )
        .unwrap(),
    ))
    .unwrap();
    assert_rejects(&late, options(), TrajectoryField::Representative, 1);
    let (text, report) = write(&late, omit()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Representative, 1)]
    );

    // Velocities on one geometry.
    let mut moving = waters();
    moving
        .get_mut(1)
        .unwrap()
        .set_velocities(Some(
            Velocities::from_vec(Quantity::new(
                vec![Vector3::new(0.0, 0.0, 0.0); 3],
                NANOMETER_PER_PICOSECOND,
            ))
            .unwrap(),
        ))
        .unwrap();
    assert_rejects(&moving, options(), TrajectoryField::Velocities, 1);
    let (text, report) = write(&moving, omit()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Velocities, 1)]
    );

    // A geometry outside the distribution: XYZ reads every frame back as a
    // member.
    let mut subset = waters();
    subset
        .set_distribution(Some(
            Distribution::new(Members::Subset(vec![0, 2]), Partition::one_state(2)).unwrap(),
        ))
        .unwrap();
    assert_rejects(&subset, options(), TrajectoryField::Membership, 1);
    let (text, report) = write(&subset, omit()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(GeometryField::Membership, 1)]
    );

    // Without a distribution, the frames read back as one state of equally
    // weighted geometries, as from a set with one.
    let mut alone = waters();
    alone.set_distribution(None).unwrap();
    assert_eq!(write(&alone, options()).unwrap().0, plain);
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
