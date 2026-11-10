//! XYZ writing of ensembles, one frame per sample.

use std::sync::Arc;

use kekule::diagnostics::{StructureField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::geometry::{PeriodicCell, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Ensemble, EnsemblePartition, Members, Order, SampleSet, TimeAxis,
    TrajectoryFrame, Velocities,
};
use kekule::traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule::traj::{TrajectoryField, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryWriter};
use kekule::units::{
    Quantity, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND,
};

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

/// The three-sample, one-state, equally weighted ensemble `kekule::xyz`
/// interprets from [`WATERS`].
fn waters() -> Ensemble {
    kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap()
        .into_ensemble()
}

/// What any XYZ file of [`waters`] leaves out: the inferred bonds.
const BONDS: WriteOmission = WriteOmission::Bonds { bonds: 2 };

fn options() -> XyzWriteOptions {
    XyzWriteOptions::default()
        .with_decimal_places(3)
        .with_comment("sample")
}

/// Writes `ensemble` and returns the text and report.
fn write(
    ensemble: &Ensemble,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), TrajectoryIoError> {
    let mut writer = XyzWriter::new(Vec::new(), ensemble.shared_topology(), options, "e.xyz")?;
    writer.write_ensemble(ensemble)?;
    let report = writer.report().clone();
    Ok((String::from_utf8(writer.finish()?).unwrap(), report))
}

/// Asserts that writing `ensemble` with `options` fails before any byte for
/// `field` at `frame`.
fn assert_rejects(
    ensemble: &Ensemble,
    options: XyzWriteOptions,
    field: TrajectoryField,
    frame: usize,
) {
    let mut writer =
        XyzWriter::new(Vec::new(), ensemble.shared_topology(), options, "e.xyz").unwrap();
    let error = writer.write_ensemble(ensemble).unwrap_err();
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
fn an_equally_weighted_ensemble_writes_one_frame_per_sample_and_reads_back() {
    let ensemble = waters();
    let (text, report) = write(&ensemble, options()).unwrap();
    assert_eq!(report.omissions(), [BONDS]);
    assert_eq!(
        text,
        "3\nsample\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n\
         3\nsample\nO 0.000 0.000 0.100\nH 0.957 0.000 0.100\nH -0.240 0.927 0.100\n\
         3\nsample\nO 0.000 0.000 0.200\nH 0.957 0.000 0.200\nH -0.240 0.927 0.200\n"
    );

    // Reading the file back gives one state of equally weighted samples.
    let read = kekule::xyz::parse_str(&text, Default::default())
        .unwrap()
        .interpret_ensemble(Default::default())
        .unwrap()
        .into_ensemble();
    assert_eq!(read.len(), ensemble.len());
    assert_eq!(read.state_count(), 1);
    assert_eq!(read.temperature(), None);
    for (written, read) in ensemble.iter().zip(&read) {
        assert_eq!(read.weight(), 1.0);
        assert_eq!(read.state().payload().energy(), None);
        let written = written.as_model_view().positions().values();
        let read = read.as_model_view().positions().values();
        for (a, b) in written.value().iter().zip(read.value().iter()) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }
}

#[test]
fn unequal_weights_are_rejected_or_omitted_with_a_report() {
    let mut ensemble = waters();
    ensemble.get_mut(2).unwrap().set_weight(3.0).unwrap();
    assert_rejects(&ensemble, options(), TrajectoryField::Weight, 2);
    let omit = options().with_ensemble_metadata(UnrepresentablePolicy::Omit);
    let (text, report) = write(&ensemble, omit).unwrap();
    assert_eq!(text, write(&waters(), options()).unwrap().0);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::Weight, 3)]
    );
}

#[test]
fn several_states_are_rejected_or_omitted_with_a_report() {
    let mut ensemble = waters();
    ensemble
        .regroup(EnsemblePartition::new(vec![0, 1, 1], vec![0, 1]).unwrap())
        .unwrap();
    assert_rejects(&ensemble, options(), TrajectoryField::State, 1);
    let omit = options().with_ensemble_metadata(UnrepresentablePolicy::Omit);
    let (_, report) = write(&ensemble, omit).unwrap();
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::State, 3)]
    );
}

#[test]
fn energies_are_rejected_or_omitted_with_a_report() {
    let mut ensemble = waters();
    ensemble
        .regroup(EnsemblePartition::new(vec![0, 1, 1], vec![0, 1]).unwrap())
        .unwrap();
    // Only the second state, of the last two samples, has an energy.
    ensemble
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-1.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    assert_rejects(&ensemble, options(), TrajectoryField::Energy, 1);
    let omit = options().with_ensemble_metadata(UnrepresentablePolicy::Omit);
    let (_, report) = write(&ensemble, omit).unwrap();
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(StructureField::Energy, 2),
            state_omission(StructureField::State, 3),
        ]
    );
}

#[test]
fn a_temperature_is_rejected_or_omitted_with_a_report() {
    let mut ensemble = waters();
    ensemble
        .set_temperature(Some(Quantity::new(298.15, KELVIN)))
        .unwrap();
    assert_rejects(&ensemble, options(), TrajectoryField::Temperature, 0);
    let omit = options().with_ensemble_metadata(UnrepresentablePolicy::Omit);
    let (_, report) = write(&ensemble, omit).unwrap();
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::Temperature, 3)]
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
    let mut ensemble = waters();
    ensemble.get_mut(1).unwrap().set_cell(Some(cell()));
    ensemble.get_mut(2).unwrap().set_cell(Some(cell()));
    assert_rejects(&ensemble, options(), TrajectoryField::PeriodicCell, 1);
    // A cell is not ensemble metadata.
    assert_rejects(
        &ensemble,
        options().with_ensemble_metadata(UnrepresentablePolicy::Omit),
        TrajectoryField::PeriodicCell,
        1,
    );
    let (text, report) =
        write(&ensemble, options().with_cell(UnrepresentablePolicy::Omit)).unwrap();
    assert_eq!(text, write(&waters(), options()).unwrap().0);
    assert_eq!(
        report.omissions(),
        [BONDS, WriteOmission::Cell { structures: 2 }]
    );

    // Frames follow the same policy, and the report counts every frame and
    // sample written.
    let topology = ensemble.shared_topology();
    let frame = TrajectoryFrame::new(ensemble.get(0).unwrap().as_model_view().positions().clone())
        .with_cell(cell());
    let trajectory =
        kekule::structure::Trajectory::from_items(Arc::clone(&topology), [frame]).unwrap();
    let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options(), "f.xyz").unwrap();
    let error = writer.write_frame(trajectory.get(0).unwrap()).unwrap_err();
    assert_eq!(error.field(), Some(TrajectoryField::PeriodicCell));
    let mut writer = XyzWriter::new(
        Vec::new(),
        topology,
        options().with_cell(UnrepresentablePolicy::Omit),
        "f.xyz",
    )
    .unwrap();
    writer.write_frame(trajectory.get(0).unwrap()).unwrap();
    assert_eq!(
        writer.report().omissions(),
        [BONDS, WriteOmission::Cell { structures: 1 }]
    );
    writer.write_ensemble(&ensemble).unwrap();
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
fn every_kind_of_sample_state_is_reported_together() {
    let mut ensemble = waters();
    ensemble
        .regroup(EnsemblePartition::state_per_sample(3))
        .unwrap();
    ensemble.get_mut(0).unwrap().set_weight(2.0).unwrap();
    ensemble
        .state_mut(0)
        .unwrap()
        .set_energy(Some(Quantity::new(0.0, KILOJOULE_PER_MOLE)))
        .unwrap();
    ensemble
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    // The first sample's energy comes before its weight.
    assert_rejects(&ensemble, options(), TrajectoryField::Energy, 0);
    let omit = options().with_ensemble_metadata(UnrepresentablePolicy::Omit);
    let (_, report) = write(&ensemble, omit).unwrap();
    assert_eq!(
        report.omissions(),
        [
            BONDS,
            state_omission(StructureField::Weight, 3),
            state_omission(StructureField::Energy, 1),
            state_omission(StructureField::State, 3),
            state_omission(StructureField::Temperature, 3),
        ]
    );
}

#[test]
fn occupancies_and_b_factors_are_always_omitted_with_a_report() {
    let mut ensemble = waters();
    let mut sample = ensemble.get_mut(0).unwrap();
    sample
        .set_occupancies(Some(vec![Some(1.0), Some(0.5), None]))
        .unwrap();
    sample
        .set_b_factors(Some(Quantity::new(
            vec![Some(10.0), None, None],
            kekule::units::SQUARE_ANGSTROM,
        )))
        .unwrap();
    ensemble
        .get_mut(1)
        .unwrap()
        .set_occupancies(Some(vec![Some(1.0); 3]))
        .unwrap();
    let (text, report) = write(&ensemble, options()).unwrap();
    assert_eq!(text, write(&waters(), options()).unwrap().0);
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

    let mut ensemble = waters();
    ensemble.properties_mut().insert(key(), value()).unwrap();
    assert_rejects(
        &ensemble,
        options(),
        TrajectoryField::CollectionProperties,
        0,
    );

    let mut ensemble = waters();
    ensemble
        .regroup(EnsemblePartition::new(vec![0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    ensemble
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key(), value())
        .unwrap();
    // Named by the state's representative.
    assert_rejects(&ensemble, options(), TrajectoryField::StateProperties, 2);

    let mut ensemble = waters();
    ensemble
        .get_mut(1)
        .unwrap()
        .properties_mut()
        .owner_mut()
        .insert(key(), value())
        .unwrap();
    assert_rejects(&ensemble, options(), TrajectoryField::Properties, 1);
}

#[test]
fn foreign_and_empty_ensembles_are_rejected() {
    let ensemble = waters();
    let other = kekule::smiles::to_topology("O").unwrap();
    let mut writer = XyzWriter::new(Vec::new(), other, options(), "e.xyz").unwrap();
    let error = writer.write_ensemble(&ensemble).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::TopologyMismatch);

    let empty = Ensemble::new(ensemble.shared_topology());
    let mut writer =
        XyzWriter::new(Vec::new(), ensemble.shared_topology(), options(), "e.xyz").unwrap();
    let error = writer.write_ensemble(&empty).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
    assert!(writer.writer().is_empty());
}

#[test]
fn ensemble_rules_match_the_pdb_writer() {
    // The XYZ writer's sample-state rules are kekule's PDB and mmCIF ones.
    // PDB omits the fixed implicit-hydrogen counts XYZ interpretation gives.
    let mut ensemble = waters();
    ensemble.get_mut(1).unwrap().set_weight(0.5).unwrap();
    ensemble
        .set_temperature(Some(Quantity::new(310.0, KELVIN)))
        .unwrap();
    let pdb = kekule::pdb::write(
        &ensemble,
        kekule::pdb::PdbWriteOptions::default()
            .with_ensemble_metadata(UnrepresentablePolicy::Omit)
            .with_implicit_hydrogens(UnrepresentablePolicy::Omit),
    )
    .unwrap()
    .1;
    let (_, xyz) = write(
        &ensemble,
        options().with_ensemble_metadata(UnrepresentablePolicy::Omit),
    )
    .unwrap();
    let state = |report: &WriteReport| {
        report
            .omissions()
            .iter()
            .copied()
            .filter(|omission| matches!(omission, WriteOmission::StructureState { .. }))
            .collect::<Vec<_>>()
    };
    assert_eq!(state(&xyz), state(&pdb));
    assert_eq!(
        state(&xyz),
        [
            state_omission(StructureField::Weight, 3),
            state_omission(StructureField::Temperature, 3),
        ]
    );
}

/// Writes `set` and returns the text and report.
fn write_set(
    set: &SampleSet,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), TrajectoryIoError> {
    let mut writer = XyzWriter::new(Vec::new(), set.shared_topology(), options, "e.xyz")?;
    writer.write_sample_set(set)?;
    let report = writer.report().clone();
    Ok((String::from_utf8(writer.finish()?).unwrap(), report))
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

#[test]
fn a_sample_set_writes_what_the_ensemble_it_converts_from_writes() {
    let mut weighted = waters();
    weighted.get_mut(2).unwrap().set_weight(3.0).unwrap();
    let mut clustered = waters();
    clustered
        .regroup(EnsemblePartition::new(vec![0, 1, 1], vec![0, 1]).unwrap())
        .unwrap();
    clustered
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-1.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    clustered
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
    for (index, ensemble) in [waters(), weighted, clustered, boxed, occupied, labelled]
        .into_iter()
        .enumerate()
    {
        let set = SampleSet::from_ensemble(ensemble.clone());
        for policy in [UnrepresentablePolicy::Fail, UnrepresentablePolicy::Omit] {
            let old = write(
                &ensemble,
                options().with_ensemble_metadata(policy).with_cell(policy),
            );
            let new = write_set(
                &set,
                options().with_sample_metadata(policy).with_cell(policy),
            );
            assert_eq!(outcome(new), outcome(old), "ensemble {index}, {policy:?}");
        }
    }
}

#[test]
fn sample_state_xyz_cannot_hold_is_rejected_or_omitted_with_a_report() {
    let omit = options().with_sample_metadata(UnrepresentablePolicy::Omit);
    let assert_set_rejects = |set: &SampleSet, field: TrajectoryField, frame: usize| {
        let mut writer =
            XyzWriter::new(Vec::new(), set.shared_topology(), options(), "e.xyz").unwrap();
        let error = writer.write_sample_set(set).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::UnsupportedField);
        assert_eq!((error.field(), error.frame()), (Some(field), Some(frame)));
        assert!(writer.writer().is_empty());
    };
    let plain = write(&waters(), options()).unwrap().0;

    // A time order.
    let mut timed = SampleSet::from_ensemble(waters());
    timed
        .set_order(Some(Order::Time(
            TimeAxis::new(Some(Quantity::new(vec![0.0, 1.0, 2.0], PICOSECOND)), None).unwrap(),
        )))
        .unwrap();
    assert_set_rejects(&timed, TrajectoryField::Time, 0);
    let (text, report) = write_set(&timed, omit.clone()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::Time, 3)]
    );

    // Velocities on one sample.
    let mut moving = SampleSet::from_ensemble(waters());
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
    assert_set_rejects(&moving, TrajectoryField::Velocities, 1);
    let (text, report) = write_set(&moving, omit.clone()).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::Velocities, 1)]
    );

    // A sample outside the distribution: XYZ reads every frame back as a
    // member.
    let mut subset = SampleSet::from_ensemble(waters());
    subset
        .set_distribution(Some(
            Distribution::new(Members::Subset(vec![0, 2]), EnsemblePartition::one_state(2))
                .unwrap(),
        ))
        .unwrap();
    assert_set_rejects(&subset, TrajectoryField::Membership, 1);
    let (text, report) = write_set(&subset, omit).unwrap();
    assert_eq!(text, plain);
    assert_eq!(
        report.omissions(),
        [BONDS, state_omission(StructureField::Membership, 1)]
    );

    // Without a distribution, the frames read back as one state of equally
    // weighted samples, as from an ensemble of one.
    let mut alone = SampleSet::from_ensemble(waters());
    alone.set_distribution(None).unwrap();
    assert_eq!(write_set(&alone, options()).unwrap().0, plain);
    let empty = SampleSet::new(alone.shared_topology());
    let mut writer =
        XyzWriter::new(Vec::new(), alone.shared_topology(), options(), "e.xyz").unwrap();
    let error = writer.write_sample_set(&empty).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
}

#[test]
fn path_writers_report_the_cells_their_frames_omit() {
    let ensemble = waters();
    let topology = ensemble.shared_topology();
    let frames = ensemble.iter().map(|sample| {
        TrajectoryFrame::new(sample.as_model_view().positions().clone()).with_cell(cell())
    });
    let trajectory = kekule::structure::Trajectory::from_items(topology, frames).unwrap();
    let path = std::env::temp_dir().join(format!("kekule-traj-{}-cells.xyz", std::process::id()));
    let saved = kekule::traj::io::write_trajectory(
        &path,
        &trajectory,
        kekule::traj::io::TrajectoryWriteOptions::default()
            .with_xyz_options(options().with_cell(UnrepresentablePolicy::Omit)),
    );
    let _ = std::fs::remove_file(&path);
    assert_eq!(
        saved.unwrap().omissions(),
        [BONDS, WriteOmission::Cell { structures: 3 }]
    );
}
