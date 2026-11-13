//! Readers that interpret sample sets: what each format states of the
//! distribution and order, and the sets they build.

use kekule::diagnostics::{
    DiagnosticKind, FormatErrorDetail, StructureField, UnrepresentablePolicy, WriteOmission,
};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::mmcif::{self, MmcifWriteOptions};
use kekule::pdb::{self, PdbWriteOptions};
use kekule::sdf::{self, SdfDataField, SdfRecordInterpretation, SdfWriteOptions};
use kekule::structure::{
    Distribution, Members, Model, Order, Partition, Place, Positions, Sample, SampleSet, TimeAxis,
    Velocities,
};
use kekule::topology::{AtomSelection, TopologyAtomIndex};
use kekule::units::{
    Quantity, Unit, ANGSTROM, CANONICAL_VELOCITY_UNIT, DEGREE, KILOJOULE_PER_MOLE, NANOMETER,
    PICOSECOND,
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
";

/// Two waters as a two-model PDB entry.
fn pdb_waters() -> String {
    let set = kekule::xyz::parse_str(WATERS, Default::default())
        .unwrap()
        .interpret_sample_set(Default::default())
        .unwrap()
        .into_sample_set();
    pdb::write(
        &set,
        PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit),
    )
    .unwrap()
    .0
}

#[test]
fn pdb_and_mmcif_sample_sets_are_one_unweighted_state() {
    let text = pdb_waters();
    assert_eq!(crate::fingerprint(&text), 0xb1bb_deac_c50e_359a);
    let document = pdb::parse_str(&text, Default::default()).unwrap();
    let set = document.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(set.reports().len(), 2);
    assert_eq!(set.sample_set().len(), 2);
    assert_eq!(set.sample_set().order(), None);
    let distribution = set.sample_set().distribution().unwrap();
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.weights(), None);
    assert_eq!(crate::native_fingerprint(&set), 0x2aba_a9f4_5933_d9f3);
    assert_eq!(
        crate::fingerprint(format!(
            "{:?}",
            pdb::write(&set, PdbWriteOptions::default()).unwrap()
        )),
        0xaaa2_2f01_7606_146c
    );

    let text = mmcif::write([&set], MmcifWriteOptions::default())
        .unwrap()
        .0;
    assert_eq!(crate::fingerprint(&text), 0x3a95_627c_2c78_254c);
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let set = document.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(crate::native_fingerprint(&set), 0xfb08_1261_f062_f25a);
    let block = &document.blocks()[0];
    let from_block = block.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(from_block.reports(), set.reports());
    assert_eq!(
        crate::native_fingerprint(&from_block),
        0xfb08_1261_f062_f25a
    );
}

fn interpret_records(text: &str) -> Vec<SdfRecordInterpretation> {
    sdf::parse_str(text, Default::default())
        .unwrap()
        .records()
        .iter()
        .map(|record| record.interpret(Default::default()).unwrap())
        .collect()
}

/// Two methane conformers in two states, weighted 1:3, the second at
/// -2.5 kJ/mol.
fn weighted_methanes() -> SampleSet {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let models = [0.0, 1.0].map(|x| {
        Model::new(
            topology.clone(),
            Positions::new(Quantity::new(
                [kekule::math::Point3::new(x, 0.0, 0.0)],
                ANGSTROM,
            ))
            .unwrap(),
        )
        .unwrap()
    });
    let mut set = SampleSet::from_models(models).unwrap();
    set.set_distribution(Some(
        Distribution::new(Members::All, Partition::state_per_sample(2))
            .unwrap()
            .with_weights(vec![1.0, 3.0])
            .unwrap(),
    ))
    .unwrap();
    set.distribution_mut()
        .unwrap()
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(-2.5, KILOJOULE_PER_MOLE)))
        .unwrap();
    set
}

/// How a three-sample methane set is ordered.
#[derive(Debug, Clone, Copy)]
enum Axis {
    /// By steps 0, 4, and 8.
    Step,
    /// By coordinates -60, -30, and 0 in a unit.
    Coordinate(Unit),
}

/// A three-sample set of one methane carbon along `axis`.
fn methane_series(axis: Axis) -> SampleSet {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let mut set =
        SampleSet::from_samples(topology, (0..3).map(|_| Sample::new(Positions::zeros(1))))
            .unwrap();
    let order = match axis {
        Axis::Step => Order::Step((0..3).map(|index| index * 4).collect()),
        Axis::Coordinate(unit) => Order::Coordinate {
            unit,
            values: (0..3).map(|index| f64::from(index) * 30.0 - 60.0).collect(),
        },
    };
    set.set_order(Some(order)).unwrap();
    set
}

#[test]
fn sdf_sample_sets_are_the_distributions_and_orders_of_their_records() {
    // Weights, states and energies make the distribution.
    let set = weighted_methanes();
    assert_eq!(crate::native_fingerprint(&set), 0x8fe6_53ce_8dea_5590);
    let text = sdf::write(&set, SdfWriteOptions::default()).unwrap().0;
    assert_eq!(crate::fingerprint(&text), 0xa96d_0167_e29e_65c6);
    let read = sdf::sample_set(&interpret_records(&text)).unwrap();
    assert_eq!(read.order(), None);
    assert_eq!(read.distribution().unwrap().state_count(), 2);
    assert_eq!(read.distribution(), set.distribution());
    assert_eq!(crate::native_fingerprint(&read), 0xf4be_d5a7_57e0_1dcf);

    // Steps or coordinates make the order.
    for (axis, native, written, reread) in [
        (
            Axis::Step,
            0x8346_2618_fe0b_dd41,
            0x27fd_0aed_408e_6608,
            0x97db_3505_645b_a4f1,
        ),
        (
            Axis::Coordinate(DEGREE),
            0x7be2_e975_7d64_a7b5,
            0xd40b_391d_7512_5e32,
            0xd492_afa0_6f9d_cd74,
        ),
    ] {
        let series = methane_series(axis);
        assert_eq!(crate::native_fingerprint(&series), native, "{axis:?}");
        let text = sdf::write(&series, SdfWriteOptions::default()).unwrap().0;
        assert_eq!(crate::fingerprint(&text), written, "{axis:?}");
        let set = sdf::sample_set(&interpret_records(&text)).unwrap();
        assert_eq!(crate::native_fingerprint(&set), reread, "{axis:?}");
        assert_eq!(set.order(), series.order(), "{axis:?}");
        assert!(matches!(
            (axis, set.order()),
            (Axis::Step, Some(Order::Step(_)))
                | (Axis::Coordinate(_), Some(Order::Coordinate { .. }))
        ));
        // The records have no distribution field, so the set has no
        // distribution.
        assert_eq!(set.distribution(), None);
        for (left, right) in set.iter().zip(series.iter()) {
            assert_eq!(left.payload(), right.payload());
        }
    }
}

#[test]
fn sdf_records_that_do_not_form_an_order_are_refused() {
    let series = methane_series(Axis::Step);
    let text = sdf::write(&series, SdfWriteOptions::default()).unwrap().0;
    let mut records = interpret_records(&text);
    // A record without the step its neighbours have.
    let stepless = interpret_records(
        &sdf::write(
            [sdf::SdfRecordSource::model(
                &series.get(0).unwrap().to_model(),
            )],
            SdfWriteOptions::default(),
        )
        .unwrap()
        .0,
    );
    records.push(stepless[0].clone());
    let error = sdf::sample_set(&records).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.to_string(),
        "Cannot interpret the SDF input in record 4. The data field \"kekule.step\" is missing. Other records have it."
    );

    // Steps and coordinates together.
    let scan = methane_series(Axis::Coordinate(ANGSTROM));
    let mut mixed = interpret_records(&sdf::write(&scan, SdfWriteOptions::default()).unwrap().0);
    mixed.extend(interpret_records(&text));
    let error = sdf::sample_set(&mixed).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert!(matches!(
        error.detail(),
        Some(FormatErrorDetail::DataField { .. })
    ));
}

#[test]
fn sdf_records_of_different_topologies_do_not_form_a_sample_set() {
    let water = Model::new(
        kekule::smiles::to_topology("O").unwrap(),
        Positions::new(Quantity::new([Point3::new(0.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let records = [
        methane_record(0.0, &[]),
        SdfRecordInterpretation::new("", water, Vec::new()),
    ];
    let error = sdf::sample_set(&records).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.to_string(),
        "Cannot interpret the SDF input in record 2. The records do not form a sample set."
    );
    assert!(std::error::Error::source(&error).is_some());
}

/// One methane record at `x` with `fields`.
fn methane_record(x: f64, fields: &[(&str, &str)]) -> SdfRecordInterpretation {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let model = Model::new(
        topology,
        Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let fields = fields
        .iter()
        .map(|(name, value)| SdfDataField::new(*name, *value))
        .collect();
    SdfRecordInterpretation::new("", model, fields)
}

/// The set an SDF text of `set` reads back as.
fn sdf_round_trip(set: &SampleSet, options: SdfWriteOptions) -> SampleSet {
    sdf::sample_set(&interpret_records(&sdf::write(set, options).unwrap().0)).unwrap()
}

#[test]
fn sdf_sample_sets_state_only_what_the_records_state() {
    // Records without distribution fields give no distribution.
    let bare = [methane_record(0.0, &[]), methane_record(1.0, &[])];
    let set = sdf::sample_set(&bare).unwrap();
    assert_eq!(set.distribution(), None);
    assert_eq!(set.order(), None);
    assert_eq!(crate::native_fingerprint(&set), 0xdc67_2d95_20d3_270a);
    // A set without a distribution survives a round trip.
    let read = sdf_round_trip(&set, SdfWriteOptions::default());
    assert_eq!(read.distribution(), None);
    assert_eq!(read.len(), 2);

    // States without weights give a distribution without weights.
    let states = [
        methane_record(
            0.0,
            &[(sdf::STATE_FIELD, "0"), (sdf::REPRESENTATIVE_FIELD, "true")],
        ),
        methane_record(1.0, &[(sdf::STATE_FIELD, "0")]),
    ];
    let set = sdf::sample_set(&states).unwrap();
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.weights(), None);
    // One field on one record is enough for a distribution.
    let warm = [
        methane_record(0.0, &[(sdf::TEMPERATURE_FIELD, "300")]),
        methane_record(1.0, &[(sdf::TEMPERATURE_FIELD, "300")]),
    ];
    let set = sdf::sample_set(&warm).unwrap();
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.state_count(), 2);
    assert_eq!(distribution.weights(), None);
    assert!(distribution.temperature().is_some());
    // A weight on one record stores every weight.
    let weighted = [
        methane_record(0.0, &[(sdf::WEIGHT_FIELD, "2")]),
        methane_record(1.0, &[]),
    ];
    let set = sdf::sample_set(&weighted).unwrap();
    assert_eq!(
        set.distribution().unwrap().weights(),
        Some([2.0, 1.0].as_slice())
    );

    // A set with a distribution and no order survives a round trip.
    let set = weighted_methanes();
    let read = sdf_round_trip(&set, SdfWriteOptions::default());
    assert_eq!(crate::native_fingerprint(&read), 0xf4be_d5a7_57e0_1dcf);
    assert_eq!(read.order(), None);
    assert_eq!(read.distribution(), set.distribution());
}

/// A two-sample methane set ordered by coordinates in `unit`.
fn coordinate_set(unit: Unit) -> SampleSet {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let mut set = SampleSet::from_samples(
        topology,
        [0.0, 1.0].map(|x| {
            Sample::new(
                Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap(),
            )
        }),
    )
    .unwrap();
    set.set_order(Some(Order::Coordinate {
        unit,
        values: vec![-1.5, 0.25],
    }))
    .unwrap();
    set
}

#[test]
fn sdf_coordinate_units_are_the_unit_constants() {
    use kekule::units::*;
    let named = [
        DIMENSIONLESS,
        METER,
        NANOMETER,
        ANGSTROM,
        BOHR,
        CUBIC_BOHR,
        PER_NANOMETER,
        SQUARE_ANGSTROM,
        CUBIC_NANOMETER,
        SQUARE_NANOMETER,
        KILOGRAM,
        DALTON,
        SECOND,
        PICOSECOND,
        FEMTOSECOND,
        KELVIN,
        MOLAR,
        MOLE,
        COULOMB,
        ELEMENTARY_CHARGE,
        RADIAN,
        DEGREE,
        JOULE,
        KILOJOULE,
        KILOCALORIE,
        KILOJOULE_PER_MOLE,
        KILOCALORIE_PER_MOLE,
        KILOJOULE_PER_MOLE_PER_ANGSTROM,
        KILOJOULE_PER_MOLE_PER_NANOMETER,
        KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM,
        KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER,
        NANOMETER_PER_PICOSECOND,
    ];
    for unit in named {
        let set = coordinate_set(unit);
        let read = sdf_round_trip(&set, SdfWriteOptions::default());
        assert_eq!(read.order(), set.order(), "{unit}");
        let Some(Order::Coordinate {
            unit: read_unit, ..
        }) = read.order()
        else {
            panic!("{unit} gives a coordinate order");
        };
        assert_eq!(read_unit.to_string(), unit.to_string());
    }

    // A unit that is not a constant, by scale or by symbol, has no field.
    for unit in [
        Unit::new(Dimension::LENGTH, 2.0e-10, None).unwrap(),
        Unit::new(Dimension::LENGTH, 2.0e-10, Some("A")).unwrap(),
    ] {
        let set = coordinate_set(unit);
        let error = sdf::write(&set, SdfWriteOptions::default()).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Unsupported, "{error}");
        assert_eq!(
            error.detail(),
            Some(&FormatErrorDetail::UnrepresentableState(
                StructureField::OrderCoordinate
            ))
        );
        assert!(error.to_string().contains("coordinate unit"), "{error}");
        let omit = SdfWriteOptions::default().with_sample_metadata(UnrepresentablePolicy::Omit);
        let (text, report) = sdf::write(&set, omit).unwrap();
        assert!(!text.contains(sdf::COORDINATE_FIELD), "{text}");
        assert_eq!(
            report.omissions(),
            [WriteOmission::StructureState {
                field: StructureField::OrderCoordinate,
                structures: set.len(),
            }]
        );
        assert_eq!(sdf_round_trip(&set, omit).order(), None);
    }
}

#[test]
fn sdf_coordinate_units_must_match_the_first_record() {
    let records = [
        methane_record(
            0.0,
            &[
                (sdf::COORDINATE_FIELD, "1"),
                (sdf::COORDINATE_UNIT_FIELD, "A"),
            ],
        ),
        methane_record(
            1.0,
            &[
                (sdf::COORDINATE_FIELD, "2"),
                (sdf::COORDINATE_UNIT_FIELD, "nm"),
            ],
        ),
    ];
    assert_eq!(
        sdf::sample_set(&records).unwrap_err().to_string(),
        "Cannot interpret the SDF input in record 2. The data field \"kekule.coordinate.unit\" \
         differs from that of record 1."
    );
}

#[test]
fn sdf_place_fields_appear_once_per_record() {
    for (name, value) in [
        (sdf::STEP_FIELD, "4"),
        (sdf::COORDINATE_FIELD, "1.5"),
        (sdf::COORDINATE_UNIT_FIELD, "A"),
    ] {
        let mut fields = vec![
            (sdf::COORDINATE_FIELD, "1.5"),
            (sdf::COORDINATE_UNIT_FIELD, "A"),
        ];
        if name == sdf::STEP_FIELD {
            fields = vec![(sdf::STEP_FIELD, "4")];
        }
        let first = methane_record(0.0, &fields);
        fields.push((name, value));
        let records = [first, methane_record(1.0, &fields)];
        let error = sdf::sample_set(&records).unwrap_err();
        assert_eq!(error.kind(), DiagnosticKind::Inconsistent, "{error}");
        assert_eq!(error.record(), Some(2), "{error}");
        assert!(error.to_string().contains("appears twice."), "{error}");
        assert_eq!(
            error.detail(),
            Some(&FormatErrorDetail::DataField {
                name: name.to_owned()
            })
        );
    }
}

#[test]
fn xyz_sample_sets_have_no_distribution() {
    let document = kekule::xyz::parse_str(WATERS, Default::default()).unwrap();
    let set = document.interpret_sample_set(Default::default()).unwrap();
    assert_eq!(
        crate::fingerprint(format!("{:?}", set.reports())),
        0xa502_1322_7468_2c3c
    );
    assert_eq!(set.sample_set().order(), None);
    assert_eq!(set.sample_set().distribution(), None);
    assert_eq!(
        crate::native_fingerprint(set.sample_set()),
        0xcce8_08a1_364d_2b57
    );
    // With a distribution of one unweighted state it is what readers of
    // formats with states give.
    let mut set = set.into_sample_set();
    set.set_distribution(Some(
        Distribution::new(Members::All, Partition::one_state(2)).unwrap(),
    ))
    .unwrap();
    assert_eq!(crate::native_fingerprint(&set), 0x0f08_570b_3fb1_51ac);
}

/// Four frames of one water, each with a cell, velocities, a time, and a
/// step.
fn water_frames() -> SampleSet {
    let model = kekule::xyz::parse_str(
        "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n",
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_model();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(20.0, 21.0, 22.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    let mut frames = SampleSet::new(model.shared_topology());
    let axis = TimeAxis::new(
        Some(Quantity::new(Vec::new(), PICOSECOND)),
        Some(Vec::new()),
    )
    .unwrap();
    frames.set_order(Some(Order::Time(axis))).unwrap();
    for index in 0..4_u32 {
        let shift = f64::from(index) * 0.25;
        let points = model
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| Point3::new(point.x + shift, point.y, point.z - shift))
            .collect::<Vec<_>>();
        let mut sample =
            Sample::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
                .with_cell(cell);
        sample
            .set_velocities(Some(
                Velocities::new(Quantity::new(
                    [Vector3::new(1.0, 2.0, 3.0); 3],
                    CANONICAL_VELOCITY_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        let place = Place::Time {
            time: Some(Quantity::new(f64::from(index), PICOSECOND)),
            step: Some(u64::from(index) * 10),
        };
        frames.push(sample, Some(place)).unwrap();
    }
    frames
}

/// An ensemble of the frames `[3, 1, 3]` built with the calls that replace
/// `Trajectory::to_ensemble`, as the changelog states them.
fn ensemble_of_frames(
    frames: &SampleSet,
    atoms: Option<&AtomSelection>,
    weights: Option<Vec<f64>>,
) -> SampleSet {
    let mut ensemble = SampleSet::from_reader(frames, [3, 1, 3]).unwrap();
    ensemble.set_order(None).unwrap();
    ensemble.clear_velocities_and_forces();
    if let Some(atoms) = atoms {
        ensemble = ensemble.subset(atoms).unwrap();
    }
    let mut distribution = Distribution::new(Members::All, Partition::one_state(3)).unwrap();
    if let Some(weights) = weights {
        distribution = distribution.with_weights(weights).unwrap();
    }
    ensemble.set_distribution(Some(distribution)).unwrap();
    ensemble
}

/// The calls replacing `Trajectory::to_ensemble` write the PDB and mmCIF
/// bytes it wrote; the fingerprints were recorded from its output.
#[test]
fn ensembles_of_trajectory_frames_write_what_trajectory_sampling_wrote() {
    let frames = water_frames();
    let atoms = AtomSelection::from_indices(
        frames.shared_topology(),
        [TopologyAtomIndex::new(0), TopologyAtomIndex::new(1)],
    )
    .unwrap();
    let omit = UnrepresentablePolicy::Omit;
    let hydrogens = |atoms| WriteOmission::ImplicitHydrogens { atoms };
    let weights = WriteOmission::StructureState {
        field: StructureField::Weight,
        structures: 3,
    };
    let cases = [
        (
            ensemble_of_frames(&frames, None, None),
            false,
            vec![hydrogens(3)],
            0x7856_ea1c_9be2_f84b,
            0xdc86_05a8_6729_3231,
        ),
        (
            ensemble_of_frames(&frames, Some(&atoms), None),
            false,
            vec![hydrogens(2)],
            0xa036_6f06_4366_5371,
            0x3c9d_6a15_f7ef_9fb6,
        ),
        (
            ensemble_of_frames(&frames, Some(&atoms), Some(vec![2.0, 1.0, 0.5])),
            true,
            vec![weights, hydrogens(2)],
            0xa036_6f06_4366_5371,
            0x3c9d_6a15_f7ef_9fb6,
        ),
    ];
    for (ensemble, weighted, omissions, pdb_bytes, mmcif_bytes) in cases {
        let mut pdb_options = PdbWriteOptions::default().with_implicit_hydrogens(omit);
        let mut mmcif_options = MmcifWriteOptions::default().with_implicit_hydrogens(omit);
        if weighted {
            pdb_options = pdb_options.with_sample_metadata(omit);
            mmcif_options = mmcif_options.with_sample_metadata(omit);
        }
        let (text, report) = pdb::write(&ensemble, pdb_options).unwrap();
        assert_eq!(crate::fingerprint(&text), pdb_bytes, "{text}");
        assert_eq!(report.omissions(), omissions);
        let (text, report) = mmcif::write([&ensemble], mmcif_options).unwrap();
        assert_eq!(crate::fingerprint(&text), mmcif_bytes, "{text}");
        assert_eq!(report.omissions(), omissions);
    }
}
