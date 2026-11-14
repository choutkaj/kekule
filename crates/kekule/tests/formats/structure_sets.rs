//! Readers that interpret structure sets: what each format states of the
//! order, that none states a distribution, and the sets they build.

use kekule::diagnostics::{
    DiagnosticKind, FormatErrorDetail, GeometryField, UnrepresentablePolicy, WriteOmission,
};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::mmcif::{self, MmcifWriteOptions};
use kekule::pdb::{self, PdbWriteOptions};
use kekule::sdf::{self, SdfDataField, SdfRecordInterpretation, SdfWriteOptions};
use kekule::structure::{
    Distribution, Ensemble, Geometry, Order, Partition, Place, Positions, Samples, StateData,
    Structure, StructureSet, TimeAxis, Velocities,
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
        .interpret_structure_set(Default::default())
        .unwrap()
        .into_structure_set();
    pdb::write(
        &set,
        PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit),
    )
    .unwrap()
    .0
}

#[test]
fn pdb_and_mmcif_structure_sets_are_their_models() {
    let text = pdb_waters();
    let model = |z: &str| {
        [
            ("1", "O1", "0.000", "0.000", "O"),
            ("2", "H2", "0.957", "0.000", "H"),
            ("3", "H3", "-0.240", "0.927", "H"),
        ]
        .map(|(serial, name, x, y, element)| {
            format!(
                "HETATM{serial:>5}  {name:<3} UNL A   1    {x:>8}{y:>8}{z:>8}                       \
                 {element}"
            )
        })
        .join("\n")
    };
    assert_eq!(
        text,
        format!(
            "MODEL        1\n{}\nENDMDL\nMODEL        2\n{}\nENDMDL\nEND\n",
            model("0.000"),
            model("0.100")
        )
    );
    let document = pdb::parse_str(&text, Default::default()).unwrap();
    let set = document
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(set.reports().len(), 2);
    assert_eq!(set.structure_set().len(), 2);
    assert_eq!(set.structure_set().order(), None);
    assert_waters(set.structure_set());
    // The interpretation writes back the same entry, losing nothing.
    let (written, report) = pdb::write(&set, PdbWriteOptions::default()).unwrap();
    assert_eq!(written, text);
    assert!(report.is_empty());

    let text = mmcif::write([&set], MmcifWriteOptions::default())
        .unwrap()
        .0;
    assert_eq!(text.matches("data_").count(), 1, "{text}");
    let rows = text
        .lines()
        .filter(|line| line.starts_with("HETATM"))
        .map(|line| {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            // The atom name, coordinates, and model number.
            [fields[3], fields[10], fields[11], fields[12], fields[20]].join(" ")
        })
        .collect::<Vec<_>>();
    assert_eq!(
        rows,
        [
            "O1 0.000 0.000 0.000 1",
            "H2 0.957 0.000 0.000 1",
            "H3 -0.240 0.927 0.000 1",
            "O1 0.000 0.000 0.100 2",
            "H2 0.957 0.000 0.100 2",
            "H3 -0.240 0.927 0.100 2",
        ]
    );
    let document = mmcif::parse_str(&text, Default::default()).unwrap();
    let set = document
        .interpret_structure_set(Default::default())
        .unwrap();
    assert_eq!(set.reports().len(), 2);
    assert_eq!(set.structure_set().order(), None);
    assert_waters(set.structure_set());
    let block = &document.blocks()[0];
    let from_block = block.interpret_structure_set(Default::default()).unwrap();
    assert_eq!(from_block.reports(), set.reports());
    assert!(from_block
        .topology()
        .same_layout(set.structure_set().topology()));
    for (left, right) in from_block.structure_set().iter().zip(set.structure_set()) {
        assert_eq!(left.payload(), right.payload());
    }
}

/// Asserts that `set` holds the two waters of [`WATERS`]: one water topology
/// and the second geometry 0.1 A above the first.
fn assert_waters(set: &StructureSet) {
    assert_eq!(set.len(), 2);
    let elements = set
        .topology()
        .atoms()
        .map(|atom| atom.element().symbol())
        .collect::<Vec<_>>();
    assert_eq!(elements, ["O", "H", "H"]);
    assert_eq!(set.topology().bond_count(), 2);
    let expected = [[0.0, 0.0, 0.0], [0.957, 0.0, 0.0], [-0.24, 0.927, 0.0]];
    for (geometry, z) in set.iter().zip([0.0, 0.1]) {
        for (index, [x, y, _]) in expected.into_iter().enumerate() {
            let point = geometry
                .positions()
                .get(index)
                .unwrap()
                .value_in(ANGSTROM)
                .unwrap();
            for (actual, expected) in [point.x, point.y, point.z].into_iter().zip([x, y, z]) {
                assert!((actual - expected).abs() <= 1e-9, "{actual} != {expected}");
            }
        }
    }
}

fn interpret_records(text: &str) -> Vec<SdfRecordInterpretation> {
    sdf::parse_str(text, Default::default())
        .unwrap()
        .records()
        .iter()
        .map(|record| record.interpret(Default::default()).unwrap())
        .collect()
}

/// Two methane conformers at 0 and 1 A along x as an ensemble of two states,
/// weighted 1:3, the second at -2.5 kJ/mol, at 300 K.
fn weighted_methanes() -> Ensemble {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let structures = [0.0, 1.0].map(|x| {
        Structure::new(
            topology.clone(),
            Positions::new(Quantity::new(
                [kekule::math::Point3::new(x, 0.0, 0.0)],
                ANGSTROM,
            ))
            .unwrap(),
        )
        .unwrap()
    });
    let set = StructureSet::from_structures(structures).unwrap();
    let distribution = Distribution::new(Samples::All, Partition::state_per_sample(2))
        .unwrap()
        .with_log_weights(vec![0.0, 3.0_f64.ln()])
        .unwrap()
        .with_states(vec![
            StateData::new(),
            StateData::new()
                .with_energy(Quantity::new(-2.5, KILOJOULE_PER_MOLE))
                .unwrap(),
        ])
        .unwrap()
        .with_temperature(Quantity::new(300.0, kekule::units::KELVIN))
        .unwrap();
    Ensemble::new(set, distribution).unwrap()
}

/// How a three-geometry methane set is ordered.
#[derive(Debug, Clone, Copy)]
enum Axis {
    /// By steps 0, 4, and 8.
    Step,
    /// By coordinates -60, -30, and 0 in a unit.
    Coordinate(Unit),
}

/// A three-geometry set of one methane carbon along `axis`.
fn methane_series(axis: Axis) -> StructureSet {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let mut set =
        StructureSet::from_geometries(topology, (0..3).map(|_| Geometry::new(Positions::zeros(1))))
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
fn sdf_structure_sets_are_the_orders_of_their_records() {
    // An ensemble's set writes its geometries alone: the records read back
    // as an unordered set of the same positions.
    let ensemble = weighted_methanes();
    let text = sdf::write(ensemble.set(), SdfWriteOptions::default())
        .unwrap()
        .0;
    assert!(!text.contains("kekule."), "{text}");
    let records = interpret_records(&text);
    assert!(records.iter().all(|record| record.data_fields().is_empty()));
    let read = sdf::structure_set(&records).unwrap();
    assert_eq!(read.order(), None);
    assert_eq!(read.len(), 2);
    for (left, right) in read.iter().zip(ensemble.set()) {
        assert_eq!(left.positions(), right.positions());
    }

    // Steps or coordinates make the order.
    for (axis, fields) in [
        (
            Axis::Step,
            vec![
                vec![(sdf::STEP_FIELD, "0")],
                vec![(sdf::STEP_FIELD, "4")],
                vec![(sdf::STEP_FIELD, "8")],
            ],
        ),
        (
            Axis::Coordinate(DEGREE),
            ["-60", "-30", "0"]
                .map(|value| {
                    vec![
                        (sdf::COORDINATE_FIELD, value),
                        (sdf::COORDINATE_UNIT_FIELD, "deg"),
                    ]
                })
                .to_vec(),
        ),
    ] {
        let series = methane_series(axis);
        let text = sdf::write(&series, SdfWriteOptions::default()).unwrap().0;
        let records = interpret_records(&text);
        let written = records
            .iter()
            .map(|record| {
                record
                    .data_fields()
                    .iter()
                    .map(|field| (field.name().to_owned(), field.value().to_owned()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let fields = fields
            .iter()
            .map(|record| {
                record
                    .iter()
                    .map(|(name, value)| ((*name).to_owned(), (*value).to_owned()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        assert_eq!(written, fields, "{axis:?}");
        let set = sdf::structure_set(&records).unwrap();
        assert_eq!(set.order(), series.order(), "{axis:?}");
        assert!(matches!(
            (axis, set.order()),
            (Axis::Step, Some(Order::Step(_)))
                | (Axis::Coordinate(_), Some(Order::Coordinate { .. }))
        ));
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
            [sdf::SdfRecordSource::structure(
                &series.get(0).unwrap().to_structure(),
            )],
            SdfWriteOptions::default(),
        )
        .unwrap()
        .0,
    );
    records.push(stepless[0].clone());
    let error = sdf::structure_set(&records).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.to_string(),
        "Cannot interpret the SDF input in record 4. The data field \"kekule.step\" is missing. Other records have it."
    );

    // Steps and coordinates together.
    let scan = methane_series(Axis::Coordinate(ANGSTROM));
    let mut mixed = interpret_records(&sdf::write(&scan, SdfWriteOptions::default()).unwrap().0);
    mixed.extend(interpret_records(&text));
    let error = sdf::structure_set(&mixed).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert!(matches!(
        error.detail(),
        Some(FormatErrorDetail::DataField { .. })
    ));
}

#[test]
fn sdf_records_of_different_topologies_do_not_form_a_structure_set() {
    let water = Structure::new(
        kekule::smiles::to_topology("O").unwrap(),
        Positions::new(Quantity::new([Point3::new(0.0, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let records = [
        methane_record(0.0, &[]),
        SdfRecordInterpretation::new("", water, Vec::new()),
    ];
    let error = sdf::structure_set(&records).unwrap_err();
    assert_eq!(error.kind(), DiagnosticKind::Inconsistent);
    assert_eq!(
        error.to_string(),
        "Cannot interpret the SDF input in record 2. The records do not form a structure set."
    );
    assert!(std::error::Error::source(&error).is_some());
}

/// One methane record at `x` with `fields`.
fn methane_record(x: f64, fields: &[(&str, &str)]) -> SdfRecordInterpretation {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let structure = Structure::new(
        topology,
        Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], ANGSTROM)).unwrap(),
    )
    .unwrap();
    let fields = fields
        .iter()
        .map(|(name, value)| SdfDataField::new(*name, *value))
        .collect();
    SdfRecordInterpretation::new("", structure, fields)
}

/// The set an SDF text of `set` reads back as.
fn sdf_round_trip(set: &StructureSet, options: SdfWriteOptions) -> StructureSet {
    sdf::structure_set(&interpret_records(&sdf::write(set, options).unwrap().0)).unwrap()
}

#[test]
fn sdf_structure_sets_state_only_what_the_records_state() {
    // Records without place fields give an unordered set.
    let bare = [methane_record(0.0, &[]), methane_record(1.0, &[])];
    let set = sdf::structure_set(&bare).unwrap();
    assert_eq!(set.order(), None);
    assert_eq!(set.len(), 2);
    let close = |(left, right): (f64, f64), expected: (f64, f64)| {
        (left - expected.0).abs() <= 1e-12 && (right - expected.1).abs() <= 1e-12
    };
    let x = |set: &StructureSet, geometry: usize| {
        set.get(geometry)
            .unwrap()
            .positions()
            .get(0)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            .x
    };
    assert!(close((x(&set, 0), x(&set, 1)), (0.0, 1.0)));
    let read = sdf_round_trip(&set, SdfWriteOptions::default());
    assert_eq!(read.order(), None);
    assert_eq!(read.len(), 2);
    assert!(close((x(&read, 0), x(&read, 1)), (0.0, 1.0)));

    // The fields an earlier writer gave a distribution are ordinary data: the
    // records form the same unordered set, whatever the fields hold.
    let old = [
        [
            methane_record(
                0.0,
                &[("kekule.state", "0"), ("kekule.representative", "true")],
            ),
            methane_record(1.0, &[("kekule.state", "0")]),
        ],
        [
            methane_record(0.0, &[("kekule.temperature.K", "300")]),
            methane_record(1.0, &[("kekule.temperature.K", "310")]),
        ],
        [
            methane_record(0.0, &[("kekule.weight", "2")]),
            methane_record(
                1.0,
                &[("kekule.weight", "heavy"), ("kekule.member", "true")],
            ),
        ],
    ];
    for records in &old {
        let set = sdf::structure_set(records).unwrap();
        assert_eq!(set.order(), None);
        assert!(close((x(&set, 0), x(&set, 1)), (0.0, 1.0)));
        assert!(records.iter().all(|record| record
            .data_fields()
            .iter()
            .all(|field| field.name().starts_with("kekule."))));
    }
}

/// A two-geometry methane set ordered by coordinates in `unit`.
fn coordinate_set(unit: Unit) -> StructureSet {
    let topology = kekule::smiles::to_topology("C").unwrap();
    let mut set = StructureSet::from_geometries(
        topology,
        [0.0, 1.0].map(|x| {
            Geometry::new(
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
                GeometryField::OrderCoordinate
            ))
        );
        assert!(error.to_string().contains("coordinate unit"), "{error}");
        let omit = SdfWriteOptions::default().with_geometry_metadata(UnrepresentablePolicy::Omit);
        let (text, report) = sdf::write(&set, omit).unwrap();
        assert!(!text.contains(sdf::COORDINATE_FIELD), "{text}");
        assert_eq!(
            report.omissions(),
            [WriteOmission::GeometryState {
                field: GeometryField::OrderCoordinate,
                geometries: set.len(),
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
        sdf::structure_set(&records).unwrap_err().to_string(),
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
        let error = sdf::structure_set(&records).unwrap_err();
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
fn xyz_structure_sets_have_no_order() {
    let document = kekule::xyz::parse_str(WATERS, Default::default()).unwrap();
    let set = document
        .interpret_structure_set(Default::default())
        .unwrap();
    let reports = set.reports();
    assert_eq!(
        reports
            .iter()
            .map(|report| (report.frame(), report.comment()))
            .collect::<Vec<_>>(),
        [(0, "first"), (1, "second")]
    );
    assert_eq!(set.structure_set().order(), None);
    assert_waters(set.structure_set());
    // An ensemble of one unweighted state over every geometry reads the
    // frames as equally likely draws.
    let ensemble = Ensemble::uniform(set.into_structure_set());
    let distribution = ensemble.distribution();
    assert_eq!(distribution.samples(), &Samples::All);
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.log_weights(), None);
    assert_eq!(distribution.probabilities().unwrap(), [0.5, 0.5]);
}

/// Four frames of one water, each with a cell, velocities, a time, and a
/// step.
fn water_frames() -> StructureSet {
    let structure = kekule::xyz::parse_str(
        "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n",
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_structure();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(20.0, 21.0, 22.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    let mut frames = StructureSet::new(structure.shared_topology());
    let axis = TimeAxis::new(
        Some(Quantity::new(Vec::new(), PICOSECOND)),
        Some(Vec::new()),
    )
    .unwrap();
    frames.set_order(Some(Order::Time(axis))).unwrap();
    for index in 0..4_u32 {
        let shift = f64::from(index) * 0.25;
        let points = structure
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| Point3::new(point.x + shift, point.y, point.z - shift))
            .collect::<Vec<_>>();
        let mut geometry =
            Geometry::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
                .with_cell(cell);
        geometry
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
        frames.push(geometry, Some(place)).unwrap();
    }
    frames
}

/// The frames `[3, 1, 3]` as an unordered set without velocities, built with
/// the calls that replace `Trajectory::to_ensemble`, as the changelog states
/// them.
fn ensemble_of_frames(frames: &StructureSet, atoms: Option<&AtomSelection>) -> StructureSet {
    let mut ensemble = StructureSet::from_reader(frames, [3, 1, 3]).unwrap();
    ensemble.set_order(None).unwrap();
    ensemble.clear_velocities_and_forces();
    if let Some(atoms) = atoms {
        ensemble = ensemble.subset(atoms).unwrap();
    }
    ensemble
}

/// The calls replacing `Trajectory::to_ensemble` write one model per sampled
/// frame, its positions shifted as the frame's.
#[test]
fn ensembles_of_trajectory_frames_write_one_model_per_sampled_frame() {
    let frames = water_frames();
    let atoms = AtomSelection::from_indices(
        frames.shared_topology(),
        [TopologyAtomIndex::new(0), TopologyAtomIndex::new(1)],
    )
    .unwrap();
    let omit = UnrepresentablePolicy::Omit;
    let hydrogens = |atoms| WriteOmission::ImplicitHydrogens { atoms };
    let cases = [
        (ensemble_of_frames(&frames, None), 3, vec![hydrogens(3)]),
        (
            ensemble_of_frames(&frames, Some(&atoms)),
            2,
            vec![hydrogens(2)],
        ),
    ];
    // Frame i is the water shifted by 0.25 i nm along x and back along z:
    // frames 3, 1, and 3 in angstroms, as the writers round them.
    let water = [("O1", 0.0, 0.0), ("H2", 0.957, 0.0), ("H3", -0.24, 0.927)];
    let shifts = [7.5, 2.5, 7.5];
    for (ensemble, atoms, omissions) in cases {
        assert_eq!(ensemble.len(), 3);
        assert_eq!(ensemble.order(), None);
        let pdb_options = PdbWriteOptions::default().with_implicit_hydrogens(omit);
        let mmcif_options = MmcifWriteOptions::default().with_implicit_hydrogens(omit);
        let (text, report) = pdb::write(&ensemble, pdb_options).unwrap();
        assert_eq!(report.omissions(), omissions);
        assert!(text.starts_with(
            "CRYST1   20.000   21.000   22.000  90.00  90.00  90.00 P 1           1\n"
        ));
        assert_eq!(text.matches("MODEL ").count(), 3, "{text}");
        let expected = shifts
            .iter()
            .flat_map(|shift| {
                water[..atoms]
                    .iter()
                    .enumerate()
                    .map(move |(serial, (name, x, y))| {
                        format!(
                            "HETATM{:>5}  {name:<3} UNL A   1    {:>8.3}{:>8.3}{:>8.3}",
                            serial + 1,
                            x + shift,
                            y,
                            -shift
                        )
                    })
            })
            .collect::<Vec<_>>();
        let atom_rows = text
            .lines()
            .filter(|line| line.starts_with("HETATM"))
            .map(|line| line[..54].to_owned())
            .collect::<Vec<_>>();
        assert_eq!(atom_rows, expected);
        let (text, report) = mmcif::write([&ensemble], mmcif_options).unwrap();
        assert_eq!(report.omissions(), omissions);
        let rows = text
            .lines()
            .filter(|line| line.starts_with("HETATM"))
            .map(|line| {
                let fields = line.split_whitespace().collect::<Vec<_>>();
                [fields[3], fields[10], fields[11], fields[12], fields[20]].join(" ")
            })
            .collect::<Vec<_>>();
        let expected = shifts
            .iter()
            .enumerate()
            .flat_map(|(model, shift)| {
                water[..atoms].iter().map(move |(name, x, y)| {
                    format!("{name} {:.3} {y:.3} {:.3} {}", x + shift, -shift, model + 1)
                })
            })
            .collect::<Vec<_>>();
        assert_eq!(rows, expected);
    }
}
