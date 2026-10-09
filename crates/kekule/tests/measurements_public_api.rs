use kekule::{
    geometry::{PeriodicCell, Point3, Vector3},
    smiles,
    structure::{
        measure::{self, ConnectivityCheck, MeasurementError},
        Model, Positions,
    },
    topology::{AtomSelection, InstanceAtomId, SelectionError},
    units::{Quantity, Unit, ANGSTROM, DEGREE, NANOMETER, PICOSECOND},
};

fn model(points: [Point3; 4], unit: Unit) -> Model {
    let molecule = smiles::to_molecules("CCCC").unwrap().pop().unwrap();
    Model::from_molecule(
        molecule.clone(),
        &Positions::new(Quantity::new(points, unit)).unwrap(),
    )
    .unwrap()
}

fn right_angle() -> [Point3; 4] {
    [
        Point3::new(1.0, 0.0, 0.0),
        Point3::origin(),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 1.0),
    ]
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}

#[test]
fn named_measurements_have_units_and_a_defined_dihedral_sign() {
    for (unit, scale) in [(ANGSTROM, 1.0), (NANOMETER, 0.1)] {
        let model = model(
            right_angle().map(|p| Point3::new(p.x * scale, p.y * scale, p.z * scale)),
            unit,
        );
        let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(model.topology().atom_ids()).unwrap();
        close(
            measure::distance(model.as_model_view(), a, b)
                .unwrap()
                .value_in(ANGSTROM)
                .unwrap(),
            1.0,
        );
        close(
            measure::angle(model.as_model_view(), a, b, c)
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            90.0,
        );
        close(
            measure::dihedral(model.as_model_view(), a, b, c, d)
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            -90.0,
        );
        close(
            measure::distance(model.as_model_view(), a, a)
                .unwrap()
                .into_value(),
            0.0,
        );
    }
    let reflected = model(right_angle().map(|p| Point3::new(p.x, p.y, -p.z)), ANGSTROM);
    let a = reflected.topology().atom_ids();
    close(
        measure::dihedral(reflected.as_model_view(), a[0], a[1], a[2], a[3])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
}

#[test]
fn consecutive_bonds_preserve_measurements_in_both_directions() {
    // Include single, double, and triple bonds: order does not affect adjacency.
    let molecule = smiles::to_molecules("C=CC#N").unwrap().pop().unwrap();
    let model = Model::from_molecule(
        molecule.clone(),
        &Positions::new(Quantity::new(right_angle(), ANGSTROM)).unwrap(),
    )
    .unwrap();
    let ids = <[InstanceAtomId; 4]>::try_from(model.topology().atom_ids()).unwrap();
    let mut reversed = ids;
    reversed.reverse();
    for [a, b, c, d] in [ids, reversed] {
        assert_eq!(
            measure::distance_with_connectivity(
                model.as_model_view(),
                a,
                b,
                ConnectivityCheck::ConsecutiveBonds
            ),
            measure::distance(model.as_model_view(), a, b)
        );
        assert_eq!(
            measure::angle_with_connectivity(
                model.as_model_view(),
                a,
                b,
                c,
                ConnectivityCheck::ConsecutiveBonds
            ),
            measure::angle(model.as_model_view(), a, b, c)
        );
        assert_eq!(
            measure::dihedral_with_connectivity(
                model.as_model_view(),
                a,
                b,
                c,
                d,
                ConnectivityCheck::ConsecutiveBonds
            ),
            measure::dihedral(model.as_model_view(), a, b, c, d)
        );
    }
}

#[test]
fn optional_connectivity_reports_each_missing_consecutive_pair() {
    let model = model(right_angle(), NANOMETER);
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(model.topology().atom_ids()).unwrap();
    assert_eq!(
        ConnectivityCheck::default(),
        ConnectivityCheck::Unrestricted
    );
    for [x, y] in [[a, c], [a, a]] {
        assert_eq!(
            measure::distance_with_connectivity(
                model.as_model_view(),
                x,
                y,
                ConnectivityCheck::ConsecutiveBonds
            ),
            Err(MeasurementError::MissingBond { a: x, b: y })
        );
        let unrestricted = measure::distance(model.as_model_view(), x, y).unwrap();
        assert_eq!(
            measure::distance_with_connectivity(
                model.as_model_view(),
                x,
                y,
                ConnectivityCheck::default()
            )
            .unwrap(),
            unrestricted
        );
    }
    for (atoms, missing) in [([a, c, d], [a, c]), ([a, b, d], [b, d])] {
        let [x, y, z] = atoms;
        assert_eq!(
            measure::angle_with_connectivity(
                model.as_model_view(),
                x,
                y,
                z,
                ConnectivityCheck::ConsecutiveBonds
            ),
            Err(MeasurementError::MissingBond {
                a: missing[0],
                b: missing[1]
            })
        );
        let unrestricted = measure::angle(model.as_model_view(), x, y, z).unwrap();
        assert_eq!(
            measure::angle_with_connectivity(
                model.as_model_view(),
                x,
                y,
                z,
                ConnectivityCheck::default()
            )
            .unwrap(),
            unrestricted
        );
    }
    for (atoms, missing) in [
        ([a, c, b, d], [a, c]), // Both first and last pairs missing: report first.
        ([b, a, d, c], [a, d]),
        ([a, b, c, a], [c, a]),
    ] {
        let [w, x, y, z] = atoms;
        assert_eq!(
            measure::dihedral_with_connectivity(
                model.as_model_view(),
                w,
                x,
                y,
                z,
                ConnectivityCheck::ConsecutiveBonds
            ),
            Err(MeasurementError::MissingBond {
                a: missing[0],
                b: missing[1]
            })
        );
        let unrestricted = measure::dihedral(model.as_model_view(), w, x, y, z).unwrap();
        assert_eq!(
            measure::dihedral_with_connectivity(
                model.as_model_view(),
                w,
                x,
                y,
                z,
                ConnectivityCheck::default()
            )
            .unwrap(),
            unrestricted
        );
    }
}

#[test]
fn reused_definitions_do_not_create_bonds_between_instances() {
    let molecule = smiles::to_molecules("CC").unwrap().pop().unwrap();
    let mut builder = Model::builder();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    let points = right_angle();
    builder
        .add_instance(
            definition,
            &Positions::new(Quantity::new([points[0], points[1]], NANOMETER)).unwrap(),
        )
        .unwrap();
    builder
        .add_instance(
            definition,
            &Positions::new(Quantity::new([points[2], points[3]], NANOMETER)).unwrap(),
        )
        .unwrap();
    let model = builder.build().unwrap();
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(model.topology().atom_ids()).unwrap();
    let missing = Err(MeasurementError::MissingBond { a: b, b: c });
    assert_eq!(
        measure::distance_with_connectivity(
            model.as_model_view(),
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds
        ),
        missing
    );
    assert_eq!(
        measure::angle_with_connectivity(
            model.as_model_view(),
            a,
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds
        ),
        missing
    );
    assert_eq!(
        measure::dihedral_with_connectivity(
            model.as_model_view(),
            a,
            b,
            c,
            d,
            ConnectivityCheck::ConsecutiveBonds
        ),
        missing
    );
    assert!(measure::distance_with_connectivity(
        model.as_model_view(),
        b,
        c,
        ConnectivityCheck::Unrestricted
    )
    .is_ok());
    assert!(measure::angle_with_connectivity(
        model.as_model_view(),
        a,
        b,
        c,
        ConnectivityCheck::Unrestricted
    )
    .is_ok());
    assert!(measure::dihedral_with_connectivity(
        model.as_model_view(),
        a,
        b,
        c,
        d,
        ConnectivityCheck::Unrestricted
    )
    .is_ok());
}

#[test]
fn connectivity_checks_reject_invalid_ids_before_missing_bonds() {
    let model = model(right_angle(), NANOMETER);
    let [a, b, c, _] = <[InstanceAtomId; 4]>::try_from(model.topology().atom_ids()).unwrap();
    for invalid in [
        InstanceAtomId::new(a.molecule(), kekule::core::AtomId::new(99)),
        InstanceAtomId::new(kekule::topology::MoleculeInstanceId::new(99), a.atom()),
    ] {
        let error = Err(MeasurementError::InvalidAtomId(invalid));
        assert_eq!(
            measure::distance_with_connectivity(
                model.as_model_view(),
                a,
                invalid,
                ConnectivityCheck::ConsecutiveBonds
            ),
            error
        );
        assert_eq!(
            measure::angle_with_connectivity(
                model.as_model_view(),
                a,
                c,
                invalid,
                ConnectivityCheck::ConsecutiveBonds
            ),
            error
        );
        assert_eq!(
            measure::dihedral_with_connectivity(
                model.as_model_view(),
                a,
                c,
                b,
                invalid,
                ConnectivityCheck::ConsecutiveBonds
            ),
            error
        );
    }
}

#[test]
fn bonded_geometry_is_still_checked_independently() {
    let line = model(
        [
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
        NANOMETER,
    );
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(line.topology().atom_ids()).unwrap();
    close(
        measure::angle_with_connectivity(
            line.as_model_view(),
            a,
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        180.0,
    );
    assert_eq!(
        measure::dihedral_with_connectivity(
            line.as_model_view(),
            a,
            b,
            c,
            d,
            ConnectivityCheck::ConsecutiveBonds
        ),
        Err(MeasurementError::DegenerateGeometry)
    );
    // Nonconsecutive repeats are allowed when the geometry remains defined.
    close(
        measure::angle_with_connectivity(
            line.as_model_view(),
            a,
            b,
            a,
            ConnectivityCheck::ConsecutiveBonds,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        0.0,
    );
    let coincident = model([Point3::origin(); 4], NANOMETER);
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(coincident.topology().atom_ids()).unwrap();
    assert_eq!(
        measure::distance_with_connectivity(
            coincident.as_model_view(),
            a,
            b,
            ConnectivityCheck::ConsecutiveBonds
        )
        .unwrap()
        .into_value(),
        0.0
    );
    assert_eq!(
        measure::angle_with_connectivity(
            coincident.as_model_view(),
            a,
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds
        ),
        Err(MeasurementError::DegenerateGeometry)
    );
    assert_eq!(
        measure::dihedral_with_connectivity(
            coincident.as_model_view(),
            a,
            b,
            c,
            d,
            ConnectivityCheck::ConsecutiveBonds
        ),
        Err(MeasurementError::DegenerateGeometry)
    );
}

#[test]
fn tiny_finite_geometry_normalizes_without_reciprocal_overflow() {
    let tiny = model(
        right_angle().map(|p| Point3::new(p.x * 1e-310, p.y * 1e-310, p.z * 1e-310)),
        NANOMETER,
    );
    let a = tiny.topology().atom_ids();
    close(
        measure::angle(tiny.as_model_view(), a[0], a[1], a[2])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
    close(
        measure::dihedral(tiny.as_model_view(), a[0], a[1], a[2], a[3])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        -90.0,
    );
}

#[test]
fn undefined_geometry_and_numerical_overflow_reject_without_nan_results() {
    let line = model(
        [
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
        NANOMETER,
    );
    let a = line.topology().atom_ids();
    close(
        measure::angle(line.as_model_view(), a[0], a[1], a[2])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        180.0,
    );
    assert_eq!(
        measure::angle(line.as_model_view(), a[0], a[0], a[1]),
        Err(MeasurementError::DegenerateGeometry)
    );
    assert_eq!(
        measure::dihedral(line.as_model_view(), a[0], a[1], a[2], a[3]),
        Err(MeasurementError::DegenerateGeometry)
    );
    let invalid = InstanceAtomId::new(a[0].molecule(), kekule::core::AtomId::new(99));
    assert!(matches!(
        measure::distance(line.as_model_view(), invalid, a[0]),
        Err(MeasurementError::InvalidAtomId(_))
    ));
    let huge = model(
        [
            Point3::new(f64::MAX, 0.0, 0.0),
            Point3::new(-f64::MAX, 0.0, 0.0),
            Point3::origin(),
            Point3::origin(),
        ],
        NANOMETER,
    );
    let a = huge.topology().atom_ids();
    assert_eq!(
        measure::distance(huge.as_model_view(), a[0], a[1]),
        Err(MeasurementError::NumericalFailure)
    );
    let large = model(
        [
            Point3::new(1e200, 0.0, 0.0),
            Point3::origin(),
            Point3::new(0.0, 1e200, 0.0),
            Point3::new(0.0, 1e200, 1e200),
        ],
        NANOMETER,
    );
    let a = large.topology().atom_ids();
    close(
        measure::angle(large.as_model_view(), a[0], a[1], a[2])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
}

#[test]
fn spatial_selection_is_inclusive_unit_aware_and_specific_to_each_view() {
    let mut model = model(right_angle(), NANOMETER);
    let top = model.shared_topology();
    let ids = top.atom_ids();
    let all = AtomSelection::all(&top);
    let reference = AtomSelection::from_atoms(&top, [ids[1]]).unwrap();
    // Test exact inclusion in canonical units, separately from conversion:
    // 10 * the floating-point Angstrom scale can round just below 1 nm.
    let selected = measure::within(
        model.as_model_view(),
        &all,
        &reference,
        Quantity::new(1.0, NANOMETER),
    )
    .unwrap();
    assert_eq!(
        selected.atom_ids().collect::<Vec<_>>(),
        [ids[0], ids[1], ids[2]]
    );
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &all,
            &reference,
            Quantity::new(11.0, ANGSTROM)
        )
        .unwrap(),
        selected
    );
    let candidates = all.difference(&reference).unwrap();
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &candidates,
            &reference,
            Quantity::new(0.0, NANOMETER)
        )
        .unwrap()
        .indices()
        .len(),
        0
    );
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &all,
            &reference,
            Quantity::new(0.0, NANOMETER)
        )
        .unwrap(),
        reference
    );
    let empty = AtomSelection::from_atoms(&top, []).unwrap();
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &all,
            &empty,
            Quantity::new(1.0, NANOMETER)
        )
        .unwrap(),
        empty
    );
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &empty,
            &all,
            Quantity::new(1.0, NANOMETER)
        )
        .unwrap(),
        empty
    );
    for cutoff in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            measure::within(
                model.as_model_view(),
                &empty,
                &reference,
                Quantity::new(cutoff, NANOMETER)
            ),
            Err(MeasurementError::InvalidCutoff)
        );
    }
    assert!(matches!(
        measure::within(
            model.as_model_view(),
            &all,
            &reference,
            Quantity::new(1.0, PICOSECOND)
        ),
        Err(MeasurementError::Unit(_))
    ));
    // A perceived snapshot shares the layout, so its selections stay usable.
    let perceived = AtomSelection::all(&std::sync::Arc::new(top.perceived().unwrap()));
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &perceived,
            &reference,
            Quantity::new(1.0, NANOMETER)
        ),
        Ok(selected.clone())
    );
    // An independently published equal topology does not.
    let foreign = AtomSelection::all(&smiles::to_topology("CCCC").unwrap());
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &foreign,
            &reference,
            Quantity::new(1.0, NANOMETER)
        ),
        Err(MeasurementError::Selection(
            SelectionError::TopologyMismatch
        ))
    );
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &all,
            &foreign,
            Quantity::new(1.0, NANOMETER)
        ),
        Err(MeasurementError::Selection(
            SelectionError::TopologyMismatch
        ))
    );
    model
        .set_position(
            ids[0],
            Quantity::new(Point3::new(10.0, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    assert_eq!(selected.atom_ids().len(), 3); // The earlier atom set stays fixed.
    assert_eq!(
        measure::within(
            model.as_model_view(),
            &all,
            &reference,
            Quantity::new(1.0, NANOMETER)
        )
        .unwrap()
        .atom_ids()
        .collect::<Vec<_>>(),
        [ids[1], ids[2]]
    );
    model.conformation_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    close(
        measure::distance(model.as_model_view(), ids[0], ids[1])
            .unwrap()
            .value_in(NANOMETER)
            .unwrap(),
        10.0,
    );
}
