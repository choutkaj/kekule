use kekule::{
    geometry::{PeriodicCell, Point3, Vector3},
    smiles,
    structure::{
        measure::{self, MeasurementError},
        Model, Positions,
    },
    topology::{AtomSelection, InstanceAtomId, SelectionError},
    units::{Quantity, Unit, ANGSTROM, DEGREE, NANOMETER, PICOSECOND},
};

fn model(points: [Point3; 4], unit: Unit) -> Model {
    let molecule = smiles::to_molecules("CCCC").unwrap().pop().unwrap();
    Model::from_molecule(
        &molecule,
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
        let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(model.atom_ids()).unwrap();
        close(
            measure::distance(model.view(), a, b)
                .unwrap()
                .value_in(ANGSTROM)
                .unwrap(),
            1.0,
        );
        close(
            measure::angle(model.view(), a, b, c)
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            90.0,
        );
        close(
            measure::dihedral(model.view(), a, b, c, d)
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            -90.0,
        );
        close(
            measure::distance(model.view(), a, a).unwrap().into_value(),
            0.0,
        );
    }
    let reflected = model(right_angle().map(|p| Point3::new(p.x, p.y, -p.z)), ANGSTROM);
    let a = reflected.atom_ids();
    close(
        measure::dihedral(reflected.view(), a[0], a[1], a[2], a[3])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
}

#[test]
fn tiny_finite_geometry_normalizes_without_reciprocal_overflow() {
    let tiny = model(
        right_angle().map(|p| Point3::new(p.x * 1e-310, p.y * 1e-310, p.z * 1e-310)),
        NANOMETER,
    );
    let a = tiny.atom_ids();
    close(
        measure::angle(tiny.view(), a[0], a[1], a[2])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
    close(
        measure::dihedral(tiny.view(), a[0], a[1], a[2], a[3])
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
    let a = line.atom_ids();
    close(
        measure::angle(line.view(), a[0], a[1], a[2])
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        180.0,
    );
    assert_eq!(
        measure::angle(line.view(), a[0], a[0], a[1]),
        Err(MeasurementError::DegenerateGeometry)
    );
    assert_eq!(
        measure::dihedral(line.view(), a[0], a[1], a[2], a[3]),
        Err(MeasurementError::DegenerateGeometry)
    );
    let invalid = InstanceAtomId::new(a[0].molecule(), kekule::core::AtomId::new(99));
    assert!(matches!(
        measure::distance(line.view(), invalid, a[0]),
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
    let a = huge.atom_ids();
    assert_eq!(
        measure::distance(huge.view(), a[0], a[1]),
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
    let a = large.atom_ids();
    close(
        measure::angle(large.view(), a[0], a[1], a[2])
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
        model.view(),
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
            model.view(),
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
            model.view(),
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
            model.view(),
            &all,
            &reference,
            Quantity::new(0.0, NANOMETER)
        )
        .unwrap(),
        reference
    );
    let empty = AtomSelection::from_atoms(&top, []).unwrap();
    assert_eq!(
        measure::within(model.view(), &all, &empty, Quantity::new(1.0, NANOMETER)).unwrap(),
        empty
    );
    assert_eq!(
        measure::within(model.view(), &empty, &all, Quantity::new(1.0, NANOMETER)).unwrap(),
        empty
    );
    for cutoff in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            measure::within(
                model.view(),
                &empty,
                &reference,
                Quantity::new(cutoff, NANOMETER)
            ),
            Err(MeasurementError::InvalidCutoff)
        );
    }
    assert!(matches!(
        measure::within(
            model.view(),
            &all,
            &reference,
            Quantity::new(1.0, PICOSECOND)
        ),
        Err(MeasurementError::Unit(_))
    ));
    let foreign = AtomSelection::all(&std::sync::Arc::new(top.perceived().unwrap()));
    assert_eq!(
        measure::within(
            model.view(),
            &foreign,
            &reference,
            Quantity::new(1.0, NANOMETER)
        ),
        Err(MeasurementError::Selection(
            SelectionError::TopologyMismatch
        ))
    );
    assert_eq!(
        measure::within(model.view(), &all, &foreign, Quantity::new(1.0, NANOMETER)),
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
            model.view(),
            &all,
            &reference,
            Quantity::new(1.0, NANOMETER)
        )
        .unwrap()
        .atom_ids()
        .collect::<Vec<_>>(),
        [ids[1], ids[2]]
    );
    model.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    close(
        measure::distance(model.view(), ids[0], ids[1])
            .unwrap()
            .value_in(NANOMETER)
            .unwrap(),
        10.0,
    );
}
