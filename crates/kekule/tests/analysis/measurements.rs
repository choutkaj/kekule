use kekule::{
    diagnostics::TopologyMismatch,
    math::{PeriodicCell, Point3, Vector3},
    smiles,
    structure::{
        measure::{self, ConnectivityCheck},
        CoordinateError, Positions, Structure,
    },
    topology::{AtomSelection, InstanceAtomId},
    units::{Quantity, Unit, ANGSTROM, DEGREE, NANOMETER, PICOSECOND},
};

fn structure(points: [Point3; 4], unit: Unit) -> Structure {
    let molecule = smiles::to_molecules("CCCC").unwrap().pop().unwrap();
    Structure::from_molecule(
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
        let structure = structure(
            right_angle().map(|p| Point3::new(p.x * scale, p.y * scale, p.z * scale)),
            unit,
        );
        let [a, b, c, d] =
            <[InstanceAtomId; 4]>::try_from(structure.topology().atom_ids()).unwrap();
        close(
            measure::distance(
                structure.as_structure_view(),
                a,
                b,
                ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap(),
            1.0,
        );
        close(
            measure::angle(
                structure.as_structure_view(),
                a,
                b,
                c,
                ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
            90.0,
        );
        close(
            measure::dihedral(
                structure.as_structure_view(),
                a,
                b,
                c,
                d,
                ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
            -90.0,
        );
        close(
            measure::distance(
                structure.as_structure_view(),
                a,
                a,
                ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .into_value(),
            0.0,
        );
    }
    let reflected = structure(right_angle().map(|p| Point3::new(p.x, p.y, -p.z)), ANGSTROM);
    let a = reflected.topology().atom_ids();
    close(
        measure::dihedral(
            reflected.as_structure_view(),
            a[0],
            a[1],
            a[2],
            a[3],
            ConnectivityCheck::Unrestricted,
        )
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
    let structure = Structure::from_molecule(
        molecule.clone(),
        &Positions::new(Quantity::new(right_angle(), ANGSTROM)).unwrap(),
    )
    .unwrap();
    let ids = <[InstanceAtomId; 4]>::try_from(structure.topology().atom_ids()).unwrap();
    let mut reversed = ids;
    reversed.reverse();
    for [a, b, c, d] in [ids, reversed] {
        assert_eq!(
            measure::distance(
                structure.as_structure_view(),
                a,
                b,
                ConnectivityCheck::ConsecutiveBonds
            ),
            measure::distance(
                structure.as_structure_view(),
                a,
                b,
                ConnectivityCheck::Unrestricted
            )
        );
        assert_eq!(
            measure::angle(
                structure.as_structure_view(),
                a,
                b,
                c,
                ConnectivityCheck::ConsecutiveBonds
            ),
            measure::angle(
                structure.as_structure_view(),
                a,
                b,
                c,
                ConnectivityCheck::Unrestricted
            )
        );
        assert_eq!(
            measure::dihedral(
                structure.as_structure_view(),
                a,
                b,
                c,
                d,
                ConnectivityCheck::ConsecutiveBonds
            ),
            measure::dihedral(
                structure.as_structure_view(),
                a,
                b,
                c,
                d,
                ConnectivityCheck::Unrestricted
            )
        );
    }
}

#[test]
fn optional_connectivity_reports_each_missing_consecutive_pair() {
    let structure = structure(right_angle(), NANOMETER);
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(structure.topology().atom_ids()).unwrap();
    assert_eq!(
        ConnectivityCheck::default(),
        ConnectivityCheck::Unrestricted
    );
    for [x, y] in [[a, c], [a, a]] {
        assert_eq!(
            measure::distance(
                structure.as_structure_view(),
                x,
                y,
                ConnectivityCheck::ConsecutiveBonds
            ),
            Err(CoordinateError::MissingBond { a: x, b: y })
        );
        let unrestricted = measure::distance(
            structure.as_structure_view(),
            x,
            y,
            ConnectivityCheck::Unrestricted,
        )
        .unwrap();
        assert_eq!(
            measure::distance(
                structure.as_structure_view(),
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
            measure::angle(
                structure.as_structure_view(),
                x,
                y,
                z,
                ConnectivityCheck::ConsecutiveBonds
            ),
            Err(CoordinateError::MissingBond {
                a: missing[0],
                b: missing[1]
            })
        );
        let unrestricted = measure::angle(
            structure.as_structure_view(),
            x,
            y,
            z,
            ConnectivityCheck::Unrestricted,
        )
        .unwrap();
        assert_eq!(
            measure::angle(
                structure.as_structure_view(),
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
            measure::dihedral(
                structure.as_structure_view(),
                w,
                x,
                y,
                z,
                ConnectivityCheck::ConsecutiveBonds
            ),
            Err(CoordinateError::MissingBond {
                a: missing[0],
                b: missing[1]
            })
        );
        let unrestricted = measure::dihedral(
            structure.as_structure_view(),
            w,
            x,
            y,
            z,
            ConnectivityCheck::Unrestricted,
        )
        .unwrap();
        assert_eq!(
            measure::dihedral(
                structure.as_structure_view(),
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
    let mut builder = Structure::builder();
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
    let structure = builder.build().unwrap();
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(structure.topology().atom_ids()).unwrap();
    let missing = Err(CoordinateError::MissingBond { a: b, b: c });
    assert_eq!(
        measure::distance(
            structure.as_structure_view(),
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds
        ),
        missing
    );
    assert_eq!(
        measure::angle(
            structure.as_structure_view(),
            a,
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds
        ),
        missing
    );
    assert_eq!(
        measure::dihedral(
            structure.as_structure_view(),
            a,
            b,
            c,
            d,
            ConnectivityCheck::ConsecutiveBonds
        ),
        missing
    );
    assert!(measure::distance(
        structure.as_structure_view(),
        b,
        c,
        ConnectivityCheck::Unrestricted
    )
    .is_ok());
    assert!(measure::angle(
        structure.as_structure_view(),
        a,
        b,
        c,
        ConnectivityCheck::Unrestricted
    )
    .is_ok());
    assert!(measure::dihedral(
        structure.as_structure_view(),
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
    let structure = structure(right_angle(), NANOMETER);
    let [a, b, c, _] = <[InstanceAtomId; 4]>::try_from(structure.topology().atom_ids()).unwrap();
    for invalid in [
        InstanceAtomId::new(a.instance(), kekule::molecule::AtomId::new(99)),
        InstanceAtomId::new(kekule::topology::MoleculeInstanceId::new(99), a.atom()),
    ] {
        let error = Err(CoordinateError::InvalidAtomId(invalid));
        assert_eq!(
            measure::distance(
                structure.as_structure_view(),
                a,
                invalid,
                ConnectivityCheck::ConsecutiveBonds
            ),
            error
        );
        assert_eq!(
            measure::angle(
                structure.as_structure_view(),
                a,
                c,
                invalid,
                ConnectivityCheck::ConsecutiveBonds
            ),
            error
        );
        assert_eq!(
            measure::dihedral(
                structure.as_structure_view(),
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
    let line = structure(
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
        measure::angle(
            line.as_structure_view(),
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
        measure::dihedral(
            line.as_structure_view(),
            a,
            b,
            c,
            d,
            ConnectivityCheck::ConsecutiveBonds
        ),
        Err(CoordinateError::DegenerateGeometry)
    );
    // Nonconsecutive repeats are allowed when the geometry remains defined.
    close(
        measure::angle(
            line.as_structure_view(),
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
    let coincident = structure([Point3::origin(); 4], NANOMETER);
    let [a, b, c, d] = <[InstanceAtomId; 4]>::try_from(coincident.topology().atom_ids()).unwrap();
    assert_eq!(
        measure::distance(
            coincident.as_structure_view(),
            a,
            b,
            ConnectivityCheck::ConsecutiveBonds
        )
        .unwrap()
        .into_value(),
        0.0
    );
    assert_eq!(
        measure::angle(
            coincident.as_structure_view(),
            a,
            b,
            c,
            ConnectivityCheck::ConsecutiveBonds
        ),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(
        measure::dihedral(
            coincident.as_structure_view(),
            a,
            b,
            c,
            d,
            ConnectivityCheck::ConsecutiveBonds
        ),
        Err(CoordinateError::DegenerateGeometry)
    );
}

#[test]
fn tiny_finite_geometry_normalizes_without_reciprocal_overflow() {
    let tiny = structure(
        right_angle().map(|p| Point3::new(p.x * 1e-310, p.y * 1e-310, p.z * 1e-310)),
        NANOMETER,
    );
    let a = tiny.topology().atom_ids();
    close(
        measure::angle(
            tiny.as_structure_view(),
            a[0],
            a[1],
            a[2],
            ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        90.0,
    );
    close(
        measure::dihedral(
            tiny.as_structure_view(),
            a[0],
            a[1],
            a[2],
            a[3],
            ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        -90.0,
    );
}

#[test]
fn undefined_geometry_and_numerical_overflow_reject_without_nan_results() {
    let line = structure(
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
        measure::angle(
            line.as_structure_view(),
            a[0],
            a[1],
            a[2],
            ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        180.0,
    );
    assert_eq!(
        measure::angle(
            line.as_structure_view(),
            a[0],
            a[0],
            a[1],
            ConnectivityCheck::Unrestricted
        ),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(
        measure::dihedral(
            line.as_structure_view(),
            a[0],
            a[1],
            a[2],
            a[3],
            ConnectivityCheck::Unrestricted
        ),
        Err(CoordinateError::DegenerateGeometry)
    );
    let invalid = InstanceAtomId::new(a[0].instance(), kekule::molecule::AtomId::new(99));
    assert!(matches!(
        measure::distance(
            line.as_structure_view(),
            invalid,
            a[0],
            ConnectivityCheck::Unrestricted
        ),
        Err(CoordinateError::InvalidAtomId(_))
    ));
    let huge = structure(
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
        measure::distance(
            huge.as_structure_view(),
            a[0],
            a[1],
            ConnectivityCheck::Unrestricted
        ),
        Err(CoordinateError::NumericalFailure)
    );
    let large = structure(
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
        measure::angle(
            large.as_structure_view(),
            a[0],
            a[1],
            a[2],
            ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        90.0,
    );
}

#[test]
fn spatial_selection_is_inclusive_unit_aware_and_specific_to_each_view() {
    let mut structure = structure(right_angle(), NANOMETER);
    let top = structure.shared_topology();
    let ids = top.atom_ids();
    let all = AtomSelection::all(top.clone());
    let reference = AtomSelection::from_atoms(top.clone(), [ids[1]]).unwrap();
    // Test exact inclusion in canonical units, separately from conversion:
    // 10 * the floating-point Angstrom scale can round just below 1 nm.
    let selected = measure::within(
        structure.as_structure_view(),
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
            structure.as_structure_view(),
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
            structure.as_structure_view(),
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
            structure.as_structure_view(),
            &all,
            &reference,
            Quantity::new(0.0, NANOMETER)
        )
        .unwrap(),
        reference
    );
    let empty = AtomSelection::from_atoms(top.clone(), []).unwrap();
    assert_eq!(
        measure::within(
            structure.as_structure_view(),
            &all,
            &empty,
            Quantity::new(1.0, NANOMETER)
        )
        .unwrap(),
        empty
    );
    assert_eq!(
        measure::within(
            structure.as_structure_view(),
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
                structure.as_structure_view(),
                &empty,
                &reference,
                Quantity::new(cutoff, NANOMETER)
            ),
            Err(CoordinateError::InvalidCutoff)
        );
    }
    assert!(matches!(
        measure::within(
            structure.as_structure_view(),
            &all,
            &reference,
            Quantity::new(1.0, PICOSECOND)
        ),
        Err(CoordinateError::Unit(_))
    ));
    // A perceived snapshot shares the layout, so its selections stay usable.
    let perceived = AtomSelection::all(top.perceived().unwrap());
    assert_eq!(
        measure::within(
            structure.as_structure_view(),
            &perceived,
            &reference,
            Quantity::new(1.0, NANOMETER)
        ),
        Ok(selected.clone())
    );
    // An independently published equal topology does not.
    let foreign = AtomSelection::all(smiles::to_topology("CCCC").unwrap());
    assert_eq!(
        measure::within(
            structure.as_structure_view(),
            &foreign,
            &reference,
            Quantity::new(1.0, NANOMETER)
        ),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Selection
        ))
    );
    assert_eq!(
        measure::within(
            structure.as_structure_view(),
            &all,
            &foreign,
            Quantity::new(1.0, NANOMETER)
        ),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Selection
        ))
    );
    structure
        .set_position(
            ids[0],
            Quantity::new(Point3::new(10.0, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    assert_eq!(selected.atom_ids().len(), 3); // The earlier atom set stays fixed.
    assert_eq!(
        measure::within(
            structure.as_structure_view(),
            &all,
            &reference,
            Quantity::new(1.0, NANOMETER)
        )
        .unwrap()
        .atom_ids()
        .collect::<Vec<_>>(),
        [ids[1], ids[2]]
    );
    structure.geometry_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    close(
        measure::distance(
            structure.as_structure_view(),
            ids[0],
            ids[1],
            ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .value_in(NANOMETER)
        .unwrap(),
        10.0,
    );
}

/// Large selections, which `within` searches through a grid, select what
/// measuring every pair selects, in candidate order: seeded atoms, atoms on
/// a 0.25 nm lattice, and atoms exactly at the cutoff.
#[test]
fn spatial_selection_of_many_atoms_matches_every_pair() {
    let count = 400;
    let mut state = 0x853c_49e6_748f_ea9b_u64;
    let mut unit = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let molecule = smiles::to_molecules(&"C".repeat(count))
        .unwrap()
        .pop()
        .unwrap();
    for (spread, cutoff) in [
        (3.0, 0.0),
        (3.0, 0.25),
        (3.0, 0.4),
        (8.0, 1.0),
        (1e250, 0.5),
    ] {
        let mut points = Vec::new();
        for index in 0..count {
            let point = match index % 4 {
                0 if index > 0 => {
                    let previous: Point3 = points[index - 1];
                    Point3::new(previous.x, previous.y + cutoff, previous.z)
                }
                1 => Point3::new(
                    (unit() * 12.0).floor() * 0.25,
                    (unit() * 12.0).floor() * 0.25,
                    (unit() * 12.0).floor() * 0.25,
                ),
                _ => Point3::new(unit() * spread, unit() * 3.0, unit() * 3.0),
            };
            points.push(point);
        }
        let structure = Structure::from_molecule(
            molecule.clone(),
            &Positions::new(Quantity::new(&points[..], NANOMETER)).unwrap(),
        )
        .unwrap();
        let top = structure.shared_topology();
        let ids = top.atom_ids();
        let candidates = AtomSelection::from_atoms(
            top.clone(),
            ids.iter()
                .copied()
                .step_by(2)
                .chain(ids.iter().copied().skip(1).step_by(6)),
        )
        .unwrap();
        let reference =
            AtomSelection::from_atoms(top.clone(), ids.iter().copied().skip(1).step_by(3)).unwrap();
        let index = |atom: InstanceAtomId| ids.iter().position(|id| *id == atom).unwrap();
        let expected = candidates
            .atom_ids()
            .filter(|&atom| {
                reference
                    .atom_ids()
                    .any(|other| (points[index(atom)] - points[index(other)]).norm() <= cutoff)
            })
            .collect::<Vec<_>>();
        assert!(!expected.is_empty());
        let selected = measure::within(
            structure.as_structure_view(),
            &candidates,
            &reference,
            Quantity::new(cutoff, NANOMETER),
        )
        .unwrap();
        assert_eq!(
            selected.atom_ids().collect::<Vec<_>>(),
            expected,
            "{spread} {cutoff}"
        );
    }
}
