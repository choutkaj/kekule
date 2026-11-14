use super::*;
use crate::{
    math::{Matrix3, PeriodicCell},
    properties::{PropertyKey, PropertyValue},
    smiles,
    structure::Positions,
    units::{ANGSTROM, DEGREE, NANOMETER, PICOSECOND, RADIAN},
};

fn fixture(smiles: &str, points: &[Point3]) -> Structure {
    Structure::new(
        smiles::to_topology(smiles).unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn chain() -> Structure {
    fixture(
        "CCCC",
        &[
            Point3::new(0.0, 1.0, 0.0),
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
        ],
    )
}

fn ids(structure: &Structure) -> [InstanceAtomId; 4] {
    structure.topology().atom_ids().try_into().unwrap()
}
fn point(structure: &Structure, atom: InstanceAtomId) -> Point3 {
    structure.position(atom).unwrap().into_value()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}
fn close_point(a: Point3, b: Point3) {
    close((a - b).norm(), 0.0);
}

#[test]
fn distance_moves_entire_branched_fragment_and_preserves_other_instances() {
    let mut structure = fixture(
        "CC(C)C.O",
        &[
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(3.0, 4.0, 5.0),
        ],
    );
    let atoms = structure.topology().atom_ids().to_vec();
    let baseline = structure.clone();
    let edit = DistanceEdit::new(structure.shared_topology(), [atoms[0], atoms[1]]).unwrap();
    assert_eq!(
        edit.moving_atoms().atom_ids().collect::<Vec<_>>(),
        atoms[1..4]
    );
    assert_eq!(edit.atoms(), [atoms[0], atoms[1]]);
    structure
        .set_distance(atoms[0], atoms[1], Quantity::new(15.0, ANGSTROM))
        .unwrap();
    close(
        edit.measure(structure.as_structure_view())
            .unwrap()
            .value_in(NANOMETER)
            .unwrap(),
        1.5,
    );
    for &atom in &atoms[1..4] {
        close_point(
            point(&structure, atom),
            point(&baseline, atom) + Vector3::new(0.5, 0.0, 0.0),
        );
    }
    for &atom in &[atoms[0], atoms[4]] {
        assert_eq!(point(&structure, atom), point(&baseline, atom));
    }
    assert!(Arc::ptr_eq(
        &structure.shared_topology(),
        &baseline.shared_topology()
    ));
}

#[test]
fn angle_moves_fragment_rigidly_and_reversal_moves_opposite_side() {
    let baseline = chain();
    let [a, b, c, d] = ids(&baseline);
    let edit = AngleEdit::new(baseline.shared_topology(), [a, b, c]).unwrap();
    assert_eq!(edit.atoms(), [a, b, c]);
    for target in [0.0, 35.0, 120.0, 180.0] {
        let mut structure = baseline.clone();
        structure
            .set_angle(a, b, c, Quantity::new(target, DEGREE))
            .unwrap();
        close(
            edit.measure(structure.as_structure_view())
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            target,
        );
        assert_eq!(point(&structure, a), point(&baseline, a));
        assert_eq!(point(&structure, b), point(&baseline, b));
        close(
            (point(&structure, d) - point(&structure, c)).norm(),
            (point(&baseline, d) - point(&baseline, c)).norm(),
        );
        close((point(&structure, c) - point(&structure, b)).norm(), 1.0);
    }
    let mut reversed = baseline.clone();
    reversed
        .set_angle(c, b, a, Quantity::new(120.0, DEGREE))
        .unwrap();
    assert_ne!(point(&reversed, a), point(&baseline, a));
    assert_eq!(point(&reversed, c), point(&baseline, c));
    assert_eq!(point(&reversed, d), point(&baseline, d));
}

#[test]
fn dihedral_scan_uses_current_coordinates_and_matches_measurement_sign() {
    let mut structure = chain();
    let baseline = structure.clone();
    let [a, b, c, d] = ids(&structure);
    let edit = DihedralEdit::new(structure.shared_topology(), [a, b, c, d]).unwrap();
    assert_eq!(edit.atoms(), [a, b, c, d]);
    close(
        edit.measure(structure.as_structure_view())
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        45.0,
    );
    for target in [
        -180.0, -179.0, -90.0, 0.0, 120.0, 179.0, 180.0, 540.0, -720.0,
    ] {
        edit.apply(&mut structure, Quantity::new(target, DEGREE))
            .unwrap();
        let actual = measure::dihedral(
            structure.as_structure_view(),
            a,
            b,
            c,
            d,
            measure::ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .into_value();
        close(signed_angle(actual - target.to_radians()), 0.0);
        for atom in [a, b, c] {
            assert_eq!(point(&structure, atom), point(&baseline, atom));
        }
        close(
            (point(&structure, d) - point(&structure, c)).norm(),
            2.0_f64.sqrt(),
        );
        let once = structure.clone();
        edit.apply(&mut structure, Quantity::new(target, DEGREE))
            .unwrap();
        close_point(point(&structure, d), point(&once, d));
    }
    // A changed reference plane must be read at application time.
    structure
        .set_position(a, Quantity::new(Point3::new(0.0, 0.0, 1.0), NANOMETER))
        .unwrap();
    edit.apply(&mut structure, Quantity::new(0.0, RADIAN))
        .unwrap();
    close_point(point(&structure, d), Point3::new(1.0, 0.0, 2.0_f64.sqrt()));
    let mut independent = baseline.clone();
    independent
        .set_dihedral(a, b, c, d, Quantity::new(90.0, DEGREE))
        .unwrap();
    close_point(
        point(&independent, d),
        Point3::new(1.0, 0.0, 2.0_f64.sqrt()),
    );
    // Reversing all references retains the measurement sign, and moves A's side.
    let mut reversed = baseline.clone();
    reversed
        .set_dihedral(d, c, b, a, Quantity::new(-90.0, DEGREE))
        .unwrap();
    close(
        measure::dihedral(
            reversed.as_structure_view(),
            a,
            b,
            c,
            d,
            measure::ConnectivityCheck::Unrestricted,
        )
        .unwrap()
        .value_in(DEGREE)
        .unwrap(),
        -90.0,
    );
    assert_eq!(point(&reversed, d), point(&baseline, d));
}

#[test]
fn half_turn_ties_and_large_finite_targets_wrap_without_overflow() {
    assert_eq!(signed_angle(-PI), PI);
    assert_eq!(signed_angle(PI), PI);
    let mut structure = chain();
    let [a, b, c, d] = ids(&structure);
    let edit = DihedralEdit::new(structure.shared_topology(), [a, b, c, d]).unwrap();
    for target in [PI, -PI, f64::MAX, -f64::MAX] {
        edit.apply(&mut structure, Quantity::new(target, RADIAN))
            .unwrap();
        close(
            signed_angle(
                edit.measure(structure.as_structure_view())
                    .unwrap()
                    .into_value()
                    - signed_angle(target),
            ),
            0.0,
        );
    }
}

#[test]
fn rings_reject_automatic_edits_but_explicit_edits_can_deform_them() {
    let mut ring = fixture(
        "C1CCC1",
        &[
            Point3::new(0.0, 1.0, 0.0),
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
    );
    let [a, b, c, d] = ids(&ring);
    let topology = ring.shared_topology();
    assert_eq!(
        DistanceEdit::new(topology.clone(), [a, b]).unwrap_err(),
        CoordinateError::InseparableFragment
    );
    assert_eq!(
        AngleEdit::new(topology.clone(), [a, b, c]).unwrap_err(),
        CoordinateError::InseparableFragment
    );
    assert_eq!(
        DihedralEdit::new(topology.clone(), [a, b, c, d]).unwrap_err(),
        CoordinateError::InseparableFragment
    );
    let baseline = ring.clone();
    let selected = AtomSelection::from_atoms(topology.clone(), [b, c, d]).unwrap();
    let edit = DihedralEdit::with_moving_atoms(topology.clone(), [a, b, c, d], &selected).unwrap();
    edit.apply(&mut ring, Quantity::new(90.0, DEGREE)).unwrap();
    for atom in [a, b, c] {
        assert_eq!(point(&ring, atom), point(&baseline, atom));
    }
    close((point(&ring, d) - point(&ring, a)).norm(), 3.0_f64.sqrt());
    close(
        edit.measure(ring.as_structure_view())
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
}

#[test]
fn explicit_distance_and_angle_allow_disconnected_atoms_and_keep_pivots() {
    let mut structure = fixture(
        "C.C.C.C",
        &[
            Point3::new(0.0, 1.0, 0.0),
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 5.0),
        ],
    );
    let [a, b, c, d] = ids(&structure);
    let topology = structure.shared_topology();
    let moving = AtomSelection::from_atoms(topology.clone(), [b, c]).unwrap();
    let baseline = structure.clone();
    let angle = AngleEdit::with_moving_atoms(topology.clone(), [a, b, c], &moving).unwrap();
    angle
        .apply(&mut structure, Quantity::new(45.0, DEGREE))
        .unwrap();
    for atom in [a, b, d] {
        assert_eq!(point(&structure, atom), point(&baseline, atom));
    }
    close(
        angle
            .measure(structure.as_structure_view())
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        45.0,
    );
    let distance = DistanceEdit::with_moving_atoms(topology.clone(), [a, c], &moving).unwrap();
    distance
        .apply(&mut structure, Quantity::new(2.0, NANOMETER))
        .unwrap();
    close(
        distance
            .measure(structure.as_structure_view())
            .unwrap()
            .into_value(),
        2.0,
    );
    assert_ne!(point(&structure, b), point(&baseline, b));
    assert_eq!(point(&structure, d), point(&baseline, d));
}

#[test]
fn preparation_validates_references_connectivity_and_selections() {
    let structure = chain();
    let [a, b, c, d] = ids(&structure);
    let topology = structure.shared_topology();
    assert_eq!(
        DistanceEdit::new(topology.clone(), [a, a]).unwrap_err(),
        CoordinateError::RepeatedAtom(a)
    );
    assert_eq!(
        AngleEdit::new(topology.clone(), [a, c, d]).unwrap_err(),
        CoordinateError::MissingBond { a, b: c }
    );
    assert_eq!(
        DihedralEdit::new(topology.clone(), [a, b, d, c]).unwrap_err(),
        CoordinateError::MissingBond { a: b, b: d }
    );
    let larger = fixture("CCCCC", &[Point3::origin(); 5]);
    let invalid = larger.topology().atom_ids()[4];
    assert_eq!(
        DistanceEdit::new(topology.clone(), [a, invalid]).unwrap_err(),
        CoordinateError::InvalidAtomId(invalid)
    );
    for selected in [vec![], vec![a, d], vec![b, c]] {
        let moving = AtomSelection::from_atoms(topology.clone(), selected).unwrap();
        assert_eq!(
            DihedralEdit::with_moving_atoms(topology.clone(), [a, b, c, d], &moving).unwrap_err(),
            CoordinateError::InvalidMovingSelection
        );
    }
    let different = chain();
    let foreign = AtomSelection::all(different.shared_topology());
    assert_eq!(
        DistanceEdit::with_moving_atoms(topology.clone(), [a, b], &foreign).unwrap_err(),
        CoordinateError::TopologyMismatch(TopologyMismatch::Selection)
    );
    // Consecutive bonds need not be single or chemically rotatable.
    let unsaturated = fixture("C=CC#N", structure.positions().values().value());
    assert!(DihedralEdit::new(unsaturated.shared_topology(), ids(&unsaturated)).is_ok());
}

#[test]
fn prepared_edits_reject_independently_equal_topologies_atomically() {
    let structure = chain();
    let [a, b, c, d] = ids(&structure);
    let edit = DihedralEdit::new(structure.shared_topology(), [a, b, c, d]).unwrap();
    let mut other = chain();
    let baseline = other.clone();
    assert_eq!(
        edit.measure(other.as_structure_view()),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Structure
        ))
    );
    assert_eq!(
        edit.apply(&mut other, Quantity::new(30.0, DEGREE)),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Structure
        ))
    );
    assert_eq!(other, baseline);
    let mut clone = structure.clone();
    edit.apply(&mut clone, Quantity::new(30.0, DEGREE)).unwrap();
}

#[test]
fn invalid_targets_and_units_leave_structure_unchanged() {
    let mut structure = chain();
    let [a, b, c, d] = ids(&structure);
    let baseline = structure.clone();
    for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            structure.set_distance(a, b, Quantity::new(target, NANOMETER)),
            Err(CoordinateError::InvalidTarget)
        );
        assert_eq!(structure, baseline);
    }
    for target in [-1.0, 181.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            structure.set_angle(a, b, c, Quantity::new(target, DEGREE)),
            Err(CoordinateError::InvalidTarget)
        );
        assert_eq!(structure, baseline);
    }
    assert_eq!(
        structure.set_dihedral(a, b, c, d, Quantity::new(f64::NAN, DEGREE)),
        Err(CoordinateError::InvalidTarget)
    );
    assert!(matches!(
        structure.set_angle(a, b, c, Quantity::new(1.0, NANOMETER)),
        Err(CoordinateError::Unit(_))
    ));
    assert!(matches!(
        structure.set_distance(a, b, Quantity::new(1.0, PICOSECOND)),
        Err(CoordinateError::Unit(_))
    ));
    assert_eq!(structure, baseline);
}

#[test]
fn degenerate_geometry_is_rejected_without_inventing_axes() {
    let mut structure = chain();
    let [a, b, c, d] = ids(&structure);
    structure
        .set_position(a, Quantity::new(Point3::new(-1.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    let straight = structure.clone();
    structure
        .set_angle(a, b, c, Quantity::new(180.0, DEGREE))
        .unwrap();
    assert_eq!(structure, straight);
    assert_eq!(
        structure.set_angle(a, b, c, Quantity::new(90.0, DEGREE)),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(
        structure.set_dihedral(a, b, c, d, Quantity::new(0.0, RADIAN)),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(structure, straight);
    structure
        .set_position(a, structure.position(b).unwrap())
        .unwrap();
    let coincident = structure.clone();
    assert_eq!(
        structure.set_distance(a, b, Quantity::new(1.0, NANOMETER)),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(
        structure.set_angle(a, b, c, Quantity::new(90.0, DEGREE)),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(structure, coincident);
}

#[test]
fn numerical_failures_in_late_atoms_do_not_publish_earlier_atoms() {
    let mut structure = chain();
    let [a, b, c, d] = ids(&structure);
    structure
        .set_position(d, Quantity::new(Point3::new(f64::MAX, 0.0, 0.0), NANOMETER))
        .unwrap();
    let baseline = structure.clone();
    assert_eq!(
        structure.set_distance(b, c, Quantity::new(f64::MAX, NANOMETER)),
        Err(CoordinateError::NumericalFailure)
    );
    assert_eq!(structure, baseline);
    let all = AtomSelection::all(structure.shared_topology());
    assert_eq!(
        structure.translate(
            &all,
            Quantity::new(Vector3::new(f64::MAX, 0.0, 0.0), NANOMETER)
        ),
        Err(CoordinateError::NumericalFailure)
    );
    assert_eq!(structure, baseline);
    let transform =
        RigidTransform::from_canonical(Matrix3::identity(), Vector3::new(f64::MAX, 0.0, 0.0))
            .unwrap();
    assert_eq!(
        structure.transform_atoms(&all, &transform),
        Err(CoordinateError::NumericalFailure)
    );
    assert_eq!(structure, baseline);
    structure
        .set_position(
            a,
            Quantity::new(Point3::new(-f64::MAX, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    let baseline = structure.clone();
    assert_eq!(
        structure.rotate(
            &all,
            Quantity::new(Point3::new(f64::MAX, 0.0, 0.0), NANOMETER),
            Vector3::new(0.0, 0.0, 1.0),
            Quantity::new(90.0, DEGREE)
        ),
        Err(CoordinateError::NumericalFailure)
    );
    assert_eq!(structure, baseline);
}

#[test]
fn finite_but_unachievable_distance_returns_numerical_failure() {
    let mut structure = fixture(
        "CC",
        &[
            Point3::new(1.0e20, 0.0, 0.0),
            Point3::new(1.0e20 + 1.0e6, 0.0, 0.0),
        ],
    );
    let [a, b] = structure.topology().atom_ids().try_into().unwrap();
    let baseline = structure.clone();
    assert_eq!(
        structure.set_distance(a, b, Quantity::new(1.0, NANOMETER)),
        Err(CoordinateError::NumericalFailure)
    );
    assert_eq!(structure, baseline);
}

#[test]
fn rigid_operations_use_units_right_handed_rotation_and_exact_selection() {
    let mut structure = chain();
    let [a, b, c, d] = ids(&structure);
    let baseline = structure.clone();
    let selected = AtomSelection::from_atoms(structure.shared_topology(), [c, d]).unwrap();
    structure
        .translate(
            &selected,
            Quantity::new(Vector3::new(10.0, 20.0, 30.0), ANGSTROM),
        )
        .unwrap();
    close_point(point(&structure, c), Point3::new(2.0, 2.0, 3.0));
    structure
        .rotate(
            &selected,
            Quantity::new(Point3::new(10.0, 20.0, 30.0), ANGSTROM),
            Vector3::new(0.0, 0.0, 8.0),
            Quantity::new(90.0, DEGREE),
        )
        .unwrap();
    close_point(point(&structure, c), Point3::new(1.0, 3.0, 3.0));
    close_point(point(&structure, d), Point3::new(0.0, 3.0, 4.0));
    for atom in [a, b] {
        assert_eq!(point(&structure, atom), point(&baseline, atom));
    }
    let transform =
        RigidTransform::from_canonical(Matrix3::identity(), Vector3::new(1.0, 0.0, 0.0)).unwrap();
    structure.transform_atoms(&selected, &transform).unwrap();
    close_point(point(&structure, c), Point3::new(2.0, 3.0, 3.0));
    // Both huge and subnormal finite axes are accepted.
    for scale in [f64::MAX, f64::from_bits(1)] {
        let mut rotated = baseline.clone();
        rotated
            .rotate(
                &selected,
                Quantity::new(Point3::origin(), NANOMETER),
                Vector3::new(0.0, 0.0, scale),
                Quantity::new(90.0, DEGREE),
            )
            .unwrap();
        close_point(point(&rotated, c), Point3::new(0.0, 1.0, 0.0));
    }
}

#[test]
fn empty_selections_still_validate_binding_and_motion_arguments() {
    let mut structure = chain();
    let baseline = structure.clone();
    let empty = AtomSelection::empty(structure.shared_topology());
    let zero = Quantity::new(Point3::origin(), NANOMETER);
    let axis = Vector3::new(0.0, 0.0, 1.0);
    structure
        .translate(
            &empty,
            Quantity::new(Vector3::new(1.0, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    structure
        .rotate(&empty, zero, axis, Quantity::new(90.0, DEGREE))
        .unwrap();
    structure
        .transform_atoms(&empty, &RigidTransform::identity())
        .unwrap();
    assert_eq!(structure, baseline);
    assert_eq!(
        structure.rotate(&empty, zero, Vector3::zero(), Quantity::new(0.0, DEGREE)),
        Err(CoordinateError::DegenerateGeometry)
    );
    assert_eq!(
        structure.rotate(&empty, zero, axis, Quantity::new(f64::NAN, DEGREE)),
        Err(CoordinateError::InvalidTarget)
    );
    assert!(matches!(
        structure.rotate(&empty, zero, axis, Quantity::new(1.0, NANOMETER)),
        Err(CoordinateError::Unit(_))
    ));
    assert!(matches!(
        structure.translate(&empty, Quantity::new(Vector3::zero(), PICOSECOND)),
        Err(CoordinateError::Unit(_))
    ));
    assert_eq!(
        structure.translate(
            &empty,
            Quantity::new(Vector3::new(f64::INFINITY, 0.0, 0.0), NANOMETER)
        ),
        Err(CoordinateError::NonFiniteInput)
    );
    assert_eq!(
        structure.rotate(
            &empty,
            Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), NANOMETER),
            axis,
            Quantity::new(0.0, DEGREE)
        ),
        Err(CoordinateError::NonFiniteInput)
    );
    assert_eq!(
        structure.rotate(
            &empty,
            zero,
            Vector3::new(f64::NAN, 0.0, 0.0),
            Quantity::new(0.0, DEGREE)
        ),
        Err(CoordinateError::NonFiniteInput)
    );
    let foreign = AtomSelection::empty(chain().shared_topology());
    assert_eq!(
        structure.translate(&foreign, Quantity::new(Vector3::zero(), NANOMETER)),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Selection
        ))
    );
    assert_eq!(
        structure.rotate(&foreign, zero, axis, Quantity::new(0.0, DEGREE)),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Selection
        ))
    );
    assert_eq!(
        structure.transform_atoms(&foreign, &RigidTransform::identity()),
        Err(CoordinateError::TopologyMismatch(
            TopologyMismatch::Selection
        ))
    );
    assert_eq!(structure, baseline);
}

#[test]
fn all_edits_preserve_shared_topology_cell_properties_and_chemistry() {
    let mut structure = chain();
    let [a, b, c, d] = ids(&structure);
    structure.geometry_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(5.0, 6.0, 7.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    structure
        .geometry_mut()
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("test.note").unwrap(),
            PropertyValue::String("preserved".into()),
        )
        .unwrap();
    let a_index = structure.topology().atom_index(a).unwrap();
    structure
        .geometry_mut()
        .set_occupancy(a_index, Some(0.75))
        .unwrap();
    let baseline = structure.clone();
    let all = AtomSelection::all(structure.shared_topology());
    structure
        .set_distance(b, c, Quantity::new(2.0, NANOMETER))
        .unwrap();
    structure
        .set_angle(a, b, c, Quantity::new(100.0, DEGREE))
        .unwrap();
    structure
        .set_dihedral(a, b, c, d, Quantity::new(-45.0, DEGREE))
        .unwrap();
    structure
        .translate(&all, Quantity::new(Vector3::new(1.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    structure
        .rotate(
            &all,
            Quantity::new(Point3::origin(), NANOMETER),
            Vector3::new(0.0, 0.0, 1.0),
            Quantity::new(30.0, DEGREE),
        )
        .unwrap();
    structure
        .transform_atoms(&all, &RigidTransform::identity())
        .unwrap();
    assert!(Arc::ptr_eq(
        &structure.shared_topology(),
        &baseline.shared_topology()
    ));
    assert_eq!(structure.cell(), baseline.cell());
    assert_eq!(structure.properties(), baseline.properties());
}

#[test]
fn geometry_edits_work_in_oblique_frames_and_at_different_length_scales() {
    for scale in [1.0e-6, 1.0, 1.0e6] {
        let mut structure = chain();
        let [a, b, c, d] = ids(&structure);
        let points = structure
            .positions()
            .values()
            .value()
            .iter()
            .map(|p| Point3::new(p.x * scale, p.y * scale, p.z * scale))
            .collect::<Vec<_>>();
        structure
            .geometry_mut()
            .set_positions(Quantity::new(points, NANOMETER))
            .unwrap();
        let all = AtomSelection::all(structure.shared_topology());
        structure
            .rotate(
                &all,
                Quantity::new(Point3::origin(), NANOMETER),
                Vector3::new(1.0, 2.0, 3.0),
                Quantity::new(37.0, DEGREE),
            )
            .unwrap();
        structure
            .translate(
                &all,
                Quantity::new(Vector3::new(7.0, -3.0, 1.0) * scale, NANOMETER),
            )
            .unwrap();
        let before = structure.clone();
        structure
            .set_distance(b, c, Quantity::new(1.75 * scale, NANOMETER))
            .unwrap();
        close(
            measure::distance(
                structure.as_structure_view(),
                b,
                c,
                measure::ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .into_value()
                / scale,
            1.75,
        );
        structure
            .set_angle(a, b, c, Quantity::new(47.0, DEGREE))
            .unwrap();
        close(
            measure::angle(
                structure.as_structure_view(),
                a,
                b,
                c,
                measure::ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
            47.0,
        );
        structure
            .set_dihedral(a, b, c, d, Quantity::new(-139.0, DEGREE))
            .unwrap();
        close(
            measure::dihedral(
                structure.as_structure_view(),
                a,
                b,
                c,
                d,
                measure::ConnectivityCheck::Unrestricted,
            )
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
            -139.0,
        );
        assert_eq!(point(&structure, a), point(&before, a));
        assert_eq!(point(&structure, b), point(&before, b));
        close(
            (point(&structure, c) - point(&structure, d)).norm() / scale,
            2.0_f64.sqrt(),
        );
    }
}

#[test]
fn automatic_fragments_may_contain_a_ring_on_the_moving_side() {
    let structure = fixture("CCC1CC1", &[Point3::origin(); 5]);
    let atoms = structure.topology().atom_ids();
    let edit = AngleEdit::new(structure.shared_topology(), [atoms[0], atoms[1], atoms[2]]).unwrap();
    assert_eq!(
        edit.moving_atoms().atom_ids().collect::<Vec<_>>(),
        atoms[2..]
    );
}

#[test]
fn prepared_selection_is_a_snapshot_and_repeated_definitions_stay_independent() {
    let molecule = smiles::to_molecules("CCCC").unwrap().pop().unwrap();
    let mut builder = crate::topology::TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let baseline = chain();
    let points = baseline.positions().values().value().repeat(2);
    let mut structure = Structure::new(
        builder.build().unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap();
    let topology = structure.shared_topology();
    let [a, b, c, d] = structure.topology().atom_ids()[..4].try_into().unwrap();
    let mut moving = AtomSelection::from_atoms(topology.clone(), [d]).unwrap();
    let edit = DihedralEdit::with_moving_atoms(topology.clone(), [a, b, c, d], &moving).unwrap();
    moving.clear();
    assert_eq!(edit.moving_atoms().len(), 1);
    let before = structure.clone();
    edit.apply(&mut structure, Quantity::new(-90.0, DEGREE))
        .unwrap();
    structure
        .set_distance(b, c, Quantity::new(2.0, NANOMETER))
        .unwrap();
    for &atom in &structure.topology().atom_ids()[4..] {
        assert_eq!(point(&structure, atom), point(&before, atom));
    }
}

#[test]
fn transform_atoms_uses_rotation_and_translation_together() {
    let mut structure = chain();
    let baseline = structure.clone();
    let [a, b, c, d] = ids(&structure);
    let selected = AtomSelection::from_atoms(structure.shared_topology(), [c, d]).unwrap();
    let transform = RigidTransform::from_canonical(
        Matrix3::from_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ),
        Vector3::new(2.0, 3.0, 4.0),
    )
    .unwrap();
    structure.transform_atoms(&selected, &transform).unwrap();
    close_point(point(&structure, c), Point3::new(2.0, 4.0, 4.0));
    close_point(point(&structure, d), Point3::new(1.0, 4.0, 5.0));
    for atom in [a, b] {
        assert_eq!(point(&structure, atom), point(&baseline, atom));
    }
}
