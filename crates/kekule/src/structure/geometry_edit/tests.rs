use super::*;
use crate::{
    geometry::{Matrix3, PeriodicCell},
    properties::{PropertyKey, PropertyValue},
    smiles,
    structure::Positions,
    units::{ANGSTROM, DEGREE, NANOMETER, PICOSECOND, RADIAN},
};

fn fixture(smiles: &str, points: &[Point3]) -> Model {
    Model::new(
        smiles::to_topology(smiles).unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn chain() -> Model {
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

fn ids(model: &Model) -> [InstanceAtomId; 4] {
    model.topology().atom_ids().try_into().unwrap()
}
fn point(model: &Model, atom: InstanceAtomId) -> Point3 {
    model.position(atom).unwrap().into_value()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-10, "{a} != {b}");
}
fn close_point(a: Point3, b: Point3) {
    close((a - b).norm(), 0.0);
}

#[test]
fn distance_moves_entire_branched_fragment_and_preserves_other_instances() {
    let mut model = fixture(
        "CC(C)C.O",
        &[
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 1.0),
            Point3::new(3.0, 4.0, 5.0),
        ],
    );
    let atoms = model.topology().atom_ids().to_vec();
    let baseline = model.clone();
    let edit = DistanceEdit::new(&model.shared_topology(), [atoms[0], atoms[1]]).unwrap();
    assert_eq!(
        edit.moving_atoms().atom_ids().collect::<Vec<_>>(),
        atoms[1..4]
    );
    assert_eq!(edit.atoms(), [atoms[0], atoms[1]]);
    model
        .set_distance(atoms[0], atoms[1], Quantity::new(15.0, ANGSTROM))
        .unwrap();
    close(
        edit.measure(model.as_model_view())
            .unwrap()
            .value_in(NANOMETER)
            .unwrap(),
        1.5,
    );
    for &atom in &atoms[1..4] {
        close_point(
            point(&model, atom),
            point(&baseline, atom) + Vector3::new(0.5, 0.0, 0.0),
        );
    }
    for &atom in &[atoms[0], atoms[4]] {
        assert_eq!(point(&model, atom), point(&baseline, atom));
    }
    assert!(Arc::ptr_eq(
        &model.shared_topology(),
        &baseline.shared_topology()
    ));
}

#[test]
fn angle_moves_fragment_rigidly_and_reversal_moves_opposite_side() {
    let baseline = chain();
    let [a, b, c, d] = ids(&baseline);
    let edit = AngleEdit::new(&baseline.shared_topology(), [a, b, c]).unwrap();
    assert_eq!(edit.atoms(), [a, b, c]);
    for target in [0.0, 35.0, 120.0, 180.0] {
        let mut model = baseline.clone();
        model
            .set_angle(a, b, c, Quantity::new(target, DEGREE))
            .unwrap();
        close(
            edit.measure(model.as_model_view())
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            target,
        );
        assert_eq!(point(&model, a), point(&baseline, a));
        assert_eq!(point(&model, b), point(&baseline, b));
        close(
            (point(&model, d) - point(&model, c)).norm(),
            (point(&baseline, d) - point(&baseline, c)).norm(),
        );
        close((point(&model, c) - point(&model, b)).norm(), 1.0);
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
    let mut model = chain();
    let baseline = model.clone();
    let [a, b, c, d] = ids(&model);
    let edit = DihedralEdit::new(&model.shared_topology(), [a, b, c, d]).unwrap();
    assert_eq!(edit.atoms(), [a, b, c, d]);
    close(
        edit.measure(model.as_model_view())
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        45.0,
    );
    for target in [
        -180.0, -179.0, -90.0, 0.0, 120.0, 179.0, 180.0, 540.0, -720.0,
    ] {
        edit.apply(&mut model, Quantity::new(target, DEGREE))
            .unwrap();
        let actual = measure::dihedral(model.as_model_view(), a, b, c, d)
            .unwrap()
            .into_value();
        close(signed_angle(actual - target.to_radians()), 0.0);
        for atom in [a, b, c] {
            assert_eq!(point(&model, atom), point(&baseline, atom));
        }
        close((point(&model, d) - point(&model, c)).norm(), 2.0_f64.sqrt());
        let once = model.clone();
        edit.apply(&mut model, Quantity::new(target, DEGREE))
            .unwrap();
        close_point(point(&model, d), point(&once, d));
    }
    // A changed reference plane must be read at application time.
    model
        .set_position(a, Quantity::new(Point3::new(0.0, 0.0, 1.0), NANOMETER))
        .unwrap();
    edit.apply(&mut model, Quantity::new(0.0, RADIAN)).unwrap();
    close_point(point(&model, d), Point3::new(1.0, 0.0, 2.0_f64.sqrt()));
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
        measure::dihedral(reversed.as_model_view(), a, b, c, d)
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
    let mut model = chain();
    let [a, b, c, d] = ids(&model);
    let edit = DihedralEdit::new(&model.shared_topology(), [a, b, c, d]).unwrap();
    for target in [PI, -PI, f64::MAX, -f64::MAX] {
        edit.apply(&mut model, Quantity::new(target, RADIAN))
            .unwrap();
        close(
            signed_angle(
                edit.measure(model.as_model_view()).unwrap().into_value() - signed_angle(target),
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
        DistanceEdit::new(&topology, [a, b]).unwrap_err(),
        GeometryEditError::InseparableFragment
    );
    assert_eq!(
        AngleEdit::new(&topology, [a, b, c]).unwrap_err(),
        GeometryEditError::InseparableFragment
    );
    assert_eq!(
        DihedralEdit::new(&topology, [a, b, c, d]).unwrap_err(),
        GeometryEditError::InseparableFragment
    );
    let baseline = ring.clone();
    let selected = AtomSelection::from_atoms(&topology, [b, c, d]).unwrap();
    let edit = DihedralEdit::with_moving_atoms(&topology, [a, b, c, d], &selected).unwrap();
    edit.apply(&mut ring, Quantity::new(90.0, DEGREE)).unwrap();
    for atom in [a, b, c] {
        assert_eq!(point(&ring, atom), point(&baseline, atom));
    }
    close((point(&ring, d) - point(&ring, a)).norm(), 3.0_f64.sqrt());
    close(
        edit.measure(ring.as_model_view())
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        90.0,
    );
}

#[test]
fn explicit_distance_and_angle_allow_disconnected_atoms_and_keep_pivots() {
    let mut model = fixture(
        "C.C.C.C",
        &[
            Point3::new(0.0, 1.0, 0.0),
            Point3::origin(),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(3.0, 4.0, 5.0),
        ],
    );
    let [a, b, c, d] = ids(&model);
    let topology = model.shared_topology();
    let moving = AtomSelection::from_atoms(&topology, [b, c]).unwrap();
    let baseline = model.clone();
    let angle = AngleEdit::with_moving_atoms(&topology, [a, b, c], &moving).unwrap();
    angle
        .apply(&mut model, Quantity::new(45.0, DEGREE))
        .unwrap();
    for atom in [a, b, d] {
        assert_eq!(point(&model, atom), point(&baseline, atom));
    }
    close(
        angle
            .measure(model.as_model_view())
            .unwrap()
            .value_in(DEGREE)
            .unwrap(),
        45.0,
    );
    let distance = DistanceEdit::with_moving_atoms(&topology, [a, c], &moving).unwrap();
    distance
        .apply(&mut model, Quantity::new(2.0, NANOMETER))
        .unwrap();
    close(
        distance
            .measure(model.as_model_view())
            .unwrap()
            .into_value(),
        2.0,
    );
    assert_ne!(point(&model, b), point(&baseline, b));
    assert_eq!(point(&model, d), point(&baseline, d));
}

#[test]
fn preparation_validates_references_connectivity_and_selections() {
    let model = chain();
    let [a, b, c, d] = ids(&model);
    let topology = model.shared_topology();
    assert_eq!(
        DistanceEdit::new(&topology, [a, a]).unwrap_err(),
        GeometryEditError::RepeatedAtom(a)
    );
    assert_eq!(
        AngleEdit::new(&topology, [a, c, d]).unwrap_err(),
        GeometryEditError::MissingBond { a, b: c }
    );
    assert_eq!(
        DihedralEdit::new(&topology, [a, b, d, c]).unwrap_err(),
        GeometryEditError::MissingBond { a: b, b: d }
    );
    let larger = fixture("CCCCC", &[Point3::origin(); 5]);
    let invalid = larger.topology().atom_ids()[4];
    assert_eq!(
        DistanceEdit::new(&topology, [a, invalid]).unwrap_err(),
        GeometryEditError::InvalidAtomId(invalid)
    );
    for selected in [vec![], vec![a, d], vec![b, c]] {
        let moving = AtomSelection::from_atoms(&topology, selected).unwrap();
        assert_eq!(
            DihedralEdit::with_moving_atoms(&topology, [a, b, c, d], &moving).unwrap_err(),
            GeometryEditError::InvalidMovingSelection
        );
    }
    let different = chain();
    let foreign = AtomSelection::all(&different.shared_topology());
    assert_eq!(
        DistanceEdit::with_moving_atoms(&topology, [a, b], &foreign).unwrap_err(),
        GeometryEditError::TopologyMismatch
    );
    // Consecutive bonds need not be single or chemically rotatable.
    let unsaturated = fixture("C=CC#N", model.positions().values().value());
    assert!(DihedralEdit::new(&unsaturated.shared_topology(), ids(&unsaturated)).is_ok());
}

#[test]
fn prepared_edits_reject_independently_equal_topologies_atomically() {
    let model = chain();
    let [a, b, c, d] = ids(&model);
    let edit = DihedralEdit::new(&model.shared_topology(), [a, b, c, d]).unwrap();
    let mut other = chain();
    let baseline = other.clone();
    assert_eq!(
        edit.measure(other.as_model_view()),
        Err(GeometryEditError::TopologyMismatch)
    );
    assert_eq!(
        edit.apply(&mut other, Quantity::new(30.0, DEGREE)),
        Err(GeometryEditError::TopologyMismatch)
    );
    assert_eq!(other, baseline);
    let mut clone = model.clone();
    edit.apply(&mut clone, Quantity::new(30.0, DEGREE)).unwrap();
}

#[test]
fn invalid_targets_and_units_leave_model_unchanged() {
    let mut model = chain();
    let [a, b, c, d] = ids(&model);
    let baseline = model.clone();
    for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            model.set_distance(a, b, Quantity::new(target, NANOMETER)),
            Err(GeometryEditError::InvalidTarget)
        );
        assert_eq!(model, baseline);
    }
    for target in [-1.0, 181.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            model.set_angle(a, b, c, Quantity::new(target, DEGREE)),
            Err(GeometryEditError::InvalidTarget)
        );
        assert_eq!(model, baseline);
    }
    assert_eq!(
        model.set_dihedral(a, b, c, d, Quantity::new(f64::NAN, DEGREE)),
        Err(GeometryEditError::InvalidTarget)
    );
    assert!(matches!(
        model.set_angle(a, b, c, Quantity::new(1.0, NANOMETER)),
        Err(GeometryEditError::Unit(_))
    ));
    assert!(matches!(
        model.set_distance(a, b, Quantity::new(1.0, PICOSECOND)),
        Err(GeometryEditError::Unit(_))
    ));
    assert_eq!(model, baseline);
}

#[test]
fn degenerate_geometry_is_rejected_without_inventing_axes() {
    let mut model = chain();
    let [a, b, c, d] = ids(&model);
    model
        .set_position(a, Quantity::new(Point3::new(-1.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    let straight = model.clone();
    model
        .set_angle(a, b, c, Quantity::new(180.0, DEGREE))
        .unwrap();
    assert_eq!(model, straight);
    assert_eq!(
        model.set_angle(a, b, c, Quantity::new(90.0, DEGREE)),
        Err(GeometryEditError::DegenerateGeometry)
    );
    assert_eq!(
        model.set_dihedral(a, b, c, d, Quantity::new(0.0, RADIAN)),
        Err(GeometryEditError::DegenerateGeometry)
    );
    assert_eq!(model, straight);
    model.set_position(a, model.position(b).unwrap()).unwrap();
    let coincident = model.clone();
    assert_eq!(
        model.set_distance(a, b, Quantity::new(1.0, NANOMETER)),
        Err(GeometryEditError::DegenerateGeometry)
    );
    assert_eq!(
        model.set_angle(a, b, c, Quantity::new(90.0, DEGREE)),
        Err(GeometryEditError::DegenerateGeometry)
    );
    assert_eq!(model, coincident);
}

#[test]
fn numerical_failures_in_late_atoms_do_not_publish_earlier_atoms() {
    let mut model = chain();
    let [a, b, c, d] = ids(&model);
    model
        .set_position(d, Quantity::new(Point3::new(f64::MAX, 0.0, 0.0), NANOMETER))
        .unwrap();
    let baseline = model.clone();
    assert_eq!(
        model.set_distance(b, c, Quantity::new(f64::MAX, NANOMETER)),
        Err(GeometryEditError::NumericalFailure)
    );
    assert_eq!(model, baseline);
    let all = AtomSelection::all(&model.shared_topology());
    assert_eq!(
        model.translate(
            &all,
            Quantity::new(Vector3::new(f64::MAX, 0.0, 0.0), NANOMETER)
        ),
        Err(GeometryEditError::NumericalFailure)
    );
    assert_eq!(model, baseline);
    let transform =
        RigidTransform::new(Matrix3::identity(), Vector3::new(f64::MAX, 0.0, 0.0)).unwrap();
    assert_eq!(
        model.apply_transform(&all, &transform),
        Err(GeometryEditError::NumericalFailure)
    );
    assert_eq!(model, baseline);
    model
        .set_position(
            a,
            Quantity::new(Point3::new(-f64::MAX, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    let baseline = model.clone();
    assert_eq!(
        model.rotate(
            &all,
            Quantity::new(Point3::new(f64::MAX, 0.0, 0.0), NANOMETER),
            Vector3::new(0.0, 0.0, 1.0),
            Quantity::new(90.0, DEGREE)
        ),
        Err(GeometryEditError::NumericalFailure)
    );
    assert_eq!(model, baseline);
}

#[test]
fn finite_but_unachievable_distance_returns_numerical_failure() {
    let mut model = fixture(
        "CC",
        &[
            Point3::new(1.0e20, 0.0, 0.0),
            Point3::new(1.0e20 + 1.0e6, 0.0, 0.0),
        ],
    );
    let [a, b] = model.topology().atom_ids().try_into().unwrap();
    let baseline = model.clone();
    assert_eq!(
        model.set_distance(a, b, Quantity::new(1.0, NANOMETER)),
        Err(GeometryEditError::NumericalFailure)
    );
    assert_eq!(model, baseline);
}

#[test]
fn rigid_operations_use_units_right_handed_rotation_and_exact_selection() {
    let mut model = chain();
    let [a, b, c, d] = ids(&model);
    let baseline = model.clone();
    let selected = AtomSelection::from_atoms(&model.shared_topology(), [c, d]).unwrap();
    model
        .translate(
            &selected,
            Quantity::new(Vector3::new(10.0, 20.0, 30.0), ANGSTROM),
        )
        .unwrap();
    close_point(point(&model, c), Point3::new(2.0, 2.0, 3.0));
    model
        .rotate(
            &selected,
            Quantity::new(Point3::new(10.0, 20.0, 30.0), ANGSTROM),
            Vector3::new(0.0, 0.0, 8.0),
            Quantity::new(90.0, DEGREE),
        )
        .unwrap();
    close_point(point(&model, c), Point3::new(1.0, 3.0, 3.0));
    close_point(point(&model, d), Point3::new(0.0, 3.0, 4.0));
    for atom in [a, b] {
        assert_eq!(point(&model, atom), point(&baseline, atom));
    }
    let transform = RigidTransform::new(Matrix3::identity(), Vector3::new(1.0, 0.0, 0.0)).unwrap();
    model.apply_transform(&selected, &transform).unwrap();
    close_point(point(&model, c), Point3::new(2.0, 3.0, 3.0));
    // Both huge and subnormal finite axes are accepted.
    for scale in [f64::MAX, f64::from_bits(1)] {
        let mut sample = baseline.clone();
        sample
            .rotate(
                &selected,
                Quantity::new(Point3::origin(), NANOMETER),
                Vector3::new(0.0, 0.0, scale),
                Quantity::new(90.0, DEGREE),
            )
            .unwrap();
        close_point(point(&sample, c), Point3::new(0.0, 1.0, 0.0));
    }
}

#[test]
fn empty_selections_still_validate_binding_and_motion_arguments() {
    let mut model = chain();
    let baseline = model.clone();
    let empty = AtomSelection::empty(&model.shared_topology());
    let zero = Quantity::new(Point3::origin(), NANOMETER);
    let axis = Vector3::new(0.0, 0.0, 1.0);
    model
        .translate(
            &empty,
            Quantity::new(Vector3::new(1.0, 0.0, 0.0), NANOMETER),
        )
        .unwrap();
    model
        .rotate(&empty, zero, axis, Quantity::new(90.0, DEGREE))
        .unwrap();
    model
        .apply_transform(&empty, &RigidTransform::identity())
        .unwrap();
    assert_eq!(model, baseline);
    assert_eq!(
        model.rotate(&empty, zero, Vector3::zero(), Quantity::new(0.0, DEGREE)),
        Err(GeometryEditError::DegenerateGeometry)
    );
    assert_eq!(
        model.rotate(&empty, zero, axis, Quantity::new(f64::NAN, DEGREE)),
        Err(GeometryEditError::InvalidTarget)
    );
    assert!(matches!(
        model.rotate(&empty, zero, axis, Quantity::new(1.0, NANOMETER)),
        Err(GeometryEditError::Unit(_))
    ));
    assert!(matches!(
        model.translate(&empty, Quantity::new(Vector3::zero(), PICOSECOND)),
        Err(GeometryEditError::Unit(_))
    ));
    assert_eq!(
        model.translate(
            &empty,
            Quantity::new(Vector3::new(f64::INFINITY, 0.0, 0.0), NANOMETER)
        ),
        Err(GeometryEditError::NonFiniteInput)
    );
    assert_eq!(
        model.rotate(
            &empty,
            Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), NANOMETER),
            axis,
            Quantity::new(0.0, DEGREE)
        ),
        Err(GeometryEditError::NonFiniteInput)
    );
    assert_eq!(
        model.rotate(
            &empty,
            zero,
            Vector3::new(f64::NAN, 0.0, 0.0),
            Quantity::new(0.0, DEGREE)
        ),
        Err(GeometryEditError::NonFiniteInput)
    );
    let foreign = AtomSelection::empty(&chain().shared_topology());
    assert_eq!(
        model.translate(&foreign, Quantity::new(Vector3::zero(), NANOMETER)),
        Err(GeometryEditError::TopologyMismatch)
    );
    assert_eq!(
        model.rotate(&foreign, zero, axis, Quantity::new(0.0, DEGREE)),
        Err(GeometryEditError::TopologyMismatch)
    );
    assert_eq!(
        model.apply_transform(&foreign, &RigidTransform::identity()),
        Err(GeometryEditError::TopologyMismatch)
    );
    assert_eq!(model, baseline);
}

#[test]
fn all_edits_preserve_shared_topology_cell_properties_and_chemistry() {
    let mut model = chain();
    let [a, b, c, d] = ids(&model);
    model.conformation_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(5.0, 6.0, 7.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    model
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("test.note").unwrap(),
            PropertyValue::String("preserved".into()),
        )
        .unwrap();
    let a_index = model.topology().atom_index(a).unwrap();
    model
        .conformation_mut()
        .set_occupancy(a_index, Some(0.75))
        .unwrap();
    let baseline = model.clone();
    let all = AtomSelection::all(&model.shared_topology());
    model
        .set_distance(b, c, Quantity::new(2.0, NANOMETER))
        .unwrap();
    model
        .set_angle(a, b, c, Quantity::new(100.0, DEGREE))
        .unwrap();
    model
        .set_dihedral(a, b, c, d, Quantity::new(-45.0, DEGREE))
        .unwrap();
    model
        .translate(&all, Quantity::new(Vector3::new(1.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    model
        .rotate(
            &all,
            Quantity::new(Point3::origin(), NANOMETER),
            Vector3::new(0.0, 0.0, 1.0),
            Quantity::new(30.0, DEGREE),
        )
        .unwrap();
    model
        .apply_transform(&all, &RigidTransform::identity())
        .unwrap();
    assert!(Arc::ptr_eq(
        &model.shared_topology(),
        &baseline.shared_topology()
    ));
    assert_eq!(model.cell(), baseline.cell());
    assert_eq!(model.properties(), baseline.properties());
}

#[test]
fn geometry_edits_work_in_oblique_frames_and_at_different_length_scales() {
    for scale in [1.0e-6, 1.0, 1.0e6] {
        let mut model = chain();
        let [a, b, c, d] = ids(&model);
        let points = model
            .positions()
            .values()
            .value()
            .iter()
            .map(|p| Point3::new(p.x * scale, p.y * scale, p.z * scale))
            .collect::<Vec<_>>();
        model
            .conformation_mut()
            .set_positions(Quantity::new(points, NANOMETER))
            .unwrap();
        let all = AtomSelection::all(&model.shared_topology());
        model
            .rotate(
                &all,
                Quantity::new(Point3::origin(), NANOMETER),
                Vector3::new(1.0, 2.0, 3.0),
                Quantity::new(37.0, DEGREE),
            )
            .unwrap();
        model
            .translate(
                &all,
                Quantity::new(Vector3::new(7.0, -3.0, 1.0) * scale, NANOMETER),
            )
            .unwrap();
        let before = model.clone();
        model
            .set_distance(b, c, Quantity::new(1.75 * scale, NANOMETER))
            .unwrap();
        close(
            measure::distance(model.as_model_view(), b, c)
                .unwrap()
                .into_value()
                / scale,
            1.75,
        );
        model
            .set_angle(a, b, c, Quantity::new(47.0, DEGREE))
            .unwrap();
        close(
            measure::angle(model.as_model_view(), a, b, c)
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            47.0,
        );
        model
            .set_dihedral(a, b, c, d, Quantity::new(-139.0, DEGREE))
            .unwrap();
        close(
            measure::dihedral(model.as_model_view(), a, b, c, d)
                .unwrap()
                .value_in(DEGREE)
                .unwrap(),
            -139.0,
        );
        assert_eq!(point(&model, a), point(&before, a));
        assert_eq!(point(&model, b), point(&before, b));
        close(
            (point(&model, c) - point(&model, d)).norm() / scale,
            2.0_f64.sqrt(),
        );
    }
}

#[test]
fn automatic_fragments_may_contain_a_ring_on_the_moving_side() {
    let model = fixture("CCC1CC1", &[Point3::origin(); 5]);
    let atoms = model.topology().atom_ids();
    let edit = AngleEdit::new(&model.shared_topology(), [atoms[0], atoms[1], atoms[2]]).unwrap();
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
    let mut model = Model::new(
        builder.build().unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap();
    let topology = model.shared_topology();
    let [a, b, c, d] = model.topology().atom_ids()[..4].try_into().unwrap();
    let mut moving = AtomSelection::from_atoms(&topology, [d]).unwrap();
    let edit = DihedralEdit::with_moving_atoms(&topology, [a, b, c, d], &moving).unwrap();
    moving.clear();
    assert_eq!(edit.moving_atoms().len(), 1);
    let before = model.clone();
    edit.apply(&mut model, Quantity::new(-90.0, DEGREE))
        .unwrap();
    model
        .set_distance(b, c, Quantity::new(2.0, NANOMETER))
        .unwrap();
    for &atom in &model.topology().atom_ids()[4..] {
        assert_eq!(point(&model, atom), point(&before, atom));
    }
}

#[test]
fn apply_transform_uses_rotation_and_translation_together() {
    let mut model = chain();
    let baseline = model.clone();
    let [a, b, c, d] = ids(&model);
    let selected = AtomSelection::from_atoms(&model.shared_topology(), [c, d]).unwrap();
    let transform = RigidTransform::new(
        Matrix3::from_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ),
        Vector3::new(2.0, 3.0, 4.0),
    )
    .unwrap();
    model.apply_transform(&selected, &transform).unwrap();
    close_point(point(&model, c), Point3::new(2.0, 4.0, 4.0));
    close_point(point(&model, d), Point3::new(1.0, 4.0, 5.0));
    for atom in [a, b] {
        assert_eq!(point(&model, atom), point(&baseline, atom));
    }
}
