use super::*;
use crate::math::PeriodicCell;
use crate::molecule::{Atom, BondOrder, Element};
use crate::structure::{Positions, Structure};
use crate::topology::{AtomSelection, Topology, TopologyBuilder};
use crate::units::{Quantity, NANOMETER};
use std::sync::Arc;

fn topology(atom_count: usize) -> Arc<Topology> {
    let mut graph = crate::molecule::MoleculeEditor::new();
    let mut previous = None;
    for _ in 0..atom_count {
        let atom = graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(parent) = previous {
            graph.add_bond(parent, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    let molecule = graph.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

fn structure(topology: &Arc<Topology>, points: &[Point3]) -> Structure {
    let positions = Positions::new(Quantity::new(points, NANOMETER)).unwrap();
    Structure::new(Arc::clone(topology), positions).unwrap()
}

fn structures(moving: &[Point3], reference: &[Point3]) -> (Arc<Topology>, Structure, Structure) {
    assert_eq!(moving.len(), reference.len());
    let topology = topology(moving.len());
    let moving = structure(&topology, moving);
    let reference = structure(&topology, reference);
    (topology, moving, reference)
}

fn all(topology: &Arc<Topology>) -> AtomSelection {
    AtomSelection::from_atoms(topology.clone(), topology.atom_ids().iter().copied()).unwrap()
}

fn fixture_points() -> [Point3; 5] {
    [
        Point3::new(-1.2, 0.4, 2.1),
        Point3::new(0.7, -1.5, 0.3),
        Point3::new(2.4, 0.8, -0.6),
        Point3::new(-0.3, 2.2, 1.1),
        Point3::new(1.3, 1.7, 2.8),
    ]
}

fn non_axis_transform() -> RigidTransform {
    let axis = Vector3::new(1.0, -2.0, 3.0);
    let axis = axis / axis.norm();
    let angle = 0.731_f64;
    let cosine = angle.cos();
    let sine = angle.sin();
    let one_minus_cosine = 1.0 - cosine;
    let rotation = Matrix3::from_columns(
        Vector3::new(
            cosine + axis.x * axis.x * one_minus_cosine,
            axis.y * axis.x * one_minus_cosine + axis.z * sine,
            axis.z * axis.x * one_minus_cosine - axis.y * sine,
        ),
        Vector3::new(
            axis.x * axis.y * one_minus_cosine - axis.z * sine,
            cosine + axis.y * axis.y * one_minus_cosine,
            axis.z * axis.y * one_minus_cosine + axis.x * sine,
        ),
        Vector3::new(
            axis.x * axis.z * one_minus_cosine + axis.y * sine,
            axis.y * axis.z * one_minus_cosine - axis.x * sine,
            cosine + axis.z * axis.z * one_minus_cosine,
        ),
    );
    RigidTransform::from_canonical(rotation, Vector3::new(2.3, -4.1, 0.8)).unwrap()
}

fn transform_points(points: &[Point3], transform: RigidTransform) -> Vec<Point3> {
    points
        .iter()
        .copied()
        .map(|point| transform.transform_point(point))
        .collect()
}

fn assert_close(left: f64, right: f64, tolerance: f64) {
    assert!(
        (left - right).abs() <= tolerance,
        "{left} is not within {tolerance} of {right}"
    );
}

fn assert_point_close(left: Point3, right: Point3, tolerance: f64) {
    assert_close(left.x, right.x, tolerance);
    assert_close(left.y, right.y, tolerance);
    assert_close(left.z, right.z, tolerance);
}

fn assert_transform_close(left: RigidTransform, right: RigidTransform, tolerance: f64) {
    for (left, right) in left
        .rotation()
        .columns()
        .into_iter()
        .zip(right.rotation().columns())
    {
        assert_close(left.x, right.x, tolerance);
        assert_close(left.y, right.y, tolerance);
        assert_close(left.z, right.z, tolerance);
    }
    let left = left.translation().into_value();
    let right = right.translation().into_value();
    assert_close(left.x, right.x, tolerance);
    assert_close(left.y, right.y, tolerance);
    assert_close(left.z, right.z, tolerance);
}

#[test]
fn identity_fit_reports_canonical_units_and_preserves_inputs() {
    let points = fixture_points();
    let (topology, moving, reference) = structures(&points, &points);
    let selection = all(&topology);
    let moving_before = moving.positions().values().value().to_vec();
    let reference_before = reference.positions().values().value().to_vec();
    let selection_before = selection.clone();

    let result = crate::alignment::fit(
        moving.as_structure_view(),
        reference.as_structure_view(),
        &selection,
        Default::default(),
    )
    .unwrap();

    assert_transform_close(result.transform(), RigidTransform::identity(), 1.0e-12);
    assert_close(result.rmsd().into_value(), 0.0, 1.0e-14);
    assert_eq!(result.rmsd().unit(), CANONICAL_LENGTH_UNIT);
    assert_eq!(result.selected_atom_count(), points.len());
    assert_eq!(
        moving.positions().values().into_value(),
        moving_before.as_slice()
    );
    assert_eq!(
        reference.positions().values().into_value(),
        reference_before.as_slice()
    );
    assert_eq!(selection, selection_before);
}

#[test]
fn translation_fit_maps_moving_to_reference_not_the_inverse() {
    let moving = fixture_points();
    let expected =
        RigidTransform::from_canonical(Matrix3::identity(), Vector3::new(3.2, -1.7, 4.4)).unwrap();
    let reference = transform_points(&moving, expected);
    let (topology, moving_structure, reference_structure) = structures(&moving, &reference);

    let result = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &all(&topology),
        Default::default(),
    )
    .unwrap();

    assert_transform_close(result.transform(), expected, 1.0e-12);
    for (moving, reference) in moving.into_iter().zip(reference) {
        assert_point_close(
            result.transform().transform_point(moving),
            reference,
            1.0e-11,
        );
        assert!((result.transform().transform_point(reference).x - moving.x).abs() > 1.0);
    }
}

#[test]
fn asymmetric_non_axis_rotation_and_translation_are_recovered() {
    let moving = fixture_points();
    let expected = non_axis_transform();
    let reference = transform_points(&moving, expected);
    let (topology, moving, reference_structure) = structures(&moving, &reference);

    let result = crate::alignment::fit(
        moving.as_structure_view(),
        reference_structure.as_structure_view(),
        &all(&topology),
        Default::default(),
    )
    .unwrap();

    assert_transform_close(result.transform(), expected, 2.0e-12);
    assert!(result.rmsd().into_value() < 2.0e-12);
}

#[test]
fn subset_fit_ignores_unselected_atoms() {
    let moving = fixture_points();
    let expected = non_axis_transform();
    let mut reference = transform_points(&moving, expected);
    reference[3] = Point3::new(100.0, -80.0, 40.0);
    reference[4] = Point3::new(-70.0, 60.0, -30.0);
    let (topology, moving_structure, reference_structure) = structures(&moving, &reference);
    let selection =
        AtomSelection::from_atoms(topology.clone(), topology.atom_ids()[..3].iter().copied())
            .unwrap();

    let result = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        Default::default(),
    )
    .unwrap();

    assert_eq!(result.selected_atom_count(), 3);
    for index in 0..3 {
        assert_point_close(
            result.transform().transform_point(moving[index]),
            reference[index],
            3.0e-11,
        );
    }
}

#[test]
fn noisy_fit_has_expected_nonzero_improved_rmsd() {
    let moving = fixture_points();
    let expected = non_axis_transform();
    let mut reference = transform_points(&moving, expected);
    reference[0].x += 0.08;
    reference[1].y -= 0.05;
    reference[4].z += 0.03;
    let (topology, moving_structure, reference_structure) = structures(&moving, &reference);
    let selection = all(&topology);
    let result = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        Default::default(),
    )
    .unwrap();
    let unaligned = moving
        .iter()
        .zip(&reference)
        .map(|(moving, reference)| (*moving - *reference).norm_squared())
        .sum::<f64>()
        / moving.len() as f64;
    let unaligned = unaligned.sqrt();

    assert!(result.rmsd().into_value() > 0.01);
    assert!(result.rmsd().into_value() < 0.08);
    assert!(result.rmsd().into_value() < unaligned);
}

#[test]
fn explicit_weights_change_fit_and_common_scale_does_not() {
    let moving = fixture_points();
    let expected = non_axis_transform();
    let mut reference = transform_points(&moving, expected);
    reference[4].x += 2.5;
    reference[4].y -= 1.5;
    let (topology, moving_structure, reference_structure) = structures(&moving, &reference);
    let selection = all(&topology);
    let uniform = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        Default::default(),
    )
    .unwrap();
    let weights = [10.0, 10.0, 10.0, 10.0, 0.1];
    let scaled_weights = [30.0, 30.0, 30.0, 30.0, 0.3];
    let weighted = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        AlignmentOptions {
            weighting: Weighting::Explicit(&weights),
            ..AlignmentOptions::default()
        },
    )
    .unwrap();
    let scaled = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        AlignmentOptions {
            weighting: Weighting::Explicit(&scaled_weights),
            ..AlignmentOptions::default()
        },
    )
    .unwrap();

    let uniform_error = (uniform.transform().transform_point(moving[0]) - reference[0]).norm();
    let weighted_error = (weighted.transform().transform_point(moving[0]) - reference[0]).norm();
    assert!(weighted_error < uniform_error);
    assert_transform_close(weighted.transform(), scaled.transform(), 2.0e-12);
    assert_close(
        weighted.rmsd().into_value(),
        scaled.rmsd().into_value(),
        2.0e-12,
    );
}

#[test]
fn explicit_weights_follow_selection_indices_order_not_semantic_id_insertion_order() {
    let moving = fixture_points();
    let expected = non_axis_transform();
    let mut reference = transform_points(&moving, expected);
    reference[0].x += 1.7;
    reference[4].x -= 1.9;
    let (topology, moving_structure, reference_structure) = structures(&moving, &reference);
    let atom_ids = topology.atom_ids();
    let selection = AtomSelection::from_atoms(
        topology.clone(),
        [atom_ids[4], atom_ids[1], atom_ids[3], atom_ids[0]],
    )
    .unwrap();

    assert_eq!(
        selection
            .indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>(),
        [0, 1, 3, 4]
    );

    // These weights make dense index 0 dominant in selection.indices()
    // order. Reversing the endpoint weights instead makes the first
    // inserted semantic ID (dense index 4) dominant.
    let selection_order_weights = [100.0, 1.0, 1.0, 0.1];
    let insertion_order_interpretation = [0.1, 1.0, 1.0, 100.0];
    let selection_order_fit = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        AlignmentOptions {
            weighting: Weighting::Explicit(&selection_order_weights),
            ..AlignmentOptions::default()
        },
    )
    .unwrap();
    let insertion_order_fit = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &selection,
        AlignmentOptions {
            weighting: Weighting::Explicit(&insertion_order_interpretation),
            ..AlignmentOptions::default()
        },
    )
    .unwrap();

    let selected_dense_zero_error =
        (selection_order_fit.transform().transform_point(moving[0]) - reference[0]).norm();
    let insertion_dense_zero_error =
        (insertion_order_fit.transform().transform_point(moving[0]) - reference[0]).norm();
    let selected_dense_four_error =
        (selection_order_fit.transform().transform_point(moving[4]) - reference[4]).norm();
    let insertion_dense_four_error =
        (insertion_order_fit.transform().transform_point(moving[4]) - reference[4]).norm();
    assert!(selected_dense_zero_error < insertion_dense_zero_error * 0.1);
    assert!(insertion_dense_four_error < selected_dense_four_error * 0.1);
}

#[test]
fn weighted_noisy_fit_matches_high_precision_svd_golden() {
    let moving = [
        Point3::new(-1.4, 0.2, 2.1),
        Point3::new(0.3, -1.7, 0.5),
        Point3::new(2.2, 0.9, -0.8),
        Point3::new(-0.6, 2.4, 1.3),
        Point3::new(1.5, 1.1, 2.7),
        Point3::new(-2.0, -0.9, -1.2),
    ];
    let reference = [
        Point3::new(1.873307692308, -4.135076923077, -0.201230769231),
        Point3::new(3.672153846154, -2.542538461538, 1.554384615385),
        Point3::new(1.272692307692, -0.149923076923, 2.197230769231),
        Point3::new(-0.144923076923, -2.663769230769, -0.328307692308),
        Point3::new(0.199692307692, -3.570923076923, 2.290230769231),
        Point3::new(4.000538461538, -1.492384615385, -1.161153846154),
    ];
    let weights = [0.7, 3.25, 1.1, 5.5, 0.4, 2.75];
    let (topology, moving_structure, reference_structure) = structures(&moving, &reference);

    // Independent reference: mpmath 1.3.0 at 100 decimal digits, using
    // weighted centroids/cross-covariance and its general SVD, followed by
    // V * diag(1, 1, det(V * U^T)) * U^T. The 3e-12 transform and 5e-14
    // RMSD tolerances cover only f64 rounding relative to that reference.
    let expected = RigidTransform::from_canonical(
        Matrix3::from_columns(
            Vector3::new(
                -0.23369675064066617,
                0.30951178758237726,
                0.9217311332962317,
            ),
            Vector3::new(
                -0.9253543656919893,
                0.22023658458981056,
                -0.30856951356702656,
            ),
            Vector3::new(-0.2985048184448125, -0.925039620857278, 0.2349395096816475),
        ),
        Vector3::new(2.333984523603927, -1.795397432011596, 0.6560333319022371),
    )
    .unwrap();
    let expected_rmsd = 0.029397741136029633_f64;

    let result = crate::alignment::fit(
        moving_structure.as_structure_view(),
        reference_structure.as_structure_view(),
        &all(&topology),
        AlignmentOptions {
            weighting: Weighting::Explicit(&weights),
            ..AlignmentOptions::default()
        },
    )
    .unwrap();

    assert_transform_close(result.transform(), expected, 3.0e-12);
    assert_close(result.rmsd().into_value(), expected_rmsd, 5.0e-14);
    assert_close(result.transform().rotation().determinant(), 1.0, 1.0e-12);
}

#[test]
fn mirrored_coordinates_keep_a_proper_nonzero_residual() {
    let moving = fixture_points();
    let reference = moving.map(|point| Point3::new(-point.x, point.y, point.z));
    let (topology, moving, reference) = structures(&moving, &reference);

    let result = crate::alignment::fit(
        moving.as_structure_view(),
        reference.as_structure_view(),
        &all(&topology),
        Default::default(),
    )
    .unwrap();

    assert_close(result.transform().rotation().determinant(), 1.0, 1.0e-12);
    assert!(result.rmsd().into_value() > 0.1);
    RigidTransform::new(
        result.transform().rotation(),
        result.transform().translation(),
    )
    .unwrap();
}

#[test]
fn planar_non_collinear_selection_succeeds() {
    let moving = [
        Point3::new(-1.0, -0.5, 0.0),
        Point3::new(2.0, -0.5, 0.0),
        Point3::new(0.2, 1.7, 0.0),
        Point3::new(-0.4, 0.8, 0.0),
    ];
    let expected = non_axis_transform();
    let reference = transform_points(&moving, expected);
    let (topology, moving, reference) = structures(&moving, &reference);

    let result = crate::alignment::fit(
        moving.as_structure_view(),
        reference.as_structure_view(),
        &all(&topology),
        Default::default(),
    )
    .unwrap();

    assert_transform_close(result.transform(), expected, 3.0e-12);
    assert!(result.rmsd().into_value() < 2.0e-12);
}

#[test]
fn rank_two_point_sets_with_rank_one_cross_covariance_are_rejected() {
    let moving = [
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
    ];
    let reference = [
        Point3::new(1.0, 1.0, 0.0),
        Point3::new(-1.0, 1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
        Point3::new(0.0, -1.0, 0.0),
    ];
    let (topology, moving, reference) = structures(&moving, &reference);

    assert_eq!(
        crate::alignment::fit(
            moving.as_structure_view(),
            reference.as_structure_view(),
            &all(&topology),
            Default::default()
        ),
        Err(AlignmentError::DegenerateGeometry {
            geometry: AlignmentGeometry::CrossCovariance,
        })
    );
}

#[test]
fn fewer_than_three_selected_atoms_are_structured_failures() {
    let points = fixture_points();
    let (topology, moving, reference) = structures(&points, &points);
    for count in 0..3 {
        let selection = AtomSelection::from_atoms(
            topology.clone(),
            topology.atom_ids()[..count].iter().copied(),
        )
        .unwrap();
        assert_eq!(
            crate::alignment::fit(
                moving.as_structure_view(),
                reference.as_structure_view(),
                &selection,
                Default::default()
            ),
            Err(AlignmentError::InsufficientSelectedAtoms {
                selected: count,
                minimum: 3,
            })
        );
    }
}

#[test]
fn coincident_and_collinear_geometry_are_rejected() {
    let valid = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ];
    let coincident = [Point3::new(2.0, -1.0, 0.5); 4];
    let collinear = [
        Point3::new(-2.0, 0.0, 0.0),
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
    ];
    let (topology, moving, reference) = structures(&coincident, &valid);
    assert_eq!(
        crate::alignment::fit(
            moving.as_structure_view(),
            reference.as_structure_view(),
            &all(&topology),
            Default::default()
        ),
        Err(AlignmentError::DegenerateGeometry {
            geometry: AlignmentGeometry::Moving,
        })
    );
    let (topology, moving, reference) = structures(&valid, &collinear);
    assert_eq!(
        crate::alignment::fit(
            moving.as_structure_view(),
            reference.as_structure_view(),
            &all(&topology),
            Default::default()
        ),
        Err(AlignmentError::DegenerateGeometry {
            geometry: AlignmentGeometry::Reference,
        })
    );
}

#[test]
fn near_collinear_rank_test_is_scale_relative_and_thresholded() {
    let points = |height: f64, scale: f64| {
        [
            Point3::new(-scale, 0.0, 0.0),
            Point3::new(scale, 0.0, 0.0),
            Point3::new(0.0, height * scale, 0.0),
            Point3::new(0.5 * scale, -height * scale, 0.0),
        ]
    };
    for scale in [1.0, 1.0e8] {
        let accepted = points(1.0e-5, scale);
        let (topology, moving, reference) = structures(&accepted, &accepted);
        assert!(crate::alignment::fit(
            moving.as_structure_view(),
            reference.as_structure_view(),
            &all(&topology),
            Default::default()
        )
        .is_ok());

        let rejected = points(1.0e-7, scale);
        let (topology, moving, reference) = structures(&rejected, &rejected);
        assert!(matches!(
            crate::alignment::fit(
                moving.as_structure_view(),
                reference.as_structure_view(),
                &all(&topology),
                Default::default()
            ),
            Err(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Moving
            })
        ));
    }
}

#[test]
fn exact_shared_topology_is_required_for_views_and_selection() {
    let points = fixture_points();
    let topology_a = topology(points.len());
    let topology_b = topology(points.len());
    assert!(topology_a.same_layout(&topology_b));
    assert!(!Arc::ptr_eq(&topology_a, &topology_b));
    let moving_a = structure(&topology_a, &points);
    let reference_a = structure(&topology_a, &points);
    let reference_b = structure(&topology_b, &points);
    assert_eq!(
        crate::alignment::fit(
            moving_a.as_structure_view(),
            reference_b.as_structure_view(),
            &all(&topology_a),
            Default::default()
        ),
        Err(AlignmentError::TopologyMismatch(
            TopologyMismatch::Reference
        ))
    );
    assert_eq!(
        crate::alignment::fit(
            moving_a.as_structure_view(),
            reference_a.as_structure_view(),
            &all(&topology_b),
            Default::default()
        ),
        Err(AlignmentError::TopologyMismatch(
            TopologyMismatch::Selection
        ))
    );
}

#[test]
fn invalid_explicit_weights_are_structured_failures() {
    let points = fixture_points();
    let (topology, moving, reference) = structures(&points, &points);
    let selection = all(&topology);
    let fit = |weights: &[f64]| {
        crate::alignment::fit(
            moving.as_structure_view(),
            reference.as_structure_view(),
            &selection,
            AlignmentOptions {
                weighting: Weighting::Explicit(weights),
                ..AlignmentOptions::default()
            },
        )
    };
    assert_eq!(
        fit(&[1.0]),
        Err(AlignmentError::WeightCountMismatch(LengthMismatch::new(
            Countable::Weight,
            points.len(),
            1
        )))
    );
    for (weight, expected) in [
        (
            0.0,
            AlignmentError::NonPositiveWeight { selection_index: 2 },
        ),
        (
            -1.0,
            AlignmentError::NonPositiveWeight { selection_index: 2 },
        ),
        (
            f64::NAN,
            AlignmentError::NonFiniteWeight { selection_index: 2 },
        ),
        (
            f64::INFINITY,
            AlignmentError::NonFiniteWeight { selection_index: 2 },
        ),
    ] {
        let weights = [1.0, 1.0, weight, 1.0, 1.0];
        assert_eq!(fit(&weights), Err(expected));
    }
}

#[test]
fn periodic_policy_fits_stored_coordinates_by_default_and_can_reject() {
    let moving_points = fixture_points();
    let expected = non_axis_transform();
    let reference_points = transform_points(&moving_points, expected);
    let (topology, mut moving, mut reference) = structures(&moving_points, &reference_points);
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(3.0, 4.0, 5.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    moving.geometry_mut().set_cell(Some(cell));
    reference.geometry_mut().set_cell(Some(cell));
    let selection = all(&topology);
    let moving_before = moving.positions().values().value().to_vec();

    assert_eq!(
        crate::alignment::fit(
            moving.as_structure_view(),
            reference.as_structure_view(),
            &selection,
            AlignmentOptions {
                periodic_policy: PeriodicPolicy::RejectPeriodic,
                ..Default::default()
            }
        ),
        Err(AlignmentError::PeriodicCoordinates {
            moving: true,
            reference: true,
        })
    );
    let result = crate::alignment::fit(
        moving.as_structure_view(),
        reference.as_structure_view(),
        &selection,
        Default::default(),
    )
    .unwrap();
    assert_transform_close(result.transform(), expected, 3.0e-12);
    assert_eq!(moving.cell(), Some(&cell));
    assert_eq!(reference.cell(), Some(&cell));
    assert_eq!(
        moving.positions().values().into_value(),
        moving_before.as_slice()
    );
}

#[test]
fn large_coordinate_offsets_retain_small_internal_geometry() {
    let offset = Vector3::new(1.0e12, -2.0e12, 3.0e12);
    let moving = fixture_points().map(|point| point + offset);
    let rotation = non_axis_transform().rotation();
    let expected =
        RigidTransform::from_canonical(rotation, Vector3::new(-4.0e11, 7.0e11, 2.0e11)).unwrap();
    let reference = transform_points(&moving, expected);
    let (topology, moving, reference_structure) = structures(&moving, &reference);

    let result = crate::alignment::fit(
        moving.as_structure_view(),
        reference_structure.as_structure_view(),
        &all(&topology),
        Default::default(),
    )
    .unwrap();

    assert_close(result.transform().rotation().determinant(), 1.0, 1.0e-12);
    for (moving, reference) in moving.positions().values().value().iter().zip(reference) {
        assert_point_close(
            result.transform().transform_point(*moving),
            reference,
            2.0e-3,
        );
    }
    assert!(result.rmsd().into_value() < 1.0e-3);
}
