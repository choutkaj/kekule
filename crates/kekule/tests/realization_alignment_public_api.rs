//! Superposition and RMSD over realization collections (ensembles and
//! trajectories), against collection items or independent references.

use std::sync::Arc;

use kekule::alignment::{
    kabsch, AlignedRmsdOptions, AlignmentError, AlignmentGeometry, AlignmentOptions,
    AtomCorrespondence, AtomCorrespondenceError, CorrespondenceSide, PeriodicPolicy, Reference,
    Weighting,
};
use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::geometry::{Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Ensemble, EnsembleMember, Forces, Model, Positions, Realization, Realizations, Trajectory,
    TrajectoryFrame, Velocities,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex};
use kekule::units::{
    Quantity, ANGSTROM, CANONICAL_FORCE_UNIT, CANONICAL_LENGTH_UNIT, CANONICAL_VELOCITY_UNIT,
    DIMENSIONLESS, NANOMETER, PICOSECOND,
};

fn chain(atoms: usize) -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    let mut previous = None;
    for _ in 0..atoms {
        let atom = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(previous) = previous {
            editor.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    Arc::new(Topology::from_molecule(editor.finish().unwrap()).unwrap())
}

fn selection(topology: &Arc<Topology>, indices: &[usize]) -> AtomSelection {
    AtomSelection::from_atoms(
        topology,
        indices.iter().map(|index| topology.atom_ids()[*index]),
    )
    .unwrap()
}

fn frame(points: &[Point3]) -> TrajectoryFrame {
    TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
}

fn quarter_turn() -> RigidTransform {
    RigidTransform::new(
        Matrix3::from_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ),
        Vector3::new(4.0, -2.0, 1.0),
    )
    .unwrap()
}

fn transformed(points: &[Point3], transform: RigidTransform) -> Vec<Point3> {
    points
        .iter()
        .map(|point| transform.transform_point(*point))
        .collect()
}

fn spread() -> [Point3; 4] {
    [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 1.5, 0.0),
        Point3::new(0.5, 0.4, 1.2),
    ]
}

fn collinear() -> [Point3; 4] {
    [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(3.0, 0.0, 0.0),
    ]
}

fn payloads<P: Realization + Clone>(collection: &Realizations<P>) -> Vec<P> {
    collection
        .iter()
        .map(|item| item.payload().clone())
        .collect()
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "expected {expected}, received {actual}"
    );
}

fn assert_vector_close(actual: Vector3, expected: Vector3, tolerance: f64) {
    assert!(
        (actual - expected).norm() <= tolerance,
        "expected {expected:?}, received {actual:?}"
    );
}

fn periodic_cell(edge: f64) -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(edge, edge, edge), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn direct_rmsd_does_not_fit_and_explicit_weights_follow_pair_order() {
    let topology = chain(2);
    let reference = [Point3::origin(), Point3::origin()];
    let moving = [Point3::new(1.0, 0.0, 0.0), Point3::new(3.0, 0.0, 0.0)];
    let trajectory =
        Trajectory::from_items(Arc::clone(&topology), [frame(&reference), frame(&moving)]).unwrap();
    let all = AtomSelection::all(&topology);

    let uniform = trajectory.rmsd(0, &all).unwrap();
    assert_eq!(uniform.unit(), CANONICAL_LENGTH_UNIT);
    assert_close(uniform.value()[0], 0.0, 1.0e-14);
    assert_close(uniform.value()[1], 5.0_f64.sqrt(), 1.0e-14);

    let weighted = trajectory
        .rmsd_with_options(
            0,
            &all,
            AlignmentOptions {
                weighting: Weighting::Explicit(&[3.0, 1.0]),
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
    assert_close(weighted.value()[1], 3.0_f64.sqrt(), 1.0e-14);
    let scaled = trajectory
        .rmsd_with_options(
            0,
            &all,
            AlignmentOptions {
                weighting: Weighting::Explicit(&[30.0, 10.0]),
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
    assert_eq!(scaled, weighted);
}

#[test]
fn split_and_fused_fit_measure_workflows_agree_for_distinct_selections() {
    let topology = chain(4);
    let moving = spread();
    let mut reference = transformed(&moving, quarter_turn());
    reference[3].x += 2.0;
    let trajectory =
        Trajectory::from_items(Arc::clone(&topology), [frame(&reference), frame(&moving)]).unwrap();
    let fit = selection(&topology, &[0, 1, 2]);
    let measured = selection(&topology, &[3]);

    let fused = trajectory.aligned_rmsd(0, &fit, &measured).unwrap();
    assert_close(fused.value()[0], 0.0, 1.0e-12);
    assert_close(fused.value()[1], 2.0, 2.0e-12);

    let mut split = trajectory.clone();
    let report = split.superpose(0, &fit).unwrap();
    assert_eq!(report.reference_index(), Some(0));
    assert_eq!(report.len(), 2);
    assert_close(
        report.alignment(1).unwrap().rmsd().into_value(),
        0.0,
        1.0e-12,
    );
    let split = split.rmsd(0, &measured).unwrap();
    for (split, fused) in split.value().iter().zip(fused.value()) {
        assert_close(*split, *fused, 2.0e-12);
    }
}

#[test]
fn superposition_rotates_complete_geometric_state_and_preserves_metadata() {
    let topology = chain(4);
    let moving = spread();
    let transform = quarter_turn();
    let reference = transformed(&moving, transform);
    let vectors = [
        Vector3::new(10.0, 0.0, 0.0),
        Vector3::new(0.0, 20.0, 0.0),
        Vector3::new(0.0, 0.0, 30.0),
    ];
    let moving_cell = PeriodicCell::new(Quantity::new(vectors, ANGSTROM), [true; 3]).unwrap();
    let reference_cell = PeriodicCell::new(
        Quantity::new(
            vectors.map(|vector| transform.transform_vector(vector)),
            ANGSTROM,
        ),
        [true; 3],
    )
    .unwrap();

    let mut reference_frame = frame(&reference);
    reference_frame
        .conformation_mut()
        .set_cell(Some(reference_cell));
    let mut moving_frame = frame(&moving);
    moving_frame.conformation_mut().set_cell(Some(moving_cell));
    moving_frame
        .set_velocities(Some(
            Velocities::new(Quantity::new(
                vec![Vector3::new(1.0, 0.0, 0.0); 4],
                CANONICAL_VELOCITY_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
    moving_frame
        .set_forces(Some(
            Forces::new(Quantity::new(
                vec![Vector3::new(0.0, 1.0, 0.0); 4],
                CANONICAL_FORCE_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
    moving_frame
        .set_time(Some(Quantity::new(2.5, PICOSECOND)))
        .unwrap();
    moving_frame.set_step(Some(25));
    let score = PropertyKey::new("score").unwrap();
    let label = PropertyKey::new("label").unwrap();
    let mut conformation = moving_frame.conformation_mut();
    let mut properties = conformation.properties_mut();
    properties
        .atoms_mut()
        .set_value(
            score.clone(),
            TopologyAtomIndex::new(0),
            Some(PropertyValue::Real {
                value: 0.7,
                unit: DIMENSIONLESS,
            }),
        )
        .unwrap();
    properties
        .owner_mut()
        .insert(label.clone(), PropertyValue::String("moving".into()))
        .unwrap();
    let mut trajectory =
        Trajectory::from_items(Arc::clone(&topology), [reference_frame, moving_frame]).unwrap();
    let expected_properties = trajectory.get(1).unwrap().properties().clone();

    let report = trajectory
        .superpose(0, &AtomSelection::all(&topology))
        .unwrap();
    assert_eq!(report.alignments().len(), 2);
    let aligned = trajectory.get(1).unwrap();
    for (actual, expected) in aligned
        .positions()
        .values()
        .value()
        .iter()
        .zip(reference.iter())
    {
        assert!((*actual - *expected).norm() <= 3.0e-12);
    }
    let cell = aligned.cell().copied().unwrap();
    assert_eq!(cell.periodic_axes(), reference_cell.periodic_axes());
    for (actual, expected) in cell
        .vectors()
        .into_value()
        .into_iter()
        .zip(reference_cell.vectors().into_value())
    {
        assert_vector_close(actual, expected, 1.0e-12);
    }
    for velocity in aligned.velocities().unwrap().values().value().iter() {
        assert_vector_close(*velocity, Vector3::new(0.0, 1.0, 0.0), 2.0e-12);
    }
    for force in aligned.forces().unwrap().values().value().iter() {
        assert_vector_close(*force, Vector3::new(-1.0, 0.0, 0.0), 2.0e-12);
    }
    assert_eq!(aligned.time(), Some(Quantity::new(2.5, PICOSECOND)));
    assert_eq!(aligned.step(), Some(25));
    assert_eq!(aligned.properties(), &expected_properties);
    assert_eq!(
        aligned.properties().owner().get(&label),
        Some(&PropertyValue::String("moving".into()))
    );
}

#[test]
fn superposition_failure_leaves_the_complete_collection_unchanged() {
    let topology = chain(4);
    let reference = spread();
    let moving = transformed(&reference, quarter_turn());
    let mut trajectory = Trajectory::from_items(
        Arc::clone(&topology),
        [frame(&reference), frame(&moving), frame(&collinear())],
    )
    .unwrap();
    let before = payloads(&trajectory);
    let all = AtomSelection::all(&topology);

    assert_eq!(
        trajectory.superpose(3, &all),
        Err(AlignmentError::ReferenceOutOfRange { index: 3, len: 3 })
    );
    // The last item fails after earlier items fitted; nothing is published.
    assert_eq!(
        trajectory.superpose(0, &all),
        Err(AlignmentError::Item {
            index: 2,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Moving,
            }),
        })
    );
    assert_eq!(
        trajectory.superpose(0, &selection(&topology, &[0, 1])),
        Err(AlignmentError::InsufficientSelectedAtoms {
            selected: 2,
            minimum: 3,
        })
    );
    assert_eq!(payloads(&trajectory), before);
}

#[test]
fn rmsd_reports_reference_selection_weight_and_periodic_failures() {
    let topology = chain(4);
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
    ];
    let mut periodic = frame(&points);
    periodic
        .conformation_mut()
        .set_cell(Some(periodic_cell(1.0)));
    let trajectory = Trajectory::from_items(Arc::clone(&topology), [periodic]).unwrap();
    let all = AtomSelection::all(&topology);

    assert_eq!(
        trajectory.rmsd(1, &all),
        Err(AlignmentError::ReferenceOutOfRange { index: 1, len: 1 })
    );
    // Collection-wide failures are reported once, not per item.
    assert_eq!(
        trajectory.rmsd(0, &selection(&topology, &[])),
        Err(AlignmentError::EmptySelection)
    );
    assert_eq!(
        trajectory.rmsd(0, &AtomSelection::all(&chain(4))),
        Err(AlignmentError::SelectionTopologyMismatch)
    );
    let with_weights = |weights: &[f64]| {
        trajectory.rmsd_with_options(
            0,
            &all,
            AlignmentOptions {
                weighting: Weighting::Explicit(weights),
                ..AlignmentOptions::default()
            },
        )
    };
    assert_eq!(
        with_weights(&[1.0]),
        Err(AlignmentError::WeightCountMismatch {
            expected: 4,
            actual: 1,
        })
    );
    assert_eq!(
        with_weights(&[1.0, 1.0, f64::NAN, 1.0]),
        Err(AlignmentError::NonFiniteWeight { selection_index: 2 })
    );
    assert_eq!(
        with_weights(&[1.0, 1.0, 0.0, 1.0]),
        Err(AlignmentError::NonPositiveWeight { selection_index: 2 })
    );
    assert_eq!(
        trajectory.rmsd_with_options(
            0,
            &all,
            AlignmentOptions {
                periodic_policy: PeriodicPolicy::RejectPeriodic,
                ..AlignmentOptions::default()
            },
        ),
        Err(AlignmentError::Item {
            index: 0,
            source: Box::new(AlignmentError::PeriodicCoordinates {
                moving: true,
                reference: true,
            }),
        })
    );
    assert_eq!(trajectory.rmsd(0, &all).unwrap().value(), &[0.0]);
}

#[test]
fn aligned_rmsd_reports_the_failing_item_without_mutating_input() {
    let topology = chain(4);
    let trajectory = Trajectory::from_items(
        Arc::clone(&topology),
        [frame(&spread()), frame(&collinear())],
    )
    .unwrap();
    let before = payloads(&trajectory);
    let all = AtomSelection::all(&topology);
    assert_eq!(
        trajectory.aligned_rmsd(0, &all, &all),
        Err(AlignmentError::Item {
            index: 1,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Moving,
            }),
        })
    );
    assert_eq!(payloads(&trajectory), before);
}

#[test]
fn split_and_fused_workflows_keep_collection_and_item_properties() {
    let topology = chain(3);
    let reference = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let moving = reference.map(|point| Point3::new(point.x + 4.0, point.y - 2.0, point.z + 1.0));
    let key = PropertyKey::new("simulation").unwrap();
    let mut annotated = frame(&moving);
    annotated
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(key.clone(), PropertyValue::Int(7))
        .unwrap();
    let mut trajectory =
        Trajectory::from_items(Arc::clone(&topology), [frame(&reference), annotated]).unwrap();
    trajectory
        .properties_mut()
        .insert(key, PropertyValue::String("run_1".into()))
        .unwrap();
    let owner_properties = trajectory.properties().clone();
    let item_properties = trajectory.get(1).unwrap().properties().clone();
    let all = AtomSelection::all(&topology);
    let reject = AlignmentOptions {
        periodic_policy: PeriodicPolicy::RejectPeriodic,
        ..AlignmentOptions::default()
    };

    assert!(
        trajectory
            .rmsd_with_options(0, &all, reject)
            .unwrap()
            .value()[1]
            > 0.4
    );
    let fused = trajectory
        .aligned_rmsd_with_options(
            0,
            &all,
            &all,
            AlignedRmsdOptions {
                fit: reject,
                ..AlignedRmsdOptions::default()
            },
        )
        .unwrap();
    assert!(fused.value()[1] < 1.0e-12);

    trajectory.superpose(0, &all).unwrap();
    assert!(Arc::ptr_eq(&topology, &trajectory.shared_topology()));
    assert_eq!(trajectory.properties(), &owner_properties);
    assert_eq!(trajectory.get(1).unwrap().properties(), &item_properties);
    assert!(trajectory.rmsd(0, &all).unwrap().value()[1] < 1.0e-12);
}

#[test]
fn superposition_uses_stored_periodic_coordinates_and_explicit_weights() {
    let topology = chain(3);
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.1, 0.0, 0.0),
        Point3::new(0.0, 0.1, 0.0),
    ];
    let mut reference = frame(&points);
    let mut moving = frame(&points.map(|p| Point3::new(p.x + 0.3, p.y - 0.2, p.z + 0.1)));
    reference
        .conformation_mut()
        .set_cell(Some(periodic_cell(1.0)));
    moving.conformation_mut().set_cell(Some(periodic_cell(1.0)));
    let original = Trajectory::from_items(Arc::clone(&topology), [reference, moving]).unwrap();
    let all = AtomSelection::all(&topology);
    assert!(original.rmsd(0, &all).unwrap().value()[1] > 0.3);

    let mut weighted = original.clone();
    weighted
        .superpose_with_options(
            0,
            &all,
            AlignmentOptions {
                weighting: Weighting::Explicit(&[1.0, 2.0, 3.0]),
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
    assert!(weighted.rmsd(0, &all).unwrap().value()[1] < 1.0e-12);

    let mut rejected = original.clone();
    assert!(rejected
        .superpose_with_options(
            0,
            &all,
            AlignmentOptions {
                periodic_policy: PeriodicPolicy::RejectPeriodic,
                ..AlignmentOptions::default()
            },
        )
        .is_err());
    assert_eq!(payloads(&rejected), payloads(&original));
    rejected.superpose(0, &all).unwrap();
    for (actual, expected) in rejected.iter().zip(weighted.iter()) {
        for (a, b) in actual
            .positions()
            .values()
            .value()
            .iter()
            .zip(expected.positions().values().value().iter())
        {
            assert!((*a - *b).norm() < 1.0e-12);
        }
        // Uniform and weighted exact fits agree up to roundoff.
        let (actual, expected) = (actual.cell().unwrap(), expected.cell().unwrap());
        assert_eq!(actual.periodic_axes(), expected.periodic_axes());
        for (a, b) in actual
            .vectors()
            .into_value()
            .into_iter()
            .zip(expected.vectors().into_value())
        {
            assert_vector_close(a, b, 1.0e-12);
        }
    }
}

#[test]
fn ensembles_share_the_collection_contract_and_keep_weights() {
    let topology = chain(4);
    let reference = spread();
    let mut member = EnsembleMember::new(
        Positions::new(Quantity::new(
            transformed(&reference, quarter_turn()),
            NANOMETER,
        ))
        .unwrap(),
    );
    member.set_weight(Some(0.25)).unwrap();
    let mut ensemble = Ensemble::from_items(Arc::clone(&topology), [member]).unwrap();
    let model = Model::new(
        Arc::clone(&topology),
        Positions::new(Quantity::new(reference, NANOMETER)).unwrap(),
    )
    .unwrap();
    let all = AtomSelection::all(&topology);
    assert!(ensemble.rmsd(&model, &all).unwrap().value()[0] > 1.0);
    let report = ensemble.superpose(&model, &all).unwrap();
    assert_eq!(report.reference_index(), None);
    assert!(ensemble.rmsd(&model, &all).unwrap().value()[0] < 1.0e-12);
    assert_eq!(ensemble.get(0).unwrap().weight(), Some(0.25));

    // The free kernel and collection items share one view contract.
    let item = ensemble.get(0).unwrap();
    let fit = kabsch(item.as_model_view(), model.as_model_view(), &all).unwrap();
    assert!(fit.rmsd().into_value() < 1.0e-12);
    let reference: Reference<'_> = (&model).into();
    assert!(matches!(reference, Reference::View(_)));
    assert!(matches!(Reference::from(0), Reference::Index(0)));
}

fn moved(point: Point3) -> Point3 {
    Point3::new(-point.y + 2.0, point.x - 3.0, point.z + 1.0)
}

fn chain_model(points: &[Point3]) -> Model {
    Model::new(
        chain(points.len()),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

/// A moving trajectory, an independent reference whose atoms are permuted and
/// rigidly moved, a fit correspondence, and a measured correspondence.
fn correspondence_fixture() -> (Trajectory, Model, AtomCorrespondence, AtomCorrespondence) {
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 3.0, 0.0),
        Point3::new(0.0, 0.0, 4.0),
        Point3::new(2.0, 2.0, 2.0),
    ];
    let reference = chain_model(&[
        moved(points[2]),
        moved(points[0]),
        moved(points[3]),
        moved(points[1]),
        moved(points[4]),
        Point3::new(99.0, 99.0, 99.0),
    ]);
    let topology = chain(points.len());
    let reference_ids = reference.topology().atom_ids();
    let fit = AtomCorrespondence::from_pairs(
        &topology,
        &reference.shared_topology(),
        [(3, 2), (0, 1), (2, 0), (1, 3)].map(|(a, b)| (topology.atom_ids()[a], reference_ids[b])),
    )
    .unwrap();
    let measurement = AtomCorrespondence::from_pairs(
        &topology,
        &reference.shared_topology(),
        [(topology.atom_ids()[4], reference_ids[4])],
    )
    .unwrap();
    let frames = (0..2).map(|index| {
        let mut points = points;
        points[4].x += index as f64;
        let mut frame = frame(&points);
        frame
            .set_time(Some(Quantity::new(index as f64 + 2.0, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index as u64 + 20));
        frame
            .set_velocities(Some(
                Velocities::new(Quantity::new(
                    vec![Vector3::new(1.0, 0.0, 0.0); 5],
                    CANONICAL_VELOCITY_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        frame
            .set_forces(Some(
                Forces::new(Quantity::new(
                    vec![Vector3::new(0.0, 1.0, 0.0); 5],
                    CANONICAL_FORCE_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        let mut conformation = frame.conformation_mut();
        conformation.set_cell(Some(periodic_cell(20.0)));
        let mut properties = conformation.properties_mut();
        properties
            .owner_mut()
            .insert(
                PropertyKey::new("label").unwrap(),
                PropertyValue::String(format!("frame {index}")),
            )
            .unwrap();
        properties
            .atoms_mut()
            .set_value(
                PropertyKey::new("tag").unwrap(),
                TopologyAtomIndex::new(4),
                Some(PropertyValue::String("moving site".into())),
            )
            .unwrap();
        frame
    });
    let mut trajectory = Trajectory::from_items(topology, frames).unwrap();
    trajectory
        .properties_mut()
        .insert(
            PropertyKey::new("source").unwrap(),
            PropertyValue::String("moving system".into()),
        )
        .unwrap();
    (trajectory, reference, fit, measurement)
}

#[test]
fn correspondence_superposition_keeps_identity_and_transforms_every_vector_field() {
    let (trajectory, reference, fit, _) = correspondence_fixture();
    let reference_before = reference.clone();
    let weights = [8.0, 1.0, 4.0, 2.0];
    let mut aligned = trajectory.clone();
    let report = aligned
        .superpose_with_options(
            &reference,
            &fit,
            AlignmentOptions {
                weighting: Weighting::Explicit(&weights),
                ..AlignmentOptions::default()
            },
        )
        .unwrap();
    assert_eq!(report.reference_index(), None);
    assert!(Arc::ptr_eq(
        &trajectory.shared_topology(),
        &aligned.shared_topology()
    ));
    assert_eq!(aligned.properties(), trajectory.properties());
    for (actual, source) in aligned.iter().zip(trajectory.iter()) {
        assert_eq!(actual.properties(), source.properties());
        assert_eq!(actual.time(), source.time());
        assert_eq!(actual.step(), source.step());
        for (p, q) in actual
            .positions()
            .values()
            .value()
            .iter()
            .zip(source.positions().values().value().iter())
        {
            assert!((*p - moved(*q)).norm() < 1.0e-12);
        }
        for v in actual.velocities().unwrap().values().value().iter() {
            assert_vector_close(*v, Vector3::new(0.0, 1.0, 0.0), 1.0e-12);
        }
        for f in actual.forces().unwrap().values().value().iter() {
            assert_vector_close(*f, Vector3::new(-1.0, 0.0, 0.0), 1.0e-12);
        }
        let cell = actual.cell().unwrap();
        assert_vector_close(
            cell.vectors().value()[0],
            Vector3::new(0.0, 20.0, 0.0),
            1.0e-12,
        );
        assert_eq!(cell.periodic_axes(), [true; 3]);
    }
    assert_eq!(reference, reference_before);
}

#[test]
fn correspondence_fit_and_measurement_are_independent_and_weighted_in_pair_order() {
    let (trajectory, reference, fit, measurement) = correspondence_fixture();
    let fused = trajectory
        .aligned_rmsd(&reference, &fit, &measurement)
        .unwrap();
    assert!(fused.value()[0].abs() < 1.0e-12);
    assert!((fused.value()[1] - 1.0).abs() < 1.0e-12);
    let mut loaded = trajectory.clone();
    loaded.superpose(&reference, &fit).unwrap();
    let split = loaded.rmsd(&reference, &measurement).unwrap();
    for (a, b) in fused.value().iter().zip(split.value()) {
        assert!((a - b).abs() < 1.0e-12);
    }

    let top = trajectory.shared_topology();
    let reference_ids = reference.topology().atom_ids();
    let measured = AtomCorrespondence::from_pairs(
        &top,
        &reference.shared_topology(),
        [
            (top.atom_ids()[4], reference_ids[4]),
            (top.atom_ids()[0], reference_ids[1]),
        ],
    )
    .unwrap();
    let weighted = trajectory
        .aligned_rmsd_with_options(
            &reference,
            &fit,
            &measured,
            AlignedRmsdOptions {
                measurement_weighting: Weighting::Explicit(&[3.0, 1.0]),
                ..AlignedRmsdOptions::default()
            },
        )
        .unwrap();
    assert!((weighted.value()[1] - 0.75_f64.sqrt()).abs() < 1.0e-12);
    let direct = trajectory.rmsd(&reference, &measured).unwrap();
    for (frame, actual) in trajectory.iter().zip(direct.value()) {
        let expected = measured
            .index_pairs()
            .iter()
            .map(|&(a, b)| {
                (frame.positions().values().value()[a.index()]
                    - reference.positions().values().value()[b.index()])
                .norm_squared()
            })
            .sum::<f64>()
            / 2.0;
        assert!((actual - expected.sqrt()).abs() < 1.0e-12);
    }
    assert!(direct.value()[0] > 1.0);
}

#[test]
fn correspondence_failures_preserve_whole_collections() {
    let (mut trajectory, reference, fit, measurement) = correspondence_fixture();
    let atoms = trajectory.topology().atom_count();
    trajectory
        .push(TrajectoryFrame::new(Positions::zeros(atoms)))
        .unwrap();
    let before = payloads(&trajectory);
    assert!(matches!(
        trajectory.superpose(&reference, &fit),
        Err(AlignmentError::Item { index: 2, source })
            if matches!(*source, AlignmentError::DegenerateGeometry { .. })
    ));
    assert_eq!(payloads(&trajectory), before);
    assert!(matches!(
        trajectory.aligned_rmsd(&reference, &fit, &measurement),
        Err(AlignmentError::Item { index: 2, .. })
    ));
    let wrong = chain_model(reference.positions().values().value());
    assert_eq!(
        trajectory.rmsd(&wrong, &measurement),
        Err(AlignmentError::Correspondence(
            AtomCorrespondenceError::TopologyMismatch {
                side: CorrespondenceSide::Reference
            }
        ))
    );
    let mut empty = Trajectory::new(wrong.shared_topology());
    assert_eq!(
        empty.superpose(&reference, &fit),
        Err(AlignmentError::Correspondence(
            AtomCorrespondenceError::TopologyMismatch {
                side: CorrespondenceSide::Moving
            }
        ))
    );
}

#[test]
fn correspondence_measurement_validates_weights_empty_pairs_and_periodic_policy() {
    let (trajectory, reference, fit, measurement) = correspondence_fixture();
    let empty = AtomCorrespondence::from_pairs(
        &trajectory.shared_topology(),
        &reference.shared_topology(),
        [],
    )
    .unwrap();
    assert_eq!(
        trajectory.rmsd(&reference, &empty),
        Err(AlignmentError::EmptySelection)
    );
    assert!(matches!(
        trajectory.rmsd_with_options(
            &reference,
            &measurement,
            AlignmentOptions {
                weighting: Weighting::Explicit(&[]),
                ..AlignmentOptions::default()
            },
        ),
        Err(AlignmentError::WeightCountMismatch { .. })
    ));
    let reject = AlignmentOptions {
        periodic_policy: PeriodicPolicy::RejectPeriodic,
        ..AlignmentOptions::default()
    };
    assert!(matches!(
        trajectory.rmsd_with_options(&reference, &measurement, reject),
        Err(AlignmentError::Item { index: 0, source })
            if matches!(*source, AlignmentError::PeriodicCoordinates { moving: true, reference: false })
    ));
    assert!(matches!(
        trajectory.aligned_rmsd_with_options(
            &reference,
            &fit,
            &measurement,
            AlignedRmsdOptions {
                fit: reject,
                ..AlignedRmsdOptions::default()
            },
        ),
        Err(AlignmentError::Item { index: 0, source })
            if matches!(*source, AlignmentError::PeriodicCoordinates { .. })
    ));
}
