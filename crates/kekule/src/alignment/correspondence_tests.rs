use super::*;
use crate::core::{Atom, AtomId, BondOrder, Element, MoleculeEditor};
use crate::structure::{Model, Positions};
use crate::topology::{InstanceAtomId, Topology};
use crate::units::{ANGSTROM, NANOMETER};

fn model(points: &[Point3]) -> Model {
    let mut editor = MoleculeEditor::new();
    let mut previous = None;
    for _ in points {
        let atom = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(previous) = previous {
            editor.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    Model::new(
        Topology::from_molecule(&editor.finish().unwrap()).unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn points() -> [Point3; 4] {
    [
        Point3::new(0., 0., 0.),
        Point3::new(2., 0., 0.),
        Point3::new(0., 3., 0.),
        Point3::new(0., 0., 4.),
    ]
}

fn transformed(p: Point3) -> Point3 {
    Point3::new(-p.y + 2., p.x - 3., p.z + 1.)
}

fn pairs(moving: &Model, reference: &Model, order: &[(usize, usize)]) -> AtomCorrespondence {
    AtomCorrespondence::from_pairs(
        &moving.shared_topology(),
        &reference.shared_topology(),
        order
            .iter()
            .map(|&(a, b)| (moving.atom_ids()[a], reference.atom_ids()[b])),
    )
    .unwrap()
}

#[test]
fn independent_layouts_align_without_changing_identity_or_coordinates() {
    let moving = model(&points());
    let reference = model(&points().map(transformed));
    let before = (moving.clone(), reference.clone());
    let correspondence = AtomCorrespondence::from_same_layout(
        &moving.shared_topology(),
        &reference.shared_topology(),
    )
    .unwrap();
    assert_eq!(correspondence.len(), 4);
    assert!(!correspondence.is_empty());
    let selection = AtomSelection::all(&moving.shared_topology());
    assert_eq!(
        kabsch(moving.view(), reference.view(), &selection),
        Err(AlignmentError::TopologyMismatch)
    );
    let fit = kabsch_with_correspondence(moving.view(), reference.view(), &correspondence).unwrap();
    assert!(fit.rmsd().value_in(ANGSTROM).unwrap() < 1e-10);
    for p in points() {
        assert!((fit.transform().transform_point(p) - transformed(p)).norm() < 1e-12);
    }
    assert_eq!((moving, reference), before);
}

#[test]
fn explicit_pairs_preserve_order_and_fit_reordered_partial_models_in_mixed_units() {
    let source = points();
    let moving = model(&source);
    let target = [
        transformed(source[2]),
        transformed(source[0]),
        Point3::new(100., 100., 100.),
        transformed(source[3]),
        transformed(source[1]),
    ];
    let mut reference = model(&target);
    // Keep the same physical positions while exercising a different input unit.
    reference = Model::new(
        reference.shared_topology(),
        Positions::new(Quantity::new(
            target.map(|p| Point3::new(p.x * 10., p.y * 10., p.z * 10.)),
            ANGSTROM,
        ))
        .unwrap(),
    )
    .unwrap();
    assert!(matches!(
        AtomCorrespondence::from_same_layout(
            &moving.shared_topology(),
            &reference.shared_topology()
        ),
        Err(AtomCorrespondenceError::LayoutMismatch)
    ));
    let order = [(3, 3), (0, 1), (2, 0), (1, 4)];
    let correspondence = pairs(&moving, &reference, &order);
    assert_eq!(
        correspondence.atom_pairs().collect::<Vec<_>>(),
        order.map(|(a, b)| (moving.atom_ids()[a], reference.atom_ids()[b]))
    );
    let fit = kabsch_with_correspondence(moving.view(), reference.view(), &correspondence).unwrap();
    assert_eq!(fit.selected_atom_count(), 4);
    for p in source {
        assert!((fit.transform().transform_point(p) - transformed(p)).norm() < 1e-12);
    }
}

#[test]
fn correspondence_weights_follow_pairs_and_share_existing_numerical_kernel() {
    let moving = model(&points());
    let mut target = points().map(transformed);
    target[2].x += 0.4;
    let reference = model(&target);
    let order = [(3, 3), (0, 0), (2, 2), (1, 1)];
    let correspondence = pairs(&moving, &reference, &order);
    let weights = [8., 1., 4., 2.];
    let options = KabschOptions {
        weighting: AlignmentWeighting::Explicit(&weights),
        ..Default::default()
    };
    let fit = kabsch_with_correspondence_and_options(
        moving.view(),
        reference.view(),
        &correspondence,
        options,
    )
    .unwrap();
    let reference_on_shared_topology =
        Model::new(moving.shared_topology(), reference.positions().clone()).unwrap();
    let reordered_weights = [1., 2., 4., 8.];
    let expected = kabsch_with_options(
        moving.view(),
        reference_on_shared_topology.view(),
        &AtomSelection::all(&moving.shared_topology()),
        KabschOptions {
            weighting: AlignmentWeighting::Explicit(&reordered_weights),
            ..Default::default()
        },
    )
    .unwrap();
    assert!((fit.rmsd().into_value() - expected.rmsd().into_value()).abs() < 1e-12);
    for p in points() {
        assert!(
            (fit.transform().transform_point(p) - expected.transform().transform_point(p)).norm()
                < 1e-12
        );
    }
    assert!(fit.rmsd().into_value() > 0.01);
    for weights in [&[1., 2.][..], &[1., 0., 1., 1.], &[1., f64::NAN, 1., 1.]] {
        assert!(kabsch_with_correspondence_and_options(
            moving.view(),
            reference.view(),
            &correspondence,
            KabschOptions {
                weighting: AlignmentWeighting::Explicit(weights),
                ..Default::default()
            }
        )
        .is_err());
    }
}

#[test]
fn correspondence_rejects_invalid_and_duplicate_atoms_on_each_side() {
    let moving = model(&points());
    let reference = model(&points());
    let a = moving.atom_ids()[0];
    let b = reference.atom_ids()[0];
    let invalid = InstanceAtomId::new(a.molecule(), AtomId::new(99));
    for (side, input) in [
        (CorrespondenceSide::Moving, [(invalid, b)]),
        (CorrespondenceSide::Reference, [(a, invalid)]),
    ] {
        assert_eq!(
            AtomCorrespondence::from_pairs(
                &moving.shared_topology(),
                &reference.shared_topology(),
                input
            )
            .unwrap_err(),
            AtomCorrespondenceError::InvalidAtom {
                side,
                pair_index: 0,
                atom: invalid
            }
        );
    }
    for (side, repeated) in [
        (CorrespondenceSide::Moving, (a, reference.atom_ids()[1])),
        (CorrespondenceSide::Reference, (moving.atom_ids()[1], b)),
    ] {
        assert!(
            matches!(AtomCorrespondence::from_pairs(&moving.shared_topology(), &reference.shared_topology(), [(a,b), repeated]), Err(AtomCorrespondenceError::DuplicateAtom { side: actual, pair_index: 1, .. }) if actual == side)
        );
    }
}

#[test]
fn correspondence_rejects_new_snapshots_even_when_layouts_still_match() {
    let mut moving = model(&points());
    let mut reference = model(&points());
    let correspondence = AtomCorrespondence::from_same_layout(
        &moving.shared_topology(),
        &reference.shared_topology(),
    )
    .unwrap();
    let original_moving = moving.clone();
    moving.perceive().unwrap();
    assert!(moving
        .topology()
        .same_layout(correspondence.moving_topology()));
    assert!(matches!(
        kabsch_with_correspondence(moving.view(), reference.view(), &correspondence),
        Err(AlignmentError::Correspondence(
            AtomCorrespondenceError::TopologyMismatch {
                side: CorrespondenceSide::Moving
            }
        ))
    ));
    reference.perceive().unwrap();
    assert!(reference
        .topology()
        .same_layout(correspondence.reference_topology()));
    assert!(matches!(
        kabsch_with_correspondence(original_moving.view(), reference.view(), &correspondence),
        Err(AlignmentError::Correspondence(
            AtomCorrespondenceError::TopologyMismatch {
                side: CorrespondenceSide::Reference
            }
        ))
    ));
}

#[test]
fn correspondence_retains_alignment_geometry_and_periodic_rejections() {
    let mut moving = model(&points());
    let reference = model(&[Point3::new(0., 0., 0.); 4]);
    let correspondence = pairs(&moving, &reference, &[(0, 0), (1, 1), (2, 2), (3, 3)]);
    assert!(matches!(
        kabsch_with_correspondence(moving.view(), reference.view(), &correspondence),
        Err(AlignmentError::DegenerateGeometry { .. })
    ));
    let short = pairs(&moving, &reference, &[(0, 0), (1, 1)]);
    assert!(matches!(
        kabsch_with_correspondence(moving.view(), reference.view(), &short),
        Err(AlignmentError::InsufficientSelectedAtoms { selected: 2, .. })
    ));
    let empty = pairs(&moving, &reference, &[]);
    assert!(empty.is_empty());
    assert!(matches!(
        kabsch_with_correspondence(moving.view(), reference.view(), &empty),
        Err(AlignmentError::InsufficientSelectedAtoms { selected: 0, .. })
    ));
    moving.set_cell(Some(
        crate::geometry::PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(10., 0., 0.),
                    Vector3::new(0., 10., 0.),
                    Vector3::new(0., 0., 10.),
                ],
                NANOMETER,
            ),
            [true; 3],
        )
        .unwrap(),
    ));
    assert!(matches!(
        kabsch_with_correspondence_and_options(
            moving.view(),
            reference.view(),
            &correspondence,
            KabschOptions {
                periodic_policy: PeriodicAlignmentPolicy::RejectPeriodic,
                ..Default::default()
            }
        ),
        Err(AlignmentError::PeriodicCoordinates { moving: true, .. })
    ));
}
