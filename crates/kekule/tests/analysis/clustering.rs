//! Partitioning ensembles into clusters and collapsing them.

use std::sync::Arc;

use kekule::alignment::{kabsch, AlignmentError, AlignmentGeometry, SymmetricSelection};
use kekule::clustering::{
    ClusterOrder, ClusteringError, ClusteringOptions, EnsemblePartition, PartitionError,
};
use kekule::geometry::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Ensemble, EnsembleError, EnsembleMember, Positions, StructureError};
use kekule::topology::{AtomSelection, Topology};
use kekule::units::{Quantity, UnitError, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

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

/// A chiral four-point shape, or its mirror image, optionally moved rigidly.
fn shape(mirror: bool, moved: bool) -> Vec<Point3> {
    let z = if mirror { -1.2 } else { 1.2 };
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(0.0, 1.5, 0.0),
        Point3::new(0.5, 0.4, z),
    ];
    points
        .into_iter()
        .map(|point| {
            if moved {
                quarter_turn().transform_point(point)
            } else {
                point
            }
        })
        .collect()
}

fn member(points: Vec<Point3>, weight: f64) -> EnsembleMember {
    EnsembleMember::new(
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
        weight,
    )
    .unwrap()
}

/// Shape A, A moved, its mirror B, and B moved, weighted 1, 3, 2, 1.
fn two_shapes() -> (Arc<Topology>, Ensemble) {
    let topology = kekule::smiles::to_topology("CCCC").unwrap();
    let members = [
        member(shape(false, false), 1.0),
        member(shape(false, true), 3.0),
        member(shape(true, false), 2.0),
        member(shape(true, true), 1.0),
    ];
    let ensemble = Ensemble::from_items(Arc::clone(&topology), members).unwrap();
    (topology, ensemble)
}

/// The fitted RMSD between the two shapes.
fn shape_distance(ensemble: &Ensemble, atoms: &AtomSelection) -> f64 {
    kabsch(
        ensemble.get(0).unwrap().as_model_view(),
        ensemble.get(2).unwrap().as_model_view(),
        atoms,
    )
    .unwrap()
    .rmsd()
    .into_value()
}

fn partition(assignments: &[usize], representatives: &[usize]) -> EnsemblePartition {
    EnsemblePartition::new(assignments.to_vec(), representatives.to_vec()).unwrap()
}

#[test]
fn radius_clustering_groups_members_around_the_first_visited_leaders() {
    let (topology, mut ensemble) = two_shapes();
    let atoms = AtomSelection::all(&topology);
    let distance = shape_distance(&ensemble, &atoms);
    let radius = |factor: f64| Quantity::new(distance * factor, NANOMETER);

    // By descending weight the moved A (3) leads, then B (2) founds its own
    // cluster; the remaining members join the cluster of their shape.
    let clusters = ensemble.cluster_by_radius(radius(0.5), &atoms).unwrap();
    assert_eq!(clusters, partition(&[0, 0, 1, 1], &[1, 2]));
    assert_eq!(clusters.members(1).collect::<Vec<_>>(), [2, 3]);
    // A radius reaching the other shape leaves one cluster.
    assert_eq!(
        ensemble.cluster_by_radius(radius(2.0), &atoms).unwrap(),
        partition(&[0, 0, 0, 0], &[1])
    );

    // By ascending energy the lowest member of each shape leads.
    for (index, energy) in [0.0, 5.0, 1.0, 3.0].into_iter().enumerate() {
        ensemble
            .get_mut(index)
            .unwrap()
            .set_energy(Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
            .unwrap();
    }
    let by_energy = ClusteringOptions {
        order: ClusterOrder::Energy,
        ..ClusteringOptions::default()
    };
    assert_eq!(
        ensemble
            .cluster_by_radius_with_options(radius(0.5), &atoms, by_energy)
            .unwrap(),
        partition(&[0, 0, 1, 1], &[0, 2])
    );
    ensemble.get_mut(3).unwrap().set_energy(None).unwrap();
    assert_eq!(
        ensemble.cluster_by_radius_with_options(radius(0.5), &atoms, by_energy),
        Err(ClusteringError::Ensemble(Box::new(
            EnsembleError::MissingEnergy { index: 3 }
        )))
    );

    let empty = Ensemble::new(Arc::clone(&topology));
    assert_eq!(
        empty.cluster_by_radius(radius(0.5), &atoms).unwrap(),
        partition(&[], &[])
    );
}

#[test]
fn radius_clustering_rejects_invalid_radii_and_unfittable_members() {
    let (topology, ensemble) = two_shapes();
    let atoms = AtomSelection::all(&topology);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            ensemble.cluster_by_radius(Quantity::new(value, NANOMETER), &atoms),
            Err(ClusteringError::InvalidRadius)
        );
    }
    assert_eq!(
        ensemble.cluster_by_radius(Quantity::new(1.0, KELVIN), &atoms),
        Err(ClusteringError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KELVIN,
                to: NANOMETER,
            }
        )))
    );
    let radius = Quantity::new(0.1, NANOMETER);
    let two_atoms =
        AtomSelection::from_atoms(&topology, topology.atom_ids()[..2].iter().copied()).unwrap();
    assert_eq!(
        ensemble.cluster_by_radius(radius, &two_atoms),
        Err(ClusteringError::Alignment(Box::new(
            AlignmentError::InsufficientSelectedAtoms {
                selected: 2,
                minimum: 3,
            }
        )))
    );

    // The collinear member is placed second and cannot be fitted.
    let collinear = (0..4)
        .map(|index| Point3::new(f64::from(index), 0.0, 0.0))
        .collect();
    let degenerate = Ensemble::from_items(
        Arc::clone(&topology),
        [member(shape(false, false), 2.0), member(collinear, 1.0)],
    )
    .unwrap();
    assert_eq!(
        degenerate.cluster_by_radius(radius, &atoms),
        Err(ClusteringError::Alignment(Box::new(AlignmentError::Item {
            index: 1,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Moving,
            }),
        })))
    );
}

#[test]
fn symmetric_radius_clustering_merges_rotamers() {
    // Acetate and a copy with its oxygens exchanged, lifted out of the plane
    // so no rotation exchanges them.
    let topology = kekule::smiles::to_topology("CC(=O)[O-]").unwrap();
    let acetate = |o2: f64, o3: f64| {
        vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.15, 0.0, 0.0),
            Point3::new(0.22, o2, 0.05),
            Point3::new(0.22, o3, 0.05),
        ]
    };
    let ensemble = Ensemble::from_items(
        Arc::clone(&topology),
        [
            member(acetate(0.11, -0.11), 1.0),
            member(acetate(-0.11, 0.11), 1.0),
        ],
    )
    .unwrap();
    let all = AtomSelection::all(&topology);
    let radius = Quantity::new(0.01, NANOMETER);
    assert_eq!(
        ensemble.cluster_by_radius(radius, &all).unwrap(),
        partition(&[0, 1], &[0, 1])
    );
    let symmetric = SymmetricSelection::new(&all).unwrap();
    let rotamers = ensemble.cluster_by_radius(radius, &symmetric).unwrap();
    assert_eq!(rotamers, partition(&[0, 0], &[0]));
    let collapsed = ensemble.collapse(&rotamers).unwrap();
    assert_eq!(collapsed.len(), 1);
    assert_eq!(collapsed.get(0).unwrap().weight(), 2.0);
}

#[test]
fn collapse_keeps_representatives_with_summed_cluster_weights() {
    let (_, mut ensemble) = two_shapes();
    let key = PropertyKey::new("source").unwrap();
    ensemble
        .properties_mut()
        .insert(key.clone(), PropertyValue::Int(7))
        .unwrap();
    ensemble
        .get_mut(2)
        .unwrap()
        .set_energy(Some(Quantity::new(7.0, KILOJOULE_PER_MOLE)))
        .unwrap();

    let collapsed = ensemble
        .collapse(&partition(&[0, 0, 1, 1], &[1, 2]))
        .unwrap();
    assert_eq!(collapsed.len(), 2);
    assert!(collapsed.topology().shares_layout(ensemble.topology()));
    assert_eq!(
        collapsed
            .iter()
            .map(|member| member.weight())
            .collect::<Vec<_>>(),
        [4.0, 3.0]
    );
    for (cluster, representative) in [(0, 1), (1, 2)] {
        assert_eq!(
            collapsed.get(cluster).unwrap().positions(),
            ensemble.get(representative).unwrap().positions()
        );
        assert_eq!(
            collapsed.get(cluster).unwrap().energy(),
            ensemble.get(representative).unwrap().energy()
        );
    }
    assert_eq!(
        collapsed.properties().get(&key),
        Some(&PropertyValue::Int(7))
    );

    assert_eq!(
        ensemble.collapse(&partition(&[0, 0, 0], &[0])).unwrap_err(),
        ClusteringError::PartitionMismatch {
            expected: 4,
            actual: 3,
        }
    );
    let topology = ensemble.shared_topology();
    let huge = Ensemble::from_items(
        topology,
        [
            member(shape(false, false), f64::MAX),
            member(shape(false, true), f64::MAX),
        ],
    )
    .unwrap();
    assert_eq!(
        huge.collapse(&partition(&[0, 0], &[0])).unwrap_err(),
        ClusteringError::Ensemble(Box::new(EnsembleError::Structure(
            StructureError::InvalidWeight
        )))
    );
}

#[test]
fn partitions_assign_every_member_to_a_cluster_containing_its_representative() {
    let valid = partition(&[0, 1, 0], &[2, 1]);
    assert_eq!(valid.member_count(), 3);
    assert_eq!(valid.cluster_count(), 2);
    assert_eq!(valid.assignments(), [0, 1, 0]);
    assert_eq!(valid.representatives(), [2, 1]);
    assert_eq!(valid.members(0).collect::<Vec<_>>(), [0, 2]);
    assert_eq!(valid.members(5).count(), 0);

    assert_eq!(
        EnsemblePartition::new(vec![0, 2], vec![0, 1]),
        Err(PartitionError::ClusterOutOfRange {
            member: 1,
            cluster: 2,
            clusters: 2,
        })
    );
    assert_eq!(
        EnsemblePartition::new(vec![0, 1], vec![0, 3]),
        Err(PartitionError::RepresentativeOutOfRange {
            cluster: 1,
            member: 3,
            members: 2,
        })
    );
    // A cluster whose representative belongs elsewhere would be empty.
    assert_eq!(
        EnsemblePartition::new(vec![0, 0], vec![0, 1]),
        Err(PartitionError::RepresentativeNotInCluster {
            cluster: 1,
            member: 1,
        })
    );
}
