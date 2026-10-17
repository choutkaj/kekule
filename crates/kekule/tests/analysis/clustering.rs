//! Partitioning ensembles into clusters, hierarchical clustering on distance
//! matrices, and collapsing clusters to representatives.

use std::sync::Arc;

use kekule::alignment::{
    self, AlignmentError, AlignmentGeometry, AtomCorrespondence, FitAtoms, SymmetricSelection,
};
use kekule::clustering::{
    ClusterOrder, ClusteringError, ClusteringOptions, Cut, Dendrogram, DistanceMatrix,
    EnsemblePartition, Linkage, PartitionError,
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
    alignment::fit(
        ensemble.get(0).unwrap().as_model_view(),
        ensemble.get(2).unwrap().as_model_view(),
        atoms,
        Default::default(),
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
    let clusters = ensemble
        .cluster_by_radius(radius(0.5), &atoms, Default::default())
        .unwrap();
    assert_eq!(clusters, partition(&[0, 0, 1, 1], &[1, 2]));
    assert_eq!(clusters.members(1).collect::<Vec<_>>(), [2, 3]);
    // A radius reaching the other shape leaves one cluster.
    assert_eq!(
        ensemble
            .cluster_by_radius(radius(2.0), &atoms, Default::default())
            .unwrap(),
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
    let by_energy = ClusteringOptions::default().with_order(ClusterOrder::Energy);
    assert_eq!(
        ensemble
            .cluster_by_radius(radius(0.5), &atoms, by_energy)
            .unwrap(),
        partition(&[0, 0, 1, 1], &[0, 2])
    );
    ensemble.get_mut(3).unwrap().set_energy(None).unwrap();
    assert_eq!(
        ensemble.cluster_by_radius(radius(0.5), &atoms, by_energy),
        Err(ClusteringError::Ensemble(Box::new(
            EnsembleError::MissingEnergy { index: 3 }
        )))
    );

    let empty = Ensemble::new(Arc::clone(&topology));
    assert_eq!(
        empty
            .cluster_by_radius(radius(0.5), &atoms, Default::default())
            .unwrap(),
        partition(&[], &[])
    );
}

#[test]
fn radius_clustering_rejects_invalid_radii_and_unfittable_members() {
    let (topology, ensemble) = two_shapes();
    let atoms = AtomSelection::all(&topology);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            ensemble.cluster_by_radius(Quantity::new(value, NANOMETER), &atoms, Default::default()),
            Err(ClusteringError::InvalidRadius)
        );
    }
    assert_eq!(
        ensemble.cluster_by_radius(Quantity::new(1.0, KELVIN), &atoms, Default::default()),
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
        ensemble.cluster_by_radius(radius, &two_atoms, Default::default()),
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
        degenerate.cluster_by_radius(radius, &atoms, Default::default()),
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
        ensemble
            .cluster_by_radius(radius, &all, Default::default())
            .unwrap(),
        partition(&[0, 1], &[0, 1])
    );
    let symmetric = SymmetricSelection::new(&all).unwrap();
    let rotamers = ensemble
        .cluster_by_radius(radius, &symmetric, Default::default())
        .unwrap();
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

/// Condensed distances between points on a line, in nanometers.
fn line_distances(points: &[f64]) -> DistanceMatrix {
    let mut values = Vec::new();
    for (index, first) in points.iter().enumerate() {
        for second in &points[index + 1..] {
            values.push((second - first).abs());
        }
    }
    DistanceMatrix::from_condensed(points.len(), Quantity::new(values, NANOMETER)).unwrap()
}

/// Points at 0, 1, 4, 10, and 11.5 nm: two tight pairs and a point nearer
/// the first pair.
fn five_points() -> DistanceMatrix {
    line_distances(&[0.0, 1.0, 4.0, 10.0, 11.5])
}

/// Three members whose summed distances exceed `f64::MAX` nm: 2.3e308 from
/// member 0, 2.1e308 from member 1, and 2.0e308 from member 2.
fn huge_distances() -> DistanceMatrix {
    DistanceMatrix::from_condensed(3, Quantity::new(vec![1.2e308, 1.1e308, 0.9e308], NANOMETER))
        .unwrap()
}

/// Each merge as (clusters, height in nm, size).
fn merges(dendrogram: &Dendrogram) -> Vec<([usize; 2], f64, usize)> {
    dendrogram
        .merges()
        .iter()
        .map(|merge| (merge.clusters(), merge.height().into_value(), merge.size()))
        .collect()
}

#[test]
fn distance_matrices_store_each_unordered_pair_once() {
    let matrix =
        DistanceMatrix::from_condensed(3, Quantity::new(vec![1.0, 2.0, 3.0], NANOMETER)).unwrap();
    assert_eq!(matrix.len(), 3);
    assert_eq!(
        matrix.condensed(),
        Quantity::new(&[1.0, 2.0, 3.0][..], NANOMETER)
    );
    for (first, second, expected) in [
        (0, 1, 1.0),
        (1, 0, 1.0),
        (0, 2, 2.0),
        (2, 1, 3.0),
        (2, 2, 0.0),
    ] {
        assert_eq!(
            matrix.get(first, second),
            Some(Quantity::new(expected, NANOMETER))
        );
    }
    assert_eq!(matrix.get(3, 0), None);

    assert_eq!(
        DistanceMatrix::from_condensed(3, Quantity::new(vec![1.0, 2.0], NANOMETER)),
        Err(ClusteringError::CondensedLength {
            members: 3,
            expected: 3,
            actual: 2,
        })
    );
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            DistanceMatrix::from_condensed(3, Quantity::new(vec![1.0, invalid, 3.0], NANOMETER)),
            Err(ClusteringError::InvalidDistance { index: 1 })
        );
    }
    assert_eq!(
        DistanceMatrix::from_condensed(2, Quantity::new(vec![1.0], KELVIN)),
        Err(ClusteringError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KELVIN,
                to: NANOMETER,
            }
        )))
    );
    assert_eq!(
        DistanceMatrix::from_condensed(usize::MAX, Quantity::new(Vec::new(), NANOMETER)),
        Err(ClusteringError::MatrixTooLarge {
            members: usize::MAX
        })
    );
    for members in [0, 1] {
        let trivial =
            DistanceMatrix::from_condensed(members, Quantity::new(Vec::new(), NANOMETER)).unwrap();
        assert_eq!(trivial.dendrogram(Linkage::Average).merges().len(), 0);
    }
}

#[test]
fn dendrograms_merge_the_closest_clusters_under_each_linkage() {
    // Members are clusters 0-4; merges create 5-8. Heights derived by hand:
    // single takes the closest pair between clusters, complete the
    // farthest, and average the mean over all pairs.
    let matrix = five_points();
    assert_eq!(
        merges(&matrix.dendrogram(Linkage::Single)),
        [
            ([0, 1], 1.0, 2),
            ([3, 4], 1.5, 2),
            ([2, 5], 3.0, 3),
            ([6, 7], 6.0, 5)
        ]
    );
    assert_eq!(
        merges(&matrix.dendrogram(Linkage::Complete)),
        [
            ([0, 1], 1.0, 2),
            ([3, 4], 1.5, 2),
            ([2, 5], 4.0, 3),
            ([6, 7], 11.5, 5)
        ]
    );
    let average = matrix.dendrogram(Linkage::Average);
    assert_eq!(average.member_count(), 5);
    assert_eq!(
        merges(&average)[..3],
        [([0, 1], 1.0, 2), ([3, 4], 1.5, 2), ([2, 5], 3.5, 3)]
    );
    // The mean of the six distances between {0, 1, 4} and {10, 11.5} is
    // 54.5 / 6 = 109 / 12.
    let last = average.merges()[3];
    assert_eq!((last.clusters(), last.size()), ([6, 7], 5));
    assert!((last.height().into_value() - 109.0 / 12.0).abs() <= 1.0e-12);

    // Equal distances merge deterministically, keeping the earlier pair.
    let even = line_distances(&[0.0, 1.0, 2.0]).dendrogram(Linkage::Single);
    assert_eq!(merges(&even), [([0, 1], 1.0, 2), ([2, 3], 1.0, 3)]);
}

#[test]
fn average_linkage_heights_stay_finite_for_the_largest_finite_distances() {
    // The mean of equal distances is that distance, however close to
    // `f64::MAX`.
    let matrix =
        DistanceMatrix::from_condensed(3, Quantity::new(vec![1.0e308; 3], NANOMETER)).unwrap();
    let dendrogram = matrix.dendrogram(Linkage::Average);
    assert_eq!(
        merges(&dendrogram),
        [([0, 1], 1.0e308, 2), ([2, 3], 1.0e308, 3)]
    );
    assert_eq!(
        dendrogram
            .cut(Cut::Height(Quantity::new(1.1e308, NANOMETER)))
            .unwrap(),
        [0, 0, 0]
    );
}

#[test]
fn dendrogram_cuts_apply_merges_by_height_or_cluster_count() {
    let dendrogram = five_points().dendrogram(Linkage::Average);
    let count = |clusters| dendrogram.cut(Cut::Count(clusters));
    assert_eq!(count(1).unwrap(), [0, 0, 0, 0, 0]);
    assert_eq!(count(2).unwrap(), [0, 0, 0, 1, 1]);
    assert_eq!(count(3).unwrap(), [0, 0, 1, 2, 2]);
    assert_eq!(count(5).unwrap(), [0, 1, 2, 3, 4]);
    for requested in [0, 6] {
        assert_eq!(
            count(requested),
            Err(ClusteringError::InvalidClusterCount {
                requested,
                members: 5,
            })
        );
    }
    // Height cuts are inclusive.
    let height = |value| dendrogram.cut(Cut::Height(Quantity::new(value, NANOMETER)));
    assert_eq!(height(1.5).unwrap(), [0, 0, 1, 2, 2]);
    assert_eq!(height(0.5).unwrap(), [0, 1, 2, 3, 4]);
    assert_eq!(height(100.0).unwrap(), [0, 0, 0, 0, 0]);
    for invalid in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(height(invalid), Err(ClusteringError::InvalidHeight));
    }

    let empty = DistanceMatrix::from_condensed(0, Quantity::new(Vec::new(), NANOMETER))
        .unwrap()
        .dendrogram(Linkage::Single);
    assert_eq!(empty.cut(Cut::Count(0)).unwrap(), Vec::<usize>::new());
    assert_eq!(
        empty.cut(Cut::Count(1)),
        Err(ClusteringError::InvalidClusterCount {
            requested: 1,
            members: 0,
        })
    );
}

#[test]
fn pairwise_rmsd_fits_every_unordered_pair_once() {
    let (topology, ensemble) = two_shapes();
    let atoms = AtomSelection::all(&topology);
    let matrix = ensemble
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(matrix.len(), 4);
    for first in 0..4 {
        for second in first + 1..4 {
            let fit = alignment::fit(
                ensemble.get(first).unwrap().as_model_view(),
                ensemble.get(second).unwrap().as_model_view(),
                &atoms,
                Default::default(),
            )
            .unwrap();
            assert_eq!(matrix.get(first, second), Some(fit.rmsd()));
            assert_eq!(matrix.get(second, first), Some(fit.rmsd()));
        }
    }

    let empty = Ensemble::new(Arc::clone(&topology))
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert!(empty.is_empty());

    let collinear = (0..4)
        .map(|index| Point3::new(f64::from(index), 0.0, 0.0))
        .collect();
    let degenerate = Ensemble::from_items(
        Arc::clone(&topology),
        [member(shape(false, false), 1.0), member(collinear, 1.0)],
    )
    .unwrap();
    assert_eq!(
        degenerate.pairwise_fitted_rmsd(&atoms, Default::default()),
        Err(ClusteringError::PairFit {
            first: 0,
            second: 1,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })
    );
}

#[test]
fn hierarchical_clusters_are_represented_by_weighted_medoids() {
    let topology = kekule::smiles::to_topology("CCCC").unwrap();
    let ensemble = |weights: [f64; 5]| {
        Ensemble::from_items(
            Arc::clone(&topology),
            weights.map(|weight| member(shape(false, false), weight)),
        )
        .unwrap()
    };
    let distances = five_points();

    // Uniformly weighted, member 1 is closest to the rest of {0, 1, 2}
    // (summed distances 5, 4, 7), and the tie in {3, 4} keeps member 3.
    let uniform = ensemble([1.0; 5])
        .cluster_hierarchical(&distances, Linkage::Average, Cut::Count(2))
        .unwrap();
    assert_eq!(uniform.assignments(), [0, 0, 0, 1, 1]);
    assert_eq!(uniform.representatives(), [1, 3]);
    // Heavy members pull the medoid: weighted sums 4.1, 3.1, 0.7 and 0.3,
    // 0.15 after scaling by the largest weight.
    let weighted = ensemble([1.0, 1.0, 10.0, 1.0, 2.0])
        .cluster_hierarchical(&distances, Linkage::Average, Cut::Count(2))
        .unwrap();
    assert_eq!(weighted.representatives(), [2, 4]);

    assert_eq!(
        ensemble([1.0; 5]).cluster_hierarchical(
            &line_distances(&[0.0, 1.0]),
            Linkage::Average,
            Cut::Count(1)
        ),
        Err(ClusteringError::MatrixMismatch {
            expected: 5,
            actual: 2,
        })
    );

    // Summed distances beyond `f64::MAX` still compare: 2.3e308, 2.1e308,
    // and 2.0e308 nm make member 2 the medoid.
    let three = Ensemble::from_items(
        Arc::clone(&topology),
        [1.0; 3].map(|weight| member(shape(false, false), weight)),
    )
    .unwrap();
    assert_eq!(
        three
            .cluster_hierarchical(&huge_distances(), Linkage::Complete, Cut::Count(1))
            .unwrap(),
        partition(&[0, 0, 0], &[2])
    );

    // Fitted RMSDs separate the two shapes under a complete-linkage cut.
    let (topology, shapes) = two_shapes();
    let atoms = AtomSelection::all(&topology);
    let distance = shape_distance(&shapes, &atoms);
    let matrix = shapes
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    let clusters = shapes
        .cluster_hierarchical(
            &matrix,
            Linkage::Complete,
            Cut::Height(Quantity::new(distance / 2.0, NANOMETER)),
        )
        .unwrap();
    assert_eq!(clusters.assignments(), [0, 0, 1, 1]);
}

#[test]
fn member_distances_do_not_depend_on_member_order() {
    let (topology, ensemble) = two_shapes();
    let reversed = ensemble.select([3, 2, 1, 0]).unwrap();
    let all = AtomSelection::all(&topology);
    let symmetric = SymmetricSelection::new(&all).unwrap();
    for atoms in [FitAtoms::from(&all), FitAtoms::from(&symmetric)] {
        let forward = ensemble
            .pairwise_fitted_rmsd(atoms, Default::default())
            .unwrap();
        let backward = reversed
            .pairwise_fitted_rmsd(atoms, Default::default())
            .unwrap();
        for first in 0..4 {
            for second in 0..4 {
                let a = forward.get(first, second).unwrap().into_value();
                let b = backward.get(3 - first, 3 - second).unwrap().into_value();
                // Fitting either member onto the other gives the same RMSD up
                // to rounding.
                assert!((a - b).abs() <= 1.0e-12, "{first} {second}: {a} != {b}");
            }
        }
    }

    // Exchanging two atoms makes each fit depend on which member moves.
    let ids = topology.atom_ids();
    let exchanged = AtomCorrespondence::from_pairs(
        &topology,
        &topology,
        [
            (ids[0], ids[1]),
            (ids[1], ids[0]),
            (ids[2], ids[2]),
            (ids[3], ids[3]),
        ],
    )
    .unwrap();
    assert_eq!(
        ensemble.pairwise_fitted_rmsd(&exchanged, Default::default()),
        Err(ClusteringError::DirectedCorrespondence)
    );
    assert_eq!(
        ensemble.cluster_by_radius(
            Quantity::new(0.1, NANOMETER),
            &exchanged,
            Default::default()
        ),
        Err(ClusteringError::DirectedCorrespondence)
    );
    // The pairing is rejected before members are examined.
    assert_eq!(
        Ensemble::new(Arc::clone(&topology)).pairwise_fitted_rmsd(&exchanged, Default::default()),
        Err(ClusteringError::DirectedCorrespondence)
    );
}

/// Five members of one shape with the given weights.
fn weighted_members(weights: [f64; 5]) -> Ensemble {
    Ensemble::from_items(
        kekule::smiles::to_topology("CCCC").unwrap(),
        weights.map(|weight| member(shape(false, false), weight)),
    )
    .unwrap()
}

#[test]
fn k_medoids_minimizes_the_weighted_distance_to_medoids() {
    let distances = five_points();
    // Points 0, 1, 4, 10, 11.5 nm. BUILD picks 4 (summed distance 20.5),
    // then 10; swapping 4 for 1 lowers the total from 8.5 to the optimum
    // 5.5, which {1, 11.5} only ties.
    let uniform = weighted_members([1.0; 5]);
    let two = uniform.cluster_k_medoids(&distances, 2).unwrap();
    assert_eq!(two.assignments(), [0, 0, 0, 1, 1]);
    assert_eq!(two.representatives(), [1, 3]);

    // Weighting 4 by 10 and 11.5 by 2 makes {4, 11.5} optimal at 8.5.
    let weighted = weighted_members([1.0, 1.0, 10.0, 1.0, 2.0]);
    let two = weighted.cluster_k_medoids(&distances, 2).unwrap();
    assert_eq!(two.assignments(), [0, 0, 0, 1, 1]);
    assert_eq!(two.representatives(), [2, 4]);

    // One cluster gathers around the weighted medoid; one per member keeps
    // every member alone.
    assert_eq!(
        uniform.cluster_k_medoids(&distances, 1).unwrap(),
        partition(&[0, 0, 0, 0, 0], &[2])
    );
    assert_eq!(
        uniform.cluster_k_medoids(&distances, 5).unwrap(),
        partition(&[0, 1, 2, 3, 4], &[0, 1, 2, 3, 4])
    );

    // Summed distances beyond `f64::MAX` still compare: member 2 is the
    // medoid at 2.0e308 nm.
    let three = Ensemble::from_items(
        kekule::smiles::to_topology("CCCC").unwrap(),
        [1.0; 3].map(|weight| member(shape(false, false), weight)),
    )
    .unwrap();
    assert_eq!(
        three.cluster_k_medoids(&huge_distances(), 1).unwrap(),
        partition(&[0, 0, 0], &[2])
    );
}

#[test]
fn k_medoids_keeps_each_medoid_in_its_own_cluster() {
    // Identical members: every choice ties, so the lowest indices become
    // medoids, and the free member joins the lower one.
    let topology = kekule::smiles::to_topology("CCCC").unwrap();
    let identical = Ensemble::from_items(
        topology,
        [1.0; 3].map(|weight| member(shape(false, false), weight)),
    )
    .unwrap();
    let zeros = DistanceMatrix::from_condensed(3, Quantity::new(vec![0.0; 3], NANOMETER)).unwrap();
    assert_eq!(
        identical.cluster_k_medoids(&zeros, 2).unwrap(),
        partition(&[0, 1, 0], &[0, 1])
    );
}

#[test]
fn k_medoids_rejects_invalid_cluster_counts_and_matrices() {
    let ensemble = weighted_members([1.0; 5]);
    let distances = five_points();
    for requested in [0, 6] {
        assert_eq!(
            ensemble.cluster_k_medoids(&distances, requested),
            Err(ClusteringError::InvalidClusterCount {
                requested,
                members: 5,
            })
        );
    }
    assert_eq!(
        ensemble.cluster_k_medoids(&line_distances(&[0.0, 1.0]), 1),
        Err(ClusteringError::MatrixMismatch {
            expected: 5,
            actual: 2,
        })
    );
    let topology = kekule::smiles::to_topology("CCCC").unwrap();
    let empty = Ensemble::new(topology);
    let none = DistanceMatrix::from_condensed(0, Quantity::new(Vec::new(), NANOMETER)).unwrap();
    assert_eq!(
        empty.cluster_k_medoids(&none, 0).unwrap(),
        partition(&[], &[])
    );

    // Fitted RMSDs separate the two shapes into two medoid clusters.
    let (topology, shapes) = two_shapes();
    let matrix = shapes
        .pairwise_fitted_rmsd(&AtomSelection::all(&topology), Default::default())
        .unwrap();
    assert_eq!(
        shapes.cluster_k_medoids(&matrix, 2).unwrap().assignments(),
        [0, 0, 1, 1]
    );
}
