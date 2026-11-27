//! Clustering structure sets by radius, hierarchically, and around medoids,
//! weighted or not, and the populations, free energies, and collapse of the
//! clusterings.

use std::f64::consts::E;
use std::sync::Arc;

use kekule::alignment::{
    self, AlignmentError, AlignmentGeometry, AtomCorrespondence, FitAtoms, SymmetricSelection,
};
use kekule::clustering::{Clustering, ClusteringError, Cut, Dendrogram, DistanceMatrix, Linkage};
use kekule::diagnostics::{Countable, IndexOutOfRange};
use kekule::math::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::statistics::{Weights, WeightsError};
use kekule::structure::{Geometry, Positions, StructureSet};
use kekule::topology::{AtomSelection, Topology};
use kekule::units::{Quantity, UnitError, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

fn quarter_turn() -> RigidTransform {
    RigidTransform::new(
        Matrix3::from_columns(
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 1.0),
        ),
        Quantity::new(Vector3::new(4.0, -2.0, 1.0), kekule::units::NANOMETER),
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

fn geometry(points: Vec<Point3>) -> Geometry {
    Geometry::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
}

/// A set of geometries at `members` and weights of their linear weights.
fn conformers(
    topology: impl Into<Arc<Topology>>,
    members: impl IntoIterator<Item = (Vec<Point3>, f64)>,
) -> (StructureSet, Weights) {
    let (points, weights): (Vec<_>, Vec<_>) = members.into_iter().unzip();
    let set = StructureSet::from_geometries(topology, points.into_iter().map(geometry)).unwrap();
    let weights = given(&weights);
    (set, weights)
}

/// Given weights of these linear weights.
fn given(weights: &[f64]) -> Weights {
    Weights::given(
        weights.iter().map(|weight| Some(weight.ln())).collect(),
        None,
    )
    .unwrap()
}

/// Asserts that each value is within a few ulps of its expected value.
fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert!(
            (actual - expected).abs() <= 4.0 * f64::EPSILON * expected.abs().max(1.0),
            "{actual} != {expected}"
        );
    }
}

/// The log populations of `clustering` under `weights`, every cluster
/// populated.
fn log_populations(clustering: &Clustering, weights: &Weights) -> Vec<f64> {
    clustering
        .log_populations(weights)
        .unwrap()
        .into_iter()
        .map(Option::unwrap)
        .collect()
}

/// Shape A, A moved, its mirror B, and B moved, weighted 1, 3, 2, 1.
fn two_shapes() -> (Arc<Topology>, StructureSet, Weights) {
    let topology = kekule::smiles::to_topology("CCCC").unwrap();
    let members = [
        (shape(false, false), 1.0),
        (shape(false, true), 3.0),
        (shape(true, false), 2.0),
        (shape(true, true), 1.0),
    ];
    let (set, weights) = conformers(Arc::clone(&topology), members);
    (topology, set, weights)
}

/// The fitted RMSD between the two shapes.
fn shape_distance(set: &StructureSet, atoms: &AtomSelection) -> f64 {
    alignment::fit(
        set.get(0).unwrap().as_structure_view(),
        set.get(2).unwrap().as_structure_view(),
        atoms,
        Default::default(),
    )
    .unwrap()
    .rmsd()
    .into_value()
}

/// A clustering of every geometry.
fn clustering(assignments: &[usize], representatives: &[usize]) -> Clustering {
    Clustering::new(
        assignments.iter().copied().map(Some).collect(),
        representatives.to_vec(),
    )
    .unwrap()
}

#[test]
fn radius_clustering_groups_members_around_the_heaviest_leaders() {
    let (topology, set, weights) = two_shapes();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    let distance = shape_distance(&set, &atoms);
    let radius = |factor: f64| Quantity::new(distance * factor, NANOMETER);

    // By descending weight the moved A (3) leads, then B (2) founds its own
    // cluster; the remaining members join the cluster of their shape.
    let clusters = set
        .cluster_by_radius(radius(0.5), &atoms, Some(&weights), Default::default())
        .unwrap();
    assert_eq!(clusters, clustering(&[0, 0, 1, 1], &[1, 2]));
    assert_eq!(clusters.members(1).collect::<Vec<_>>(), [2, 3]);
    // The clusters' populations sum their samples' weights: 1 + 3 and 2 + 1
    // of 7.
    assert_close(
        &log_populations(&clusters, &weights),
        &[(4.0_f64 / 7.0).ln(), (3.0_f64 / 7.0).ln()],
    );
    // Unweighted, the set visits its geometries in geometry order: A leads,
    // then B founds the second cluster.
    assert_eq!(
        set.cluster_by_radius(radius(0.5), &atoms, None, Default::default())
            .unwrap(),
        clustering(&[0, 0, 1, 1], &[0, 2])
    );
    // A radius reaching the other shape leaves one cluster.
    assert_eq!(
        set.cluster_by_radius(radius(2.0), &atoms, Some(&weights), Default::default())
            .unwrap(),
        clustering(&[0, 0, 0, 0], &[1])
    );
    // Geometries that are not samples are in no cluster.
    let some = weights.restricted([0, 2, 3]).unwrap();
    let partial = set
        .cluster_by_radius(radius(0.5), &atoms, Some(&some), Default::default())
        .unwrap();
    assert_eq!(
        partial,
        Clustering::new(vec![Some(1), None, Some(0), Some(0)], vec![2, 0]).unwrap()
    );
    // Weights of another set are refused.
    assert_eq!(
        set.cluster_by_radius(
            radius(0.5),
            &atoms,
            Some(&Weights::equal(3)),
            Default::default()
        ),
        Err(ClusteringError::Weights(WeightsError::CoverageMismatch {
            covered: 3,
            geometries: 4,
        }))
    );

    let empty = StructureSet::new(Arc::clone(&topology));
    assert_eq!(
        empty
            .cluster_by_radius(radius(0.5), &atoms, None, Default::default())
            .unwrap(),
        clustering(&[], &[])
    );
}

#[test]
fn radius_clustering_rejects_invalid_radii_and_unfittable_members() {
    let (topology, set, weights) = two_shapes();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    let weights = Some(&weights);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            set.cluster_by_radius(
                Quantity::new(value, NANOMETER),
                &atoms,
                weights,
                Default::default()
            ),
            Err(ClusteringError::InvalidRadius)
        );
    }
    assert_eq!(
        set.cluster_by_radius(
            Quantity::new(1.0, KELVIN),
            &atoms,
            weights,
            Default::default()
        ),
        Err(ClusteringError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KELVIN,
                to: NANOMETER,
            }
        )))
    );
    let radius = Quantity::new(0.1, NANOMETER);
    let two_atoms = AtomSelection::from_atoms(
        Arc::clone(&topology),
        topology.atom_ids()[..2].iter().copied(),
    )
    .unwrap();
    assert_eq!(
        set.cluster_by_radius(radius, &two_atoms, weights, Default::default()),
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
    let (degenerate, weights) = conformers(
        Arc::clone(&topology),
        [(shape(false, false), 2.0), (collinear, 1.0)],
    );
    assert_eq!(
        degenerate.cluster_by_radius(radius, &atoms, Some(&weights), Default::default()),
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
    let (set, weights) = conformers(
        Arc::clone(&topology),
        [(acetate(0.11, -0.11), 1.0), (acetate(-0.11, 0.11), 1.0)],
    );
    let all = AtomSelection::all(Arc::clone(&topology));
    let radius = Quantity::new(0.01, NANOMETER);
    assert_eq!(
        set.cluster_by_radius(radius, &all, Some(&weights), Default::default())
            .unwrap(),
        clustering(&[0, 1], &[0, 1])
    );
    let symmetric = SymmetricSelection::new(&all).unwrap();
    let rotamers = set
        .cluster_by_radius(radius, &symmetric, Some(&weights), Default::default())
        .unwrap();
    assert_eq!(rotamers, clustering(&[0, 0], &[0]));
    let (representatives, collapsed) = rotamers.collapse(&weights).unwrap();
    assert_eq!(representatives, [0]);
    // The one cluster holds the whole population.
    assert_eq!(collapsed.log_weight(0), Some(0.0));
}

#[test]
fn a_heavier_sample_in_no_cluster_leaves_the_cluster_populations_exact() {
    // Geometry 2 outweighs the clustered ones by e^(10^16) and is in no
    // cluster, so the clusters weigh 1 and e among themselves.
    let clustering = Clustering::new(vec![Some(0), Some(1), None], vec![0, 1]).unwrap();
    let weights = Weights::given(vec![Some(0.0), Some(1.0), Some(1e16)], None).unwrap();
    let expected = [1.0 / (1.0 + E), E / (1.0 + E)];
    let probabilities = log_populations(&clustering, &weights)
        .into_iter()
        .map(f64::exp)
        .collect::<Vec<_>>();
    assert_close(&probabilities, &expected);
    let (representatives, collapsed) = clustering.collapse(&weights).unwrap();
    assert_eq!(representatives, [0, 1]);
    assert_close(&collapsed.probabilities().unwrap(), &expected);
}

#[test]
fn collapsing_clusters_keeps_representatives_with_their_populations() {
    let (_, mut set, weights) = two_shapes();
    let key = PropertyKey::new("source").unwrap();
    set.properties_mut()
        .insert(key.clone(), PropertyValue::Int(7))
        .unwrap();

    let clusters = clustering(&[0, 0, 1, 1], &[1, 2]);
    let (representatives, collapsed) = clusters.collapse(&weights).unwrap();
    assert_eq!(representatives, [1, 2]);
    // Each log weight is its cluster's log population: the weights 1 + 3 and
    // 2 + 1 of 7.
    assert_eq!(collapsed.len(), 2);
    assert_close(
        &[
            collapsed.log_weight(0).unwrap(),
            collapsed.log_weight(1).unwrap(),
        ],
        &[(4.0_f64 / 7.0).ln(), (3.0_f64 / 7.0).ln()],
    );
    // Selecting the representatives keeps their geometries and the set's
    // values.
    let kept = set.select(representatives.iter().copied()).unwrap();
    assert!(kept.topology().shares_layout(set.topology()));
    for (cluster, representative) in [(0, 1), (1, 2)] {
        assert_eq!(
            kept.get(cluster).unwrap().positions(),
            set.get(representative).unwrap().positions()
        );
    }
    assert_eq!(kept.properties().get(&key), Some(&PropertyValue::Int(7)));
    // Collapsed weights over the kept set give the same cluster shares.
    assert_close(&collapsed.probabilities().unwrap(), &[4.0 / 7.0, 3.0 / 7.0]);

    // Weights of another set are refused.
    assert_eq!(
        clusters.log_populations(&Weights::equal(3)),
        Err(ClusteringError::Weights(WeightsError::CoverageMismatch {
            covered: 3,
            geometries: 4,
        }))
    );
    // Weights of e^1000, beyond `f64::MAX`, sum without overflow: the first
    // cluster weighs 2 e^1000 against e^999.
    let huge = Weights::given(vec![Some(1000.0), Some(1000.0), Some(999.0)], None).unwrap();
    let (_, collapsed) = clustering(&[0, 0, 1], &[0, 2]).collapse(&huge).unwrap();
    let total = (2.0 + (-1.0_f64).exp()).ln();
    assert_close(
        &[
            collapsed.log_weight(0).unwrap(),
            collapsed.log_weight(1).unwrap(),
        ],
        &[2.0_f64.ln() - total, -1.0 - total],
    );
}

#[test]
fn cluster_free_energies_are_those_of_their_boltzmann_populations() {
    let (_, mut set, _) = two_shapes();
    let energy = PropertyKey::new("energy").unwrap();
    set.columns_mut()
        .insert(
            energy.clone(),
            PropertyColumn::Real {
                unit: KILOJOULE_PER_MOLE,
                values: [1.0, 3.0, 7.0, 9.0].map(Some).to_vec(),
            },
        )
        .unwrap();
    let weights = Weights::boltzmann(&set, &energy, Quantity::new(300.0, KELVIN)).unwrap();
    let clusters = clustering(&[0, 0, 1, 1], &[1, 2]);
    // Each cluster's free energy is -RT ln(e^(-a/RT) + e^(-b/RT)), relative
    // to the lowest.
    let rt: f64 = 8.314_462_618_153_24e-3 * 300.0;
    let free = |a: f64, b: f64| -rt * ((-a / rt).exp() + (-b / rt).exp()).ln();
    let energies = clusters.free_energies(&weights).unwrap();
    assert_eq!(energies.unit(), KILOJOULE_PER_MOLE);
    let energies = energies.into_value();
    assert_eq!(energies[0], Some(0.0));
    let expected = free(7.0, 9.0) - free(1.0, 3.0);
    assert!((energies[1].unwrap() - expected).abs() < 1e-12);
    // The collapsed weights keep the temperature, so the free energies of
    // the representatives are the clusters'.
    let (representatives, collapsed) = clusters.collapse(&weights).unwrap();
    let alone = Clustering::new(vec![Some(0), Some(1)], vec![0, 1]).unwrap();
    assert_eq!(collapsed.temperature(), weights.temperature());
    let again = alone.free_energies(&collapsed).unwrap().into_value();
    assert_eq!(representatives, [1, 2]);
    assert!((again[1].unwrap() - expected).abs() < 1e-12);
    // Weights without a temperature have no free energies.
    assert_eq!(
        clusters.free_energies(&Weights::equal(4)),
        Err(ClusteringError::Weights(WeightsError::MissingTemperature))
    );
    assert_eq!(
        ClusteringError::Weights(WeightsError::MissingTemperature).to_string(),
        "The weights have no temperature."
    );
}

#[test]
fn clusters_rank_cover_and_select_by_population() {
    // Clusters of weights 1, 3 + 1, and 2; geometry 4 in none.
    let clusters = Clustering::new(
        vec![Some(0), Some(1), Some(2), Some(1), None],
        vec![0, 1, 2],
    )
    .unwrap();
    let weights = given(&[1.0, 3.0, 2.0, 1.0, 100.0]);
    // Geometry 4 is in no cluster, so it is left out of the total.
    assert_close(
        &log_populations(&clusters, &weights),
        &[1.0 / 7.0, 4.0 / 7.0, 2.0 / 7.0].map(f64::ln),
    );
    assert_eq!(clusters.ranked_by_population(&weights).unwrap(), [1, 2, 0]);
    assert_eq!(clusters.covering(&weights, 0.5).unwrap(), [1]);
    assert_eq!(clusters.covering(&weights, 0.6).unwrap(), [1, 2]);
    assert_eq!(clusters.covering(&weights, 1.0).unwrap(), [1, 2, 0]);
    for fraction in [0.0, 1.5, f64::NAN] {
        assert_eq!(
            clusters.covering(&weights, fraction),
            Err(ClusteringError::Weights(WeightsError::InvalidFraction))
        );
    }
    assert_eq!(clusters.geometries_of(&[1, 0]).unwrap(), [1, 3, 0]);
    assert_eq!(
        clusters.geometries_of(&[3]),
        Err(ClusteringError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Cluster,
            3,
            3
        )))
    );
    // A cluster without samples has no population and ranks last.
    let some = weights.restricted([1, 2, 3]).unwrap();
    assert_eq!(
        clusters.log_populations(&some).unwrap()[0],
        None,
        "cluster 0 has no samples"
    );
    assert_eq!(clusters.ranked_by_population(&some).unwrap(), [1, 2, 0]);

    // Selecting renumbers the kept clusters in cluster order; a cluster
    // whose representative is not selected is represented by its first
    // selected geometry, and a duplicate is in its cluster twice.
    let selected = clusters.select([3, 4, 2, 3]).unwrap();
    assert_eq!(
        selected,
        Clustering::new(vec![Some(0), None, Some(1), Some(0)], vec![0, 2]).unwrap()
    );
    assert_eq!(
        clusters.select([1, 1]).unwrap(),
        Clustering::new(vec![Some(0), Some(0)], vec![0]).unwrap()
    );
    assert_eq!(
        clusters.select([5]),
        Err(ClusteringError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            5,
            5
        )))
    );
    // Selected weights and clustering agree with the selected set.
    let picked = [3, 2, 3];
    assert_close(
        &log_populations(
            &clusters.select(picked).unwrap(),
            &weights.select(picked).unwrap(),
        ),
        &[2.0 / 4.0, 2.0 / 4.0].map(f64::ln),
    );
}

#[test]
fn clusterings_contain_their_representatives() {
    let valid = Clustering::new(vec![Some(0), Some(1), Some(0)], vec![2, 1]).unwrap();
    assert_eq!(valid.len(), 3);
    assert_eq!(valid.cluster_count(), 2);
    assert_eq!(valid.cluster_of(2), Some(0));
    assert_eq!(valid.cluster_of(5), None);
    assert_eq!(valid.representative(0), Some(2));
    assert_eq!(valid.representatives(), [2, 1]);
    assert_eq!(valid.members(0).collect::<Vec<_>>(), [0, 2]);
    assert_eq!(valid.members(5).count(), 0);
    assert_eq!(valid.sizes(), [2, 1]);
    assert!(Clustering::new(Vec::new(), Vec::new()).unwrap().is_empty());

    for (cluster_of, representatives, error, message) in [
        (
            vec![Some(0), Some(2)],
            vec![0, 1],
            ClusteringError::ClusterOutOfRange {
                geometry: 1,
                cluster: 2,
                clusters: 2,
            },
            "Geometry 1 is in cluster 2, which has no representative.",
        ),
        (
            vec![Some(0), Some(1)],
            vec![0, 3],
            ClusteringError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 3, 2)),
            "Geometry 3 does not exist. There are 2 geometries.",
        ),
        // A cluster whose representative belongs elsewhere would be empty.
        (
            vec![Some(0), Some(0)],
            vec![0, 1],
            ClusteringError::RepresentativeNotMember {
                cluster: 1,
                geometry: 1,
            },
            "Geometry 1 represents cluster 1 but is not in it.",
        ),
        (
            vec![Some(0), None],
            vec![1],
            ClusteringError::RepresentativeNotMember {
                cluster: 0,
                geometry: 1,
            },
            "Geometry 1 represents cluster 0 but is not in it.",
        ),
    ] {
        let actual = Clustering::new(cluster_of, representatives).unwrap_err();
        assert_eq!(actual, error);
        assert_eq!(actual.to_string(), message);
    }
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
    assert_eq!(matrix.member_count(), 3);
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
fn pairwise_rmsd_fits_every_unordered_pair_of_samples_once() {
    let (topology, set, weights) = two_shapes();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, None, Default::default())
        .unwrap();
    assert_eq!(matrix.member_count(), 4);
    let fit = |first: usize, second: usize| {
        alignment::fit(
            set.get(first).unwrap().as_structure_view(),
            set.get(second).unwrap().as_structure_view(),
            &atoms,
            Default::default(),
        )
        .unwrap()
        .rmsd()
    };
    for first in 0..4 {
        for second in first + 1..4 {
            assert_eq!(matrix.get(first, second), Some(fit(first, second)));
            assert_eq!(matrix.get(second, first), Some(fit(first, second)));
        }
    }
    // With weights, the members are the samples in geometry order.
    let some = weights.restricted([1, 3]).unwrap();
    let samples = set
        .pairwise_fitted_rmsd(&atoms, Some(&some), Default::default())
        .unwrap();
    assert_eq!(samples.member_count(), 2);
    assert_eq!(samples.get(0, 1), Some(fit(1, 3)));

    let empty = StructureSet::new(Arc::clone(&topology))
        .pairwise_fitted_rmsd(&atoms, None, Default::default())
        .unwrap();
    assert!(empty.is_empty());

    let collinear = (0..4)
        .map(|index| Point3::new(f64::from(index), 0.0, 0.0))
        .collect();
    let (degenerate, _) = conformers(
        Arc::clone(&topology),
        [(shape(false, false), 1.0), (collinear, 1.0)],
    );
    assert_eq!(
        degenerate.pairwise_fitted_rmsd(&atoms, None, Default::default()),
        Err(ClusteringError::PairFit {
            first: 0,
            second: 1,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })
    );
}

/// Five members of one shape and weights of these linear weights.
fn weighted_members(weights: [f64; 5]) -> (StructureSet, Weights) {
    conformers(
        kekule::smiles::to_topology("CCCC").unwrap(),
        weights.map(|weight| (shape(false, false), weight)),
    )
}

#[test]
fn hierarchical_clusters_are_represented_by_weighted_medoids() {
    let distances = five_points();

    // Uniformly weighted, member 1 is closest to the rest of {0, 1, 2}
    // (summed distances 5, 4, 7), and the tie in {3, 4} keeps member 3.
    let (set, uniform) = weighted_members([1.0; 5]);
    let clusters = set
        .cluster_hierarchical(&distances, Some(&uniform), Linkage::Average, Cut::Count(2))
        .unwrap();
    assert_eq!(clusters, clustering(&[0, 0, 0, 1, 1], &[1, 3]));
    // Unweighted is the same as equal weights.
    assert_eq!(
        set.cluster_hierarchical(&distances, None, Linkage::Average, Cut::Count(2))
            .unwrap(),
        clusters
    );
    // Heavy members pull the medoid: weighted sums 4.1, 3.1, 0.7 and 0.3,
    // 0.15 after scaling by the largest weight.
    let (set, heavy) = weighted_members([1.0, 1.0, 10.0, 1.0, 2.0]);
    let weighted = set
        .cluster_hierarchical(&distances, Some(&heavy), Linkage::Average, Cut::Count(2))
        .unwrap();
    assert_eq!(weighted.representatives(), [2, 4]);
    // The dendrogram cut alone, without the distances, gives the same
    // clusters, each represented by its heaviest member, the first of equal
    // weights; here the heaviest members are the medoids.
    let dendrogram = distances.dendrogram(Linkage::Average);
    assert_eq!(
        dendrogram.clustering(Cut::Count(2), Some(&heavy)).unwrap(),
        weighted
    );
    assert_eq!(
        dendrogram.clustering(Cut::Count(2), None).unwrap(),
        clustering(&[0, 0, 0, 1, 1], &[0, 3])
    );

    assert_eq!(
        set.cluster_hierarchical(
            &line_distances(&[0.0, 1.0]),
            Some(&heavy),
            Linkage::Average,
            Cut::Count(1)
        ),
        Err(ClusteringError::MatrixMismatch {
            expected: 5,
            actual: 2,
        })
    );
    // A matrix over the samples clusters the samples: geometry 2 is none.
    let some = heavy.restricted([0, 1, 3, 4]).unwrap();
    let partial = set
        .cluster_hierarchical(
            &line_distances(&[0.0, 1.0, 10.0, 11.5]),
            Some(&some),
            Linkage::Average,
            Cut::Count(2),
        )
        .unwrap();
    assert_eq!(
        partial,
        Clustering::new(vec![Some(0), Some(0), None, Some(1), Some(1)], vec![0, 4]).unwrap()
    );

    // Summed distances beyond `f64::MAX` still compare: 2.3e308, 2.1e308,
    // and 2.0e308 nm make member 2 the medoid.
    let (three, _) = conformers(
        kekule::smiles::to_topology("CCCC").unwrap(),
        [1.0; 3].map(|weight| (shape(false, false), weight)),
    );
    assert_eq!(
        three
            .cluster_hierarchical(&huge_distances(), None, Linkage::Complete, Cut::Count(1))
            .unwrap(),
        clustering(&[0, 0, 0], &[2])
    );

    // Fitted RMSDs separate the two shapes under a complete-linkage cut.
    let (topology, shapes, weights) = two_shapes();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    let distance = shape_distance(&shapes, &atoms);
    let matrix = shapes
        .pairwise_fitted_rmsd(&atoms, Some(&weights), Default::default())
        .unwrap();
    let clusters = shapes
        .cluster_hierarchical(
            &matrix,
            Some(&weights),
            Linkage::Complete,
            Cut::Height(Quantity::new(distance / 2.0, NANOMETER)),
        )
        .unwrap();
    assert_eq!(clusters, clustering(&[0, 0, 1, 1], &[1, 2]));
}

#[test]
fn member_distances_do_not_depend_on_member_order() {
    let (topology, set, _) = two_shapes();
    let reversed = set.select([3, 2, 1, 0]).unwrap();
    let all = AtomSelection::all(Arc::clone(&topology));
    let symmetric = SymmetricSelection::new(&all).unwrap();
    for atoms in [FitAtoms::from(&all), FitAtoms::from(&symmetric)] {
        let forward = set
            .pairwise_fitted_rmsd(atoms, None, Default::default())
            .unwrap();
        let backward = reversed
            .pairwise_fitted_rmsd(atoms, None, Default::default())
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
        Arc::clone(&topology),
        Arc::clone(&topology),
        [
            (ids[0], ids[1]),
            (ids[1], ids[0]),
            (ids[2], ids[2]),
            (ids[3], ids[3]),
        ],
    )
    .unwrap();
    assert_eq!(
        set.pairwise_fitted_rmsd(&exchanged, None, Default::default()),
        Err(ClusteringError::DirectedCorrespondence)
    );
    assert_eq!(
        set.cluster_by_radius(
            Quantity::new(0.1, NANOMETER),
            &exchanged,
            None,
            Default::default()
        ),
        Err(ClusteringError::DirectedCorrespondence)
    );
    // The pairing is rejected before members are examined.
    assert_eq!(
        StructureSet::new(Arc::clone(&topology)).pairwise_fitted_rmsd(
            &exchanged,
            None,
            Default::default()
        ),
        Err(ClusteringError::DirectedCorrespondence)
    );
}

#[test]
fn k_medoids_minimizes_the_weighted_distance_to_medoids() {
    let distances = five_points();
    // Points 0, 1, 4, 10, 11.5 nm. BUILD picks 4 (summed distance 20.5),
    // then 10; swapping 4 for 1 lowers the total from 8.5 to the optimum
    // 5.5, which {1, 11.5} only ties.
    let (set, uniform) = weighted_members([1.0; 5]);
    let two = set
        .cluster_k_medoids(&distances, Some(&uniform), 2)
        .unwrap();
    assert_eq!(two, clustering(&[0, 0, 0, 1, 1], &[1, 3]));

    // Weighting 4 by 10 and 11.5 by 2 makes {4, 11.5} optimal at 8.5.
    let (set, heavy) = weighted_members([1.0, 1.0, 10.0, 1.0, 2.0]);
    let two = set.cluster_k_medoids(&distances, Some(&heavy), 2).unwrap();
    assert_eq!(two, clustering(&[0, 0, 0, 1, 1], &[2, 4]));
    assert_close(
        &log_populations(&two, &heavy),
        &[(12.0_f64 / 15.0).ln(), (3.0_f64 / 15.0).ln()],
    );
    // Unweighted, the set weighs its geometries equally, as the uniform
    // members.
    assert_eq!(
        set.cluster_k_medoids(&distances, None, 2).unwrap(),
        clustering(&[0, 0, 0, 1, 1], &[1, 3])
    );

    // One cluster gathers around the weighted medoid; one per member keeps
    // every member alone.
    assert_eq!(
        set.cluster_k_medoids(&distances, None, 1).unwrap(),
        clustering(&[0, 0, 0, 0, 0], &[2])
    );
    assert_eq!(
        set.cluster_k_medoids(&distances, None, 5).unwrap(),
        clustering(&[0, 1, 2, 3, 4], &[0, 1, 2, 3, 4])
    );

    // Summed distances beyond `f64::MAX` still compare: member 2 is the
    // medoid at 2.0e308 nm.
    let (three, _) = conformers(
        kekule::smiles::to_topology("CCCC").unwrap(),
        [1.0; 3].map(|weight| (shape(false, false), weight)),
    );
    assert_eq!(
        three.cluster_k_medoids(&huge_distances(), None, 1).unwrap(),
        clustering(&[0, 0, 0], &[2])
    );
}

#[test]
fn k_medoids_keeps_each_medoid_in_its_own_cluster() {
    // Identical members: every choice ties, so the lowest indices become
    // medoids, and the free member joins the lower one.
    let (identical, weights) = conformers(
        kekule::smiles::to_topology("CCCC").unwrap(),
        [1.0; 3].map(|weight| (shape(false, false), weight)),
    );
    let zeros = DistanceMatrix::from_condensed(3, Quantity::new(vec![0.0; 3], NANOMETER)).unwrap();
    assert_eq!(
        identical
            .cluster_k_medoids(&zeros, Some(&weights), 2)
            .unwrap(),
        clustering(&[0, 1, 0], &[0, 1])
    );
}

#[test]
fn k_medoids_rejects_invalid_cluster_counts_and_matrices() {
    let (set, weights) = weighted_members([1.0; 5]);
    let distances = five_points();
    for requested in [0, 6] {
        assert_eq!(
            set.cluster_k_medoids(&distances, Some(&weights), requested),
            Err(ClusteringError::InvalidClusterCount {
                requested,
                members: 5,
            })
        );
    }
    assert_eq!(
        set.cluster_k_medoids(&line_distances(&[0.0, 1.0]), Some(&weights), 1),
        Err(ClusteringError::MatrixMismatch {
            expected: 5,
            actual: 2,
        })
    );
    let empty = StructureSet::new(kekule::smiles::to_topology("CCCC").unwrap());
    let none = DistanceMatrix::from_condensed(0, Quantity::new(Vec::new(), NANOMETER)).unwrap();
    assert_eq!(
        empty.cluster_k_medoids(&none, None, 0).unwrap(),
        clustering(&[], &[])
    );

    // Fitted RMSDs separate the two shapes into two medoid clusters.
    let (topology, shapes, weights) = two_shapes();
    let matrix = shapes
        .pairwise_fitted_rmsd(
            &AtomSelection::all(Arc::clone(&topology)),
            Some(&weights),
            Default::default(),
        )
        .unwrap();
    let clusters = shapes
        .cluster_k_medoids(&matrix, Some(&weights), 2)
        .unwrap();
    let assignments = (0..4)
        .map(|geometry| clusters.cluster_of(geometry))
        .collect::<Vec<_>>();
    assert_eq!(assignments, [0, 0, 1, 1].map(Some));
}
