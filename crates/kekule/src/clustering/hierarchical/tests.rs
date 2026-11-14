use super::*;

/// Deterministic distances in `[0.5, 1.5)` from a linear congruential
/// generator, so no two are equal in practice.
fn random_matrix(members: usize, seed: u64) -> DistanceMatrix {
    let mut state = seed;
    let values = (0..condensed_len(members).unwrap())
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            0.5 + (state >> 11) as f64 / (1_u64 << 53) as f64
        })
        .collect();
    DistanceMatrix { members, values }
}

/// The textbook O(n^3) algorithm: recompute every cluster distance from
/// member distances and merge the closest pair.
fn naive_heights(matrix: &DistanceMatrix, linkage: Linkage) -> Vec<(f64, usize)> {
    let mut clusters = (0..matrix.member_count())
        .map(|member| vec![member])
        .collect::<Vec<_>>();
    let mut heights = Vec::new();
    while clusters.len() > 1 {
        let mut best = (0, 1, f64::INFINITY);
        for first in 0..clusters.len() {
            for second in first + 1..clusters.len() {
                let pairs = clusters[first]
                    .iter()
                    .flat_map(|a| clusters[second].iter().map(move |b| matrix.at(*a, *b)));
                let distance = match linkage {
                    Linkage::Single => pairs.fold(f64::INFINITY, f64::min),
                    Linkage::Complete => pairs.fold(0.0, f64::max),
                    Linkage::Average => {
                        let count = clusters[first].len() * clusters[second].len();
                        pairs.sum::<f64>() / count as f64
                    }
                };
                if distance < best.2 {
                    best = (first, second, distance);
                }
            }
        }
        let merged = clusters.remove(best.1);
        clusters[best.0].extend(merged);
        heights.push((best.2, clusters[best.0].len()));
    }
    heights
}

#[test]
fn condensed_pairs_invert_condensed_indices() {
    for members in 0..9 {
        let mut index = 0;
        for first in 0..members {
            for second in first + 1..members {
                assert_eq!(condensed_index(members, first, second), index);
                assert_eq!(condensed_pair(members, index), (first, second));
                index += 1;
            }
        }
        if members > 0 {
            assert_eq!(condensed_pair(members, index), (members - 1, members));
        }
    }
}

#[test]
fn objective_distances_rescale_only_matrices_whose_sums_could_overflow() {
    let ordinary = random_matrix(7, 5);
    let objective = ObjectiveDistances::new(&ordinary);
    assert_eq!(objective.scale, 1.0);
    assert_eq!(objective.at(2, 5), ordinary.at(2, 5));
    // 12 times 1.2e308 overflows until halved four times: 0.075e308
    // times 12 is 0.9e308.
    let huge = DistanceMatrix {
        members: 3,
        values: vec![1.2e308, 1.1e308, 0.9e308],
    };
    let objective = ObjectiveDistances::new(&huge);
    assert_eq!(objective.scale, 0.0625);
    assert_eq!(objective.at(1, 0), 1.2e308 * 0.0625);
}

#[test]
fn nearest_neighbor_chain_matches_the_naive_algorithm() {
    for linkage in [Linkage::Single, Linkage::Complete, Linkage::Average] {
        for (members, seed) in [(2, 1), (7, 2), (16, 3), (31, 4)] {
            let matrix = random_matrix(members, seed);
            let dendrogram = matrix.dendrogram(linkage);
            let expected = naive_heights(&matrix, linkage);
            assert_eq!(dendrogram.merges().len(), expected.len());
            for (merge, (height, size)) in dendrogram.merges().iter().zip(expected) {
                assert!(
                    (merge.height - height).abs() <= 1.0e-12,
                    "{linkage:?} {members}: {} != {height}",
                    merge.height
                );
                assert_eq!(merge.size, size, "{linkage:?} {members}");
            }
        }
    }
}
