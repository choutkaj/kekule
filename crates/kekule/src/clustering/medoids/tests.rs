use super::*;
use crate::units::{Quantity, NANOMETER};

/// Deterministic values in `[low, low + 1)` from a linear congruential
/// generator.
fn sequence(count: usize, seed: u64, low: f64) -> Vec<f64> {
    let mut state = seed;
    (0..count)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            low + (state >> 11) as f64 / (1_u64 << 53) as f64
        })
        .collect()
}

fn total(distances: ObjectiveDistances<'_>, weights: &[f64], medoids: &[usize]) -> f64 {
    Nearest::new(distances, medoids).total(weights)
}

/// Textbook PAM SWAP: recompute the full weighted total for every
/// exchange and apply the best one.
fn naive_swap(distances: ObjectiveDistances<'_>, weights: &[f64], medoids: &mut [usize]) {
    loop {
        let current = total(distances, weights, medoids);
        let mut best: Option<(usize, usize, f64)> = None;
        for candidate in 0..distances.len() {
            if medoids.contains(&candidate) {
                continue;
            }
            for slot in 0..medoids.len() {
                let mut trial = medoids.to_vec();
                trial[slot] = candidate;
                let change = total(distances, weights, &trial) - current;
                if change < -1.0e-12 * current && best.is_none_or(|(_, _, low)| change < low) {
                    best = Some((slot, candidate, change));
                }
            }
        }
        match best {
            Some((slot, candidate, _)) => medoids[slot] = candidate,
            None => return,
        }
    }
}

#[test]
fn fast_swaps_reach_the_same_medoids_as_naive_pam() {
    let mut swapped = false;
    for (members, clusters, seed) in [(6, 2, 1), (12, 3, 2), (25, 4, 3), (40, 7, 4)] {
        let pairs = members * (members - 1) / 2;
        let matrix = DistanceMatrix::from_condensed(
            members,
            Quantity::new(sequence(pairs, seed, 0.5), NANOMETER),
        )
        .unwrap();
        let distances = ObjectiveDistances::new(&matrix);
        let weights = sequence(members, seed + 100, 0.1);
        let start = build(distances, &weights, clusters);
        let mut fast = start.clone();
        swap(distances, &weights, &mut fast);
        swapped |= fast != start;
        let mut naive = start;
        naive_swap(distances, &weights, &mut naive);
        assert_eq!(fast, naive, "{members} members, {clusters} clusters");
    }
    // BUILD alone must not already be optimal everywhere.
    assert!(swapped);
}
