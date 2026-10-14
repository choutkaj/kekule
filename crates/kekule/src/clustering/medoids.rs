//! Weighted k-medoids clustering by Partitioning Around Medoids.

use crate::structure::Ensemble;

use super::hierarchical::ObjectiveDistances;
use super::{ClusteringError, DistanceMatrix, EnsemblePartition};

impl Ensemble {
    /// Partitions the members into `clusters` clusters around medoid
    /// members, minimizing the weighted total distance
    /// `sum_i w_i d(i, m(i))` from each member to its nearest medoid.
    ///
    /// This is Partitioning Around Medoids (Kaufman and Rousseeuw): a greedy
    /// BUILD adds, one at a time, the member that most lowers the total
    /// distance, then SWAP repeatedly makes the single medoid exchange that
    /// lowers it most, until none does. Swaps are evaluated as in FastPAM1
    /// (Schubert and Rousseeuw 2019), which finds the same swaps as PAM in
    /// `O(n^2)` time per swap instead of `O(k n^2)`; BUILD costs `O(k n^2)`.
    /// The result is a local optimum, not necessarily the global one.
    ///
    /// Weights are relative, so heavier members pull medoids toward
    /// themselves. Each member joins its nearest medoid, a medoid always
    /// joins itself, and equal distances favor the medoid with the lowest
    /// member index. Every choice breaks ties by lowest member index, so the
    /// result is deterministic. Clusters are numbered in order of their
    /// lowest member, and each is represented by its medoid.
    ///
    /// `clusters` ranges from one to the member count (zero only for no
    /// members). Fails when the matrix covers a different member count.
    pub fn cluster_k_medoids(
        &self,
        distances: &DistanceMatrix,
        clusters: usize,
    ) -> Result<EnsemblePartition, ClusteringError> {
        if distances.len() != self.len() {
            return Err(ClusteringError::MatrixMismatch {
                expected: self.len(),
                actual: distances.len(),
            });
        }
        let members = self.len();
        let valid = if members == 0 {
            clusters == 0
        } else {
            (1..=members).contains(&clusters)
        };
        if !valid {
            return Err(ClusteringError::InvalidClusterCount {
                requested: clusters,
                members,
            });
        }
        let weights = self.scaled_weights();
        let objective = ObjectiveDistances::new(distances);
        let mut medoids = build(objective, &weights, clusters);
        swap(objective, &weights, &mut medoids);
        Ok(partition(distances, &medoids))
    }
}

/// Greedy initialization: each step adds the non-medoid that leaves the
/// smallest weighted total distance, the lowest index on ties.
fn build(distances: ObjectiveDistances<'_>, weights: &[f64], clusters: usize) -> Vec<usize> {
    let members = distances.len();
    let mut medoids = Vec::with_capacity(clusters);
    let mut is_medoid = vec![false; members];
    let mut nearest = vec![f64::INFINITY; members];
    for _ in 0..clusters {
        let mut best = (usize::MAX, f64::INFINITY);
        for candidate in (0..members).filter(|candidate| !is_medoid[*candidate]) {
            let total = (0..members)
                .map(|other| weights[other] * nearest[other].min(distances.at(candidate, other)))
                .sum::<f64>();
            if best.0 == usize::MAX || total < best.1 {
                best = (candidate, total);
            }
        }
        let medoid = best.0;
        medoids.push(medoid);
        is_medoid[medoid] = true;
        for (other, distance) in nearest.iter_mut().enumerate() {
            *distance = distance.min(distances.at(medoid, other));
        }
    }
    medoids
}

/// Each member's nearest medoid slot and distance, and its distance to the
/// second-nearest medoid (infinite with one medoid).
struct Nearest {
    slot: Vec<usize>,
    first: Vec<f64>,
    second: Vec<f64>,
}

impl Nearest {
    fn new(distances: ObjectiveDistances<'_>, medoids: &[usize]) -> Self {
        let members = distances.len();
        let mut nearest = Self {
            slot: vec![0; members],
            first: vec![f64::INFINITY; members],
            second: vec![f64::INFINITY; members],
        };
        for member in 0..members {
            for (slot, medoid) in medoids.iter().enumerate() {
                let distance = distances.at(*medoid, member);
                if distance < nearest.first[member] {
                    nearest.second[member] = nearest.first[member];
                    nearest.first[member] = distance;
                    nearest.slot[member] = slot;
                } else if distance < nearest.second[member] {
                    nearest.second[member] = distance;
                }
            }
        }
        nearest
    }

    fn total(&self, weights: &[f64]) -> f64 {
        self.first
            .iter()
            .zip(weights)
            .map(|(distance, weight)| weight * distance)
            .sum()
    }
}

/// Applies the best improving medoid exchange until none improves, stopping
/// as well if rounding makes a predicted improvement fail to lower the
/// recomputed total, so the loop always ends.
fn swap(distances: ObjectiveDistances<'_>, weights: &[f64], medoids: &mut [usize]) {
    let mut nearest = Nearest::new(distances, medoids);
    let mut total = nearest.total(weights);
    loop {
        let Some((slot, candidate)) = best_swap(distances, weights, medoids, &nearest, total)
        else {
            return;
        };
        let previous = medoids[slot];
        medoids[slot] = candidate;
        let next = Nearest::new(distances, medoids);
        let next_total = next.total(weights);
        if next_total >= total {
            medoids[slot] = previous;
            return;
        }
        nearest = next;
        total = next_total;
    }
}

/// The exchange of a medoid slot for a non-medoid that lowers the weighted
/// total most, if it lowers it by more than rounding; ties keep the lowest
/// candidate, then the lowest slot.
///
/// For each candidate `x`, removing the medoid in slot `m` and adding `x`
/// changes member `o`'s distance from `d_n` (nearest) to `min(d(x, o), d_s)`
/// (second-nearest) when `m` is its nearest medoid, and to
/// `min(d(x, o), d_n)` otherwise. The second change does not depend on `m`,
/// so one pass over the members prices every slot at once.
fn best_swap(
    distances: ObjectiveDistances<'_>,
    weights: &[f64],
    medoids: &[usize],
    nearest: &Nearest,
    total: f64,
) -> Option<(usize, usize)> {
    let members = distances.len();
    let mut is_medoid = vec![false; members];
    for medoid in medoids {
        is_medoid[*medoid] = true;
    }
    let threshold = -1.0e-12 * total;
    let mut best: Option<(usize, usize, f64)> = None;
    let mut changes = vec![0.0; medoids.len()];
    for candidate in (0..members).filter(|candidate| !is_medoid[*candidate]) {
        changes.fill(0.0);
        let mut shared = 0.0;
        for other in 0..members {
            let to_candidate = distances.at(candidate, other);
            let first = nearest.first[other];
            let removed = to_candidate.min(nearest.second[other]) - first;
            let kept = (to_candidate - first).min(0.0);
            changes[nearest.slot[other]] += weights[other] * (removed - kept);
            shared += weights[other] * kept;
        }
        for (slot, change) in changes.iter().enumerate() {
            let change = change + shared;
            if change < threshold && best.is_none_or(|(_, _, lowest)| change < lowest) {
                best = Some((slot, candidate, change));
            }
        }
    }
    best.map(|(slot, candidate, _)| (slot, candidate))
}

/// Assigns members to their nearest medoid and numbers clusters by lowest
/// member.
fn partition(distances: &DistanceMatrix, medoids: &[usize]) -> EnsemblePartition {
    let members = distances.len();
    let mut slot_of = vec![usize::MAX; members];
    for (slot, medoid) in medoids.iter().enumerate() {
        slot_of[*medoid] = slot;
    }
    // Equal distances favor the medoid with the lowest member index.
    let mut by_member = (0..medoids.len()).collect::<Vec<_>>();
    by_member.sort_by_key(|slot| medoids[*slot]);
    for (member, slot) in slot_of.iter_mut().enumerate() {
        if *slot != usize::MAX {
            continue;
        }
        let mut best = (by_member[0], f64::INFINITY);
        for candidate in &by_member {
            let distance = distances.at(medoids[*candidate], member);
            if distance < best.1 {
                best = (*candidate, distance);
            }
        }
        *slot = best.0;
    }
    let mut labels = vec![usize::MAX; medoids.len()];
    let mut assignments = Vec::with_capacity(members);
    let mut representatives = Vec::with_capacity(medoids.len());
    for slot in slot_of {
        if labels[slot] == usize::MAX {
            labels[slot] = representatives.len();
            representatives.push(medoids[slot]);
        }
        assignments.push(labels[slot]);
    }
    EnsemblePartition {
        assignments,
        representatives,
    }
}

#[cfg(test)]
mod tests {
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
}
