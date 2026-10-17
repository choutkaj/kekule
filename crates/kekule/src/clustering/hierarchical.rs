//! Pairwise distance matrices and agglomerative hierarchical clustering.

use crate::alignment::{
    best_kabsch, check_fit_count, AlignmentOptions, FitAtoms, NormalizedWeights,
};
use crate::structure::{Ensemble, ModelView};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};

use super::{undirected, ClusteringError, EnsemblePartition};

/// Symmetric pairwise distances between the members of one ensemble, stored
/// once per unordered pair in [`CANONICAL_LENGTH_UNIT`].
///
/// The condensed order lists the pairs `(i, j)` with `i < j` row by row:
/// `(0, 1), (0, 2), ..., (0, n - 1), (1, 2), ...`, the order SciPy's
/// `squareform` uses. Build one from fitted RMSDs with
/// [`Ensemble::pairwise_fitted_rmsd`] or from any other distances with
/// [`Self::from_condensed`]. It stores `n (n - 1) / 2` values, so memory
/// grows quadratically with the member count.
#[derive(Debug, Clone, PartialEq)]
pub struct DistanceMatrix {
    members: usize,
    values: Vec<f64>,
}

impl DistanceMatrix {
    /// Validates condensed distances between `members` members: exactly
    /// `members (members - 1) / 2` finite, nonnegative lengths.
    pub fn from_condensed(
        members: usize,
        distances: Quantity<Vec<f64>>,
    ) -> Result<Self, ClusteringError> {
        let expected = condensed_len(members)?;
        let values = distances.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
        if values.len() != expected {
            return Err(ClusteringError::CondensedLength {
                members,
                expected,
                actual: values.len(),
            });
        }
        if let Some(index) = values
            .iter()
            .position(|value| !(value.is_finite() && *value >= 0.0))
        {
            return Err(ClusteringError::InvalidDistance { index });
        }
        Ok(Self { members, values })
    }

    /// The number of members.
    pub fn len(&self) -> usize {
        self.members
    }

    pub fn is_empty(&self) -> bool {
        self.members == 0
    }

    /// The distance between two members in either order; zero for a member
    /// and itself, and `None` outside the matrix.
    pub fn get(&self, first: usize, second: usize) -> Option<Quantity<f64>> {
        if first >= self.members || second >= self.members {
            return None;
        }
        Some(Quantity::new(self.at(first, second), CANONICAL_LENGTH_UNIT))
    }

    /// Every distance in condensed order.
    pub fn condensed(&self) -> Quantity<&[f64]> {
        Quantity::new(self.values.as_slice(), CANONICAL_LENGTH_UNIT)
    }

    /// Agglomerative clustering of the members under `linkage`.
    ///
    /// Repeatedly merges the two closest clusters until one remains, using
    /// the nearest-neighbor chain algorithm in `O(n^2)` time on a copy of
    /// the distances. Every supported linkage is monotone, so merges come
    /// out in nondecreasing height. Equal distances are resolved
    /// deterministically: a nearest-neighbor search keeps the previous chain
    /// member, then the lowest index, and equal heights keep merge order.
    pub fn dendrogram(&self, linkage: Linkage) -> Dendrogram {
        let joins = nearest_neighbor_chain(self, linkage);
        Dendrogram::from_joins(self.members, joins)
    }

    /// The distance between two members in canonical units.
    pub(super) fn at(&self, first: usize, second: usize) -> f64 {
        match first.cmp(&second) {
            std::cmp::Ordering::Equal => 0.0,
            std::cmp::Ordering::Less => self.values[condensed_index(self.members, first, second)],
            std::cmp::Ordering::Greater => {
                self.values[condensed_index(self.members, second, first)]
            }
        }
    }
}

/// Distances rescaled for summing weighted objectives over the members.
///
/// Every distance is multiplied by one power of two, the largest that keeps
/// `4 n` times the largest distance finite: one unless distances come within
/// that factor of `f64::MAX`. Objectives sum at most `3 n` distances
/// weighted by at most one, so they stay finite with room for rounding.
/// Scaling by a power of two is exact barring underflow, and rounding
/// commutes with it, so scaled objectives compare and tie exactly as
/// unscaled ones would without overflow. Scaled values are not lengths in
/// any unit and only ever compare with one another.
#[derive(Debug, Clone, Copy)]
pub(super) struct ObjectiveDistances<'a> {
    matrix: &'a DistanceMatrix,
    scale: f64,
}

impl<'a> ObjectiveDistances<'a> {
    pub(super) fn new(matrix: &'a DistanceMatrix) -> Self {
        let largest = matrix.values.iter().copied().fold(0.0, f64::max);
        let terms = 4.0 * matrix.members as f64;
        let mut scale = 1.0;
        while largest * scale * terms > f64::MAX {
            scale *= 0.5;
        }
        Self { matrix, scale }
    }

    /// The number of members.
    pub(super) fn len(&self) -> usize {
        self.matrix.members
    }

    /// The scaled distance between two members in either order.
    pub(super) fn at(&self, first: usize, second: usize) -> f64 {
        self.matrix.at(first, second) * self.scale
    }
}

/// How the distance between two clusters follows from member distances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum Linkage {
    /// The closest pair of members; clusters chain through close neighbors.
    Single,
    /// The farthest pair of members; a cut at height `h` leaves clusters
    /// whose members are all within `h` of one another.
    Complete,
    /// The mean over all member pairs (UPGMA), counting each member once
    /// regardless of its weight.
    Average,
}

/// Where to cut a [`Dendrogram`] into flat clusters.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum Cut {
    /// Apply every merge at most this finite, nonnegative height.
    Height(Quantity<f64>),
    /// Apply merges in order until this many clusters remain, from one to
    /// the member count (zero only for no members).
    Count(usize),
}

/// The merge tree of an agglomerative clustering.
///
/// Members are clusters `0..member_count`; merge `k` creates cluster
/// `member_count + k`, the numbering SciPy's linkage matrices use. Merges
/// are in nondecreasing height.
#[derive(Debug, Clone, PartialEq)]
pub struct Dendrogram {
    members: usize,
    merges: Vec<Merge>,
}

/// One merge of two clusters in a [`Dendrogram`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Merge {
    clusters: [usize; 2],
    height: f64,
    size: usize,
    /// One member of each merged cluster.
    members: [usize; 2],
}

impl Merge {
    /// The two merged clusters, lower number first.
    pub const fn clusters(&self) -> [usize; 2] {
        self.clusters
    }

    /// The linkage distance between them in [`CANONICAL_LENGTH_UNIT`].
    pub fn height(&self) -> Quantity<f64> {
        Quantity::new(self.height, CANONICAL_LENGTH_UNIT)
    }

    /// The number of members in the merged cluster.
    pub const fn size(&self) -> usize {
        self.size
    }
}

impl Dendrogram {
    pub fn member_count(&self) -> usize {
        self.members
    }

    /// Every merge, in nondecreasing height.
    pub fn merges(&self) -> &[Merge] {
        &self.merges
    }

    /// The flat cluster of every member after `cut`. Clusters are numbered
    /// in order of their lowest member.
    pub fn cut(&self, cut: Cut) -> Result<Vec<usize>, ClusteringError> {
        let applied = match cut {
            Cut::Height(height) => {
                let height = height.value_in(CANONICAL_LENGTH_UNIT)?;
                if !(height.is_finite() && height >= 0.0) {
                    return Err(ClusteringError::InvalidHeight);
                }
                self.merges
                    .iter()
                    .take_while(|merge| merge.height <= height)
                    .count()
            }
            Cut::Count(count) => {
                let valid = if self.members == 0 {
                    count == 0
                } else {
                    (1..=self.members).contains(&count)
                };
                if !valid {
                    return Err(ClusteringError::InvalidClusterCount {
                        requested: count,
                        members: self.members,
                    });
                }
                self.members.saturating_sub(count)
            }
        };
        let mut sets = DisjointSets::new(self.members);
        for merge in &self.merges[..applied] {
            sets.union(merge.members[0], merge.members[1]);
        }
        let mut labels = vec![usize::MAX; self.members];
        let mut assignments = Vec::with_capacity(self.members);
        let mut clusters = 0;
        for member in 0..self.members {
            let root = sets.find(member);
            if labels[root] == usize::MAX {
                labels[root] = clusters;
                clusters += 1;
            }
            assignments.push(labels[root]);
        }
        Ok(assignments)
    }

    /// Numbers merges found in any order, given as joined members and
    /// heights, after a stable sort by height.
    fn from_joins(members: usize, mut joins: Vec<(usize, usize, f64)>) -> Self {
        joins.sort_by(|a, b| a.2.total_cmp(&b.2));
        let mut sets = DisjointSets::new(members);
        let mut cluster = (0..members).collect::<Vec<_>>();
        let mut size = vec![1; members];
        let merges = joins
            .into_iter()
            .enumerate()
            .map(|(index, (first, second, height))| {
                let (a, b) = (sets.find(first), sets.find(second));
                let mut clusters = [cluster[a], cluster[b]];
                clusters.sort_unstable();
                let root = sets.union(a, b);
                cluster[root] = members + index;
                size[root] = size[a] + size[b];
                Merge {
                    clusters,
                    height,
                    size: size[root],
                    members: [first, second],
                }
            })
            .collect();
        Self { members, merges }
    }
}

impl Ensemble {
    /// Fitted RMSD between every pair of members, after a
    /// [`fit`](crate::alignment::fit) over `atoms` with the weights and
    /// periodic policy of `options`.
    ///
    /// `atoms` is a selection or a
    /// [`SymmetricSelection`](crate::alignment::SymmetricSelection) of this
    /// ensemble's layout, so each distance is the same whichever member
    /// moves; an [`AtomCorrespondence`](crate::alignment::AtomCorrespondence)
    /// is rejected. Each unordered pair is fitted once, lower member moving,
    /// so the cost is `n (n - 1) / 2` fits (times symmetry permutations). A
    /// failed fit names both members.
    pub fn pairwise_fitted_rmsd<'a>(
        &self,
        atoms: impl Into<FitAtoms<'a>>,
        options: AlignmentOptions<'_>,
    ) -> Result<DistanceMatrix, ClusteringError> {
        let atoms = undirected(atoms.into())?;
        let views = self
            .iter()
            .map(|member| member.as_model_view())
            .collect::<Vec<ModelView<'_>>>();
        let members = views.len();
        let expected = condensed_len(members)?;
        let Some(first) = views.first() else {
            return Ok(DistanceMatrix {
                members,
                values: Vec::new(),
            });
        };
        let sets = atoms.pair_sets(*first, *first)?;
        check_fit_count(atoms.len())?;
        let weights = NormalizedWeights::new(options.weighting(), atoms.len())?;
        let mut values = Vec::with_capacity(expected);
        for (index, moving) in views.iter().enumerate() {
            for (offset, reference) in views[index + 1..].iter().enumerate() {
                let fit = best_kabsch(
                    *moving,
                    *reference,
                    &sets,
                    weights,
                    options.periodic_policy(),
                )
                .map_err(|source| ClusteringError::PairFit {
                    first: index,
                    second: index + 1 + offset,
                    source: Box::new(source),
                })?;
                values.push(fit.rmsd().into_value());
            }
        }
        Ok(DistanceMatrix { members, values })
    }

    /// Hierarchical clustering of the members on `distances`, cut into flat
    /// clusters, each represented by its weighted medoid.
    ///
    /// The dendrogram depends on distances only; see
    /// [`DistanceMatrix::dendrogram`]. The representative of a cluster is
    /// the member minimizing the weight-summed distance to the cluster's
    /// members, `sum_j w_j d(i, j)`, with ties to the lowest index. Clusters
    /// are numbered in order of their lowest member. Fails when the matrix
    /// covers a different member count.
    pub fn cluster_hierarchical(
        &self,
        distances: &DistanceMatrix,
        linkage: Linkage,
        cut: Cut,
    ) -> Result<EnsemblePartition, ClusteringError> {
        if distances.len() != self.len() {
            return Err(ClusteringError::MatrixMismatch {
                expected: self.len(),
                actual: distances.len(),
            });
        }
        let assignments = distances.dendrogram(linkage).cut(cut)?;
        let representatives = self.weighted_medoids(distances, &assignments);
        Ok(EnsemblePartition {
            assignments,
            representatives,
        })
    }

    /// The weighted medoid of every cluster of `assignments`, whose clusters
    /// are numbered densely from zero.
    fn weighted_medoids(&self, distances: &DistanceMatrix, assignments: &[usize]) -> Vec<usize> {
        let weights = self.scaled_weights();
        let distances = ObjectiveDistances::new(distances);
        let clusters = assignments.iter().max().map_or(0, |last| last + 1);
        let mut members = vec![Vec::new(); clusters];
        for (member, cluster) in assignments.iter().enumerate() {
            members[*cluster].push(member);
        }
        members
            .iter()
            .map(|cluster| {
                let mut best = (cluster[0], f64::INFINITY);
                for candidate in cluster {
                    let cost = cluster
                        .iter()
                        .map(|other| weights[*other] * distances.at(*candidate, *other))
                        .sum::<f64>();
                    if cost < best.1 {
                        best = (*candidate, cost);
                    }
                }
                best.0
            })
            .collect()
    }
}

/// The number of unordered pairs of `members` members.
fn condensed_len(members: usize) -> Result<usize, ClusteringError> {
    members
        .checked_mul(members.saturating_sub(1))
        .map(|pairs| pairs / 2)
        .ok_or(ClusteringError::MatrixTooLarge { members })
}

/// The condensed position of the pair `(first, second)` with `first < second`.
fn condensed_index(members: usize, first: usize, second: usize) -> usize {
    first * members - first * (first + 1) / 2 + (second - first - 1)
}

/// Joins of the nearest-neighbor chain algorithm (Müllner 2011, sec. 3.2),
/// as `(member, member, height)` in the order they are found.
///
/// Slot `s` always holds the cluster containing member `s`; a merge keeps
/// the lower slot and retires the higher one.
fn nearest_neighbor_chain(
    distances: &DistanceMatrix,
    linkage: Linkage,
) -> Vec<(usize, usize, f64)> {
    let members = distances.members;
    let mut working = distances.values.clone();
    let at = |working: &[f64], first: usize, second: usize| {
        let (low, high) = if first < second {
            (first, second)
        } else {
            (second, first)
        };
        working[condensed_index(members, low, high)]
    };
    let mut active = vec![true; members];
    let mut size = vec![1_usize; members];
    let mut joins = Vec::with_capacity(members.saturating_sub(1));
    let mut chain: Vec<usize> = Vec::new();
    for _ in 1..members {
        if chain.is_empty() {
            chain.push(
                active
                    .iter()
                    .position(|alive| *alive)
                    .expect("two clusters remain"),
            );
        }
        let (a, b) = loop {
            let current = *chain.last().expect("nonempty chain");
            let previous = chain.len().checked_sub(2).map(|index| chain[index]);
            // Keeping the previous chain member on ties guarantees progress.
            let mut nearest = previous.map(|slot| (slot, at(&working, current, slot)));
            for (slot, alive) in active.iter().enumerate() {
                if !*alive || slot == current {
                    continue;
                }
                let distance = at(&working, current, slot);
                if nearest.is_none_or(|(_, best)| distance < best) {
                    nearest = Some((slot, distance));
                }
            }
            let (next, _) = nearest.expect("another active cluster");
            if Some(next) == previous {
                chain.truncate(chain.len() - 2);
                break (current, next);
            }
            chain.push(next);
        };
        let (low, high) = if a < b { (a, b) } else { (b, a) };
        joins.push((low, high, at(&working, low, high)));
        for (slot, alive) in active.iter().enumerate() {
            if !*alive || slot == low || slot == high {
                continue;
            }
            let (to_low, to_high) = (at(&working, low, slot), at(&working, high, slot));
            let merged = match linkage {
                Linkage::Single => to_low.min(to_high),
                Linkage::Complete => to_low.max(to_high),
                Linkage::Average => {
                    // The size-weighted mean, stepped from one input towards
                    // the other so no intermediate exceeds the larger input;
                    // summing scaled inputs first overflows near `f64::MAX`.
                    // Clamping keeps rounding from leaving the inputs' range.
                    let (low_size, high_size) = (size[low] as f64, size[high] as f64);
                    let fraction = high_size / (low_size + high_size);
                    (to_low + (to_high - to_low) * fraction)
                        .clamp(to_low.min(to_high), to_low.max(to_high))
                }
            };
            let (first, second) = if low < slot { (low, slot) } else { (slot, low) };
            working[condensed_index(members, first, second)] = merged;
        }
        size[low] += size[high];
        active[high] = false;
    }
    joins
}

/// Union-find over members with path halving.
struct DisjointSets {
    parent: Vec<usize>,
}

impl DisjointSets {
    fn new(members: usize) -> Self {
        Self {
            parent: (0..members).collect(),
        }
    }

    fn find(&mut self, mut member: usize) -> usize {
        while self.parent[member] != member {
            self.parent[member] = self.parent[self.parent[member]];
            member = self.parent[member];
        }
        member
    }

    /// Joins two sets and returns the root of the union, the lower root.
    fn union(&mut self, first: usize, second: usize) -> usize {
        let (a, b) = (self.find(first), self.find(second));
        let (root, child) = if a < b { (a, b) } else { (b, a) };
        self.parent[child] = root;
        root
    }
}

#[cfg(test)]
mod tests {
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
        let mut clusters = (0..matrix.len())
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
}
