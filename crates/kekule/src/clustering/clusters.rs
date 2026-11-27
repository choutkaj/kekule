//! A grouping of a structure set's geometries into clusters: [`Clustering`].

use crate::diagnostics::{Countable, IndexOutOfRange};
use crate::statistics::{WeightSource, Weights, WeightsError};
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT};

use super::ClusteringError;

/// An optional cluster per geometry of a structure set, and one
/// representative geometry per cluster.
///
/// Like [`Weights`], it is a reading beside its set, indexed by geometry: a
/// geometry without a cluster is in none. Clusters are numbered from zero,
/// and each contains its representative, so none is empty. The clustering
/// methods of a structure set return one
/// ([`StructureSet::cluster_by_radius`](crate::structure::StructureSet::cluster_by_radius)),
/// and a host builds its own with [`Self::new`].
///
/// A cluster's population is the exact sum of its samples' weights, so it
/// is the sum of its parts' populations; the populations of the clusters
/// are normalized among themselves ([`Self::log_populations`]).
///
/// ```
/// use kekule::clustering::Clustering;
/// use kekule::statistics::Weights;
///
/// // Geometries 0 and 2 in cluster 0, represented by 2; geometry 1 in none;
/// // geometry 3 alone in cluster 1.
/// let clustering = Clustering::new(vec![Some(0), None, Some(0), Some(1)], vec![2, 3])?;
/// assert_eq!(clustering.members(0).collect::<Vec<_>>(), [0, 2]);
/// assert_eq!(clustering.sizes(), [2, 1]);
///
/// let weights = Weights::equal(4);
/// let populations = clustering.log_populations(&weights)?;
/// assert!((populations[0].unwrap().exp() - 2.0 / 3.0).abs() < 1e-15);
/// let (representatives, collapsed) = clustering.collapse(&weights)?;
/// assert_eq!(representatives, [2, 3]);
/// assert_eq!(collapsed.log_weight(0), populations[0]);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Clustering {
    cluster_of: Vec<Option<usize>>,
    representatives: Vec<usize>,
}

impl Clustering {
    /// Validates a clustering built elsewhere: `cluster_of[g]` is the
    /// cluster of geometry `g`, if any, and `representatives[c]` is a
    /// geometry of cluster `c`.
    ///
    /// # Errors
    ///
    /// [`ClusteringError::ClusterOutOfRange`] for a geometry in a cluster
    /// without a representative, [`ClusteringError::IndexOutOfRange`] for a
    /// representative outside the geometries, and
    /// [`ClusteringError::RepresentativeNotMember`] for a representative of
    /// another cluster or of none.
    pub fn new(
        cluster_of: Vec<Option<usize>>,
        representatives: Vec<usize>,
    ) -> Result<Self, ClusteringError> {
        let clusters = representatives.len();
        for (geometry, cluster) in cluster_of.iter().enumerate() {
            if let Some(cluster) = *cluster {
                if cluster >= clusters {
                    return Err(ClusteringError::ClusterOutOfRange {
                        geometry,
                        cluster,
                        clusters,
                    });
                }
            }
        }
        for (cluster, &geometry) in representatives.iter().enumerate() {
            match cluster_of.get(geometry) {
                None => {
                    return Err(ClusteringError::IndexOutOfRange(IndexOutOfRange::new(
                        Countable::Geometry,
                        geometry,
                        cluster_of.len(),
                    )))
                }
                Some(assigned) if *assigned != Some(cluster) => {
                    return Err(ClusteringError::RepresentativeNotMember { cluster, geometry })
                }
                Some(_) => {}
            }
        }
        Ok(Self {
            cluster_of,
            representatives,
        })
    }

    /// A clustering whose invariants the caller has established.
    pub(crate) fn from_parts_unchecked(
        cluster_of: Vec<Option<usize>>,
        representatives: Vec<usize>,
    ) -> Self {
        debug_assert!(Self::new(cluster_of.clone(), representatives.clone()).is_ok());
        Self {
            cluster_of,
            representatives,
        }
    }

    /// The number of geometries the clustering covers.
    pub fn len(&self) -> usize {
        self.cluster_of.len()
    }

    /// Whether the clustering covers no geometry.
    pub fn is_empty(&self) -> bool {
        self.cluster_of.is_empty()
    }

    /// The number of clusters.
    pub fn cluster_count(&self) -> usize {
        self.representatives.len()
    }

    /// The cluster of geometry `geometry`, or `None` when it is in none or
    /// lies outside the clustering.
    pub fn cluster_of(&self, geometry: usize) -> Option<usize> {
        self.cluster_of.get(geometry).copied().flatten()
    }

    /// The geometries of cluster `cluster`, by increasing index; none for a
    /// cluster that does not exist.
    pub fn members(&self, cluster: usize) -> impl Iterator<Item = usize> + '_ {
        self.cluster_of
            .iter()
            .enumerate()
            .filter(move |(_, assigned)| **assigned == Some(cluster))
            .map(|(geometry, _)| geometry)
    }

    /// The representative geometry of cluster `cluster`, or `None` when it
    /// does not exist.
    pub fn representative(&self, cluster: usize) -> Option<usize> {
        self.representatives.get(cluster).copied()
    }

    /// The representative geometry of every cluster, in cluster order.
    pub fn representatives(&self) -> &[usize] {
        &self.representatives
    }

    /// The number of geometries of every cluster, in cluster order.
    pub fn sizes(&self) -> Vec<usize> {
        let mut sizes = vec![0; self.cluster_count()];
        for cluster in self.cluster_of.iter().flatten() {
            sizes[*cluster] += 1;
        }
        sizes
    }

    /// The natural logarithm of each cluster's probability among the
    /// clusters, in cluster order: the log-sum-exp of its samples' log
    /// weights, normalized so that the clusters' probabilities sum to one.
    /// Samples in no cluster are left out of the total. `None` for a cluster
    /// without samples.
    ///
    /// # Errors
    ///
    /// [`ClusteringError::Weights`] with
    /// [`WeightsError::CoverageMismatch`] when the weights cover another
    /// number of geometries.
    pub fn log_populations(&self, weights: &Weights) -> Result<Vec<Option<f64>>, ClusteringError> {
        weights.check_len(self.len())?;
        let clustered = || {
            self.cluster_of
                .iter()
                .zip(weights.log_weights())
                .filter_map(|(cluster, log_weight)| cluster.zip(*log_weight))
        };
        // Shifted by the largest log weight in a cluster, normalizations keep
        // their digits under a large common offset, and a heavier sample in
        // no cluster takes none.
        let shift = clustered()
            .map(|(_, log_weight)| log_weight)
            .reduce(f64::max)
            .unwrap_or(0.0);
        let mut values = vec![Vec::new(); self.cluster_count()];
        for (cluster, log_weight) in clustered() {
            values[cluster].push(log_weight - shift);
        }
        let sums = values
            .into_iter()
            .map(|values| {
                (!values.is_empty()).then(|| crate::statistics::log_sum_exp(values.into_iter()))
            })
            .collect::<Vec<_>>();
        let total = crate::statistics::log_sum_exp(sums.iter().flatten().copied());
        Ok(sums
            .into_iter()
            .map(|sum| sum.map(|sum| sum - total))
            .collect())
    }

    /// The free energy of each cluster, `−RT ln(p / p_max)` at the weights'
    /// temperature in [`CANONICAL_ENERGY_UNIT`], relative to the most
    /// populated cluster, which is at 0; `None` for a cluster without
    /// samples.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::log_populations`], and
    /// [`ClusteringError::Weights`] with [`WeightsError::MissingTemperature`]
    /// for weights without a temperature.
    pub fn free_energies(
        &self,
        weights: &Weights,
    ) -> Result<Quantity<Vec<Option<f64>>>, ClusteringError> {
        let populations = self.log_populations(weights)?;
        let temperature = weights
            .temperature()
            .ok_or(WeightsError::MissingTemperature)?;
        let thermal_energy = *crate::statistics::canonical(temperature)?
            .thermal_energy()
            .value();
        let highest = populations.iter().flatten().copied().reduce(f64::max);
        Ok(Quantity::new(
            populations
                .iter()
                .map(|population| {
                    population
                        .zip(highest)
                        .map(|(population, highest)| (highest - population) * thermal_energy)
                })
                .collect(),
            CANONICAL_ENERGY_UNIT,
        ))
    }

    /// The clusters by descending population; equal populations keep
    /// cluster order, and clusters without samples come last.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::log_populations`].
    pub fn ranked_by_population(&self, weights: &Weights) -> Result<Vec<usize>, ClusteringError> {
        let populations = self.log_populations(weights)?;
        Ok(ranked(&populations))
    }

    /// The fewest clusters, by descending population, whose populations
    /// together reach `fraction` of the total, in
    /// [`Self::ranked_by_population`] order: for example the clusters that
    /// hold 95 % of a Boltzmann population. `fraction` lies in `(0, 1]`; `1`
    /// takes every cluster.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::log_populations`], and
    /// [`ClusteringError::Weights`] with [`WeightsError::InvalidFraction`] for
    /// a fraction outside `(0, 1]`.
    pub fn covering(
        &self,
        weights: &Weights,
        fraction: f64,
    ) -> Result<Vec<usize>, ClusteringError> {
        crate::statistics::check_fraction(fraction)?;
        let populations = self.log_populations(weights)?;
        let mut ranked = ranked(&populations);
        let shares = ranked
            .iter()
            .map(|&cluster| populations[cluster].map_or(0.0, f64::exp))
            .collect::<Vec<_>>();
        ranked.truncate(crate::statistics::covered_count(fraction, &shares));
        Ok(ranked)
    }

    /// The clustering of a set of the geometries `indices`, in that order,
    /// as [`StructureSet::select`](crate::structure::StructureSet::select)
    /// copies them. Every cluster with a listed geometry stays, in cluster
    /// order, represented by the first copy of its representative or, when
    /// that is not listed, by its first listed geometry. A geometry listed
    /// twice is in its cluster twice.
    ///
    /// # Errors
    ///
    /// [`ClusteringError::IndexOutOfRange`] for the first index outside the
    /// clustering.
    pub fn select(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, ClusteringError> {
        let indices = indices
            .into_iter()
            .map(|geometry| self.check(geometry))
            .collect::<Result<Vec<_>, _>>()?;
        let mut renumbered = vec![None; self.cluster_count()];
        for geometry in &indices {
            if let Some(cluster) = self.cluster_of[*geometry] {
                renumbered[cluster] = Some(0);
            }
        }
        let mut count = 0;
        for cluster in renumbered.iter_mut().flatten() {
            *cluster = count;
            count += 1;
        }
        let cluster_of = indices
            .iter()
            .map(|geometry| self.cluster_of[*geometry].and_then(|cluster| renumbered[cluster]))
            .collect::<Vec<_>>();
        let mut representatives = vec![None; count];
        for (copy, geometry) in indices.iter().enumerate() {
            if let Some(old) = self.cluster_of[*geometry] {
                if self.representatives[old] == *geometry {
                    representatives[cluster_of[copy].expect("a listed cluster stays")]
                        .get_or_insert(copy);
                }
            }
        }
        for (copy, cluster) in cluster_of.iter().enumerate() {
            if let Some(cluster) = cluster {
                representatives[*cluster].get_or_insert(copy);
            }
        }
        let representatives = representatives
            .into_iter()
            .map(|copy| copy.expect("every kept cluster has a listed geometry"))
            .collect();
        Ok(Self::from_parts_unchecked(cluster_of, representatives))
    }

    /// The geometries of the clusters `clusters`, in that order, each
    /// cluster's by increasing index: what selecting whole clusters keeps.
    /// A cluster listed twice is listed twice.
    ///
    /// # Errors
    ///
    /// [`ClusteringError::IndexOutOfRange`] for the first cluster that does
    /// not exist.
    pub fn geometries_of(&self, clusters: &[usize]) -> Result<Vec<usize>, ClusteringError> {
        let mut geometries = Vec::new();
        for &cluster in clusters {
            if cluster >= self.cluster_count() {
                return Err(ClusteringError::IndexOutOfRange(IndexOutOfRange::new(
                    Countable::Cluster,
                    cluster,
                    self.cluster_count(),
                )));
            }
            geometries.extend(self.members(cluster));
        }
        Ok(geometries)
    }

    /// Each cluster's representative, in cluster order, with weights over
    /// just those geometries: each the cluster's population, so population
    /// shares are kept, given at the weights' temperature. A cluster without
    /// samples gives no sample. Select the representatives of the set with
    /// [`StructureSet::select`](crate::structure::StructureSet::select).
    ///
    /// It coarse-grains, losing the spread within each cluster.
    ///
    /// # Errors
    ///
    /// The errors of [`Self::log_populations`].
    pub fn collapse(&self, weights: &Weights) -> Result<(Vec<usize>, Weights), ClusteringError> {
        let populations = self.log_populations(weights)?;
        let collapsed = Weights::from_parts_unchecked(
            populations,
            WeightSource::Given {
                temperature: weights.temperature(),
            },
        );
        Ok((self.representatives.clone(), collapsed))
    }

    /// The cluster of every geometry.
    pub(crate) fn assignments(&self) -> &[Option<usize>] {
        &self.cluster_of
    }

    fn check(&self, geometry: usize) -> Result<usize, ClusteringError> {
        if geometry >= self.len() {
            return Err(ClusteringError::IndexOutOfRange(IndexOutOfRange::new(
                Countable::Geometry,
                geometry,
                self.len(),
            )));
        }
        Ok(geometry)
    }
}

/// Cluster indices by descending population, populations missing last; ties
/// keep cluster order.
fn ranked(populations: &[Option<f64>]) -> Vec<usize> {
    let mut clusters = (0..populations.len()).collect::<Vec<_>>();
    let key = |cluster: usize| populations[cluster].unwrap_or(f64::NEG_INFINITY);
    clusters.sort_by(|&a, &b| key(b).total_cmp(&key(a)));
    clusters
}
