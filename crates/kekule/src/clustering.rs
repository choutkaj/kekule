//! Partitions of an [`Ensemble`]'s samples into states of similar
//! structures.
//!
//! Clustering never changes the ensemble. It returns an
//! [`EnsemblePartition`] that assigns every sample to one state and names one
//! representative sample per state. [`Ensemble::regroup`] applies it without
//! changing any weight, so populations follow the samples, and
//! [`Ensemble::collapse`] then keeps one sample per state, carrying the
//! state's population.
//!
//! [`Ensemble::cluster_hierarchical`] cuts the [`Dendrogram`] of a
//! [`DistanceMatrix`], such as [`Ensemble::pairwise_fitted_rmsd`], at a height or
//! into a number of clusters, and represents each cluster by its weighted
//! medoid. Complete linkage cut at height `h` bounds every cluster's
//! diameter by `h`.
//!
//! [`Ensemble::cluster_k_medoids`] finds a given number of clusters directly,
//! choosing medoid samples that minimize the weight-summed distance from
//! every sample to its nearest medoid, so heavy samples attract medoids.
//!
//! [`Ensemble::cluster_by_radius`] is leader clustering on fitted RMSD. It
//! needs no matrix, so it scales to larger ensembles. With
//! a [`SymmetricSelection`](crate::alignment::SymmetricSelection) and a small
//! radius it merges rotamers, copies of one conformer that differ only by
//! exchanged equivalent atoms. Weight the conformers by energy first; the
//! merged state keeps its representative's energy and gains the population
//! of every rotamer, its degeneracy:
//!
//! ```
//! use std::error::Error;
//!
//! use kekule::alignment::SymmetricSelection;
//! use kekule::geometry::Point3;
//! use kekule::structure::{Ensemble, EnsembleSample, EnsembleState, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM, KELVIN, KILOJOULE_PER_MOLE};
//!
//! // Acetate, then a rotamer with its oxygens exchanged, then a distinct
//! // conformer 5 kJ/mol higher, each its own state.
//! let topology = kekule::smiles::to_topology("CC(=O)[O-]")?;
//! let mut ensemble = Ensemble::new(topology.clone());
//! for (o2, o3, energy) in [(1.1, -1.1, 0.0), (-1.1, 1.1, 0.0), (1.1, -0.3, 5.0)] {
//!     let points = [
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.2, o2, 0.5),
//!         Point3::new(2.2, o3, 0.5),
//!     ];
//!     let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//!     let state = EnsembleState::new().with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))?;
//!     ensemble.push_state(EnsembleSample::new(positions, 1.0)?, state)?;
//! }
//! ensemble.assign_boltzmann_populations(Quantity::new(298.15, KELVIN))?;
//!
//! let atoms = SymmetricSelection::new(&AtomSelection::all(topology.clone()))?;
//! let radius = Quantity::new(0.125, ANGSTROM);
//! let rotamers = ensemble.cluster_by_radius(radius, &atoms, Default::default())?;
//! assert_eq!(rotamers.assignments(), [0, 0, 1]);
//!
//! ensemble.regroup(rotamers)?;
//! assert_eq!(ensemble.state_count(), 2);
//! assert!(ensemble.state(0).unwrap().energy().is_some());
//! let conformers = ensemble.collapse()?;
//! assert_eq!(conformers.len(), 2);
//! assert_eq!(conformers.get(0).unwrap().weight(), ensemble.state(0).unwrap().population());
//! # Ok::<(), Box<dyn Error>>(())
//! ```

use std::fmt;

mod hierarchical;
mod medoids;
pub use hierarchical::{Cut, Dendrogram, DistanceMatrix, Linkage, Merge};

use crate::alignment::{
    best_kabsch, check_fit_count, AlignmentError, AlignmentOptions, FitAtoms, NormalizedWeights,
};
use crate::structure::{Ensemble, EnsembleError, EnsemblePartition, ModelView, StructureError};
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

/// The order in which leader clustering visits members.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum ClusterOrder {
    /// Descending sample weight, so the heaviest sample of a cluster leads
    /// it.
    #[default]
    Weight,
    /// Ascending energy of each sample's state, so a sample of the
    /// lowest-energy state leads its cluster; samples of one state keep
    /// sample order.
    Energy,
}

/// Options for [`Ensemble::cluster_by_radius`].
#[derive(Debug, Clone, Copy, Default)]
pub struct ClusteringOptions<'a> {
    alignment: AlignmentOptions<'a>,
    order: ClusterOrder,
}

impl<'a> ClusteringOptions<'a> {
    /// Fit weighting and periodic policy of every pairwise fit.
    pub const fn alignment(&self) -> AlignmentOptions<'a> {
        self.alignment
    }

    /// Sets [`Self::alignment`].
    #[must_use]
    pub const fn with_alignment(mut self, alignment: AlignmentOptions<'a>) -> Self {
        self.alignment = alignment;
        self
    }

    /// The order in which members are visited.
    pub const fn order(&self) -> ClusterOrder {
        self.order
    }

    /// Sets [`Self::order`].
    #[must_use]
    pub const fn with_order(mut self, order: ClusterOrder) -> Self {
        self.order = order;
        self
    }
}

impl Ensemble {
    /// Leader clustering on fitted RMSD within `radius`, a finite,
    /// nonnegative length.
    ///
    /// Members are visited in the [`ClusterOrder`] of `options`, by default
    /// descending weight; ties keep member order. Each member joins the
    /// cluster whose representative is nearest by RMSD after a
    /// [`fit`](crate::alignment::fit) over `atoms` with the fit weights of
    /// `options`, by default uniform, if
    /// that RMSD is at most `radius`; ties keep the earlier cluster.
    /// Otherwise it founds a new cluster and represents it. Every member is
    /// therefore within `radius` of its representative, and representatives
    /// are farther than `radius` from one another. Clusters are numbered in
    /// order of founding.
    ///
    /// `atoms` is a selection or a
    /// [`SymmetricSelection`](crate::alignment::SymmetricSelection) of this
    /// ensemble's layout, so each distance is the same whichever member
    /// moves; an [`AtomCorrespondence`](crate::alignment::AtomCorrespondence)
    /// is rejected. Each member is fitted to at most every representative, so
    /// the cost grows with members times clusters (times symmetry
    /// permutations). An empty ensemble gives an empty partition.
    ///
    /// # Errors
    ///
    /// - [`ClusteringError::Unit`] when `radius` is not a length, and
    ///   [`ClusteringError::InvalidRadius`] when it is negative or not finite.
    /// - [`ClusteringError::DirectedCorrespondence`] when `atoms` is an
    ///   [`AtomCorrespondence`](crate::alignment::AtomCorrespondence).
    /// - [`ClusteringError::Alignment`] when `atoms` does not belong to this
    ///   ensemble's layout, has fewer than three pairs, or the explicit
    ///   weights of `options` are invalid.
    ///   A failed fit is wrapped in [`AlignmentError::Item`] naming the
    ///   member being placed.
    /// - [`ClusteringError::Ensemble`] wrapping
    ///   [`EnsembleError::MissingEnergy`] when ordering by
    ///   [`ClusterOrder::Energy`] and a state has no energy.
    pub fn cluster_by_radius<'a>(
        &self,
        radius: Quantity<f64>,
        atoms: impl Into<FitAtoms<'a>>,
        options: ClusteringOptions<'_>,
    ) -> Result<EnsemblePartition, ClusteringError> {
        let atoms = undirected(atoms.into())?;
        let radius = radius.value_in(CANONICAL_LENGTH_UNIT)?;
        if !(radius.is_finite() && radius >= 0.0) {
            return Err(ClusteringError::InvalidRadius);
        }
        let order = match options.order {
            ClusterOrder::Weight => self.ranked_by_weight(),
            ClusterOrder::Energy => {
                let energies = self.state_energies()?;
                let states = self.partition().assignments();
                let mut order = (0..self.len()).collect::<Vec<_>>();
                order.sort_by(|&a, &b| energies[states[a]].total_cmp(&energies[states[b]]));
                order
            }
        };
        let views = self
            .iter()
            .map(|member| member.as_model_view())
            .collect::<Vec<ModelView<'_>>>();
        let Some(first) = views.first() else {
            return Ok(EnsemblePartition::state_per_sample(0));
        };
        let sets = atoms.pair_sets(*first, *first)?;
        check_fit_count(atoms.pair_count())?;
        let weights = NormalizedWeights::new(options.alignment.weighting(), atoms.pair_count())?;
        let policy = options.alignment.periodic_policy();

        let mut assignments = vec![0; views.len()];
        let mut representatives: Vec<usize> = Vec::new();
        for member in order {
            let mut nearest: Option<(usize, f64)> = None;
            for (cluster, representative) in representatives.iter().enumerate() {
                let rmsd = best_kabsch(
                    views[member],
                    views[*representative],
                    &sets,
                    weights,
                    policy,
                )
                .map_err(|source| AlignmentError::Item {
                    index: member,
                    source: Box::new(source),
                })?
                .rmsd()
                .into_value();
                if rmsd <= radius && nearest.is_none_or(|(_, best)| rmsd < best) {
                    nearest = Some((cluster, rmsd));
                }
            }
            assignments[member] = match nearest {
                Some((cluster, _)) => cluster,
                None => {
                    representatives.push(member);
                    representatives.len() - 1
                }
            };
        }
        Ok(EnsemblePartition::from_parts_unchecked(
            assignments,
            representatives,
        ))
    }
}

/// Clustering distances must not depend on which member moves. Fitting
/// member `a` onto `b` pairs `a`'s moving atoms with `b`'s reference atoms,
/// so a correspondence pairing different atoms of one layout gives a
/// different RMSD in each direction, and results would follow member order.
/// Selections pair atoms with themselves and symmetry permutations form a
/// group, so both give the same RMSD either way.
fn undirected(atoms: FitAtoms<'_>) -> Result<FitAtoms<'_>, ClusteringError> {
    match atoms {
        FitAtoms::Correspondence(_) => Err(ClusteringError::DirectedCorrespondence),
        atoms => Ok(atoms),
    }
}

/// Failure of an ensemble clustering operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ClusteringError {
    /// A clustering radius is not finite and nonnegative.
    InvalidRadius,
    /// A distance matrix covers a different member count than the ensemble.
    MatrixMismatch {
        /// The ensemble's member count.
        expected: usize,
        /// The matrix's member count.
        actual: usize,
    },
    /// Condensed distances have the wrong length for their member count.
    CondensedLength {
        /// The stated member count.
        members: usize,
        /// The condensed length for that member count.
        expected: usize,
        /// The number of distances supplied.
        actual: usize,
    },
    /// A distance is not finite and nonnegative; `index` is condensed.
    InvalidDistance {
        /// The condensed index of the invalid distance.
        index: usize,
    },
    /// The pairwise distances of this many members are not addressable.
    MatrixTooLarge {
        /// The requested member count.
        members: usize,
    },
    /// Fitting one pair of members failed.
    PairFit {
        /// The first member of the pair.
        first: usize,
        /// The second member of the pair.
        second: usize,
        /// The fit failure.
        source: Box<AlignmentError>,
    },
    /// A dendrogram cut height is not finite and nonnegative.
    InvalidHeight,
    /// A requested cluster count is zero or exceeds the member count.
    InvalidClusterCount {
        /// The requested cluster count.
        requested: usize,
        /// The number of members.
        members: usize,
    },
    /// Distances between members were requested over an
    /// [`AtomCorrespondence`](crate::alignment::AtomCorrespondence), whose
    /// fits depend on which member moves.
    DirectedCorrespondence,
    /// Fitting members or validating the fit atoms failed; see the source.
    Alignment(Box<AlignmentError>),
    /// Reading or building the ensemble failed; see the source.
    Ensemble(Box<EnsembleError>),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
}

impl fmt::Display for ClusteringError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidRadius => {
                formatter.write_str("clustering radius must be finite and nonnegative")
            }
            Self::MatrixMismatch { expected, actual } => write!(
                formatter,
                "distance matrix covers {actual} members, but the ensemble has {expected}"
            ),
            Self::CondensedLength {
                members,
                expected,
                actual,
            } => write!(
                formatter,
                "{members} members require {expected} condensed distances, but received {actual}"
            ),
            Self::InvalidDistance { index } => write!(
                formatter,
                "condensed distance {index} must be finite and nonnegative"
            ),
            Self::MatrixTooLarge { members } => write!(
                formatter,
                "pairwise distances of {members} members are not addressable"
            ),
            Self::PairFit {
                first,
                second,
                source,
            } => write!(
                formatter,
                "cannot fit members {first} and {second}: {source}"
            ),
            Self::InvalidHeight => {
                formatter.write_str("dendrogram cut height must be finite and nonnegative")
            }
            Self::InvalidClusterCount { requested, members } => write!(
                formatter,
                "cannot cut {members} members into {requested} clusters"
            ),
            Self::DirectedCorrespondence => formatter.write_str(
                "member distances require a selection, not a directed atom correspondence",
            ),
            Self::Alignment(_) => formatter.write_str("cannot fit members"),
            Self::Ensemble(_) => formatter.write_str("invalid ensemble"),
            Self::Unit(_) => formatter.write_str("invalid quantity unit"),
        }
    }
}

impl std::error::Error for ClusteringError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Alignment(error) | Self::PairFit { source: error, .. } => Some(error.as_ref()),
            Self::Ensemble(error) => Some(error.as_ref()),
            Self::Unit(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<AlignmentError> for ClusteringError {
    fn from(error: AlignmentError) -> Self {
        Self::Alignment(Box::new(error))
    }
}

impl From<EnsembleError> for ClusteringError {
    fn from(error: EnsembleError) -> Self {
        Self::Ensemble(Box::new(error))
    }
}

impl From<StructureError> for ClusteringError {
    fn from(error: StructureError) -> Self {
        Self::Ensemble(Box::new(EnsembleError::Structure(error)))
    }
}

impl From<UnitError> for ClusteringError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
