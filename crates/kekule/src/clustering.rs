//! Partitions of the samples of an [`Ensemble`], or of the geometries of a
//! [`StructureSet`], into states of similar structures.
//!
//! An ensemble clusters the samples of its distribution, leaving the other
//! geometries out; a structure set clusters every geometry. Results are
//! indexed by member position, the `i`th member being the `i`th sample
//! [`Distribution::samples`](crate::structure::Distribution::samples) lists,
//! or geometry `i` of a set: the [`DistanceMatrix`] of
//! [`Ensemble::pairwise_fitted_rmsd`] and each partition, whose
//! assignments and representatives are member positions. Errors name
//! geometry indices. [`pairwise_fitted_rmsd`] fits the geometries of any
//! [`StructureSetRead`](crate::structure::StructureSetRead), such as a set
//! read from a file, through the same loop, and an ensemble or a set of
//! either storage clusters a matrix.
//!
//! Clustering never changes the ensemble or the set. It returns a
//! [`Partition`] that assigns every member to one state and names one
//! representative member per state.
//! [`DistributionMut::regroup`](crate::structure::DistributionMut::regroup)
//! installs it without changing any weight, so populations follow the
//! samples, and [`Ensemble::collapse`] then keeps one sample per state,
//! carrying the state's population. A partition of a set's geometries
//! becomes a new distribution of every geometry in
//! [`Ensemble::new`] with `Distribution::new(Samples::All, partition)`.
//!
//! [`Ensemble::cluster_hierarchical`] cuts the [`Dendrogram`] of a
//! [`DistanceMatrix`], such as [`Ensemble::pairwise_fitted_rmsd`], at a
//! height or into a number of clusters, and represents each cluster by its
//! weighted medoid. Complete linkage cut at height `h` bounds every
//! cluster's diameter by `h`.
//!
//! [`Ensemble::pairwise_fitted_rmsd_with`] and
//! [`Ensemble::cluster_by_radius_with`] take a
//! [`Run`](crate::parallel::Run): fits run on its executor in contiguous
//! ranges whose results combine in a fixed order, so every executor and
//! thread count returns exactly the single-threaded result, including which
//! error is reported. Too few fits for two ranges of a few milliseconds each
//! run on the calling thread.
//!
//! [`Ensemble::cluster_k_medoids`] finds a given number of clusters
//! directly, choosing medoid members that minimize the weight-summed
//! distance from every member to its nearest medoid, so heavy members
//! attract medoids.
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
//! use kekule::math::Point3;
//! use kekule::structure::{
//!     Distribution, Ensemble, Geometry, Partition, Positions, Samples, StateData, StructureSet,
//! };
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM, KELVIN, KILOJOULE_PER_MOLE};
//!
//! // Acetate, then a rotamer with its oxygens exchanged, then a distinct
//! // conformer 5 kJ/mol higher, each its own state.
//! let topology = kekule::smiles::to_topology("CC(=O)[O-]")?;
//! let none = Distribution::new(Samples::All, Partition::state_per_sample(0))?;
//! let mut conformers = Ensemble::new(StructureSet::new(topology.clone()), none)?;
//! for (o2, o3, energy) in [(1.1, -1.1, 0.0), (-1.1, 1.1, 0.0), (1.1, -0.3, 5.0)] {
//!     let points = [
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.2, o2, 0.5),
//!         Point3::new(2.2, o3, 0.5),
//!     ];
//!     let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//!     let state = StateData::new().with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))?;
//!     conformers.push_state(Geometry::new(positions), None, state, None)?;
//! }
//! conformers
//!     .distribution_mut()
//!     .assign_boltzmann_populations(Quantity::new(298.15, KELVIN))?;
//!
//! let atoms = SymmetricSelection::new(&AtomSelection::all(topology.clone()))?;
//! let radius = Quantity::new(0.125, ANGSTROM);
//! let rotamers = conformers.cluster_by_radius(radius, &atoms, Default::default())?;
//! assert_eq!(rotamers.assignments(), [0, 0, 1]);
//!
//! conformers.distribution_mut().regroup(rotamers)?;
//! assert_eq!(conformers.states().len(), 2);
//! assert!(conformers.state(0).unwrap().energy().is_some());
//! let collapsed = conformers.collapse();
//! assert_eq!(collapsed.set().len(), 2);
//! assert_eq!(
//!     collapsed.distribution().log_weight(0),
//!     Some(conformers.state(0).unwrap().log_population()),
//! );
//! # Ok::<(), Box<dyn Error>>(())
//! ```

mod driver;
mod hierarchical;
mod medoids;
mod run;

pub use driver::{ClusterOrder, ClusteringError, ClusteringOptions};
pub use hierarchical::{
    pairwise_fitted_rmsd, pairwise_fitted_rmsd_with, Cut, Dendrogram, DistanceMatrix, Linkage,
    Merge,
};

#[cfg(doc)]
use crate::structure::{Ensemble, Partition, StructureSet};
