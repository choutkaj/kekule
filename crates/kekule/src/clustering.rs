//! Groupings of a structure set's geometries into clusters of similar
//! structures: [`Clustering`].
//!
//! Every method clusters the members of a set: the samples of
//! [`Weights`](crate::statistics::Weights), by increasing index, or every
//! geometry without weights. A [`DistanceMatrix`] is indexed by member
//! position, the `i`th sample or geometry `i`; a [`Clustering`] is indexed
//! by geometry, and geometries that are not samples are in no cluster.
//! Errors name geometry indices. Clustering never changes the set or the
//! weights, and a set over a file clusters as one in memory does: each
//! member is lent or read once and only its fitted atoms are kept.
//!
//! A clustering's populations are the exact sums of its samples' weights
//! ([`Clustering::log_populations`]), so they follow the samples whatever
//! the method, and [`Clustering::collapse`] keeps one representative per
//! cluster, weighing the cluster's population.
//!
//! [`StructureSet::cluster_hierarchical`] cuts the [`Dendrogram`] of a
//! [`DistanceMatrix`], such as [`StructureSet::pairwise_fitted_rmsd`], at a
//! height or into a number of clusters, and represents each cluster by its
//! weighted medoid. Complete linkage cut at height `h` bounds every
//! cluster's diameter by `h`. [`Dendrogram::clustering`] cuts a kept
//! dendrogram again without the distances, each cluster represented by its
//! heaviest member.
//!
//! [`StructureSet::pairwise_fitted_rmsd_with`] and
//! [`StructureSet::cluster_by_radius_with`] take a
//! [`Run`](crate::parallel::Run): fits run on its executor in contiguous
//! ranges whose results combine in a fixed order, so every executor and
//! thread count returns exactly the single-threaded result, including which
//! error is reported. Too few fits for two ranges of a few milliseconds each
//! run on the calling thread.
//!
//! [`StructureSet::cluster_k_medoids`] finds a given number of clusters
//! directly, choosing medoid members that minimize the weight-summed
//! distance from every member to its nearest medoid, so heavy members
//! attract medoids.
//!
//! [`StructureSet::cluster_by_radius`] is leader clustering on fitted RMSD.
//! It needs no matrix, so it scales to larger sets. With a
//! [`SymmetricSelection`](crate::alignment::SymmetricSelection) and a small
//! radius it merges rotamers, copies of one conformer that differ only by
//! exchanged equivalent atoms. Weighted by energy, the merged cluster gains
//! the population of every rotamer, its degeneracy, so its free energy
//! difference to any other cluster is `-RT ln(p_a / p_b)`:
//!
//! ```
//! use std::error::Error;
//!
//! use kekule::alignment::SymmetricSelection;
//! use kekule::math::Point3;
//! use kekule::properties::{PropertyColumn, PropertyKey};
//! use kekule::statistics::Weights;
//! use kekule::structure::{Geometry, Positions, StructureSet};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM, KELVIN, KILOJOULE_PER_MOLE};
//!
//! // Acetate, then a rotamer with its oxygens exchanged, then a distinct
//! // conformer 5 kJ/mol higher.
//! let topology = kekule::smiles::to_topology("CC(=O)[O-]")?;
//! let mut conformers = StructureSet::new(topology.clone());
//! for (o2, o3) in [(1.1, -1.1), (-1.1, 1.1), (1.1, -0.3)] {
//!     let points = [
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.2, o2, 0.5),
//!         Point3::new(2.2, o3, 0.5),
//!     ];
//!     let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//!     conformers.push(Geometry::new(positions), None)?;
//! }
//! let energy = PropertyKey::new("energy")?;
//! let values = vec![Some(0.0), Some(0.0), Some(5.0)];
//! conformers.columns_mut().insert(
//!     energy.clone(),
//!     PropertyColumn::Real { unit: KILOJOULE_PER_MOLE, values },
//! )?;
//! let weights = Weights::boltzmann(&conformers, &energy, Quantity::new(298.15, KELVIN))?;
//!
//! let atoms = SymmetricSelection::new(&AtomSelection::all(topology))?;
//! let radius = Quantity::new(0.125, ANGSTROM);
//! let rotamers =
//!     conformers.cluster_by_radius(radius, &atoms, Some(&weights), Default::default())?;
//! assert_eq!(rotamers.sizes(), [2, 1]);
//!
//! // Two rotamers at 0 kJ/mol lie RT ln 2 lower than one would.
//! let rt = 8.314_462_618e-3 * 298.15;
//! let free_energies = rotamers.free_energies(&weights)?.into_unit(KILOJOULE_PER_MOLE)?;
//! assert_eq!(free_energies.value()[0], Some(0.0));
//! assert!((free_energies.value()[1].unwrap() - 5.0 - rt * 2.0_f64.ln()).abs() < 1e-9);
//!
//! let (representatives, collapsed) = rotamers.collapse(&weights)?;
//! assert_eq!(representatives, [0, 2]);
//! let populations = rotamers.log_populations(&weights)?;
//! assert_eq!(collapsed.log_weight(0), populations[0]);
//! # Ok::<(), Box<dyn Error>>(())
//! ```

mod clusters;
mod driver;
mod hierarchical;
mod medoids;
mod run;

pub use clusters::Clustering;
pub use driver::ClusteringError;
pub use hierarchical::{
    pairwise_fitted_rmsd, pairwise_fitted_rmsd_with, Cut, Dendrogram, DistanceMatrix, Linkage,
    Merge,
};

#[cfg(doc)]
use crate::structure::StructureSet;
