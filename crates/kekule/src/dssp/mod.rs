//! DSSP 4-compatible protein secondary-structure assignment.
//!
//! [`assign`] is a read-only analysis of one
//! [`Structure`](crate::structure::Structure) coordinate snapshot that returns a
//! [`DsspResult`]: per-residue [`DsspResidue`] values (summary
//! [`DsspSecondaryStructure`], backbone torsions, hydrogen-bond and beta
//! partners, helix flags), aggregate [`DsspStatistics`], and a [`DsspReport`]
//! of skipped residues and work counts. [`DsspOptions`] and [`DsspLimits`]
//! configure the run; failures are [`DsspError`].
//!
//! [`DsspOptions::with_threads`] spreads one assignment over several scoped
//! threads, and [`assign_with`] runs it on a host's own thread pool through a
//! [`parallel::Executor`](crate::parallel::Executor). Every executor and thread count returns the same result
//! as the default single thread, including the order of tied hydrogen bonds
//! and every identifier.
//!
//! A [`DsspPlan`] prepares what assignment reads from a topology once, so the
//! frames of a trajectory, or the geometries of an ensemble, are each assigned
//! from their backbone coordinates alone, with the same results.
//!
//! Assignment does not parse structures, choose coordinate models or alternate
//! locations, repair hierarchy data, or install derived state in the structure.
//! Geometry uses single-precision coordinates and arithmetic, matching DSSP 4;
//! structure coordinates are unchanged and result values are exposed as `f64`.
//!
//! # Example
//!
//! ```
//! use kekule::dssp::{self, DsspChainBreak, DsspOptions, DsspSecondaryStructure};
//! use kekule::pdb;
//!
//! // A glycine dipeptide: too short for any secondary structure.
//! let entry = concat!(
//!     "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n",
//!     "ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C\n",
//!     "ATOM      3  C   GLY A   1       2.009   1.420   0.000  1.00  0.00           C\n",
//!     "ATOM      4  O   GLY A   1       1.251   2.390   0.000  1.00  0.00           O\n",
//!     "ATOM      5  N   GLY A   2       3.332   1.536   0.000  1.00  0.00           N\n",
//!     "ATOM      6  CA  GLY A   2       3.988   2.836   0.000  1.00  0.00           C\n",
//!     "ATOM      7  C   GLY A   2       5.488   2.736   0.000  1.00  0.00           C\n",
//!     "ATOM      8  O   GLY A   2       6.100   1.680   0.000  1.00  0.00           O\n",
//!     "END\n",
//! );
//! let structure = pdb::parse_str(entry, Default::default())?
//!     .interpret(Default::default())?
//!     .into_structure();
//!
//! let result = dssp::assign(&structure, DsspOptions::default())?;
//! let residues = result.residues().collect::<Vec<_>>();
//! assert_eq!(residues.len(), 2);
//! assert_eq!(residues[0].chain_break(), DsspChainBreak::NewChain);
//! assert_eq!(residues[1].chain_break(), DsspChainBreak::None);
//! assert_eq!(residues[0].secondary_structure(), DsspSecondaryStructure::Loop);
//! assert_eq!(residues[0].source().author_sequence_id(), Some("1"));
//! assert!(residues[0].phi().is_none());
//! assert_eq!(result.statistics().chain_segments(), 1);
//! assert_eq!(
//!     result.statistics().secondary_structure_count(DsspSecondaryStructure::Loop),
//!     2
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod assignment;
mod kernel;

pub use assignment::{
    assign, assign_with, DsspBackboneAtom, DsspBetaPartner, DsspBridgeDirection, DsspChainBreak,
    DsspCodeError, DsspError, DsspHelixPosition, DsspHydrogenBond, DsspLimits, DsspOptions,
    DsspPlan, DsspReport, DsspResidue, DsspResidueSource, DsspResult, DsspSecondaryStructure,
    DsspSkipReason, DsspSkippedResidue, DsspStatistics,
};
