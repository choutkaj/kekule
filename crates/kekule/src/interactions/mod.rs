//! Non-covalent interactions: hydrogen bonds, salt bridges, π stacking,
//! cation–π, halogen bonds, metal coordination, and hydrophobic contacts.
//!
//! Work splits by what it depends on. [`InteractionFeatures`] perceives once
//! per topology what can interact: donors and acceptors, charged groups,
//! aromatic rings, halogen-bond donors, metals and their ligands, and
//! hydrophobic atoms. A search then measures one structure, or two, with the
//! [`InteractionCriteria`] and reports each [`Interaction`] with its type,
//! the atoms of both ends, the points to draw it between, and its distance,
//! angles, and offset. Hydrophobic contacts are reduced as PLIP reduces
//! them, to the closest per atom and residue.
//!
//! - [`find_interactions_within`] searches one structure between two atom
//!   sets, such as a ligand and its pocket, or marked atoms and every
//!   visible atom. The sets may overlap; an interaction between two atoms
//!   of both sets is reported once.
//! - [`find_interactions_between`] searches two structures, such as one
//!   object's chain against another object's ligand: each
//!   [`InteractionSide`] brings its own features, structure, and atoms.
//!
//! Features are a separate derived result, not perception state. Both
//! searches read only the features of the selected atoms: the side with fewer
//! atoms lists its features, the other reads only atoms within reach of them,
//! far ones are discarded by bounding box, and the rest are paired on the
//! uniform point grid that group contact occupancy also uses, so a ligand against a
//! structure of 100,000 atoms takes about a millisecond. Coordinates are
//! used as stored; cells are ignored, so image periodic structures first.
//!
//! [`HydrogenBondOccupancyAccumulator`] reduces the hydrogen bonds of
//! weighted observations, such as the geometries of a set or trajectory frames, to
//! occupancies, and merges partial results as the [`crate::statistics`]
//! accumulators do.
//!
//! # Example
//!
//! ```
//! use kekule::math::Point3;
//! use kekule::interactions::{
//!     find_interactions_within, InteractionCriteria, InteractionFeatures, InteractionKind,
//! };
//! use kekule::structure::{Structure, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM};
//!
//! // Two waters 2.9 Å apart, the first donating a hydrogen to the second
//! // in heavy-atom mode (hydrogens implicit).
//! let topology = kekule::smiles::to_topology("O.O")?;
//! let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(2.9, 0.0, 0.0)];
//! let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//! let structure = Structure::new(topology.clone(), positions)?;
//!
//! let features = InteractionFeatures::new(topology.clone());
//! let first = AtomSelection::from_atoms(topology.clone(), [topology.atom_ids()[0]])?;
//! let second = AtomSelection::from_atoms(topology.clone(), [topology.atom_ids()[1]])?;
//! let found = find_interactions_within(
//!     &features,
//!     &structure,
//!     &first,
//!     &second,
//!     &InteractionCriteria::default(),
//! )?;
//! let bond = &found.interactions()[0];
//! assert_eq!(bond.kind(), InteractionKind::HydrogenBond);
//! assert!((bond.distance().value_in(ANGSTROM)? - 2.9).abs() < 1e-9);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod criteria;
mod features;
mod occupancy;
mod records;
mod search;

#[cfg(test)]
mod tests;

pub use criteria::{
    CationPiCriteria, HalogenBondCriteria, HydrogenBondCriteria, HydrophobicCriteria,
    InteractionCriteria, InteractionCriterion, MetalCoordinationCriteria, PiStackingCriteria,
    SaltBridgeCriteria,
};
pub use features::{
    AromaticRing, ChargeSign, ChargedGroup, ChargedGroupKind, HalogenBondDonor, HydrogenBondDonor,
    InteractionFeatures,
};
pub use occupancy::{
    HydrogenBondKey, HydrogenBondOccupancy, HydrogenBondOccupancyAccumulator, InteractionAtom,
};
pub use records::{
    Interaction, InteractionAngle, InteractionEnd, InteractionError, InteractionKind,
    InteractionRole, Interactions,
};
pub use search::{find_interactions_between, find_interactions_within, InteractionSide};
