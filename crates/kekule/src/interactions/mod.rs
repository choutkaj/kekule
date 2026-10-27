//! Non-covalent interactions: hydrogen bonds, salt bridges, π stacking,
//! cation–π, halogen bonds, metal coordination, and hydrophobic contacts.
//!
//! Work splits by what it depends on. [`InteractionFeatures`] perceives once
//! per topology what can interact: donors and acceptors, charged groups,
//! aromatic rings, halogen-bond donors, metals and their ligands, and
//! hydrophobic atoms. A search then measures one structure, or two, with the
//! [`InteractionCriteria`] and reports each [`Interaction`] with its type,
//! the atoms of both ends, the points to draw it between, and its distance,
//! angles, and offset.
//!
//! - [`find_interactions_within`] searches one structure between two atom
//!   sets, such as a ligand and its pocket, or marked atoms and every
//!   visible atom. The sets may overlap; an interaction between two atoms
//!   of both sets is reported once.
//! - [`find_interactions_between`] searches two structures, such as one
//!   object's chain against another object's ligand: each
//!   [`InteractionSide`] brings its own features, structure, and atoms.
//!
//! Both read only the features of the selected atoms, discard the far ones by
//! bounding box, and pair the rest on a uniform grid, so a ligand against a
//! structure of 100,000 atoms takes about a millisecond. Coordinates are
//! used as stored; cells are ignored, so image periodic structures first.
//!
//! [`HydrogenBondOccupancyAccumulator`] reduces the hydrogen bonds of
//! weighted observations, such as ensemble samples or trajectory frames, to
//! occupancies, and merges partial results as the [`crate::statistics`]
//! accumulators do.
//!
//! # Example
//!
//! ```
//! use kekule::geometry::Point3;
//! use kekule::interactions::{
//!     find_interactions_within, InteractionCriteria, InteractionFeatures, InteractionKind,
//! };
//! use kekule::structure::{Model, Positions};
//! use kekule::topology::AtomSelection;
//! use kekule::units::{Quantity, ANGSTROM};
//!
//! // Two waters 2.9 Å apart, the first donating a hydrogen to the second
//! // in heavy-atom mode (hydrogens implicit).
//! let topology = kekule::smiles::to_topology("O.O")?;
//! let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(2.9, 0.0, 0.0)];
//! let positions = Positions::new(Quantity::new(points, ANGSTROM))?;
//! let model = Model::new(topology.clone(), positions)?;
//!
//! let features = InteractionFeatures::new(topology.clone());
//! let first = AtomSelection::from_atoms(topology.clone(), [topology.atom_ids()[0]])?;
//! let second = AtomSelection::from_atoms(topology.clone(), [topology.atom_ids()[1]])?;
//! let found = find_interactions_within(
//!     &features,
//!     &model,
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
mod search;

use std::fmt;
use std::sync::Arc;

use crate::geometry::{Point3, Vector3};
use crate::topology::{SelectionError, Topology, TopologyAtomIndex};
use crate::units::{Quantity, CANONICAL_ANGLE_UNIT, CANONICAL_LENGTH_UNIT};

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
pub use search::{find_interactions_between, find_interactions_within, InteractionSide};

/// Invalid interaction inputs or a failed occupancy observation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum InteractionError {
    /// A structure or selection does not share its features' topology
    /// layout.
    TopologyMismatch,
    /// An atom selection is invalid; see the source.
    Selection(SelectionError),
    /// A criterion is not a finite, nonnegative length, or an angle between
    /// 0 and 180°.
    InvalidCriterion(InteractionCriterion),
    /// An aromatic ring's atoms have no plane.
    DegenerateRing,
    /// An observation weight is not finite and positive.
    InvalidWeight {
        /// The caller's index of the observation.
        index: usize,
    },
    /// The total observed weight left the finite numerical range.
    WeightOverflow {
        /// The caller's index of the observation.
        index: usize,
    },
    /// Observed interactions come from other topology layouts, another kind
    /// of search, or other hydrogen-bond criteria than earlier ones.
    ObservationMismatch {
        /// The caller's index of the observation.
        index: usize,
    },
    /// Two accumulators observed other topology layouts, another kind of
    /// search, or other hydrogen-bond criteria, and cannot be merged.
    MergeMismatch,
    /// The merged total weight would leave the finite numerical range.
    MergeWeightOverflow,
    /// The occupancy was finished before any observation.
    NoObservations,
}

impl fmt::Display for InteractionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                f.write_str("structure or selection does not share the features' topology layout")
            }
            Self::Selection(_) => f.write_str("invalid atom selection"),
            Self::InvalidCriterion(criterion) => write!(f, "invalid {criterion}"),
            Self::DegenerateRing => f.write_str("aromatic ring has no plane"),
            Self::InvalidWeight { index } => {
                write!(f, "observation {index} weight is not finite and positive")
            }
            Self::WeightOverflow { index } => {
                write!(f, "total weight overflowed at observation {index}")
            }
            Self::ObservationMismatch { index } => write!(
                f,
                "observation {index} comes from another search than earlier observations"
            ),
            Self::MergeMismatch => f.write_str("accumulators observed different searches"),
            Self::MergeWeightOverflow => f.write_str("merged total weight overflowed"),
            Self::NoObservations => f.write_str("no observations"),
        }
    }
}

impl std::error::Error for InteractionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Selection(error) => Some(error),
            _ => None,
        }
    }
}

impl From<SelectionError> for InteractionError {
    fn from(error: SelectionError) -> Self {
        Self::Selection(error)
    }
}

/// The type of an interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum InteractionKind {
    /// A hydrogen bond from a donor to an acceptor.
    HydrogenBond,
    /// A salt bridge between an anion and a cation.
    SaltBridge,
    /// Parallel π stacking of two aromatic rings.
    PiStackingParallel,
    /// T-shaped π stacking of two aromatic rings.
    PiStackingTShaped,
    /// A cation over an aromatic ring.
    CationPi,
    /// A halogen bond from a carbon-bound halogen to an acceptor.
    HalogenBond,
    /// A metal and an atom coordinating it.
    MetalCoordination,
    /// Two hydrophobic atoms in contact.
    HydrophobicContact,
}

impl InteractionKind {
    /// Every kind, in declaration order.
    pub const ALL: [Self; 8] = [
        Self::HydrogenBond,
        Self::SaltBridge,
        Self::PiStackingParallel,
        Self::PiStackingTShaped,
        Self::CationPi,
        Self::HalogenBond,
        Self::MetalCoordination,
        Self::HydrophobicContact,
    ];
}

impl fmt::Display for InteractionKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::HydrogenBond => "hydrogen bond",
            Self::SaltBridge => "salt bridge",
            Self::PiStackingParallel => "parallel pi stacking",
            Self::PiStackingTShaped => "T-shaped pi stacking",
            Self::CationPi => "cation-pi",
            Self::HalogenBond => "halogen bond",
            Self::MetalCoordination => "metal coordination",
            Self::HydrophobicContact => "hydrophobic contact",
        })
    }
}

/// The part one end plays in an interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum InteractionRole {
    /// A hydrogen-bond donor.
    Donor,
    /// A hydrogen-bond or halogen-bond acceptor.
    Acceptor,
    /// A negatively charged group.
    Anion,
    /// A positively charged group.
    Cation,
    /// An aromatic ring.
    Ring,
    /// A halogen-bond donor's halogen.
    HalogenDonor,
    /// A metal.
    Metal,
    /// An atom coordinating a metal.
    MetalLigand,
    /// A hydrophobic atom.
    Hydrophobic,
}

/// The angle an [`Interaction`] reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum InteractionAngle {
    /// D–H…A of a hydrogen bond with an explicit hydrogen.
    DonorHydrogenAcceptor,
    /// X–D…A at one heavy neighbor X of a hydrogen-bond donor, without
    /// explicit hydrogens.
    NeighborDonorAcceptor,
    /// D…A–Y at one heavy neighbor Y of a hydrogen-bond acceptor.
    DonorAcceptorNeighbor,
    /// C–X…A of a halogen bond.
    CarbonHalogenAcceptor,
    /// X…A–Y at one heavy neighbor Y of a halogen-bond acceptor.
    HalogenAcceptorNeighbor,
    /// The angle between two rings' planes, from 0° to 90°.
    RingPlanes,
    /// The angle between a ring's normal and the line from its center to a
    /// cation's charge center, from 0° to 90°.
    RingNormalCation,
}

/// One end of an [`Interaction`].
#[derive(Debug, Clone, PartialEq)]
pub struct InteractionEnd {
    role: InteractionRole,
    atoms: Vec<TopologyAtomIndex>,
    hydrogen: Option<TopologyAtomIndex>,
    center: Point3,
    normal: Option<Vector3>,
}

impl InteractionEnd {
    /// The part this end plays.
    pub const fn role(&self) -> InteractionRole {
        self.role
    }

    /// The end's atoms in its own side's topology: one atom, a ring's atoms
    /// in order around it, or a charged group's atoms, charged atoms first.
    pub fn atoms(&self) -> &[TopologyAtomIndex] {
        &self.atoms
    }

    /// The explicit hydrogen of a hydrogen-bond donor, when it has one.
    pub const fn hydrogen(&self) -> Option<TopologyAtomIndex> {
        self.hydrogen
    }

    /// The point to draw the interaction to: the atom's position, a ring's
    /// center, or a charged group's charge center, in the canonical length
    /// unit (nm).
    pub fn center(&self) -> Quantity<Point3> {
        Quantity::new(self.center, CANONICAL_LENGTH_UNIT)
    }

    /// A ring's unit normal; `None` for other ends.
    pub const fn normal(&self) -> Option<Vector3> {
        self.normal
    }
}

/// One non-covalent interaction between an end of the first side, or set,
/// and an end of the second.
#[derive(Debug, Clone, PartialEq)]
pub struct Interaction {
    kind: InteractionKind,
    first: InteractionEnd,
    second: InteractionEnd,
    distance: f64,
    hydrogen_distance: Option<f64>,
    offset: Option<f64>,
    angles: Vec<(InteractionAngle, f64)>,
}

impl Interaction {
    /// The interaction's type.
    pub const fn kind(&self) -> InteractionKind {
        self.kind
    }

    /// The end from the first side or set.
    pub const fn first(&self) -> &InteractionEnd {
        &self.first
    }

    /// The end from the second side or set.
    pub const fn second(&self) -> &InteractionEnd {
        &self.second
    }

    /// The distance between the ends' centers in the canonical length unit
    /// (nm): donor to acceptor, halogen to acceptor, metal to ligand, atom to
    /// atom, or between charge and ring centers.
    pub fn distance(&self) -> Quantity<f64> {
        Quantity::new(self.distance, CANONICAL_LENGTH_UNIT)
    }

    /// The hydrogen–acceptor distance of a hydrogen bond with an explicit
    /// hydrogen.
    pub fn hydrogen_distance(&self) -> Option<Quantity<f64>> {
        self.hydrogen_distance
            .map(|distance| Quantity::new(distance, CANONICAL_LENGTH_UNIT))
    }

    /// The offset of π stacking or cation–π; see [`PiStackingCriteria`] and
    /// [`CationPiCriteria`].
    pub fn offset(&self) -> Option<Quantity<f64>> {
        self.offset
            .map(|offset| Quantity::new(offset, CANONICAL_LENGTH_UNIT))
    }

    /// The angles the interaction's criteria test, one per neighbor where
    /// they concern neighbors, in the canonical angle unit (rad).
    pub fn angles(&self) -> impl ExactSizeIterator<Item = (InteractionAngle, Quantity<f64>)> + '_ {
        self.angles
            .iter()
            .map(|&(angle, value)| (angle, Quantity::new(value, CANONICAL_ANGLE_UNIT)))
    }
}

/// The interactions one search found, sorted by kind, then by the first and
/// second ends' atoms.
#[derive(Debug, Clone)]
pub struct Interactions {
    first: Arc<Topology>,
    second: Arc<Topology>,
    within: bool,
    hydrogen_bond: HydrogenBondCriteria,
    min_bond_separation: usize,
    interactions: Vec<Interaction>,
}

impl Interactions {
    /// The interactions found.
    pub fn interactions(&self) -> &[Interaction] {
        &self.interactions
    }

    /// The interactions of one kind.
    pub fn of_kind(&self, kind: InteractionKind) -> impl Iterator<Item = &Interaction> + '_ {
        self.interactions
            .iter()
            .filter(move |interaction| interaction.kind == kind)
    }

    /// Whether one structure was searched ([`find_interactions_within`]).
    pub const fn is_within(&self) -> bool {
        self.within
    }

    /// The topology of the first side's atom indices.
    pub fn first_topology(&self) -> &Topology {
        &self.first
    }

    /// The topology of the second side's atom indices; the first's when one
    /// structure was searched.
    pub fn second_topology(&self) -> &Topology {
        &self.second
    }
}

#[cfg(test)]
mod tests;
