//! Hydrogen-bond occupancy over weighted observations.

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::topology::{InstanceAtomId, Topology, TopologyAtomIndex};

use super::criteria::HydrogenBond;
use super::{InteractionError, InteractionKind, InteractionRole, Interactions};

/// An atom of one side of a search: `first` is false for an atom of the
/// second side of [`super::find_interactions_between`]. Atoms of
/// [`super::find_interactions_within`] all belong to the first side, since
/// one topology holds both sets.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct InteractionAtom {
    second: bool,
    atom: TopologyAtomIndex,
}

impl InteractionAtom {
    /// Whether the atom belongs to the first side's topology.
    pub const fn is_first(&self) -> bool {
        !self.second
    }

    /// The atom's dense index in its side's topology.
    pub const fn atom(&self) -> TopologyAtomIndex {
        self.atom
    }
}

/// A hydrogen bond by its heavy atoms, so that rotating or exchanging
/// hydrogens, as of a hydroxyl or ammonium, keep one bond.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct HydrogenBondKey {
    donor: InteractionAtom,
    acceptor: InteractionAtom,
}

impl HydrogenBondKey {
    /// The donor atom.
    pub const fn donor(&self) -> InteractionAtom {
        self.donor
    }

    /// The acceptor atom.
    pub const fn acceptor(&self) -> InteractionAtom {
        self.acceptor
    }
}

/// What every observation of one accumulator must share.
#[derive(Debug, Clone)]
struct Scope {
    first: Arc<Topology>,
    second: Arc<Topology>,
    within: bool,
    criteria: HydrogenBond,
}

impl Scope {
    fn of(interactions: &Interactions) -> Option<Self> {
        Some(Self {
            first: Arc::clone(&interactions.first),
            second: Arc::clone(&interactions.second),
            within: interactions.within,
            criteria: interactions.hydrogen_bond.canonical().ok()?,
        })
    }

    fn matches(&self, other: &Self) -> bool {
        self.within == other.within
            && self.first.shares_layout(&other.first)
            && self.second.shares_layout(&other.second)
            && self.criteria == other.criteria
    }
}

/// Weighted hydrogen-bond occupancy: for each hydrogen bond, the fraction of
/// the total observed weight in which it was present.
///
/// An observation is the [`Interactions`] of one structure, such as an
/// ensemble sample or a trajectory frame, with an explicit, finite, positive
/// weight. Every observation must come from the same kind of search
/// ([`super::find_interactions_within`] or
/// [`super::find_interactions_between`]) over the same topology layouts with
/// the same hydrogen-bond criteria. Bonds are keyed by their donor and
/// acceptor heavy atoms ([`HydrogenBondKey`]) and counted once per
/// observation. State grows with the distinct bonds observed, not with the
/// observation count.
///
/// Accumulators merge, as [`crate::statistics`] accumulators do: a host
/// observes ranges of structures into separate accumulators, in parallel,
/// and merges them in range order; the result is the same whichever threads
/// observed the ranges.
///
/// ```
/// use kekule::interactions::{
///     find_interactions_within, HydrogenBondOccupancyAccumulator, InteractionCriteria,
///     InteractionFeatures,
/// };
/// use kekule::geometry::Point3;
/// use kekule::structure::{Model, Positions};
/// use kekule::topology::AtomSelection;
/// use kekule::units::{Quantity, ANGSTROM};
///
/// let topology = kekule::smiles::to_topology("O.O")?;
/// let features = InteractionFeatures::new(topology.clone());
/// let all = AtomSelection::all(topology.clone());
/// let criteria = InteractionCriteria::default();
/// let mut accumulator = HydrogenBondOccupancyAccumulator::new();
/// // Bonded at 2.9 Å with weight 3, apart at 5 Å with weight 1.
/// for (index, (x, weight)) in [(2.9, 3.0), (5.0, 1.0)].into_iter().enumerate() {
///     let points = [Point3::new(0.0, 0.0, 0.0), Point3::new(x, 0.0, 0.0)];
///     let model = Model::new(topology.clone(), Positions::new(Quantity::new(points, ANGSTROM))?)?;
///     let found = find_interactions_within(&features, &model, &all, &all, &criteria)?;
///     accumulator.observe(index, &found, weight)?;
/// }
/// let occupancy = accumulator.finish()?;
/// let (_, fraction) = occupancy.bonds().next().unwrap();
/// assert_eq!(fraction, 0.75);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default)]
pub struct HydrogenBondOccupancyAccumulator {
    scope: Option<Scope>,
    weights: BTreeMap<HydrogenBondKey, f64>,
    total_weight: f64,
    observations: usize,
}

impl HydrogenBondOccupancyAccumulator {
    /// An empty accumulator; its first observation sets the search every
    /// later one must match.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds the hydrogen bonds of `interactions` with `weight`; `index` is
    /// the caller's index of the observation, reported in errors. Failure
    /// changes nothing.
    ///
    /// # Errors
    ///
    /// [`InteractionError::InvalidWeight`] when `weight` is not finite and
    /// positive, [`InteractionError::WeightOverflow`] when the total weight
    /// would leave the finite range, and
    /// [`InteractionError::ObservationMismatch`] when `interactions` come
    /// from another kind of search, other topology layouts, or other
    /// hydrogen-bond criteria than earlier observations.
    pub fn observe(
        &mut self,
        index: usize,
        interactions: &Interactions,
        weight: f64,
    ) -> Result<(), InteractionError> {
        if !weight.is_finite() || weight <= 0.0 {
            return Err(InteractionError::InvalidWeight { index });
        }
        let total = self.total_weight + weight;
        if !total.is_finite() {
            return Err(InteractionError::WeightOverflow { index });
        }
        let scope =
            Scope::of(interactions).ok_or(InteractionError::ObservationMismatch { index })?;
        if self
            .scope
            .as_ref()
            .is_some_and(|known| !known.matches(&scope))
        {
            return Err(InteractionError::ObservationMismatch { index });
        }
        let mut keys = interactions
            .of_kind(InteractionKind::HydrogenBond)
            .map(|bond| {
                let atom = |end: &super::InteractionEnd, second: bool| InteractionAtom {
                    second: second && !interactions.within,
                    atom: end.atoms[0],
                };
                let (first, second) = (atom(&bond.first, false), atom(&bond.second, true));
                if bond.first.role == InteractionRole::Donor {
                    HydrogenBondKey {
                        donor: first,
                        acceptor: second,
                    }
                } else {
                    HydrogenBondKey {
                        donor: second,
                        acceptor: first,
                    }
                }
            })
            .collect::<Vec<_>>();
        keys.sort_unstable();
        keys.dedup();
        // A bond's weight sums a subsequence of the same weights as the
        // total, so it stays finite.
        for key in keys {
            *self.weights.entry(key).or_insert(0.0) += weight;
        }
        self.scope.get_or_insert(scope);
        self.total_weight = total;
        self.observations += 1;
        Ok(())
    }

    /// The number of observations so far.
    pub fn observation_count(&self) -> usize {
        self.observations
    }

    /// The total weight observed so far.
    pub fn total_weight(&self) -> f64 {
        self.total_weight
    }

    /// Merges the observations of `other` as if observed after this
    /// accumulator's own; failure changes nothing, and `other` is left as it
    /// was. Weights add, so merging ranges in order gives the same result
    /// whichever threads observed them, equal to observing everything in one
    /// accumulator up to rounding, and exactly for integer weights totalling
    /// less than 2⁵³.
    ///
    /// # Errors
    ///
    /// [`InteractionError::MergeMismatch`] when the accumulators observed
    /// different searches, and [`InteractionError::MergeWeightOverflow`] when
    /// the total weight would leave the finite range.
    pub fn merge(&mut self, other: &Self) -> Result<(), InteractionError> {
        if let (Some(mine), Some(theirs)) = (&self.scope, &other.scope) {
            if !mine.matches(theirs) {
                return Err(InteractionError::MergeMismatch);
            }
        }
        let total = self.total_weight + other.total_weight;
        if !total.is_finite() {
            return Err(InteractionError::MergeWeightOverflow);
        }
        for (key, weight) in &other.weights {
            *self.weights.entry(*key).or_insert(0.0) += weight;
        }
        if self.scope.is_none() {
            self.scope.clone_from(&other.scope);
        }
        self.total_weight = total;
        self.observations += other.observations;
        Ok(())
    }

    /// Finishes with the bonds in key order.
    ///
    /// # Errors
    ///
    /// [`InteractionError::NoObservations`] when nothing was observed.
    pub fn finish(self) -> Result<HydrogenBondOccupancy, InteractionError> {
        let scope = self.scope.ok_or(InteractionError::NoObservations)?;
        let (keys, weights) = self.weights.into_iter().unzip();
        Ok(HydrogenBondOccupancy {
            first: scope.first,
            second: scope.second,
            keys,
            weights,
            total_weight: self.total_weight,
            observations: self.observations,
        })
    }
}

/// The hydrogen bonds observed, with their weights and occupancies.
#[derive(Debug, Clone)]
pub struct HydrogenBondOccupancy {
    first: Arc<Topology>,
    second: Arc<Topology>,
    keys: Vec<HydrogenBondKey>,
    weights: Vec<f64>,
    total_weight: f64,
    observations: usize,
}

impl HydrogenBondOccupancy {
    /// The topology of first-side atoms.
    pub fn first_topology(&self) -> &Topology {
        &self.first
    }

    /// The topology of second-side atoms.
    pub fn second_topology(&self) -> &Topology {
        &self.second
    }

    /// The bonds observed at least once, in key order.
    pub fn keys(&self) -> &[HydrogenBondKey] {
        &self.keys
    }

    /// The total weight of the observations with each bond, in key order.
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    /// The number of observations.
    pub fn observation_count(&self) -> usize {
        self.observations
    }

    /// The total weight of the observations.
    pub fn total_weight(&self) -> f64 {
        self.total_weight
    }

    /// The semantic identity of a bond's atom in its side's topology, for
    /// saving results across topology publications.
    pub fn atom_id(&self, atom: InteractionAtom) -> InstanceAtomId {
        let topology = if atom.second {
            &self.second
        } else {
            &self.first
        };
        topology
            .atom_id(atom.atom)
            .expect("observed atoms belong to their side's topology")
    }

    /// Each bond with its occupancy: the fraction of the total weight in
    /// which it was present, between zero and one.
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = (HydrogenBondKey, f64)> + '_ {
        self.keys
            .iter()
            .copied()
            .zip(self.weights.iter().map(|weight| weight / self.total_weight))
    }
}
