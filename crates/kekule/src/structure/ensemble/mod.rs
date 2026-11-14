//! A structure set read as draws from a distribution: [`Ensemble`].

use crate::topology::AtomSelection;

use super::{
    Geometry, GeometryMut, GeometryStorage, InMemory, Place, StructureSet, StructureSetError,
};

mod distribution;
mod partition;
mod states;
#[cfg(test)]
mod tests;

pub(crate) use distribution::log_sum_exp;
pub use distribution::{Distribution, DistributionMut, Samples, StateData};
pub use partition::{Partition, PartitionError};
pub use states::{StateIter, StateView};

/// A [`StructureSet`] read as draws from a [`Distribution`]: which of its
/// geometries are samples, what each weighs, and how they group into states.
///
/// A distribution is inferred, not recorded: reweighting a biased run,
/// clustering into states, and Boltzmann populations from energies are
/// interpretations, and one set often carries several. So the distribution
/// belongs to the ensemble and never to the set, and several readings of one
/// set are several ensembles sharing every coordinate.
///
/// The ensemble lends its set ([`Self::set`]) but never lends it mutably:
/// everything that changes which geometries exist goes through the ensemble,
/// which remaps the samples, so the distribution always fits its set.
/// [`Self::into_parts`] returns both without copying, and [`Self::new`]
/// pairs them again.
///
/// Statistics and clustering over an ensemble, such as [`Self::rmsf`] and
/// [`Self::cluster_by_radius`], observe its samples with their
/// probabilities and leave the other geometries out.
///
/// ```
/// use kekule::structure::{Ensemble, Order, Place, Positions, Geometry, StructureSet, TimeAxis};
/// use kekule::units::{Quantity, PICOSECOND};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let mut set = StructureSet::new(topology.clone());
/// set.set_order(Some(Order::Time(TimeAxis::new(
///     Some(Quantity::new(Vec::new(), PICOSECOND)),
///     None,
/// )?)))?;
/// for picoseconds in [0.0, 2.0, 4.0] {
///     let place = Place::Time { time: Some(Quantity::new(picoseconds, PICOSECOND)), step: None };
///     set.push(Geometry::new(Positions::zeros(topology.atom_count())), Some(place))?;
/// }
/// // Reweighted frames: the last weighs e¹⁰⁰ times the first, beyond what
/// // linear weights near one could tell apart from zero.
/// let ensemble = Ensemble::weighted(set, vec![0.0, 50.0, 100.0])?;
/// let probabilities = ensemble.distribution().probabilities()?;
/// assert!((probabilities[2] - 1.0).abs() < 1e-15);
/// assert_eq!(ensemble.distribution().log_weight(0), Some(0.0));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct Ensemble<S = InMemory> {
    set: StructureSet<S>,
    distribution: Distribution,
}

impl<S: GeometryStorage> Ensemble<S> {
    /// `distribution` over the geometries of `set`. A subset listing every
    /// geometry becomes [`Samples::All`].
    ///
    /// # Errors
    ///
    /// [`StructureSetError::IndexOutOfRange`] for a sample outside the set,
    /// and [`StructureSetError::CoverageMismatch`] when a distribution of
    /// every geometry covers another number of geometries.
    pub fn new(
        set: StructureSet<S>,
        mut distribution: Distribution,
    ) -> Result<Self, StructureSetError> {
        distribution.bind(set.len())?;
        Ok(Self { set, distribution })
    }

    /// Every geometry of `set` as an equally weighted sample of one state,
    /// represented by the first: the frames of an equilibrium trajectory or
    /// the models of an NMR bundle. An empty set has no states.
    pub fn uniform(set: StructureSet<S>) -> Self {
        let distribution = Distribution::new(Samples::All, Partition::one_state(set.len()))
            .expect("a distribution of every geometry has no subset to check");
        Self { set, distribution }
    }

    /// Every geometry of `set` as a sample of one state, represented by the
    /// first, with one natural-logarithm weight per geometry, in geometry
    /// order: the frames of a reweighted biased run.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::LogWeightCountMismatch`] for another number of
    /// log weights than geometries, and
    /// [`StructureSetError::NonFiniteLogWeight`] for the first one that is not
    /// finite.
    pub fn weighted(
        set: StructureSet<S>,
        log_weights: Vec<f64>,
    ) -> Result<Self, StructureSetError> {
        let mut ensemble = Self::uniform(set);
        ensemble
            .distribution_mut()
            .set_log_weights(Some(log_weights))?;
        Ok(ensemble)
    }

    /// The structure set whose geometries the distribution reads.
    pub fn set(&self) -> &StructureSet<S> {
        &self.set
    }

    /// The distribution over the set's geometries.
    pub fn distribution(&self) -> &Distribution {
        &self.distribution
    }

    /// The distribution for reweighting, regrouping, and editing states,
    /// through a guard that keeps its samples ([`DistributionMut`]).
    pub fn distribution_mut(&mut self) -> DistributionMut<'_> {
        DistributionMut::new(&mut self.distribution)
    }

    /// The set and the distribution, without copying.
    pub fn into_parts(self) -> (StructureSet<S>, Distribution) {
        (self.set, self.distribution)
    }

    /// One state by index, with its representative, samples, and
    /// population, or `None` out of range. The state of an ensemble in
    /// memory lends its representative geometry
    /// ([`StateView::representative_sample`]).
    pub fn state(&self, index: usize) -> Option<StateView<'_, S>> {
        (index < self.distribution.state_count()).then(|| StateView::new(self, index))
    }

    /// The states in order.
    pub fn states(&self) -> StateIter<'_, S> {
        StateIter::new(self)
    }

    /// `distribution`, already fitted to `set`, over it.
    pub(crate) fn from_parts_unchecked(set: StructureSet<S>, distribution: Distribution) -> Self {
        debug_assert!(distribution.clone().bind(set.len()).is_ok());
        Self { set, distribution }
    }

    /// The set, for changing what its storage holds; the caller then
    /// re-expresses the distribution.
    pub(crate) fn parts_mut(&mut self) -> (&mut StructureSet<S>, &mut Distribution) {
        (&mut self.set, &mut self.distribution)
    }
}

impl Ensemble {
    /// A dimension-preserving editor for one geometry, as
    /// [`StructureSet::get_mut`] gives it; the geometry stays a sample or
    /// not.
    pub fn get_mut(&mut self, index: usize) -> Option<GeometryMut<'_>> {
        self.set.get_mut(index)
    }

    /// Validates and replaces one geometry, as [`StructureSet::replace`]
    /// does, returning the previous one. The geometry keeps its place and
    /// its sample. Failure changes nothing.
    pub fn replace(
        &mut self,
        index: usize,
        geometry: Geometry,
    ) -> Result<Geometry, StructureSetError> {
        self.set.replace(index, geometry)
    }

    /// Validates and appends one geometry at `place`, as
    /// [`StructureSet::push`] does, outside the distribution, returning its
    /// index. A distribution of every geometry becomes one of the geometries
    /// before this one. Failure changes nothing.
    pub fn push(
        &mut self,
        geometry: Geometry,
        place: Option<Place>,
    ) -> Result<usize, StructureSetError> {
        let index = self.set.push(geometry, place)?;
        self.distribution
            .extend(index + 1)
            .expect("the distribution describes the geometries before the new one");
        Ok(index)
    }

    /// Validates and appends one geometry at `place`, as
    /// [`StructureSet::push`] does, as a sample of an existing state,
    /// returning its index.
    ///
    /// `log_weight` is the new sample's finite natural-logarithm weight. It
    /// is required when the distribution has weights; in one without, it
    /// gives the other samples log weight `0.0`. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::StateOutOfRange`],
    /// [`StructureSetError::LogWeightRequired`],
    /// [`StructureSetError::NonFiniteLogWeight`], and the errors of
    /// [`StructureSet::push`].
    pub fn push_sample(
        &mut self,
        geometry: Geometry,
        place: Option<Place>,
        state: usize,
        log_weight: Option<f64>,
    ) -> Result<usize, StructureSetError> {
        let index = self.set.len();
        self.distribution.check_state(state)?;
        self.distribution.check_new_log_weight(log_weight, index)?;
        self.set.push(geometry, place)?;
        self.distribution.push_sample(index, state, log_weight);
        Ok(index)
    }

    /// Validates and appends one geometry at `place`, as
    /// [`StructureSet::push`] does, as the representative of a new state,
    /// returning the state's index. `log_weight` is as in
    /// [`Self::push_sample`]. Failure changes nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::LogWeightRequired`],
    /// [`StructureSetError::NonFiniteLogWeight`], and the errors of
    /// [`StructureSet::push`].
    pub fn push_state(
        &mut self,
        geometry: Geometry,
        place: Option<Place>,
        state: StateData,
        log_weight: Option<f64>,
    ) -> Result<usize, StructureSetError> {
        let index = self.set.len();
        self.distribution.check_new_log_weight(log_weight, index)?;
        self.set.push(geometry, place)?;
        Ok(self.distribution.push_state(index, state, log_weight))
    }

    /// Removes and returns one geometry, with its place and its sample,
    /// shifting later geometries and samples down.
    ///
    /// Removing the only sample of a state removes the state, shifting later
    /// states down. The representative of a state with other samples cannot
    /// be removed ([`StructureSetError::RepresentativeRemoval`]); choose
    /// another representative with [`DistributionMut::regroup`] first.
    /// Failure changes nothing.
    pub fn remove(&mut self, index: usize) -> Result<Geometry, StructureSetError> {
        if index < self.set.len() {
            self.distribution.check_remove(index)?;
        }
        let geometry = self.set.remove(index)?;
        self.distribution.remove(index, self.set.len());
        Ok(geometry)
    }

    /// Copies geometries in the requested order, as [`StructureSet::select`]
    /// does, with the distribution restricted to them.
    ///
    /// The selected samples stay samples with their log weights, and the
    /// temperature is kept. Every state with a selected sample becomes a
    /// state of the copy, in state order, represented by the first copy of
    /// its representative or, when that is not selected, by its first
    /// selected sample. A state whose samples are all selected keeps its
    /// energy and properties; one selected only in part has neither, as a
    /// split state in [`DistributionMut::regroup`]. Use
    /// [`Self::select_states`] to keep whole states.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::IndexOutOfRange`] for the first index outside
    /// the set.
    pub fn select(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, StructureSetError> {
        let indices = indices.into_iter().collect::<Vec<_>>();
        let set = self.set.select(indices.iter().copied())?;
        Ok(Self {
            set,
            distribution: self.distribution.select(&indices),
        })
    }

    /// Copies the samples of whole states in the requested order, each state
    /// with all of its samples in geometry order (sharing their geometries),
    /// its representative, energy, and properties. Duplicates and reordering
    /// are allowed. Log weights, the temperature, the geometries' order
    /// values, and set properties are kept; geometries that are not samples
    /// are left out. [`Distribution::covering_population`],
    /// [`Distribution::within_energy_window`], and
    /// [`Distribution::ranked_by_population`] give states to extract.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::StateOutOfRange`] for a state outside the
    /// distribution.
    pub fn select_states(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, StructureSetError> {
        let indices = indices
            .into_iter()
            .map(|index| self.distribution.check_state(index).map(|()| index))
            .collect::<Result<Vec<_>, _>>()?;
        let (geometries, distribution) = self.distribution.select_states(&indices);
        Ok(Self {
            set: self.set.select(geometries)?,
            distribution,
        })
    }

    /// A new ensemble of one geometry per state, in state order: each
    /// state's representative, sharing its geometry, as the only sample of
    /// its state, with its state's population as its weight, so population
    /// shares are kept.
    ///
    /// It coarse-grains, losing the spread within each state, and drops the
    /// order. Each state keeps its energy and properties; the topology, the
    /// set properties, and the temperature are kept, and the source is
    /// unchanged.
    pub fn collapse(&self) -> Self {
        let (geometries, distribution) = self.distribution.collapse();
        let mut set = self
            .set
            .select(geometries)
            .expect("representatives are geometries of the set");
        set.set_order(None)
            .expect("an unordered set takes no order values");
        Self { set, distribution }
    }

    /// Projects every geometry onto one induced topology subset, as
    /// [`StructureSet::subset`] does. Set, state, and geometry properties
    /// are dropped; the samples, log weights, states, energies, and
    /// temperature still describe the geometries and are kept.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, StructureSetError> {
        Ok(Self {
            set: self.set.subset(selection)?,
            distribution: self.distribution.without_state_properties(),
        })
    }
}
