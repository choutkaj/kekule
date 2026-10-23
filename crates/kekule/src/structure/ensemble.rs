use crate::diagnostics::Cancellation;
use crate::structure::DenseArrayError;
use std::fmt;
use std::sync::Arc;

use crate::geometry::{PeriodicCell, RigidTransform};
use crate::molecule::Molecule;
use crate::properties::{OwnerProperties, StructureProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{
    AtomSelection, Topology, TopologyBuildError, TopologyBuilder, TopologyPerceptionError,
};
use crate::units::{Quantity, UnitError, CANONICAL_ENERGY_UNIT, CANONICAL_TEMPERATURE_UNIT};

use super::partition::EnsemblePartition;
use super::scalars::AtomScalars;
use super::state::{
    is_bound, site_mutators, site_readers, structure_mutators, structure_readers, subset_rows,
    transform_cell, transform_positions, validate_binding,
};
use super::{AsModelView, Model, ModelView, Positions, StructureError, StructureKind};

/// A molecule as the distribution of its conformations: weighted samples of
/// one shared topology, grouped into states.
///
/// The ensemble has two levels.
///
/// - An [`EnsembleSample`] is one structure with a finite, positive,
///   relative *weight*: the share of probability it carries. Only ratios of
///   weights are meaningful, so selecting, duplicating, or removing samples
///   needs no renormalization.
/// - An [`EnsembleState`] is a distinct conformation: a group of samples, one
///   of them its representative. Its *population* is the total weight of its
///   samples, derived and never stored. A state may carry an *energy*, on one
///   scale for the whole ensemble, such as a conformer's minimized energy or a
///   computed free energy, never a snapshot's instantaneous energy. Samples
///   carry none.
///
/// The ensemble owns its grouping as an [`EnsemblePartition`]: every sample
/// is in exactly one state, and every state has at least one sample, its
/// representative among them. It may record the *temperature* its
/// populations describe.
///
/// A conformer set is one state per conformer; an NMR bundle or the frames
/// of a trajectory are one state of many samples, and clustering such an
/// ensemble ([`Self::regroup`] with a partition from
/// [`crate::clustering`]) makes several. Energies never change weights
/// implicitly; [`Self::assign_boltzmann_populations`] derives state
/// populations from them explicitly, scaling the weights of each state's
/// samples together.
///
/// Samples are stored flat, in a stable order: regrouping, weighting, and
/// editing states never reorder them, so sample indices held elsewhere stay
/// valid. Sample order carries no physical meaning; use
/// [`super::Trajectory`] when order is temporal. Statistics over an ensemble
/// must use its weights.
///
/// Samples are validated against the shared topology on insertion, so they
/// can be read as [`super::ModelView`]s without further checks. Superposition
/// and RMSD are inherent methods; see [`crate::alignment`].
///
/// ```
/// use kekule::structure::{Ensemble, EnsemblePartition, EnsembleSample, Positions};
/// let topology = kekule::smiles::to_topology("O")?;
/// let atoms = topology.atom_count();
/// let samples = [
///     EnsembleSample::new(Positions::zeros(atoms), 3.0)?,
///     EnsembleSample::new(Positions::zeros(atoms), 1.0)?,
/// ];
/// let mut ensemble =
///     Ensemble::from_samples(topology, samples, EnsemblePartition::state_per_sample(2))?;
/// ensemble.normalize_weights()?;
/// assert_eq!(ensemble.get(0).unwrap().weight(), 0.75);
/// assert_eq!(ensemble.state(1).unwrap().population(), 0.25);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// Sample editors preserve dimensions and cannot overwrite a whole sample:
///
/// ```compile_fail,E0594
/// use kekule::structure::{Ensemble, EnsembleSample, Positions};
/// fn overwrite(ensemble: &mut Ensemble) {
///     *ensemble.get_mut(0).unwrap() = EnsembleSample::new(Positions::zeros(1), 1.0).unwrap();
/// }
/// ```
#[derive(Debug, Clone)]
pub struct Ensemble {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    samples: Vec<EnsembleSample>,
    partition: EnsemblePartition,
    states: Vec<EnsembleState>,
    temperature: Option<Quantity<f64>>,
}

impl Ensemble {
    /// An empty ensemble of one shared topology, with no samples or states.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            samples: Vec::new(),
            partition: EnsemblePartition::state_per_sample(0),
            states: Vec::new(),
            temperature: None,
        }
    }

    /// Validates and collects each sample in order, binding its bond
    /// property rows, and groups them by `partition` into states without
    /// energies or properties. Fails on the first invalid sample, or with
    /// [`EnsembleError::PartitionMismatch`] when `partition` covers a
    /// different number of samples.
    pub fn from_samples(
        topology: impl Into<Arc<Topology>>,
        samples: impl IntoIterator<Item = EnsembleSample>,
        partition: EnsemblePartition,
    ) -> Result<Self, EnsembleError> {
        let mut ensemble = Self::new(topology);
        let mut samples = samples.into_iter().collect::<Vec<_>>();
        for sample in &mut samples {
            ensemble.bind(sample)?;
        }
        ensemble.install(samples, partition)?;
        Ok(ensemble)
    }

    /// Collects models that share one topology layout as equally weighted
    /// samples (weight `1.0`), grouped by `partition`, keeping the first
    /// model's snapshot. Positions, cells, occupancies, B-factors, and
    /// properties move without copying.
    ///
    /// The partition says what the models are:
    /// [`EnsemblePartition::one_state`] for draws of one distribution, such
    /// as the models of an NMR bundle, and
    /// [`EnsemblePartition::state_per_sample`] for distinct conformers.
    pub fn from_models(
        models: impl IntoIterator<Item = Model>,
        partition: EnsemblePartition,
    ) -> Result<Self, EnsembleError> {
        let mut models = models.into_iter();
        let first = models.next().ok_or(EnsembleError::EmptySource)?;
        let mut ensemble = Self::new(first.shared_topology());
        let mut samples = vec![EnsembleSample::from_model(first)];
        for model in models {
            if !model.topology().shares_layout(&ensemble.topology) {
                return Err(EnsembleError::TopologyMismatch);
            }
            samples.push(EnsembleSample::from_model(model));
        }
        ensemble.install(samples, partition)?;
        Ok(ensemble)
    }

    /// Builds equally weighted samples of one molecule from dense positions
    /// in molecule atom order, grouped by `partition`.
    pub fn from_molecule_positions(
        molecule: Molecule,
        positions: impl IntoIterator<Item = Positions>,
        partition: EnsemblePartition,
    ) -> Result<Self, EnsembleError> {
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(molecule)?;
        builder.add_instance(definition)?;
        Self::from_samples(
            builder.build()?,
            positions.into_iter().map(EnsembleSample::unit),
            partition,
        )
    }

    /// Samples already bound to `topology`, forming one state, for
    /// projections from another collection of the same topology.
    pub(crate) fn from_bound_samples(
        topology: Arc<Topology>,
        properties: OwnerProperties,
        samples: Vec<EnsembleSample>,
    ) -> Self {
        debug_assert!(samples.iter().all(|sample| is_bound(
            sample.atom_count(),
            &sample.properties,
            &topology
        )));
        let partition = EnsemblePartition::one_state(samples.len());
        Self {
            topology,
            properties,
            states: vec![EnsembleState::new(); partition.state_count()],
            samples,
            partition,
            temperature: None,
        }
    }

    /// Installs bound samples grouped by `partition`, with empty states.
    fn install(
        &mut self,
        samples: Vec<EnsembleSample>,
        partition: EnsemblePartition,
    ) -> Result<(), EnsembleError> {
        if partition.sample_count() != samples.len() {
            return Err(EnsembleError::PartitionMismatch {
                expected: samples.len(),
                actual: partition.sample_count(),
            });
        }
        self.states = vec![EnsembleState::new(); partition.state_count()];
        self.samples = samples;
        self.partition = partition;
        Ok(())
    }

    /// The topology shared by every sample.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to the topology; cloning the [`Arc`] copies no topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Ensemble-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable ensemble-level annotations. Sample properties are edited
    /// through [`Self::get_mut`], state properties through
    /// [`Self::state_mut`].
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    /// Number of samples.
    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Whether the ensemble has no samples, and therefore no states.
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// One topology-bound sample by stable index.
    pub fn get(&self, index: usize) -> Option<EnsembleSampleView<'_>> {
        (index < self.samples.len()).then_some(EnsembleSampleView {
            ensemble: self,
            index,
        })
    }

    /// A dimension-preserving editor for one sample. Whole-sample
    /// replacement goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<EnsembleSampleMut<'_>> {
        self.samples
            .get_mut(index)
            .map(|sample| EnsembleSampleMut { sample })
    }

    /// Topology-bound samples in stable order.
    pub fn iter(&self) -> EnsembleIter<'_> {
        EnsembleIter {
            ensemble: self,
            indices: 0..self.samples.len(),
        }
    }

    /// The grouping of the samples into states.
    pub fn partition(&self) -> &EnsemblePartition {
        &self.partition
    }

    /// Number of states.
    pub fn state_count(&self) -> usize {
        self.states.len()
    }

    /// One state by index, with its representative, samples, and
    /// population.
    pub fn state(&self, index: usize) -> Option<EnsembleStateView<'_>> {
        (index < self.states.len()).then_some(EnsembleStateView {
            ensemble: self,
            index,
        })
    }

    /// The energy and properties of one state. Its samples and
    /// representative change only through [`Self::regroup`].
    pub fn state_mut(&mut self, index: usize) -> Option<&mut EnsembleState> {
        self.states.get_mut(index)
    }

    /// The states in order.
    pub fn states(&self) -> EnsembleStateIter<'_> {
        EnsembleStateIter {
            ensemble: self,
            indices: 0..self.states.len(),
        }
    }

    /// The temperature the populations describe, in
    /// [`CANONICAL_TEMPERATURE_UNIT`], if recorded.
    ///
    /// [`Self::assign_boltzmann_populations`] records it. Changing weights
    /// by other means does not clear it; the caller states what its
    /// populations describe.
    pub const fn temperature(&self) -> Option<Quantity<f64>> {
        self.temperature
    }

    /// Records or clears the temperature the populations describe, which
    /// must be finite and positive. No weight changes.
    pub fn set_temperature(
        &mut self,
        temperature: Option<Quantity<f64>>,
    ) -> Result<(), EnsembleError> {
        self.temperature = temperature
            .map(|temperature| {
                let kelvin = checked_temperature(temperature)?;
                Ok::<_, EnsembleError>(Quantity::new(kelvin, CANONICAL_TEMPERATURE_UNIT))
            })
            .transpose()?;
        Ok(())
    }

    /// Validates and appends one sample as a new state, the sample its
    /// representative, returning the state's index.
    pub fn push_state(
        &mut self,
        mut sample: EnsembleSample,
        state: EnsembleState,
    ) -> Result<usize, EnsembleError> {
        self.bind(&mut sample)?;
        let index = self.states.len();
        let (mut assignments, mut representatives) = self.take_partition();
        assignments.push(index);
        representatives.push(self.samples.len());
        self.partition = EnsemblePartition::from_parts_unchecked(assignments, representatives);
        self.samples.push(sample);
        self.states.push(state);
        Ok(index)
    }

    /// Validates and appends one sample to an existing state.
    pub fn push_sample(
        &mut self,
        mut sample: EnsembleSample,
        state: usize,
    ) -> Result<(), EnsembleError> {
        self.check_state(state)?;
        self.bind(&mut sample)?;
        let (mut assignments, representatives) = self.take_partition();
        assignments.push(state);
        self.partition = EnsemblePartition::from_parts_unchecked(assignments, representatives);
        self.samples.push(sample);
        Ok(())
    }

    /// Validates and replaces one sample, returning the previous one. The
    /// sample stays in its state. The index is checked first; failure leaves
    /// the ensemble unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut sample: EnsembleSample,
    ) -> Result<EnsembleSample, EnsembleError> {
        self.check_index(index)?;
        self.bind(&mut sample)?;
        Ok(std::mem::replace(&mut self.samples[index], sample))
    }

    /// Removes and returns one sample, shifting later samples down.
    ///
    /// Removing the only sample of a state removes the state, shifting
    /// later states down. The representative of a state with other samples
    /// cannot be removed ([`EnsembleError::RepresentativeRemoval`]); choose
    /// another representative with [`Self::regroup`] first.
    pub fn remove(&mut self, index: usize) -> Result<EnsembleSample, EnsembleError> {
        self.check_index(index)?;
        let state = self.partition.assignments()[index];
        let representative = self.partition.representatives()[state] == index;
        let alone = self.partition.samples(state).nth(1).is_none();
        if representative && !alone {
            return Err(EnsembleError::RepresentativeRemoval {
                sample: index,
                state,
            });
        }
        let (mut assignments, mut representatives) = self.take_partition();
        assignments.remove(index);
        if alone {
            representatives.remove(state);
            self.states.remove(state);
            for assigned in &mut assignments {
                if *assigned > state {
                    *assigned -= 1;
                }
            }
        }
        for represented in &mut representatives {
            if *represented > index {
                *represented -= 1;
            }
        }
        self.partition = EnsemblePartition::from_parts_unchecked(assignments, representatives);
        Ok(self.samples.remove(index))
    }

    /// Replaces the positions of every sample transactionally, in order.
    /// Other state is kept.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), EnsembleError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.samples.len() {
            return Err(EnsembleError::SampleCountMismatch {
                expected: self.samples.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(StructureError::Array(DenseArrayError::CountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                })
                .into());
            }
        }
        for (sample, positions) in self.samples.iter_mut().zip(positions) {
            sample.positions = positions;
        }
        Ok(())
    }

    /// Copies samples in the requested order, sharing the topology. Ranges
    /// and strides can be written as `(start..end).step_by(stride)`;
    /// duplicates and reordering are allowed. Weights, the temperature, and
    /// ensemble properties are kept.
    ///
    /// Every state with a selected sample becomes a state of the copy, in
    /// state order, represented by the first copy of its representative or,
    /// when that is not selected, by its first selected sample. A state
    /// whose samples are all selected keeps its energy and properties; a
    /// state selected only in part has neither, as a split state in
    /// [`Self::regroup`]. Use [`Self::select_states`] to keep whole states.
    pub fn select(&self, indices: impl IntoIterator<Item = usize>) -> Result<Self, EnsembleError> {
        let indices = indices
            .into_iter()
            .map(|index| self.check_index(index).map(|()| index))
            .collect::<Result<Vec<_>, _>>()?;
        let old = self.partition.assignments();
        let mut selected = vec![0; self.states.len()];
        let mut chosen = vec![false; self.samples.len()];
        for index in &indices {
            if !std::mem::replace(&mut chosen[*index], true) {
                selected[old[*index]] += 1;
            }
        }
        let sizes = self.partition.sizes();
        let mut renumbered = vec![usize::MAX; self.states.len()];
        let mut states = Vec::new();
        for (state, count) in selected.iter().enumerate() {
            if *count > 0 {
                renumbered[state] = states.len();
                states.push(if *count == sizes[state] {
                    self.states[state].clone()
                } else {
                    EnsembleState::new()
                });
            }
        }
        let mut assignments = Vec::with_capacity(indices.len());
        let mut representatives = vec![None; states.len()];
        let mut first = vec![None; states.len()];
        for (position, index) in indices.iter().enumerate() {
            let state = renumbered[old[*index]];
            assignments.push(state);
            first[state].get_or_insert(position);
            if *index == self.partition.representatives()[old[*index]] {
                representatives[state].get_or_insert(position);
            }
        }
        let representatives = representatives
            .into_iter()
            .zip(first)
            .map(|(representative, first)| {
                representative
                    .or(first)
                    .expect("every renumbered state has a selected sample")
            })
            .collect();
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            samples: indices
                .iter()
                .map(|index| self.samples[*index].clone())
                .collect(),
            partition: EnsemblePartition::from_parts_unchecked(assignments, representatives),
            states,
            temperature: self.temperature,
        })
    }

    /// Copies whole states in the requested order, each with all of its
    /// samples in sample order, its representative, energy, and properties.
    /// Duplicates and reordering are allowed. Weights, the temperature, and
    /// ensemble properties are kept, so populations stay relative to one
    /// another; [`Self::covering_population`],
    /// [`Self::within_energy_window`], and [`Self::ranked_by_population`]
    /// give states to extract.
    pub fn select_states(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, EnsembleError> {
        let indices = indices
            .into_iter()
            .map(|index| self.check_state(index).map(|()| index))
            .collect::<Result<Vec<_>, _>>()?;
        let mut samples = Vec::new();
        let mut assignments = Vec::new();
        let mut representatives = Vec::with_capacity(indices.len());
        let mut states = Vec::with_capacity(indices.len());
        for (state, old) in indices.iter().enumerate() {
            let representative = self.partition.representatives()[*old];
            for sample in self.partition.samples(*old) {
                if sample == representative {
                    representatives.push(samples.len());
                }
                samples.push(self.samples[sample].clone());
                assignments.push(state);
            }
            states.push(self.states[*old].clone());
        }
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            samples,
            partition: EnsemblePartition::from_parts_unchecked(assignments, representatives),
            states,
            temperature: self.temperature,
        })
    }

    /// Publishes one induced topology subset and projects every sample onto
    /// it. Ensemble, state, and per-sample owner properties are dropped;
    /// weights, states, energies, and the temperature describe the sampled
    /// structure and are kept.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, EnsembleError> {
        let subset = crate::topology::Topology::subset(self.topology.clone(), selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let samples = self
            .samples
            .iter()
            .map(|sample| sample.project(&atoms, &bonds))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            samples,
            partition: self.partition.clone(),
            states: self
                .states
                .iter()
                .map(|state| EnsembleState {
                    energy: state.energy,
                    properties: OwnerProperties::new(),
                })
                .collect(),
            temperature: self.temperature,
        })
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every sample,
    /// state, and the ensemble properties without copying. Failure leaves
    /// the ensemble unchanged. The new snapshot shares the original layout,
    /// so selections, buffers, and prepared potentials stay usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// [`Self::perceive`], stopping partway once `cancellation` is requested.
    ///
    /// Delegates to [`Topology::perceived_cancellable`]. A cancelled call
    /// leaves the ensemble unchanged; a completed one installs exactly what
    /// [`Self::perceive`] would.
    ///
    /// # Errors
    ///
    /// As [`Topology::perceived_cancellable`].
    pub fn perceive_cancellable(
        &mut self,
        cancellation: &Cancellation,
    ) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived_cancellable(cancellation)?);
        Ok(())
    }

    /// Length-preserving mutable access to the shared topology's static
    /// annotations.
    ///
    /// A topology shared with other owners is first replaced by a snapshot of
    /// it ([`Arc::make_mut`]), so those owners are unaffected. Either way the
    /// layout is kept, so every sample, selection, buffer, and prepared
    /// potential remains usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// Consumes the ensemble, returning its samples in order and dropping
    /// the topology, the states, and the ensemble properties.
    pub fn into_samples(self) -> Vec<EnsembleSample> {
        self.samples
    }

    /// Groups the samples into the states of `partition`, changing no
    /// weight and no sample order.
    ///
    /// The distribution stays the same; only its description changes, and
    /// populations follow the samples. A new state that contains every
    /// sample of an old state keeps the energy and properties of its
    /// representative's old state, so merging whole states (near-duplicate
    /// conformers) keeps the representative's energy, and an unchanged
    /// grouping keeps everything. A state whose old state is split apart
    /// has neither. The temperature is kept.
    ///
    /// Fails with [`EnsembleError::PartitionMismatch`] when `partition`
    /// covers a different number of samples; failure changes nothing.
    pub fn regroup(&mut self, partition: EnsemblePartition) -> Result<(), EnsembleError> {
        if partition.sample_count() != self.samples.len() {
            return Err(EnsembleError::PartitionMismatch {
                expected: self.samples.len(),
                actual: partition.sample_count(),
            });
        }
        let old = self.partition.assignments();
        let new = partition.assignments();
        // The one new state every sample of an old state went to, if it
        // stayed whole.
        let mut target = vec![None; self.states.len()];
        let mut split = vec![false; self.states.len()];
        for (old, new) in old.iter().zip(new) {
            match target[*old] {
                None => target[*old] = Some(*new),
                Some(existing) if existing != *new => split[*old] = true,
                Some(_) => {}
            }
        }
        let mut previous = std::mem::take(&mut self.states);
        self.states = partition
            .representatives()
            .iter()
            .map(|representative| {
                let state = old[*representative];
                if split[state] {
                    EnsembleState::new()
                } else {
                    std::mem::take(&mut previous[state])
                }
            })
            .collect();
        self.partition = partition;
        Ok(())
    }

    /// A new ensemble of one sample per state, in state order: a copy of
    /// each state's representative carrying the state's population as its
    /// weight, so populations are kept.
    ///
    /// It coarse-grains, losing the spread within each state. Each state
    /// keeps its energy and properties; the topology, the ensemble
    /// properties, and the temperature are shared, and the source is
    /// unchanged.
    ///
    /// Fails when a population is not a finite weight.
    pub fn collapse(&self) -> Result<Self, EnsembleError> {
        let populations = self.state_populations();
        let samples = self
            .partition
            .representatives()
            .iter()
            .zip(populations)
            .map(|(representative, population)| {
                let mut sample = self.samples[*representative].clone();
                sample.weight = checked_weight(population)?;
                Ok(sample)
            })
            .collect::<Result<Vec<_>, EnsembleError>>()?;
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            partition: EnsemblePartition::state_per_sample(samples.len()),
            samples,
            states: self.states.clone(),
            temperature: self.temperature,
        })
    }

    /// Scales every weight so the weights sum to one, keeping their ratios
    /// and therefore every probability and population share.
    ///
    /// Fails on an empty ensemble, which has no probabilities, or when a
    /// scaled weight would underflow to zero. Failure changes nothing.
    pub fn normalize_weights(&mut self) -> Result<(), EnsembleError> {
        let probabilities = self.probabilities()?;
        for (sample, weight) in self.samples.iter_mut().zip(probabilities) {
            sample.weight = weight;
        }
        Ok(())
    }

    /// The probability of each sample, in sample order: its weight divided
    /// by the total weight. The ensemble is not changed.
    ///
    /// Fails on an empty ensemble or when a probability would underflow to
    /// zero.
    pub fn probabilities(&self) -> Result<Vec<f64>, EnsembleError> {
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        Ok(scaled
            .into_iter()
            .map(|weight| checked_weight(weight / total))
            .collect::<Result<Vec<_>, _>>()?)
    }

    /// The population share of each state, in state order: the total weight
    /// of its samples divided by the total weight. The shares sum to one.
    ///
    /// Fails on an empty ensemble or when a share would underflow to zero.
    pub fn state_probabilities(&self) -> Result<Vec<f64>, EnsembleError> {
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let scaled = self.scaled_state_populations();
        let total = scaled.iter().sum::<f64>();
        Ok(scaled
            .into_iter()
            .map(|population| checked_weight(population / total))
            .collect::<Result<Vec<_>, _>>()?)
    }

    /// Kish's effective sample size over samples, `(Σw)² / Σw²`: the number
    /// of equally weighted samples that would give the same statistical
    /// precision.
    ///
    /// It is the sample count for equal weights and approaches one when a
    /// single sample dominates. Fails on an empty ensemble.
    pub fn effective_sample_size(&self) -> Result<f64, EnsembleError> {
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        let squares = scaled.iter().map(|weight| weight * weight).sum::<f64>();
        Ok(total * total / squares)
    }

    /// Sample indices by descending weight; equal weights keep sample order.
    ///
    /// Select the ten most probable samples with
    /// `ensemble.select(ensemble.ranked_by_weight().into_iter().take(10))`.
    pub fn ranked_by_weight(&self) -> Vec<usize> {
        let mut indices = (0..self.samples.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| self.samples[b].weight.total_cmp(&self.samples[a].weight));
        indices
    }

    /// State indices by descending population; equal populations keep
    /// state order.
    ///
    /// Extract the ten most populated states with
    /// `ensemble.select_states(ensemble.ranked_by_population().into_iter().take(10))`.
    pub fn ranked_by_population(&self) -> Vec<usize> {
        let populations = self.scaled_state_populations();
        let mut indices = (0..populations.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| populations[b].total_cmp(&populations[a]));
        indices
    }

    /// The fewest states, by descending population, whose populations
    /// together reach `fraction` of the total, in
    /// [`Self::ranked_by_population`] order.
    ///
    /// `fraction` must lie in `(0, 1]`; `1` selects every state, including
    /// those too light to change a floating-point total. Pass the result to
    /// [`Self::select_states`] to keep, for example, 95% of a Boltzmann
    /// population. Fails on an empty ensemble.
    pub fn covering_population(&self, fraction: f64) -> Result<Vec<usize>, EnsembleError> {
        if !(fraction > 0.0 && fraction <= 1.0) {
            return Err(EnsembleError::InvalidFraction);
        }
        if self.is_empty() {
            return Err(EnsembleError::EmptySource);
        }
        let ranked = self.ranked_by_population();
        if fraction == 1.0 {
            return Ok(ranked);
        }
        let scaled = self.scaled_state_populations();
        // Sum in ranked order so the running total ends exactly at the total.
        let total = ranked.iter().map(|&index| scaled[index]).sum::<f64>();
        let target = fraction * total;
        let mut covered = 0.0;
        let count = ranked
            .iter()
            .position(|&index| {
                covered += scaled[index];
                covered >= target
            })
            .map_or(ranked.len(), |position| position + 1);
        Ok(ranked[..count].to_vec())
    }

    /// State indices by ascending energy; equal energies keep state order.
    ///
    /// Fails with [`EnsembleError::MissingEnergy`] if any state has none.
    pub fn ranked_by_energy(&self) -> Result<Vec<usize>, EnsembleError> {
        let energies = self.state_energies()?;
        let mut indices = (0..energies.len()).collect::<Vec<_>>();
        indices.sort_by(|&a, &b| {
            energies[a]
                .partial_cmp(&energies[b])
                .expect("state energies are finite")
        });
        Ok(indices)
    }

    /// Indices, in state order, of the states whose energy lies at most
    /// `window` above the lowest state energy.
    ///
    /// `window` is a finite, nonnegative energy, such as the 6 kcal/mol
    /// window common in conformer searches. Fails if any state has no
    /// energy.
    pub fn within_energy_window(&self, window: Quantity<f64>) -> Result<Vec<usize>, EnsembleError> {
        let window = window.value_in(CANONICAL_ENERGY_UNIT)?;
        if !(window.is_finite() && window >= 0.0) {
            return Err(EnsembleError::InvalidEnergyWindow);
        }
        let energies = self.state_energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        Ok(energies
            .iter()
            .enumerate()
            .filter(|(_, energy)| *energy - minimum <= window)
            .map(|(index, _)| index)
            .collect())
    }

    /// Sets each state's population in proportion to the Boltzmann factor
    /// of its energy relative to the lowest one, `exp(-(E - E_min) / RT)`,
    /// and records `temperature`.
    ///
    /// The weights of each state's samples are scaled together, so their
    /// ratios within the state are kept, and the total weight is kept. An
    /// ensemble of one state is therefore unchanged. Samples drawn from a
    /// simulation already follow their distribution; group them into a state
    /// without an energy, which this refuses, rather than weighting them a
    /// second time.
    ///
    /// `temperature` must be finite and positive. Fails if any state has no
    /// energy ([`EnsembleError::MissingEnergy`]) or a weight would underflow
    /// to zero (restrict the ensemble to an energy window first). Failure
    /// changes nothing.
    pub fn assign_boltzmann_populations(
        &mut self,
        temperature: Quantity<f64>,
    ) -> Result<(), EnsembleError> {
        let temperature = checked_temperature(temperature)?;
        let thermal_energy = MOLAR_GAS_CONSTANT * temperature;
        if !(thermal_energy.is_finite() && thermal_energy > 0.0) {
            return Err(EnsembleError::InvalidTemperature);
        }
        let energies = self.state_energies()?;
        let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
        let targets = energies
            .iter()
            .map(|energy| (-(energy - minimum) / thermal_energy).exp())
            .collect::<Vec<_>>();
        let target_total = targets.iter().sum::<f64>();
        let scaled = self.scaled_weights();
        let total = scaled.iter().sum::<f64>();
        let maximum = self.max_weight();
        let populations = self.scaled_state_populations();
        // Each state's share of the total weight becomes its target share:
        // its samples scale by that share times the total over its
        // population. A state far lighter than the heaviest sample has a
        // population that underflows when scaled by the heaviest weight, so
        // its samples are normalized within the state instead, keeping every
        // intermediate value representable.
        let mut state_maxima = vec![0.0_f64; self.states.len()];
        for (sample, state) in self.samples.iter().zip(self.partition.assignments()) {
            state_maxima[*state] = state_maxima[*state].max(sample.weight);
        }
        let mut within = vec![0.0; self.states.len()];
        for (sample, state) in self.samples.iter().zip(self.partition.assignments()) {
            within[*state] += sample.weight / state_maxima[*state];
        }
        let shares = targets
            .iter()
            .map(|target| target / target_total)
            .collect::<Vec<_>>();
        let weights = self
            .samples
            .iter()
            .zip(self.partition.assignments())
            .map(|(sample, state)| {
                let factor = shares[*state] * (total / populations[*state]);
                let weight = if factor.is_finite() {
                    sample.weight * factor
                } else {
                    (sample.weight / state_maxima[*state])
                        * (shares[*state] * total / within[*state])
                        * maximum
                };
                checked_weight(weight)
            })
            .collect::<Result<Vec<_>, _>>()?;
        for (sample, weight) in self.samples.iter_mut().zip(weights) {
            sample.weight = weight;
        }
        self.temperature = Some(Quantity::new(temperature, CANONICAL_TEMPERATURE_UNIT));
        Ok(())
    }

    /// The largest weight, or zero for no samples.
    fn max_weight(&self) -> f64 {
        self.samples
            .iter()
            .map(|sample| sample.weight)
            .fold(0.0_f64, f64::max)
    }

    /// Weights divided by the largest one, so sums cannot overflow.
    pub(crate) fn scaled_weights(&self) -> Vec<f64> {
        let maximum = self.max_weight();
        self.samples
            .iter()
            .map(|sample| sample.weight / maximum)
            .collect()
    }

    /// The summed [`Self::scaled_weights`] of each state, in sample order.
    fn scaled_state_populations(&self) -> Vec<f64> {
        let mut populations = vec![0.0; self.states.len()];
        for (weight, state) in self
            .scaled_weights()
            .into_iter()
            .zip(self.partition.assignments())
        {
            populations[*state] += weight;
        }
        populations
    }

    /// The summed weight of each state, in sample order.
    fn state_populations(&self) -> Vec<f64> {
        let mut populations = vec![0.0; self.states.len()];
        for (sample, state) in self.samples.iter().zip(self.partition.assignments()) {
            populations[*state] += sample.weight;
        }
        populations
    }

    /// Every state energy in [`CANONICAL_ENERGY_UNIT`].
    pub(crate) fn state_energies(&self) -> Result<Vec<f64>, EnsembleError> {
        self.states
            .iter()
            .enumerate()
            .map(|(state, payload)| {
                payload
                    .energy
                    .map(Quantity::into_value)
                    .ok_or(EnsembleError::MissingEnergy { state })
            })
            .collect()
    }

    /// Applies one rigid transform per sample transactionally.
    pub(crate) fn transform_samples(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), StructureError> {
        debug_assert_eq!(transforms.len(), self.samples.len());
        let samples = self
            .samples
            .iter()
            .zip(transforms)
            .map(|(sample, transform)| sample.transformed(*transform))
            .collect::<Result<Vec<_>, _>>()?;
        self.samples = samples;
        Ok(())
    }

    /// Moves the partition out as its assignments and representatives, to
    /// be rebuilt and put back.
    fn take_partition(&mut self) -> (Vec<usize>, Vec<usize>) {
        std::mem::replace(&mut self.partition, EnsemblePartition::state_per_sample(0)).into_parts()
    }

    fn check_index(&self, index: usize) -> Result<(), EnsembleError> {
        if index >= self.samples.len() {
            return Err(EnsembleError::IndexOutOfRange {
                index,
                len: self.samples.len(),
            });
        }
        Ok(())
    }

    fn check_state(&self, state: usize) -> Result<(), EnsembleError> {
        if state >= self.states.len() {
            return Err(EnsembleError::StateOutOfRange {
                state,
                states: self.states.len(),
            });
        }
        Ok(())
    }

    fn bind(&self, sample: &mut EnsembleSample) -> Result<(), EnsembleError> {
        validate_binding(sample.atom_count(), &sample.properties, &self.topology)?;
        sample.properties.bind_bond_rows(self.topology.bond_count());
        Ok(())
    }
}

impl<'a> IntoIterator for &'a Ensemble {
    type Item = EnsembleSampleView<'a>;
    type IntoIter = EnsembleIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// One ensemble sample: positions in dense topology atom order, an optional
/// periodic cell, optional crystallographic occupancies and B-factors,
/// structure properties, and a finite, positive, relative weight.
///
/// A sample has no energy; its [`EnsembleState`] may. The sample carries no
/// topology; it is bound when inserted into an [`Ensemble`].
#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleSample {
    positions: Positions,
    cell: Option<PeriodicCell>,
    occupancies: AtomScalars,
    b_factors: AtomScalars,
    properties: StructureProperties,
    weight: f64,
}

impl EnsembleSample {
    /// Creates a non-periodic sample without occupancies, B-factors, or
    /// properties, with a finite positive relative weight.
    pub fn new(positions: Positions, weight: f64) -> Result<Self, StructureError> {
        let mut sample = Self::unit(positions);
        sample.weight = checked_weight(weight)?;
        Ok(sample)
    }

    /// A sample of an equally weighted collection.
    fn unit(positions: Positions) -> Self {
        let properties = StructureProperties::new(positions.len());
        Self {
            positions,
            cell: None,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            weight: 1.0,
        }
    }

    /// An equally weighted sample holding a model's structure state; the
    /// model's topology is dropped.
    pub(crate) fn from_model(model: Model) -> Self {
        Self {
            positions: model.positions,
            cell: model.cell,
            occupancies: model.occupancies,
            b_factors: model.b_factors,
            properties: model.properties,
            weight: 1.0,
        }
    }

    /// An equally weighted sample without occupancies or B-factors, for
    /// frames reinterpreted as an ensemble.
    pub(crate) fn from_frame_state(
        positions: Positions,
        cell: Option<PeriodicCell>,
        properties: StructureProperties,
    ) -> Self {
        Self {
            positions,
            cell,
            occupancies: AtomScalars::default(),
            b_factors: AtomScalars::default(),
            properties,
            weight: 1.0,
        }
    }

    /// Returns this sample with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    /// Number of atoms, the length of every dense per-atom array of the sample.
    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// The relative statistical weight; always finite and positive.
    pub const fn weight(&self) -> f64 {
        self.weight
    }

    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), StructureError> {
        self.weight = checked_weight(weight)?;
        Ok(())
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }

    fn model_view<'a>(
        &'a self,
        topology: &'a Arc<Topology>,
        kind: StructureKind<'a>,
    ) -> ModelView<'a> {
        ModelView::bound(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &self.occupancies,
            &self.b_factors,
            &self.properties,
            kind,
        )
    }

    /// Projects every dense row through a subset correspondence.
    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, StructureError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            occupancies: self.occupancies.select(atoms),
            b_factors: self.b_factors.select(atoms),
            properties: self.properties.project(atoms, bonds)?,
            weight: self.weight,
        })
    }

    /// Positions move and the cell rotates.
    fn transformed(&self, transform: RigidTransform) -> Result<Self, StructureError> {
        Ok(Self {
            cell: transform_cell(self.cell, transform)?,
            positions: transform_positions(&self.positions, transform)?,
            occupancies: self.occupancies.clone(),
            b_factors: self.b_factors.clone(),
            properties: self.properties.clone(),
            weight: self.weight,
        })
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    structure_readers!();
    site_readers!();
    structure_mutators!();
    site_mutators!();
}

/// One state of an [`Ensemble`]: an optional energy and owner properties.
///
/// Which samples a state groups, and which of them represents it, belong to
/// the ensemble's [`EnsemblePartition`]; [`EnsembleStateView`] reads them
/// together with the state's population.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct EnsembleState {
    energy: Option<Quantity<f64>>,
    properties: OwnerProperties,
}

impl EnsembleState {
    /// A state without an energy or properties.
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns this state with a finite energy.
    pub fn with_energy(mut self, energy: Quantity<f64>) -> Result<Self, StructureError> {
        self.set_energy(Some(energy))?;
        Ok(self)
    }

    /// The state's energy in [`CANONICAL_ENERGY_UNIT`], if known.
    ///
    /// Energies of one ensemble share one scale, such as the minimized
    /// energies of conformers or computed free energies at one level of
    /// theory, never a snapshot's instantaneous energy; only differences
    /// between states are meaningful. Record the scale in the ensemble
    /// properties.
    pub const fn energy(&self) -> Option<Quantity<f64>> {
        self.energy
    }

    /// Sets or clears a finite energy, converted to the canonical energy
    /// unit. No weight changes.
    pub fn set_energy(&mut self, energy: Option<Quantity<f64>>) -> Result<(), StructureError> {
        self.energy = energy
            .map(|energy| {
                let energy = energy.into_unit(CANONICAL_ENERGY_UNIT)?;
                if !energy.value().is_finite() {
                    return Err(StructureError::NonFiniteEnergy);
                }
                Ok(energy)
            })
            .transpose()?;
        Ok(())
    }

    /// State-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    /// Mutable state-level annotations.
    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }
}

/// The molar gas constant `N_A k_B` in kJ/(mol K), exact since the 2019 SI.
const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24e-3;

fn checked_weight(weight: f64) -> Result<f64, StructureError> {
    if weight.is_finite() && weight > 0.0 {
        Ok(weight)
    } else {
        Err(StructureError::InvalidWeight)
    }
}

/// A finite, positive temperature in [`CANONICAL_TEMPERATURE_UNIT`].
fn checked_temperature(temperature: Quantity<f64>) -> Result<f64, EnsembleError> {
    let temperature = temperature.value_in(CANONICAL_TEMPERATURE_UNIT)?;
    if temperature.is_finite() && temperature > 0.0 {
        Ok(temperature)
    } else {
        Err(EnsembleError::InvalidTemperature)
    }
}

/// One topology-bound [`EnsembleSample`] borrowed from an [`Ensemble`], with
/// its state.
///
/// Reads sample state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying.
#[derive(Debug, Clone, Copy)]
pub struct EnsembleSampleView<'a> {
    ensemble: &'a Ensemble,
    index: usize,
}

impl std::ops::Deref for EnsembleSampleView<'_> {
    type Target = EnsembleSample;

    fn deref(&self) -> &Self::Target {
        &self.ensemble.samples[self.index]
    }
}

impl<'a> EnsembleSampleView<'a> {
    /// The sample's stable index in its ensemble.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The state the sample belongs to.
    pub fn state(self) -> EnsembleStateView<'a> {
        EnsembleStateView {
            ensemble: self.ensemble,
            index: self.ensemble.partition.assignments()[self.index],
        }
    }

    /// Whether the sample represents its state.
    pub fn is_representative(self) -> bool {
        self.state().representative() == self.index
    }

    /// The topology the sample is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        &self.ensemble.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        self.ensemble.shared_topology()
    }

    /// The borrowed sample with the view's full lifetime.
    pub fn payload(self) -> &'a EnsembleSample {
        &self.ensemble.samples[self.index]
    }

    /// Borrows the sample as the common kernel input, with the view's full
    /// lifetime. Its [`StructureKind::Sample`] carries the sample's weight
    /// and state, the state's energy, and the ensemble's temperature.
    pub fn as_model_view(self) -> ModelView<'a> {
        let state = self.state();
        self.payload().model_view(
            &self.ensemble.topology,
            StructureKind::Sample {
                weight: self.weight(),
                state: state.index(),
                representative: state.representative() == self.index,
                energy: state.payload().energy.as_ref(),
                temperature: self.ensemble.temperature.as_ref(),
            },
        )
    }

    /// Materializes the sample as an owned model, dropping its weight and
    /// state.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl super::model::sealed::Sealed for EnsembleSampleView<'_> {}

impl AsModelView for EnsembleSampleView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        EnsembleSampleView::as_model_view(*self)
    }
}

/// One [`EnsembleState`] borrowed from an [`Ensemble`], with its samples,
/// representative, and population.
///
/// Reads the state's energy and properties through `Deref`.
#[derive(Debug, Clone, Copy)]
pub struct EnsembleStateView<'a> {
    ensemble: &'a Ensemble,
    index: usize,
}

impl std::ops::Deref for EnsembleStateView<'_> {
    type Target = EnsembleState;

    fn deref(&self) -> &Self::Target {
        &self.ensemble.states[self.index]
    }
}

impl<'a> EnsembleStateView<'a> {
    /// The state's index in its ensemble.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The index of the state's representative sample.
    pub fn representative(self) -> usize {
        self.ensemble.partition.representatives()[self.index]
    }

    /// The state's representative sample.
    pub fn representative_sample(self) -> EnsembleSampleView<'a> {
        EnsembleSampleView {
            ensemble: self.ensemble,
            index: self.representative(),
        }
    }

    /// The indices of the state's samples, in sample order.
    pub fn samples(self) -> impl Iterator<Item = usize> + 'a {
        self.ensemble.partition.samples(self.index)
    }

    /// The population: the total weight of the state's samples, relative
    /// like the weights. [`Ensemble::state_probabilities`] gives every
    /// state's share of the total.
    pub fn population(self) -> f64 {
        self.samples()
            .map(|sample| self.ensemble.samples[sample].weight)
            .sum()
    }

    /// The borrowed state with the view's full lifetime.
    pub fn payload(self) -> &'a EnsembleState {
        &self.ensemble.states[self.index]
    }
}

/// Dimension-preserving mutable access to one sample of an [`Ensemble`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// sample through [`Ensemble::replace`].
#[derive(Debug)]
pub struct EnsembleSampleMut<'a> {
    sample: &'a mut EnsembleSample,
}

impl std::ops::Deref for EnsembleSampleMut<'_> {
    type Target = EnsembleSample;

    fn deref(&self) -> &Self::Target {
        self.sample
    }
}

impl EnsembleSampleMut<'_> {
    /// Sets a finite positive relative weight.
    pub fn set_weight(&mut self, weight: f64) -> Result<(), StructureError> {
        self.sample.set_weight(weight)
    }

    /// Applies a rigid transform transactionally: positions move and the
    /// cell rotates.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        self.sample.apply_rigid_transform(transform)
    }

    fn storage(&mut self) -> &mut EnsembleSample {
        self.sample
    }

    structure_mutators!();
    site_mutators!();
}

/// Iterator over the topology-bound samples of an [`Ensemble`].
#[derive(Debug, Clone)]
pub struct EnsembleIter<'a> {
    ensemble: &'a Ensemble,
    indices: std::ops::Range<usize>,
}

impl<'a> Iterator for EnsembleIter<'a> {
    type Item = EnsembleSampleView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let ensemble = self.ensemble;
        self.indices
            .next()
            .map(|index| EnsembleSampleView { ensemble, index })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for EnsembleIter<'_> {}

/// Iterator over the states of an [`Ensemble`].
#[derive(Debug, Clone)]
pub struct EnsembleStateIter<'a> {
    ensemble: &'a Ensemble,
    indices: std::ops::Range<usize>,
}

impl<'a> Iterator for EnsembleStateIter<'a> {
    type Item = EnsembleStateView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let ensemble = self.ensemble;
        self.indices
            .next()
            .map(|index| EnsembleStateView { ensemble, index })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for EnsembleStateIter<'_> {}

/// Failure of an [`Ensemble`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EnsembleError {
    /// A constructor received no models, or weights of an empty ensemble
    /// were normalized.
    EmptySource,
    /// A sample index is outside the ensemble.
    IndexOutOfRange {
        /// The requested sample index.
        index: usize,
        /// The number of samples.
        len: usize,
    },
    /// A state index is outside the ensemble.
    StateOutOfRange {
        /// The requested state index.
        state: usize,
        /// The number of states.
        states: usize,
    },
    /// A model belongs to a different topology layout.
    TopologyMismatch,
    /// An ensemble-wide replacement has the wrong number of samples.
    SampleCountMismatch {
        /// The ensemble's sample count.
        expected: usize,
        /// The number of replacements supplied.
        actual: usize,
    },
    /// A partition covers a different number of samples than the ensemble.
    PartitionMismatch {
        /// The ensemble's sample count.
        expected: usize,
        /// The partition's sample count.
        actual: usize,
    },
    /// The representative of a state with other samples cannot be removed.
    RepresentativeRemoval {
        /// The representative sample.
        sample: usize,
        /// Its state.
        state: usize,
    },
    /// An energy-based operation found a state without an energy.
    MissingEnergy {
        /// The index of the state without an energy.
        state: usize,
    },
    /// A temperature is not finite and positive.
    InvalidTemperature,
    /// An energy window is not finite and nonnegative.
    InvalidEnergyWindow,
    /// A population fraction is not in `(0, 1]`.
    InvalidFraction,
    /// A sample's structural state is invalid; see the source.
    Structure(StructureError),
    /// Projecting the topology onto a selection failed; see the source.
    Subset(TopologySubsetError),
    /// Building the topology for the ensemble failed; see the source.
    TopologyBuild(TopologyBuildError),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
}

impl fmt::Display for EnsembleError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySource => formatter.write_str("ensemble source is empty"),
            Self::IndexOutOfRange { index, len } => {
                write!(formatter, "index {index} is out of range for {len} samples")
            }
            Self::StateOutOfRange { state, states } => {
                write!(
                    formatter,
                    "state {state} is out of range for {states} states"
                )
            }
            Self::TopologyMismatch => {
                formatter.write_str("model belongs to a different topology layout")
            }
            Self::SampleCountMismatch { expected, actual } => write!(
                formatter,
                "ensemble requires {expected} replacement samples, but received {actual}"
            ),
            Self::PartitionMismatch { expected, actual } => write!(
                formatter,
                "partition covers {actual} samples, but the ensemble has {expected}"
            ),
            Self::RepresentativeRemoval { sample, state } => write!(
                formatter,
                "sample {sample} represents state {state}, which has other samples"
            ),
            Self::MissingEnergy { state } => write!(formatter, "state {state} has no energy"),
            Self::InvalidTemperature => {
                formatter.write_str("temperature must be finite and positive")
            }
            Self::InvalidEnergyWindow => {
                formatter.write_str("energy window must be finite and nonnegative")
            }
            Self::InvalidFraction => {
                formatter.write_str("population fraction must be greater than 0 and at most 1")
            }
            Self::Structure(_) => formatter.write_str("invalid ensemble sample"),
            Self::Subset(_) => formatter.write_str("cannot subset ensemble"),
            Self::TopologyBuild(_) => formatter.write_str("cannot build topology"),
            Self::Unit(_) => formatter.write_str("invalid quantity unit"),
        }
    }
}

impl std::error::Error for EnsembleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Subset(error) => Some(error),
            Self::TopologyBuild(error) => Some(error),
            Self::Unit(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<StructureError> for EnsembleError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<TopologySubsetError> for EnsembleError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}

impl From<TopologyBuildError> for EnsembleError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}

impl From<UnitError> for EnsembleError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
