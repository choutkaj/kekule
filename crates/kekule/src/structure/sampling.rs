use std::sync::Arc;

use crate::properties::OwnerProperties;
use crate::topology::{AtomSelection, Topology};

use super::state::subset_rows;
use super::{Ensemble, EnsembleSample, StructureError, TrajectoryError, TrajectoryFrameView};

/// Which trajectory frames become the samples of an ensemble, with their
/// weights and an optional atom subset.
///
/// Sample `i` is frame [`Self::frames`]`[i]`, so the requested order is kept
/// and hosts can map samples back to frames. Ranges and strides are written as
/// `(start..end).step_by(stride)`; reordering and duplicates are allowed.
/// Samples weigh `1.0` each unless [`Self::with_weights`] gives one finite,
/// positive weight per frame. With [`Self::with_atoms`], the ensemble has the
/// induced subset topology, as [`super::Trajectory::subset`] publishes it.
///
/// Use it with [`super::Trajectory::to_ensemble`] for a loaded trajectory, or
/// prepare a [`TrajectorySampler`] for frames that arrive one at a time, such
/// as from a streaming reader.
#[derive(Debug, Clone)]
pub struct TrajectorySampling {
    frames: Vec<usize>,
    weights: Option<Vec<f64>>,
    atoms: Option<AtomSelection>,
}

impl TrajectorySampling {
    /// Samples `frames` in this order, equally weighted, with every atom.
    pub fn new(frames: impl IntoIterator<Item = usize>) -> Self {
        Self {
            frames: frames.into_iter().collect(),
            weights: None,
            atoms: None,
        }
    }

    /// Gives sample `i` the relative weight `weights[i]`.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::WeightCountMismatch`] unless there is one
    /// weight per frame, and [`TrajectoryError::InvalidWeight`] for the first
    /// weight that is not finite and positive.
    pub fn with_weights(
        mut self,
        weights: impl IntoIterator<Item = f64>,
    ) -> Result<Self, TrajectoryError> {
        let weights = weights.into_iter().collect::<Vec<_>>();
        if weights.len() != self.frames.len() {
            return Err(TrajectoryError::WeightCountMismatch {
                expected: self.frames.len(),
                actual: weights.len(),
            });
        }
        if let Some(sample) = weights
            .iter()
            .position(|weight| !(weight.is_finite() && *weight > 0.0))
        {
            return Err(TrajectoryError::InvalidWeight { sample });
        }
        self.weights = Some(weights);
        Ok(self)
    }

    /// Keeps only the selected atoms; the selection may come from any
    /// topology snapshot sharing the trajectory's layout.
    #[must_use]
    pub fn with_atoms(mut self, atoms: AtomSelection) -> Self {
        self.atoms = Some(atoms);
        self
    }

    /// The frame of each sample, in sample order.
    pub fn frames(&self) -> &[usize] {
        &self.frames
    }

    /// The weight of each sample, or `None` when every sample weighs `1.0`.
    pub fn weights(&self) -> Option<&[f64]> {
        self.weights.as_deref()
    }

    /// The atom subset, or `None` for every atom.
    pub fn atoms(&self) -> Option<&AtomSelection> {
        self.atoms.as_ref()
    }

    fn weight(&self, sample: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |weights| weights[sample])
    }
}

/// A [`TrajectorySampling`] prepared against one topology, collecting frames
/// in any order into an ensemble of one state.
///
/// Preparation publishes the subset topology once. [`Self::add`] then copies
/// what a sample keeps from each requested frame (positions, cell, and
/// structure properties, projected onto the subset) and nothing else, so the
/// memory held is that of the finished ensemble. Frames that no sample
/// requests are ignored, which lets a caller pass every frame of a stream.
/// [`super::Trajectory::to_ensemble`] uses it for loaded trajectories, and
/// `kekule-traj` for streaming readers.
///
/// ```
/// use std::sync::Arc;
/// use kekule::{smiles, topology::Topology};
/// use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
/// use kekule::structure::{TrajectorySampler, TrajectorySampling};
///
/// let topology = Arc::new(Topology::from_molecule(smiles::to_molecules("O")?.remove(0))?);
/// let frames = (0..4).map(|_| TrajectoryFrame::new(Positions::zeros(1)));
/// let trajectory = Trajectory::from_items(Arc::clone(&topology), frames)?;
///
/// let sampling = TrajectorySampling::new([3, 1]).with_weights([2.0, 1.0])?;
/// let mut sampler = TrajectorySampler::new(topology, &sampling)?;
/// for (index, frame) in trajectory.iter().enumerate() {
///     sampler.add(index, frame)?;
/// }
/// let ensemble = sampler.finish()?;
/// assert_eq!(ensemble.len(), 2);
/// assert_eq!(ensemble.get(0).unwrap().weight(), 2.0);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct TrajectorySampler {
    source: Arc<Topology>,
    topology: Arc<Topology>,
    /// Source atom and bond rows of the subset, or `None` for every atom.
    rows: Option<(Vec<usize>, Vec<usize>)>,
    sampling: TrajectorySampling,
    /// Sample indices ordered by frame, then by sample.
    by_frame: Vec<usize>,
    samples: Vec<Option<EnsembleSample>>,
}

impl TrajectorySampler {
    /// Prepares `sampling` for frames of `topology`.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::Subset`] when the atom subset cannot be
    /// published, such as an empty selection or one of another layout.
    pub fn new(
        topology: impl Into<Arc<Topology>>,
        sampling: &TrajectorySampling,
    ) -> Result<Self, TrajectoryError> {
        let source = topology.into();
        let (topology, rows) = match &sampling.atoms {
            Some(atoms) => {
                let subset = Topology::subset(Arc::clone(&source), atoms)?;
                let rows = subset_rows(&subset);
                (subset.shared_topology(), Some(rows))
            }
            None => (Arc::clone(&source), None),
        };
        let mut by_frame = (0..sampling.frames.len()).collect::<Vec<_>>();
        by_frame.sort_by_key(|&sample| sampling.frames[sample]);
        Ok(Self {
            source,
            topology,
            rows,
            samples: vec![None; sampling.frames.len()],
            sampling: sampling.clone(),
            by_frame,
        })
    }

    /// The topology of the ensemble: the subset, or the source topology.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The requested frames without duplicates, in ascending order: the
    /// frames a caller with random access needs to read.
    pub fn frames(&self) -> impl Iterator<Item = usize> + '_ {
        let frames = &self.sampling.frames;
        let mut previous = None;
        self.by_frame.iter().filter_map(move |&sample| {
            let frame = frames[sample];
            (previous.replace(frame) != Some(frame)).then_some(frame)
        })
    }

    /// Samples frame `index` into every sample that requests it, replacing
    /// any earlier copy, and returns whether one did.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::TopologyMismatch`] when `frame` does not
    /// share the prepared topology's layout, and [`TrajectoryError::Structure`]
    /// when it cannot be projected. The sampler is unchanged on error.
    pub fn add(
        &mut self,
        index: usize,
        frame: TrajectoryFrameView<'_>,
    ) -> Result<bool, TrajectoryError> {
        let frames = &self.sampling.frames;
        let start = self
            .by_frame
            .partition_point(|&sample| frames[sample] < index);
        let end = self
            .by_frame
            .partition_point(|&sample| frames[sample] <= index);
        if start == end {
            return Ok(false);
        }
        if !frame.topology().shares_layout(&self.source) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let sample = self.project(frame)?;
        for &slot in &self.by_frame[start..end] {
            let mut sample = sample.clone();
            sample
                .set_weight(self.sampling.weight(slot))
                .expect("sampling weights are finite and positive");
            self.samples[slot] = Some(sample);
        }
        Ok(true)
    }

    /// What a sample keeps of `frame`, projected onto the subset rows.
    fn project(&self, frame: TrajectoryFrameView<'_>) -> Result<EnsembleSample, StructureError> {
        let (positions, mut properties) = match &self.rows {
            Some((atoms, bonds)) => (
                frame.positions().select_indices(atoms)?,
                frame.properties().project(atoms, bonds)?,
            ),
            None => (frame.positions().clone(), frame.properties().clone()),
        };
        properties.bind_bond_rows(self.topology.bond_count());
        Ok(EnsembleSample::from_frame_state(
            positions,
            frame.cell().copied(),
            properties,
        ))
    }

    /// The ensemble of one state without an energy, sample `i` coming from
    /// frame `frames[i]`. It has no collection properties.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryError::FrameNotSampled`] naming the smallest
    /// requested frame that was never added.
    pub fn finish(self) -> Result<Ensemble, TrajectoryError> {
        let frames = &self.sampling.frames;
        if let Some(&sample) = self
            .by_frame
            .iter()
            .find(|&&sample| self.samples[sample].is_none())
        {
            return Err(TrajectoryError::FrameNotSampled {
                frame: frames[sample],
            });
        }
        let samples = self
            .samples
            .into_iter()
            .map(|sample| sample.expect("every requested frame was sampled"))
            .collect();
        Ok(Ensemble::from_bound_samples(
            self.topology,
            OwnerProperties::new(),
            samples,
        ))
    }
}
