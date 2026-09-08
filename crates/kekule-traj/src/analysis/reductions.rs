use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use kekule::geometry::Vector3;
use kekule::topology::{
    AtomSelection, InstanceAtomId, SelectionError, Topology, TopologyAtomIndex,
};
use kekule::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

use crate::{Trajectory, TrajectoryFrameView};

/// Invalid reduction inputs or a failed frame observation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum ReductionError {
    Selection(SelectionError),
    Unit(UnitError),
    EmptySelection,
    NoFrames,
    InvalidCutoff,
    SelfPair { pair: usize },
    DuplicatePair { pair: usize },
    TopologyMismatch { frame: usize },
    NumericalFailure { frame: usize, atom: InstanceAtomId },
    FrameCountOverflow,
}

impl fmt::Display for ReductionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Selection(error) => write!(f, "reduction selection: {error}"),
            Self::Unit(error) => write!(f, "reduction unit: {error}"),
            Self::EmptySelection => f.write_str("reduction requires at least one atom or pair"),
            Self::NoFrames => f.write_str("reduction requires at least one observed frame"),
            Self::InvalidCutoff => f.write_str("contact cutoff must be finite and nonnegative"),
            Self::SelfPair { pair } => write!(f, "contact pair {pair} repeats the same atom"),
            Self::DuplicatePair { pair } => write!(
                f,
                "contact pair {pair} duplicates an earlier unordered pair"
            ),
            Self::TopologyMismatch { frame } => {
                write!(f, "frame {frame} belongs to a different topology snapshot")
            }
            Self::NumericalFailure { frame, atom } => write!(
                f,
                "frame {frame} reduction exceeds finite numerical range at atom {atom}"
            ),
            Self::FrameCountOverflow => {
                f.write_str("frame count exceeds exact reduction arithmetic")
            }
        }
    }
}

impl std::error::Error for ReductionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Selection(error) => Some(error),
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

fn next_count(count: usize) -> Result<usize, ReductionError> {
    count
        .checked_add(1)
        .filter(|count| (*count as u64) <= (1_u64 << 53))
        .ok_or(ReductionError::FrameCountOverflow)
}

#[derive(Debug, Clone, Copy, Default)]
struct Moments {
    mean: Vector3,
    squared_deviations: f64,
}

/// Online, per-atom RMS fluctuation about each atom's mean Cartesian position.
///
/// Frames have equal weight; variance uses population normalization (division
/// by N). A single frame yields zero. The input selection is fixed and must be
/// nonempty. Coordinates are used as stored, ignoring cells. Align and perform
/// any periodic preprocessing explicitly before observation.
///
/// State uses O(selected atoms) memory, independent of frame count. Each
/// observation is transactional. Source frame indices are diagnostic labels;
/// callers may sample frames, and every successful call counts once.
///
/// ```
/// use kekule::topology::AtomSelection;
/// use kekule_traj::{analysis::{RmsfAccumulator, RmsfResult}, TrajectoryReader};
/// # fn reduce(reader: &mut impl TrajectoryReader, atoms: &AtomSelection)
/// # -> Result<RmsfResult, Box<dyn std::error::Error>> {
/// let mut rmsf = RmsfAccumulator::new(atoms)?;
/// let mut buffer = reader.frame_buffer();
/// let mut index = 0;
/// while reader.read_next(&mut buffer)? {
///     // Apply any intended periodic preprocessing and fitting here.
///     rmsf.observe(index, buffer.frame_view())?;
///     index += 1;
/// }
/// let result = rmsf.finish()?;
/// for (atom, fluctuation) in result.atoms() {
///     let residue = result.selection().topology().residue_for_atom(atom)?;
///     // Atom identity, optional residue, and the unit-bearing value stay associated.
///     # let _ = (residue, fluctuation);
/// }
/// # Ok(result)
/// # }
/// ```
#[derive(Debug)]
pub struct RmsfAccumulator {
    selection: AtomSelection,
    moments: Vec<Moments>,
    scratch: Vec<Moments>,
    frames: usize,
}

impl RmsfAccumulator {
    pub fn new(selection: &AtomSelection) -> Result<Self, ReductionError> {
        let len = selection.indices().len();
        if len == 0 {
            return Err(ReductionError::EmptySelection);
        }
        Ok(Self {
            selection: selection.clone(),
            moments: vec![Moments::default(); len],
            scratch: vec![Moments::default(); len],
            frames: 0,
        })
    }

    /// Observes one frame; failure preserves all accumulated statistics.
    pub fn observe(
        &mut self,
        frame_index: usize,
        frame: TrajectoryFrameView<'_>,
    ) -> Result<(), ReductionError> {
        self.selection
            .ensure_compatible(&frame.shared_topology())
            .map_err(|_| ReductionError::TopologyMismatch { frame: frame_index })?;
        let next = next_count(self.frames)?;
        let positions = frame.positions().values();
        for ((index, previous), pending) in self
            .selection
            .indices()
            .iter()
            .zip(&self.moments)
            .zip(&mut self.scratch)
        {
            let p = positions.value()[index.index()];
            let p = Vector3::new(p.x, p.y, p.z);
            *pending = if self.frames == 0 {
                Moments {
                    mean: p,
                    squared_deviations: 0.0,
                }
            } else {
                // Welford's update, with the squared increment written as a
                // nonnegative quantity to avoid subtracting nearby moments.
                let delta = p - previous.mean;
                let length = delta.x.hypot(delta.y).hypot(delta.z);
                let increment = length * (self.frames as f64 / next as f64).sqrt();
                Moments {
                    mean: previous.mean + delta / next as f64,
                    squared_deviations: previous.squared_deviations + increment * increment,
                }
            };
            if !pending.mean.is_finite() || !pending.squared_deviations.is_finite() {
                return Err(ReductionError::NumericalFailure {
                    frame: frame_index,
                    atom: self
                        .selection
                        .topology()
                        .atom_id(*index)
                        .expect("validated selection"),
                });
            }
        }
        std::mem::swap(&mut self.moments, &mut self.scratch);
        self.frames = next;
        Ok(())
    }

    pub fn frame_count(&self) -> usize {
        self.frames
    }

    /// Finishes with values in selection order, retaining source atom identity.
    pub fn finish(self) -> Result<RmsfResult, ReductionError> {
        if self.frames == 0 {
            return Err(ReductionError::NoFrames);
        }
        let values = self
            .moments
            .into_iter()
            .map(|m| (m.squared_deviations / self.frames as f64).sqrt())
            .collect();
        Ok(RmsfResult {
            selection: self.selection,
            values: Quantity::new(values, CANONICAL_LENGTH_UNIT),
            frames: self.frames,
        })
    }
}

/// Per-atom RMSF with its original topology-bound selection and sample count.
/// Selecting one CA per residue gives CA RMSF, not an implicit residue average.
#[derive(Debug, Clone)]
pub struct RmsfResult {
    selection: AtomSelection,
    values: Quantity<Vec<f64>>,
    frames: usize,
}

impl RmsfResult {
    pub fn selection(&self) -> &AtomSelection {
        &self.selection
    }
    pub fn values(&self) -> &Quantity<Vec<f64>> {
        &self.values
    }
    pub fn frame_count(&self) -> usize {
        self.frames
    }

    /// Resolves each value to its source atom. Residues can be obtained through
    /// `result.selection().topology().residue_for_atom(atom)`.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = (InstanceAtomId, Quantity<f64>)> + '_ {
        self.selection.atom_ids().zip(
            self.values
                .value()
                .iter()
                .map(|value| Quantity::new(*value, CANONICAL_LENGTH_UNIT)),
        )
    }
}

/// Counts frames where each specified atom pair is at Cartesian distance
/// **<= cutoff**. This is a geometric contact, not a hydrogen-bond definition.
///
/// Pairs retain input order and must be nonempty, distinct, and non-self; reversed
/// duplicates are rejected. Each frame is measured anew. Cells are ignored and
/// every observed frame has equal weight. Apply periodic preprocessing first.
/// State uses O(pairs) memory; failed observations preserve every count.
/// Cutoff comparisons use canonical floating-point values without a tolerance.
#[derive(Debug)]
pub struct ContactOccupancyAccumulator {
    topology: Arc<Topology>,
    pairs: Vec<(InstanceAtomId, InstanceAtomId)>,
    indices: Vec<(TopologyAtomIndex, TopologyAtomIndex)>,
    hits: Vec<usize>,
    scratch: Vec<bool>,
    cutoff: Quantity<f64>,
    frames: usize,
}

impl ContactOccupancyAccumulator {
    pub fn new(
        topology: &Arc<Topology>,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
        cutoff: Quantity<f64>,
    ) -> Result<Self, ReductionError> {
        let cutoff = cutoff
            .into_unit(CANONICAL_LENGTH_UNIT)
            .map_err(ReductionError::Unit)?;
        if !cutoff.value().is_finite() || *cutoff.value() < 0.0 {
            return Err(ReductionError::InvalidCutoff);
        }
        let pairs = pairs.into_iter().collect::<Vec<_>>();
        if pairs.is_empty() {
            return Err(ReductionError::EmptySelection);
        }
        let mut seen = BTreeSet::new();
        let mut indices = Vec::with_capacity(pairs.len());
        for (pair, &(a, b)) in pairs.iter().enumerate() {
            let index = |atom| {
                topology.atom_index(atom).ok_or(ReductionError::Selection(
                    SelectionError::InvalidAtomId(atom),
                ))
            };
            let (a_index, b_index) = (index(a)?, index(b)?);
            if a == b {
                return Err(ReductionError::SelfPair { pair });
            }
            if !seen.insert((a.min(b), a.max(b))) {
                return Err(ReductionError::DuplicatePair { pair });
            }
            indices.push((a_index, b_index));
        }
        let count = pairs.len();
        Ok(Self {
            topology: Arc::clone(topology),
            pairs,
            indices,
            hits: vec![0; count],
            scratch: vec![false; count],
            cutoff,
            frames: 0,
        })
    }

    /// Observes one frame, with `frame_index` used only in error diagnostics.
    pub fn observe(
        &mut self,
        frame_index: usize,
        frame: TrajectoryFrameView<'_>,
    ) -> Result<(), ReductionError> {
        if !Arc::ptr_eq(&self.topology, &frame.shared_topology()) {
            return Err(ReductionError::TopologyMismatch { frame: frame_index });
        }
        let next = next_count(self.frames)?;
        let positions = frame.positions().values();
        for ((&(a, b), pair), hit) in self.indices.iter().zip(&self.pairs).zip(&mut self.scratch) {
            let delta = positions.value()[a.index()] - positions.value()[b.index()];
            let distance = delta.x.hypot(delta.y).hypot(delta.z);
            if !distance.is_finite() {
                return Err(ReductionError::NumericalFailure {
                    frame: frame_index,
                    atom: pair.0,
                });
            }
            *hit = distance <= *self.cutoff.value();
        }
        for (count, hit) in self.hits.iter_mut().zip(&self.scratch) {
            *count += usize::from(*hit);
        }
        self.frames = next;
        Ok(())
    }

    pub fn frame_count(&self) -> usize {
        self.frames
    }

    pub fn finish(self) -> Result<ContactOccupancyResult, ReductionError> {
        if self.frames == 0 {
            return Err(ReductionError::NoFrames);
        }
        Ok(ContactOccupancyResult {
            topology: self.topology,
            pairs: self.pairs,
            hits: self.hits,
            cutoff: self.cutoff,
            frames: self.frames,
        })
    }
}

/// Contact hit counts and occupancies associated with the original atom pairs.
#[derive(Debug, Clone)]
pub struct ContactOccupancyResult {
    topology: Arc<Topology>,
    pairs: Vec<(InstanceAtomId, InstanceAtomId)>,
    hits: Vec<usize>,
    cutoff: Quantity<f64>,
    frames: usize,
}

impl ContactOccupancyResult {
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    pub fn pairs(&self) -> &[(InstanceAtomId, InstanceAtomId)] {
        &self.pairs
    }
    pub fn hit_counts(&self) -> &[usize] {
        &self.hits
    }
    pub fn frame_count(&self) -> usize {
        self.frames
    }
    pub fn cutoff(&self) -> Quantity<f64> {
        self.cutoff
    }

    /// Returns (source atom pair, fraction of observed frames in contact).
    pub fn contacts(
        &self,
    ) -> impl ExactSizeIterator<Item = ((InstanceAtomId, InstanceAtomId), f64)> + '_ {
        self.pairs.iter().copied().zip(
            self.hits
                .iter()
                .map(|hits| *hits as f64 / self.frames as f64),
        )
    }
}

impl Trajectory {
    /// Per-atom population RMSF about mean stored positions. Align and perform
    /// periodic preprocessing beforehand. Uses [`RmsfAccumulator`] internally.
    pub fn rmsf(&self, selection: &AtomSelection) -> Result<RmsfResult, ReductionError> {
        selection
            .ensure_compatible(&self.shared_topology())
            .map_err(ReductionError::Selection)?;
        let mut accumulator = RmsfAccumulator::new(selection)?;
        for (index, frame) in self.frames().enumerate() {
            accumulator.observe(index, frame)?;
        }
        accumulator.finish()
    }

    /// Cartesian contact occupancies using the same implementation as streaming.
    pub fn contact_occupancy(
        &self,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
        cutoff: Quantity<f64>,
    ) -> Result<ContactOccupancyResult, ReductionError> {
        let mut accumulator =
            ContactOccupancyAccumulator::new(&self.shared_topology(), pairs, cutoff)?;
        for (index, frame) in self.frames().enumerate() {
            accumulator.observe(index, frame)?;
        }
        accumulator.finish()
    }
}
