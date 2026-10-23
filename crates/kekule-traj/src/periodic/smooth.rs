use std::fmt;
use std::ops::Range;
use std::sync::Arc;

use kekule::geometry::{Point3, Vector3};
use kekule::structure::{Positions, TrajectoryFrameView};
use kekule::topology::Topology;

use super::stream::publish_buffer;
use super::{positions, Lattice, PeriodicError};
use crate::FrameBuffer;

/// A moving average of positions over a window of frames around each frame.
///
/// Frame `i` is smoothed from itself and the neighbors `j` that the caller
/// supplies with `|j - i| <= half_width`, each counting once. Neighbors are
/// read, never changed, so the smoothed frame can be produced in any order from
/// frames already processed by earlier pipeline steps: playback, scrubbing and
/// sequential streaming all use the same call. Smoothing comes last, after
/// making molecules whole, imaging, and superposition.
///
/// # Window ends
///
/// The window is truncated, not shifted: near either end of a trajectory,
/// frame `i` averages only the frames of `i - half_width..=i + half_width`
/// that exist, so its window is asymmetric and smaller. [`Self::window`]
/// gives these indices. A frame with no supplied neighbor, such as every frame
/// with `half_width == 0`, keeps its positions.
///
/// # Periodic images
///
/// When frame `i` has a periodic cell, each neighbor atom is first placed at
/// its nearest Cartesian image to the same atom in frame `i` under frame
/// `i`'s cell, so an atom crossing the cell between frames is not averaged
/// across it. Every frame of the window must then have a cell; frame `i`'s
/// cell alone defines the images. Without a cell in frame `i`, no frame of the
/// window may have one, and positions are averaged as stored. Imaging is per
/// atom and assumes each atom moves less than half a cell within the window.
///
/// # What changes
///
/// Positions only. Frame `i` keeps its cell, time, step, and properties, and
/// also its own velocities and forces: these are not averaged, so they belong
/// to the unsmoothed frame. Drop them when that matters.
///
/// ```
/// use std::sync::Arc;
///
/// use kekule::geometry::{PeriodicCell, Point3, Vector3};
/// use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
/// use kekule::units::{Quantity, NANOMETER};
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::{periodic::FrameSmoother, FrameBuffer};
///
/// let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
/// let topology = Arc::new(Topology::from_molecule(molecule)?);
/// let lengths = Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER);
/// let cell = PeriodicCell::orthorhombic(lengths, [true; 3])?;
/// // A helium atom crosses the x boundary of a 1 nm cell.
/// let frames = [0.9, 0.0, 0.1].map(|x| {
///     let points = vec![Point3::new(x, 0.5, 0.5)];
///     let positions = Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap();
///     TrajectoryFrame::new(positions).with_cell(cell)
/// });
/// let trajectory = Trajectory::from_items(Arc::clone(&topology), frames)?;
///
/// let smoother = FrameSmoother::new(Arc::clone(&topology), 1);
/// let window = smoother.window(1, trajectory.len());
/// let mut frame = FrameBuffer::new(topology);
/// frame.copy_from(trajectory.get(1).unwrap())?;
/// let neighbors = window
///     .filter(|&index| index != 1)
///     .map(|index| (index, trajectory.get(index).unwrap()));
/// smoother.smooth(1, &mut frame, neighbors)?;
/// let x = frame.positions().get(0).unwrap().value_in(NANOMETER)?.x;
/// assert!(x.abs() < 1.0e-12); // Not 1/3, the average across the cell.
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub struct FrameSmoother {
    topology: Arc<Topology>,
    half_width: usize,
}

impl fmt::Debug for FrameSmoother {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("FrameSmoother")
            .field("atom_count", &self.topology.atom_count())
            .field("half_width", &self.half_width)
            .finish_non_exhaustive()
    }
}

impl FrameSmoother {
    /// A smoother for frames of `topology` over windows of up to
    /// `2 * half_width + 1` frames.
    pub fn new(topology: impl Into<Arc<Topology>>, half_width: usize) -> Self {
        Self {
            topology: topology.into(),
            half_width,
        }
    }

    /// The topology whose layout every frame must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The largest index distance between a frame and a neighbor.
    pub fn half_width(&self) -> usize {
        self.half_width
    }

    /// The indices of the truncated window around `frame_index` in a
    /// trajectory of `frame_count` frames, including `frame_index` itself.
    /// Empty when `frame_index` is not below `frame_count`.
    pub fn window(&self, frame_index: usize, frame_count: usize) -> Range<usize> {
        if frame_index >= frame_count {
            return frame_count..frame_count;
        }
        let start = frame_index.saturating_sub(self.half_width);
        let end = frame_index
            .saturating_add(self.half_width)
            .saturating_add(1)
            .min(frame_count);
        start..end
    }

    /// Replaces the buffered frame `frame_index` with the average of its
    /// positions and those of `neighbors`, given as `(index, frame)` pairs.
    ///
    /// Neighbor indices must be strictly increasing, lie within
    /// [`Self::half_width`] of `frame_index`, and exclude it. Indices also
    /// label failures.
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::TopologyMismatch`] for a frame that does not
    /// share the smoother's layout, [`PeriodicError::InvalidNeighbor`] for a
    /// neighbor index outside these rules, [`PeriodicError::MissingCell`]
    /// when the window mixes frames with and without cells,
    /// [`PeriodicError::LimitExceeded`] or [`PeriodicError::NumericalFailure`]
    /// for an ill-conditioned cell or out-of-range coordinates, and
    /// [`PeriodicError::Position`] or [`PeriodicError::Publication`] when the
    /// result cannot be stored. On error the buffer is unchanged.
    pub fn smooth<'a>(
        &self,
        frame_index: usize,
        frame: &mut FrameBuffer,
        neighbors: impl IntoIterator<Item = (usize, TrajectoryFrameView<'a>)>,
    ) -> Result<(), PeriodicError> {
        let positions = self.smoothed_positions(frame_index, frame.frame_view(), neighbors)?;
        publish_buffer(frame, &positions)
    }

    pub(super) fn smoothed_positions<'a>(
        &self,
        frame_index: usize,
        frame: TrajectoryFrameView<'_>,
        neighbors: impl IntoIterator<Item = (usize, TrajectoryFrameView<'a>)>,
    ) -> Result<Positions, PeriodicError> {
        if !self.topology().shares_layout(frame.topology()) {
            return Err(PeriodicError::TopologyMismatch { frame: frame_index });
        }
        let lattice = frame
            .cell()
            .map(|cell| Lattice::new(Some(*cell), frame_index))
            .transpose()?;
        let center = frame.positions().values();
        let center = center.value();
        // Sum displacements from frame `i` rather than positions, which keeps
        // precision far from the origin and makes imaging a per-term step.
        let mut sums = vec![Vector3::zero(); center.len()];
        let mut count = 1_usize;
        let mut previous = None;
        for (index, neighbor) in neighbors {
            let ordered = previous.is_none_or(|previous| index > previous);
            if !ordered || index == frame_index || index.abs_diff(frame_index) > self.half_width {
                return Err(PeriodicError::InvalidNeighbor {
                    frame: frame_index,
                    neighbor: index,
                });
            }
            previous = Some(index);
            if !self.topology().shares_layout(neighbor.topology()) {
                return Err(PeriodicError::TopologyMismatch { frame: index });
            }
            match (neighbor.cell().is_some(), lattice.is_some()) {
                (false, true) => return Err(PeriodicError::MissingCell { frame: index }),
                (true, false) => return Err(PeriodicError::MissingCell { frame: frame_index }),
                _ => {}
            }
            let points = neighbor.positions().values();
            let points = points.value();
            for ((sum, &at), &point) in sums.iter_mut().zip(center.iter()).zip(points.iter()) {
                let delta = point - at;
                *sum += match &lattice {
                    Some(lattice) => lattice.nearest_image(delta)?.0,
                    None => delta,
                };
            }
            count += 1;
        }
        let scale = 1.0 / count as f64;
        let mut points = Vec::with_capacity(center.len());
        for (&at, sum) in center.iter().zip(sums) {
            let point: Point3 = at + sum * scale;
            if !(point - Point3::origin()).is_finite() {
                return Err(PeriodicError::NumericalFailure { frame: frame_index });
            }
            points.push(point);
        }
        positions(points, frame_index)
    }
}
