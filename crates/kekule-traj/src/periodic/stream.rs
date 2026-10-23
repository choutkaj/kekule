use std::fmt;
use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::structure::{AsModelView, ModelView, Positions, TrajectoryError, TrajectoryFrameView};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex};

use super::{checked_image, positions, Lattice, MoleculePlan, PeriodicError};
use crate::FrameBuffer;

/// A reusable bond traversal for making molecules whole and imaging them.
///
/// The plan binds one topology layout. Operations are independent between frames;
/// `frame_index` is the caller's zero-based source index for error reporting.
/// All non-position state is retained. These are the same operations used by
/// [`super::make_molecules_whole`], [`super::image_molecules`], and
/// [`super::image_onto`].
pub struct MoleculeImager {
    topology: Arc<Topology>,
    plan: MoleculePlan,
}

impl fmt::Debug for MoleculeImager {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MoleculeImager")
            .field("atom_count", &self.topology.atom_count())
            .finish_non_exhaustive()
    }
}

impl MoleculeImager {
    /// Plans the bonded molecules of `topology` once for reuse across frames.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = topology.into();
        let plan = MoleculePlan::new(&topology);
        Self { topology, plan }
    }

    /// The topology whose layout every frame must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Makes every bonded molecule of the buffered frame whole.
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::TopologyMismatch`] when the frame does not
    /// share the imager's layout, and otherwise the per-frame errors of
    /// [`super::make_molecules_whole`]. On error the buffer is unchanged.
    pub fn make_whole(
        &self,
        frame_index: usize,
        frame: &mut FrameBuffer,
    ) -> Result<(), PeriodicError> {
        let positions = self.frame_positions(frame_index, frame.as_model_view(), None)?;
        publish_buffer(frame, &positions)
    }

    /// Images whole molecules of the buffered frame around the selected anchors.
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::SelectionTopologyMismatch`],
    /// [`PeriodicError::EmptyAnchors`] or [`PeriodicError::TopologyMismatch`]
    /// for mismatched inputs, and otherwise the per-frame errors of
    /// [`super::image_molecules`]. On error the buffer is unchanged.
    pub fn image(
        &self,
        frame_index: usize,
        frame: &mut FrameBuffer,
        anchors: &AtomSelection,
    ) -> Result<(), PeriodicError> {
        let anchors = self.anchor_groups(anchors)?;
        let positions = self.frame_positions(frame_index, frame.as_model_view(), Some(&anchors))?;
        publish_buffer(frame, &positions)
    }

    /// Moves each molecule of the buffered frame to its periodic image nearest
    /// the same molecule in `reference`, as [`super::image_onto`] does.
    ///
    /// The reference is any borrowed structure of the same layout, such as a
    /// trajectory frame or another buffer's [`FrameBuffer::frame_view`], from
    /// before or after this frame: frames can be imaged onto each other in any
    /// order. To image a window of frames onto frame `i`, image each of them
    /// onto a view of frame `i`.
    ///
    /// ```
    /// use std::sync::Arc;
    ///
    /// use kekule::geometry::{PeriodicCell, Point3, Vector3};
    /// use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
    /// use kekule::units::{Quantity, NANOMETER};
    /// use kekule::{smiles, topology::Topology};
    /// use kekule_traj::{periodic::MoleculeImager, FrameBuffer};
    ///
    /// let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
    /// let topology = Arc::new(Topology::from_molecule(molecule)?);
    /// let lengths = Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER);
    /// let cell = PeriodicCell::orthorhombic(lengths, [true; 3])?;
    /// // Between two frames, a helium atom crosses the x boundary of a 1 nm cell.
    /// let frames = [0.95, 0.05].map(|x| {
    ///     let points = vec![Point3::new(x, 0.5, 0.5)];
    ///     let positions = Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap();
    ///     TrajectoryFrame::new(positions).with_cell(cell)
    /// });
    /// let trajectory = Trajectory::from_items(Arc::clone(&topology), frames)?;
    ///
    /// let imager = MoleculeImager::new(Arc::clone(&topology));
    /// let mut frame = FrameBuffer::new(topology);
    /// frame.copy_from(trajectory.get(1).unwrap())?;
    /// imager.image_onto(1, &mut frame, &trajectory.get(0).unwrap())?;
    /// let x = frame.positions().get(0).unwrap().value_in(NANOMETER)?.x;
    /// assert!((x - 1.05).abs() < 1.0e-12);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::TopologyMismatch`] when the frame does not
    /// share the imager's layout, and otherwise the errors of
    /// [`super::image_onto`]. On error the buffer is unchanged.
    pub fn image_onto(
        &self,
        frame_index: usize,
        frame: &mut FrameBuffer,
        reference: &(impl AsModelView + ?Sized),
    ) -> Result<(), PeriodicError> {
        let positions = self.reference_positions(
            frame_index,
            frame.as_model_view(),
            reference.as_model_view(),
        )?;
        publish_buffer(frame, &positions)
    }

    pub(super) fn reference_positions(
        &self,
        frame_index: usize,
        source: ModelView<'_>,
        reference: ModelView<'_>,
    ) -> Result<Positions, PeriodicError> {
        if !self.topology().shares_layout(reference.topology()) {
            return Err(PeriodicError::ReferenceTopologyMismatch);
        }
        if !self.topology().shares_layout(source.topology()) {
            return Err(PeriodicError::TopologyMismatch { frame: frame_index });
        }
        let lattice = Lattice::new(source.cell().copied(), frame_index)?;
        let source = source.positions().values();
        let reference = reference.positions().values();
        let mut points = source.value().to_vec();
        self.plan
            .image_onto(&mut points, reference.value(), &lattice)?;
        positions(points, frame_index)
    }

    pub(super) fn anchor_groups(
        &self,
        anchors: &AtomSelection,
    ) -> Result<Vec<bool>, PeriodicError> {
        if !self.topology().shares_layout(anchors.topology()) {
            return Err(PeriodicError::SelectionTopologyMismatch);
        }
        if anchors.indices().is_empty() {
            return Err(PeriodicError::EmptyAnchors);
        }
        let mut selected = vec![false; self.plan.groups.len()];
        for index in anchors.indices() {
            selected[self.plan.atom_group[index.index()]] = true;
        }
        Ok(selected)
    }

    pub(super) fn frame_positions(
        &self,
        frame_index: usize,
        source: ModelView<'_>,
        anchors: Option<&[bool]>,
    ) -> Result<Positions, PeriodicError> {
        if !self.topology().shares_layout(source.topology()) {
            return Err(PeriodicError::TopologyMismatch { frame: frame_index });
        }
        let lattice = Lattice::new(source.cell().copied(), frame_index)?;
        let mut points = self
            .plan
            .make_whole(source.positions().values().value(), &lattice)?;
        if let Some(anchors) = anchors {
            self.plan.image(&mut points, &lattice, anchors)?;
        }
        positions(points, frame_index)
    }
}

/// Sequential fractional-coordinate unwrapping with memory proportional to atom count.
///
/// Retain this object across buffers and processing chunks. The first source index
/// may be any value; subsequent indices must be consecutive. Failed operations
/// change neither the caller's buffer nor the unwrapping state, so a corrected
/// frame can be retried. Use [`Self::reset`] for a new independent sequence.
/// Available times must not decrease; missing times are allowed. This follows the
/// scientific convention and sampling assumptions of [`super::unwrap_trajectory`].
pub struct TrajectoryUnwrapper {
    topology: Arc<Topology>,
    previous: Option<PreviousFrame>,
}

impl fmt::Debug for TrajectoryUnwrapper {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TrajectoryUnwrapper")
            .field("atom_count", &self.topology.atom_count())
            .field(
                "previous_frame",
                &self.previous.as_ref().map(|previous| previous.index),
            )
            .finish_non_exhaustive()
    }
}

struct PreviousFrame {
    index: usize,
    fractional: Vec<[f64; 3]>,
    periodic: [bool; 3],
    last_time: Option<f64>,
}

struct UnwrapStep {
    positions: Positions,
    previous: PreviousFrame,
}

impl TrajectoryUnwrapper {
    /// An unwrapper for frames of `topology` that has seen no frame yet.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            previous: None,
        }
    }

    /// The topology whose layout every frame must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    /// The source index of the last successfully unwrapped frame; `None`
    /// before the first frame or after [`Self::reset`].
    pub fn last_frame_index(&self) -> Option<usize> {
        self.previous.as_ref().map(|frame| frame.index)
    }
    /// Forgets the previous frame, so the next frame starts a new sequence
    /// and keeps its coordinates.
    pub fn reset(&mut self) {
        self.previous = None;
    }

    /// Transactionally unwraps a reusable buffer and advances the sequence
    /// only on success.
    ///
    /// # Errors
    ///
    /// Returns [`PeriodicError::TopologyMismatch`] when the frame does not
    /// share the unwrapper's layout, [`PeriodicError::NonSequentialFrame`]
    /// when `frame_index` does not follow [`Self::last_frame_index`], and
    /// otherwise the per-frame errors of [`super::unwrap_trajectory`].
    pub fn unwrap(
        &mut self,
        frame_index: usize,
        frame: &mut FrameBuffer,
    ) -> Result<(), PeriodicError> {
        let step = self.prepare(frame_index, frame.frame_view())?;
        publish_buffer(frame, &step.positions)?;
        self.commit(step);
        Ok(())
    }

    pub(super) fn next_positions(
        &mut self,
        frame_index: usize,
        source: TrajectoryFrameView<'_>,
    ) -> Result<Positions, PeriodicError> {
        let step = self.prepare(frame_index, source)?;
        Ok(self.commit(step))
    }

    fn commit(&mut self, step: UnwrapStep) -> Positions {
        self.previous = Some(step.previous);
        step.positions
    }

    fn prepare(
        &self,
        frame_index: usize,
        source: TrajectoryFrameView<'_>,
    ) -> Result<UnwrapStep, PeriodicError> {
        if !self.topology().shares_layout(source.topology()) {
            return Err(PeriodicError::TopologyMismatch { frame: frame_index });
        }
        if let Some(previous) = &self.previous {
            if previous.index.checked_add(1) != Some(frame_index) {
                return Err(PeriodicError::NonSequentialFrame {
                    previous: previous.index,
                    frame: frame_index,
                });
            }
        }
        let lattice = Lattice::new(source.cell().copied(), frame_index)?;
        let time = source.time().map(|time| time.into_value());
        let mut last_time = time;
        if let Some(previous) = &self.previous {
            if previous.periodic != lattice.periodic {
                return Err(PeriodicError::PeriodicAxesChanged { frame: frame_index });
            }
            if let (Some(before), Some(now)) = (previous.last_time, time) {
                if now < before {
                    return Err(PeriodicError::NonMonotonicTime { frame: frame_index });
                }
            }
            last_time = time.or(previous.last_time);
        }
        let mut fractional = Vec::with_capacity(self.topology.atom_count());
        let mut points = Vec::with_capacity(self.topology.atom_count());
        for (atom, point) in source
            .positions()
            .values()
            .value()
            .iter()
            .copied()
            .enumerate()
        {
            let mut current = lattice.fractional(point - Point3::origin())?;
            let mut images = [0.0; 3];
            if let Some(previous) = &self.previous {
                for axis in 0..3 {
                    if lattice.periodic[axis] {
                        let delta = current[axis] - previous.fractional[atom][axis];
                        let image = checked_image(delta, frame_index)?;
                        if ((delta - image).abs() - 0.5).abs() <= 32.0 * f64::EPSILON {
                            return Err(PeriodicError::AmbiguousDisplacement {
                                frame: frame_index,
                                atom: TopologyAtomIndex::new(atom as u32),
                                axis,
                            });
                        }
                        images[axis] = image;
                        current[axis] -= image;
                    }
                }
            }
            // Preserve the first frame without an inverse/forward round trip.
            points.push(point - lattice.cartesian(images));
            fractional.push(current);
        }
        Ok(UnwrapStep {
            positions: positions(points, frame_index)?,
            previous: PreviousFrame {
                index: frame_index,
                fractional,
                periodic: lattice.periodic,
                last_time,
            },
        })
    }
}

fn publish_buffer(frame: &mut FrameBuffer, positions: &Positions) -> Result<(), PeriodicError> {
    frame
        .frame_mut()
        .set_positions(positions.values())
        .map_err(|error| PeriodicError::Publication(Box::new(TrajectoryError::Structure(error))))
}
