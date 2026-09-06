use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::structure::Positions;
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex};

use super::{checked_image, positions, Lattice, MoleculePlan, PeriodicError};
use crate::{FrameBuffer, TrajectoryFrame, TrajectoryFrameView};

/// A reusable bond traversal for making molecules whole and imaging them.
///
/// The plan binds one exact topology. Operations are independent between frames;
/// `frame_index` is the caller's zero-based source index for error reporting.
/// All non-position state is retained. These are the same operations used by
/// [`crate::Trajectory::make_molecules_whole`] and [`crate::Trajectory::image_molecules`].
pub struct MoleculeImager {
    topology: Arc<Topology>,
    plan: MoleculePlan,
}

impl MoleculeImager {
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = topology.into();
        let plan = MoleculePlan::new(&topology);
        Self { topology, plan }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn make_whole(
        &self,
        frame_index: usize,
        source: TrajectoryFrameView<'_>,
    ) -> Result<TrajectoryFrame, PeriodicError> {
        Ok(source.with_positions(self.frame_positions(frame_index, source, None)?))
    }

    pub fn make_whole_in_place(
        &self,
        frame_index: usize,
        frame: &mut FrameBuffer,
    ) -> Result<(), PeriodicError> {
        let positions = self.frame_positions(frame_index, frame.frame_view(), None)?;
        frame
            .set_positions(positions.values())
            .map_err(|e| PeriodicError::Publication(Box::new(e.into())))
    }

    pub fn image(
        &self,
        frame_index: usize,
        source: TrajectoryFrameView<'_>,
        anchors: &AtomSelection,
    ) -> Result<TrajectoryFrame, PeriodicError> {
        let anchors = self.anchor_groups(anchors)?;
        Ok(source.with_positions(self.frame_positions(frame_index, source, Some(&anchors))?))
    }

    pub fn image_in_place(
        &self,
        frame_index: usize,
        frame: &mut FrameBuffer,
        anchors: &AtomSelection,
    ) -> Result<(), PeriodicError> {
        let anchors = self.anchor_groups(anchors)?;
        let positions = self.frame_positions(frame_index, frame.frame_view(), Some(&anchors))?;
        frame
            .set_positions(positions.values())
            .map_err(|e| PeriodicError::Publication(Box::new(e.into())))
    }

    pub(super) fn anchor_groups(
        &self,
        anchors: &AtomSelection,
    ) -> Result<Vec<bool>, PeriodicError> {
        if !std::ptr::eq(self.topology(), anchors.topology()) {
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
        source: TrajectoryFrameView<'_>,
        anchors: Option<&[bool]>,
    ) -> Result<Positions, PeriodicError> {
        if !std::ptr::eq(self.topology(), source.topology()) {
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
/// scientific convention and sampling assumptions of [`crate::Trajectory::unwrap`].
pub struct TrajectoryUnwrapper {
    topology: Arc<Topology>,
    previous: Option<PreviousFrame>,
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
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            previous: None,
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    pub fn last_frame_index(&self) -> Option<usize> {
        self.previous.as_ref().map(|frame| frame.index)
    }
    pub fn reset(&mut self) {
        self.previous = None;
    }

    /// Returns the next unwrapped frame and advances state only on success.
    pub fn unwrap_frame(
        &mut self,
        frame_index: usize,
        source: TrajectoryFrameView<'_>,
    ) -> Result<TrajectoryFrame, PeriodicError> {
        let step = self.prepare(frame_index, source)?;
        Ok(source.with_positions(self.commit(step)))
    }

    /// Transactionally unwraps a reusable buffer and advances the sequence.
    pub fn unwrap_in_place(
        &mut self,
        frame_index: usize,
        frame: &mut FrameBuffer,
    ) -> Result<(), PeriodicError> {
        let step = self.prepare(frame_index, frame.frame_view())?;
        frame
            .set_positions(step.positions.values())
            .map_err(|e| PeriodicError::Publication(Box::new(e.into())))?;
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
        if !std::ptr::eq(self.topology(), source.topology()) {
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
        let time = source.time().map(|time| time.to_value());
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
