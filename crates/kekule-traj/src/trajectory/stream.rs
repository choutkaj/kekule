use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::structure::Positions;
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::Quantity;

use super::buffer::FrameBuffer;
use super::collection::Trajectory;
use super::frame::{Forces, TrajectoryFrame, TrajectoryFrameView, Velocities};
use super::{validate_atom_count, TrajectoryError};

/// Sequential reader that publishes complete frames into reusable storage.
///
/// `read_next` returns `Ok(false)` only at clean end-of-stream. Implementations
/// validate a frame completely before replacing the destination buffer.
pub trait TrajectoryReader {
    fn topology(&self) -> &Topology;

    fn shared_topology(&self) -> Arc<Topology>;

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError>;
}

/// Trajectory reader with random frame access.
///
/// Implementations may build an index eagerly; callers should inspect
/// format-specific metadata when index construction cost matters.
pub trait SeekableTrajectoryReader: TrajectoryReader {
    fn frame_count(&self) -> Option<u64>;

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError>;
}

/// Streaming writer for topology-bound frame views.
///
/// Writers reject topology mismatches and unsupported frame state rather than
/// silently dropping data.
pub trait TrajectoryWriter {
    fn topology(&self) -> &Topology;

    fn shared_topology(&self) -> Arc<Topology>;

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError>;
}

/// Sequential and seekable reader over an in-memory trajectory.
pub struct MemoryTrajectoryReader<'a> {
    trajectory: &'a Trajectory,
    cursor: usize,
}

impl<'a> MemoryTrajectoryReader<'a> {
    pub fn new(trajectory: &'a Trajectory) -> Self {
        Self {
            trajectory,
            cursor: 0,
        }
    }
}

impl TrajectoryReader for MemoryTrajectoryReader<'_> {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        let Some(frame) = self.trajectory.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.copy_from(frame.view(&self.trajectory.topology)?)?;
        self.cursor += 1;
        Ok(true)
    }
}

impl SeekableTrajectoryReader for MemoryTrajectoryReader<'_> {
    fn frame_count(&self) -> Option<u64> {
        u64::try_from(self.trajectory.len()).ok()
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        let index =
            usize::try_from(index).map_err(|_| TrajectoryError::FrameIndexOutOfRange(index))?;
        let frame = self
            .trajectory
            .frames
            .get(index)
            .ok_or(TrajectoryError::FrameIndexOutOfRange(index as u64))?;
        destination.copy_from(frame.view(&self.trajectory.topology)?)?;
        self.cursor = index.saturating_add(1);
        Ok(())
    }
}

/// In-memory writer that preserves every core frame field.
pub struct MemoryTrajectoryWriter {
    trajectory: Trajectory,
}

impl MemoryTrajectoryWriter {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            trajectory: Trajectory::new(topology),
        }
    }

    pub fn to_trajectory(self) -> Trajectory {
        self.trajectory
    }
}

impl TrajectoryWriter for MemoryTrajectoryWriter {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !Arc::ptr_eq(&self.trajectory.topology, frame.topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let positions = Positions::new(frame.positions.values())?;
        let mut owned = TrajectoryFrame::new(positions, self.trajectory.topology.bond_count());
        owned.cell = frame.cell.copied();
        owned.properties = frame.properties.clone();
        owned.velocities = frame.velocities.map(Velocities::new).transpose()?;
        owned.forces = frame.forces.map(Forces::new).transpose()?;
        owned.time = frame.time;
        owned.step = frame.step;
        self.trajectory.push(owned)
    }
}

/// Proof that a coordinate-only source order exactly matches one topology.
#[derive(Debug, Clone)]
pub struct AtomOrderAssertion {
    topology: Arc<Topology>,
    kind: AtomOrderAssertionKind,
}

impl PartialEq for AtomOrderAssertion {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.kind == other.kind
    }
}

impl Eq for AtomOrderAssertion {}

/// Evidence represented by an [`AtomOrderAssertion`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum AtomOrderAssertionKind {
    SemanticOrder,
    DeclaredTopologyOrder,
}

impl AtomOrderAssertion {
    /// Proves that an explicit semantic atom sequence is the topology's exact
    /// authoritative dense order.
    pub fn from_semantic_order(
        topology: &Arc<Topology>,
        atom_order: &[InstanceAtomId],
    ) -> Result<Self, TrajectoryError> {
        if topology.atom_ids() != atom_order {
            return Err(TrajectoryError::AtomOrderMismatch);
        }
        Ok(Self {
            topology: Arc::clone(topology),
            kind: AtomOrderAssertionKind::SemanticOrder,
        })
    }

    /// Records the caller's explicit assertion that a topology-free file uses
    /// this topology's authoritative dense atom order.
    ///
    /// This is evidence supplied by the caller, not an inference from atom
    /// count. Format readers must still validate all stronger file metadata.
    pub fn assert_file_uses_topology_order(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            kind: AtomOrderAssertionKind::DeclaredTopologyOrder,
        }
    }

    pub fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
        Arc::ptr_eq(&self.topology, topology)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub const fn kind(&self) -> AtomOrderAssertionKind {
        self.kind
    }

    /// Backward-compatible spelling for [`Self::from_semantic_order`].
    pub fn new(
        topology: &Arc<Topology>,
        atom_order: &[InstanceAtomId],
    ) -> Result<Self, TrajectoryError> {
        Self::from_semantic_order(topology, atom_order)
    }
}

/// Reference reader for a topology-free coordinate source.
pub struct CoordinateFrameReader {
    topology: Arc<Topology>,
    frames: Vec<Vec<Point3>>,
    cursor: usize,
}

impl CoordinateFrameReader {
    pub fn new(
        topology: Arc<Topology>,
        assertion: AtomOrderAssertion,
        frames: impl IntoIterator<Item = Quantity<Vec<Point3>>>,
    ) -> Result<Self, TrajectoryError> {
        if !assertion.is_compatible(&topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let frames = frames
            .into_iter()
            .map(|frame| {
                let positions = Positions::new(frame).map_err(TrajectoryError::Position)?;
                validate_atom_count(topology.atom_count(), positions.len())?;
                Ok::<_, TrajectoryError>(positions.values().value().to_vec())
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology,
            frames,
            cursor: 0,
        })
    }
}

impl TrajectoryReader for CoordinateFrameReader {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        if !Arc::ptr_eq(&self.topology, &destination.topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let Some(frame) = self.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.set_positions(Quantity::new(
            frame.as_slice(),
            kekule::units::CANONICAL_LENGTH_UNIT,
        ))?;
        destination.reset_dynamic_state();
        self.cursor += 1;
        Ok(true)
    }
}
