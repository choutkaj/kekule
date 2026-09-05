use std::sync::Arc;

use kekule::properties::{Properties, PropertyError, PropertyKey, PropertyValue};
use kekule::topology::{AtomSelection, Topology, TopologyPerceptionError};

use super::frame::TrajectoryFrame;
use super::{TrajectoryError, TrajectoryFrameView, TrajectorySliceError};

/// A deliberately loaded finite in-memory trajectory.
///
/// All frames share one exact [`Topology`] allocation and are stored in stable
/// temporal order. [`Self::push`] validates a complete frame before mutation.
/// Use streaming readers and [`super::FrameBuffer`] for trajectories that
/// should not be loaded completely into memory.
#[derive(Debug, Clone)]
pub struct Trajectory {
    pub(super) topology: Arc<Topology>,
    properties: Properties,
    pub(super) frames: Vec<TrajectoryFrame>,
}

impl Trajectory {
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: Properties::new(),
            frames: Vec::new(),
        }
    }

    pub fn from_frames(
        topology: impl Into<Arc<Topology>>,
        frames: impl IntoIterator<Item = TrajectoryFrame>,
    ) -> Result<Self, TrajectoryError> {
        let mut trajectory = Self::new(topology);
        for frame in frames {
            trajectory.push(frame)?;
        }
        Ok(trajectory)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once for the collection, independent
    /// of frame count. All frame positions, cells, velocities, forces, time,
    /// step, and properties, and collection properties are retained without
    /// copying frames. Failure leaves the entire trajectory unchanged. Other
    /// owners, readers, buffers, selections, and prepared calculations retain
    /// their original topology bindings.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    pub const fn properties(&self) -> &Properties {
        &self.properties
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        self.properties.insert(key, value)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.properties.remove(key)
    }

    pub fn clear_properties(&mut self) {
        self.properties.clear_owner();
    }

    /// Constructs one topology subset and applies its dense mapping to every frame.
    pub fn slice(&self, selection: &AtomSelection) -> Result<Self, TrajectorySliceError> {
        let subset = self.topology.subset(selection)?;
        let atom_indices = subset
            .correspondence()
            .source_atom_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let bond_indices = subset
            .correspondence()
            .source_bond_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let mut target = Self::new(subset.shared_topology());
        for frame in &self.frames {
            target.push(TrajectoryFrame {
                positions: frame.positions.select_indices(&atom_indices)?,
                cell: frame.cell,
                properties: frame
                    .properties
                    .project_realization(&atom_indices, &bond_indices)?,
                velocities: frame
                    .velocities
                    .as_ref()
                    .map(|values| values.select_indices(&atom_indices))
                    .transpose()?,
                forces: frame
                    .forces
                    .as_ref()
                    .map(|values| values.select_indices(&atom_indices))
                    .transpose()?,
                time: frame.time,
                step: frame.step,
            })?;
        }
        Ok(target)
    }

    pub fn len(&self) -> usize {
        self.frames.len()
    }

    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Returns one topology-bound frame view by stable trajectory index.
    pub fn frame(&self, index: usize) -> Option<TrajectoryFrameView<'_>> {
        self.frames.get(index).map(|frame| {
            frame
                .view(&self.topology)
                .expect("trajectory validates frame topology on insertion")
        })
    }

    #[cfg(test)]
    pub(crate) fn frame_payload(&self, index: usize) -> Option<&TrajectoryFrame> {
        self.frames.get(index)
    }

    /// Iterates topology-bound frame views in stable trajectory order.
    pub fn frames(&self) -> impl ExactSizeIterator<Item = TrajectoryFrameView<'_>> {
        self.frames.iter().map(|frame| {
            frame
                .view(&self.topology)
                .expect("trajectory validates frame topology on insertion")
        })
    }

    pub fn push(&mut self, mut frame: TrajectoryFrame) -> Result<(), TrajectoryError> {
        frame.prepare(&self.topology)?;
        self.frames.push(frame);
        Ok(())
    }

    pub(crate) fn replace_frames(
        &mut self,
        mut frames: Vec<TrajectoryFrame>,
    ) -> Result<(), TrajectoryError> {
        for frame in &mut frames {
            frame.prepare(&self.topology)?;
        }
        self.frames = frames;
        Ok(())
    }

    /// Explicitly validates that all present times are non-decreasing.
    pub fn validate_monotonic_time(&self, require_all: bool) -> Result<(), TrajectoryError> {
        let mut previous = None;
        for (index, frame) in self.frames.iter().enumerate() {
            let Some(time) = frame.time else {
                if require_all {
                    return Err(TrajectoryError::MissingTime { frame: index });
                }
                continue;
            };
            let value = time.to_value();
            if previous.is_some_and(|previous| value < previous) {
                return Err(TrajectoryError::NonMonotonicTime { frame: index });
            }
            previous = Some(value);
        }
        Ok(())
    }
}
