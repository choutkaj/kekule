use std::sync::Arc;

use kekule::properties::{Properties, PropertyError, PropertyKey, PropertyValue};
use kekule::topology::{AtomSelection, Topology};

use super::frame::TrajectoryFrame;
use super::{TrajectoryError, TrajectoryFrameView, TrajectorySliceError};

/// Deliberately loaded finite in-memory trajectory.
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

    pub fn push(&mut self, frame: TrajectoryFrame) -> Result<(), TrajectoryError> {
        frame.validate(&self.topology)?;
        self.frames.push(frame);
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
