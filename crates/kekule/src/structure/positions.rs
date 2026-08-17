use std::fmt;
use std::sync::Arc;

use crate::geometry::Point3;
use crate::topology::{InstanceAtomId, Topology, TopologyAtomIndex, TopologyMapping};
use crate::units::{Quantity, UnitError, MODEL_LENGTH_UNIT};

use super::{remap, TopologyRemapError};

/// One complete finite Cartesian array in one topology's dense atom order.
#[derive(Debug, Clone)]
pub struct Positions {
    topology: Arc<Topology>,
    values: Vec<Point3>,
}

impl PartialEq for Positions {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.values == other.values
    }
}

impl Positions {
    pub fn new<T>(topology: &Arc<Topology>, positions: Quantity<T>) -> Result<Self, PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
        let source = positions.value().as_ref();
        validate_position_count(topology, source.len())?;
        let values = source
            .iter()
            .copied()
            .enumerate()
            .map(|(index, point)| {
                let point = Point3::new(point.x * factor, point.y * factor, point.z * factor);
                if !point.is_finite() {
                    return Err(PositionError::NonFinitePosition {
                        atom: topology.atom_ids()[index],
                    });
                }
                Ok(point)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: Arc::clone(topology),
            values,
        })
    }

    pub fn zeros(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            values: vec![Point3::origin(); topology.atom_count()],
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
        Arc::ptr_eq(&self.topology, topology)
    }

    pub(super) fn topology_arc(&self) -> &Arc<Topology> {
        &self.topology
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn values(&self) -> Quantity<&[Point3]> {
        Quantity::new(self.values.as_slice(), MODEL_LENGTH_UNIT)
    }

    /// Copies complete positions through a checked topology lineage mapping.
    pub fn remap_to(
        &self,
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TopologyRemapError> {
        if !self.is_compatible(source) {
            return Err(TopologyRemapError::SourceTopologyMismatch);
        }
        let values = remap::dense_atom_values(&self.values, source, target, mapping)?;
        Ok(Self {
            topology: Arc::clone(target),
            values,
        })
    }

    /// Copies complete positions through a checked topology lineage mapping
    /// while retaining this array's allocation.
    ///
    /// Validation completes before any destination position changes.
    pub fn copy_remapped_from(
        &mut self,
        source: &Self,
        source_topology: &Arc<Topology>,
        target_topology: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<(), TopologyRemapError> {
        if !source.is_compatible(source_topology) {
            return Err(TopologyRemapError::SourceTopologyMismatch);
        }
        if !self.is_compatible(target_topology) {
            return Err(TopologyRemapError::TargetTopologyMismatch);
        }
        if source.values.len() != source_topology.atom_count() {
            return Err(TopologyRemapError::SourceAtomCountMismatch {
                expected: source_topology.atom_count(),
                actual: source.values.len(),
            });
        }
        remap::validate_complete_atom_mapping(source_topology, target_topology, mapping)?;
        for (source_index, target_index) in mapping.atom_index_pairs() {
            self.values[target_index.index()] = source.values[source_index.index()];
        }
        Ok(())
    }

    pub fn position_at(&self, index: TopologyAtomIndex) -> Result<Quantity<Point3>, PositionError> {
        self.values
            .get(index.index())
            .copied()
            .map(|point| Quantity::new(point, MODEL_LENGTH_UNIT))
            .ok_or(PositionError::InvalidAtomIndex(index))
    }

    pub fn position(
        &self,
        topology: &Arc<Topology>,
        atom: InstanceAtomId,
    ) -> Result<Quantity<Point3>, PositionError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(PositionError::InvalidAtomId(atom))?;
        self.position_at(index)
    }

    pub fn set_position(
        &mut self,
        topology: &Arc<Topology>,
        atom: InstanceAtomId,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        self.ensure_compatible(topology)?;
        let index = topology
            .atom_index(atom)
            .ok_or(PositionError::InvalidAtomId(atom))?;
        let point = position.into_unit(MODEL_LENGTH_UNIT)?.into_value();
        if !point.is_finite() {
            return Err(PositionError::NonFinitePosition { atom });
        }
        self.values[index.index()] = point;
        Ok(())
    }

    /// Replaces all positions transactionally while reusing the current
    /// allocation when the capacity permits.
    pub fn set_all<T>(
        &mut self,
        topology: &Arc<Topology>,
        positions: Quantity<T>,
    ) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = self.validate_replacement(topology, &positions)?;
        self.copy_from_validated(positions.value().as_ref(), factor);
        Ok(())
    }

    /// Validates a complete replacement without changing this array.
    ///
    /// This supports external topology-bound containers that must validate
    /// several coupled fields before publishing any of them transactionally.
    pub fn validate_all<T>(
        &self,
        topology: &Arc<Topology>,
        positions: &Quantity<T>,
    ) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.validate_replacement(topology, positions).map(drop)
    }

    fn validate_replacement<T>(
        &self,
        topology: &Arc<Topology>,
        positions: &Quantity<T>,
    ) -> Result<f64, PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.ensure_compatible(topology)?;
        let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
        let source = positions.value().as_ref();
        validate_position_count(topology, source.len())?;
        for (index, point) in source.iter().copied().enumerate() {
            let converted = Point3::new(point.x * factor, point.y * factor, point.z * factor);
            if !converted.is_finite() {
                return Err(PositionError::NonFinitePosition {
                    atom: topology.atom_ids()[index],
                });
            }
        }
        Ok(factor)
    }

    fn copy_from_validated(&mut self, source: &[Point3], factor: f64) {
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Point3::new(source.x * factor, source.y * factor, source.z * factor);
        }
    }

    fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), PositionError> {
        if !self.is_compatible(topology) {
            return Err(PositionError::TopologyMismatch);
        }
        Ok(())
    }
}

fn validate_position_count(topology: &Topology, actual: usize) -> Result<(), PositionError> {
    if actual != topology.atom_count() {
        return Err(PositionError::PositionCountMismatch {
            expected: topology.atom_count(),
            actual,
        });
    }
    crate::core::checked_fixed_id_collection_len(0, actual)
        .map_err(|_| PositionError::CapacityOverflow)?;
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PositionError {
    TopologyMismatch,
    InvalidAtomId(InstanceAtomId),
    InvalidAtomIndex(TopologyAtomIndex),
    PositionCountMismatch { expected: usize, actual: usize },
    NonFinitePosition { atom: InstanceAtomId },
    CapacityOverflow,
    Unit(UnitError),
}

impl fmt::Display for PositionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("positions belong to a different topology")
            }
            Self::InvalidAtomId(atom) => write!(formatter, "invalid topology atom: {atom}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid {index}"),
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "topology requires {expected} positions, but received {actual}"
            ),
            Self::NonFinitePosition { atom } => {
                write!(formatter, "position for atom {atom} is not finite")
            }
            Self::CapacityOverflow => {
                formatter.write_str("position count exceeds fixed-width topology capacity")
            }
            Self::Unit(error) => write!(formatter, "invalid position unit: {error}"),
        }
    }
}

impl std::error::Error for PositionError {}

impl From<UnitError> for PositionError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
