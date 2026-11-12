//! A codec's decoded frame state as a [`Sample`].

use std::sync::Arc;

use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::properties::{OwnerProperties, StructureProperties};
use crate::structure::{DenseArrayError, Forces, Positions, StructureError, Velocities};
use crate::topology::Topology;
use crate::units::Quantity;

use super::{Sample, SampleStructure};

/// The state of a frame a codec decoded, in the units it decoded them in:
/// positions, cell, velocities, forces, and frame properties.
pub(crate) struct DecodedSample {
    pub(crate) positions: Quantity<Vec<Point3>>,
    pub(crate) cell: Option<PeriodicCell>,
    pub(crate) velocities: Option<Quantity<Vec<Vector3>>>,
    pub(crate) forces: Option<Quantity<Vec<Vector3>>>,
    pub(crate) owner: Option<OwnerProperties>,
}

impl Sample {
    /// A decoded frame of `topology`'s atoms as a sample, keeping the
    /// decoded arrays and converting them to canonical units in place.
    pub(crate) fn from_decoded(
        frame: DecodedSample,
        topology: &Topology,
    ) -> Result<Self, StructureError> {
        let atoms = topology.atom_count();
        let count = |actual: usize| {
            if actual == atoms {
                Ok(())
            } else {
                Err(StructureError::Array(DenseArrayError::CountMismatch {
                    expected: atoms,
                    actual,
                }))
            }
        };
        count(frame.positions.value().len())?;
        let positions = Positions::from_vec(frame.positions)?;
        let velocities = frame
            .velocities
            .map(|values| {
                count(values.value().len())?;
                Ok::<_, StructureError>(Velocities::from_vec(values)?)
            })
            .transpose()?;
        let forces = frame
            .forces
            .map(|values| {
                count(values.value().len())?;
                Ok::<_, StructureError>(Forces::from_vec(values)?)
            })
            .transpose()?;
        let mut properties = StructureProperties::new(atoms);
        properties.bind_bond_rows(topology.bond_count());
        if let Some(owner) = frame.owner {
            *properties.owner_mut() = owner;
        }
        let mut structure = SampleStructure::from_frame_state(positions, frame.cell, properties);
        structure.velocities = velocities;
        structure.forces = forces;
        Ok(Self::from_structure(Arc::new(structure)))
    }
}
