//! A codec's decoded frame state as a [`Geometry`].

use crate::diagnostics::{Countable, LengthMismatch};
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::properties::{GeometryProperties, OwnerProperties};
use crate::structure::{DenseArrayError, Forces, GeometryError, Positions, Velocities};
use crate::topology::Topology;
use crate::units::Quantity;

use crate::structure::{Geometry, GeometryData};

/// The state of a frame a codec decoded, in the units it decoded them in:
/// positions, cell, velocities, forces, and frame properties.
pub(crate) struct DecodedGeometry {
    pub(crate) positions: Quantity<Vec<Point3>>,
    pub(crate) cell: Option<PeriodicCell>,
    pub(crate) velocities: Option<Quantity<Vec<Vector3>>>,
    pub(crate) forces: Option<Quantity<Vec<Vector3>>>,
    pub(crate) owner: Option<OwnerProperties>,
}

impl Geometry {
    /// A decoded frame of `topology`'s atoms as a geometry, keeping the
    /// decoded arrays and converting them to canonical units in place.
    pub(crate) fn from_decoded(
        frame: DecodedGeometry,
        topology: &Topology,
    ) -> Result<Self, GeometryError> {
        let atoms = topology.atom_count();
        let count = |actual: usize| {
            if actual == atoms {
                Ok(())
            } else {
                Err(GeometryError::Array(Box::new(
                    DenseArrayError::CountMismatch(LengthMismatch::new(
                        Countable::Row,
                        atoms,
                        actual,
                    )),
                )))
            }
        };
        count(frame.positions.value().len())?;
        let positions = Positions::from_vec(frame.positions)?;
        let velocities = frame
            .velocities
            .map(|values| {
                count(values.value().len())?;
                Ok::<_, GeometryError>(Velocities::from_vec(values)?)
            })
            .transpose()?;
        let forces = frame
            .forces
            .map(|values| {
                count(values.value().len())?;
                Ok::<_, GeometryError>(Forces::from_vec(values)?)
            })
            .transpose()?;
        let mut properties = GeometryProperties::new(atoms);
        properties.bind_bond_rows(topology.bond_count());
        if let Some(owner) = frame.owner {
            *properties.owner_mut() = owner;
        }
        let mut data = GeometryData::from_frame_state(positions, frame.cell, properties);
        data.velocities = velocities;
        data.forces = forces;
        Ok(Self::from_data(data))
    }
}
