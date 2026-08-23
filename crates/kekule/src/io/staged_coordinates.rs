use std::fmt;

use crate::chemistry::AtomPositionSource;
use crate::core::{AtomId, Molecule};
use crate::geometry::Point3;
use crate::structure::{PositionError, Positions};
use crate::units::{Quantity, Unit, UnitError, ANGSTROM};

/// Sparse atom-local coordinates used only while a format graph is staged.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct StagedCoordinates {
    positions: Vec<Option<Point3>>,
    unit: Unit,
}

impl StagedCoordinates {
    pub(crate) fn with_atom_capacity(
        atom_capacity: usize,
        unit: Unit,
    ) -> Result<Self, StagedCoordinateError> {
        unit.conversion_factor_to(ANGSTROM)?;
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(atom_capacity)
            .map_err(|_| StagedCoordinateError::CapacityOverflow)?;
        positions.resize(atom_capacity, None);
        Ok(Self { positions, unit })
    }

    pub(crate) const fn unit(&self) -> Unit {
        self.unit
    }

    pub(crate) fn set_position(
        &mut self,
        atom: AtomId,
        point: Quantity<Point3>,
    ) -> Result<(), StagedCoordinateError> {
        let point = point.to_unit(self.unit)?.to_value();
        if !point.is_finite() {
            return Err(StagedCoordinateError::NonFinitePosition { atom });
        }
        let required_len = atom
            .index()
            .checked_add(1)
            .ok_or(StagedCoordinateError::CapacityOverflow)?;
        if self.positions.len() < required_len {
            self.positions
                .try_reserve(required_len - self.positions.len())
                .map_err(|_| StagedCoordinateError::CapacityOverflow)?;
            self.positions.resize(required_len, None);
        }
        self.positions[atom.index()] = Some(point);
        Ok(())
    }

    pub(crate) fn position(&self, atom: AtomId) -> Option<Quantity<Point3>> {
        self.position_value(atom)
            .map(|point| Quantity::new(point, self.unit))
    }

    pub(crate) fn to_positions(
        &self,
        molecule: &Molecule,
    ) -> Result<Positions, StagedCoordinateError> {
        let values = molecule
            .atom_ids()
            .map(|atom| {
                self.position_value(atom)
                    .ok_or(StagedCoordinateError::MissingPosition { atom })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Positions::new(Quantity::new(values, self.unit))?)
    }
}

impl AtomPositionSource for StagedCoordinates {
    fn position_value(&self, atom: AtomId) -> Option<Point3> {
        self.positions.get(atom.index()).copied().flatten()
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum StagedCoordinateError {
    CapacityOverflow,
    MissingPosition { atom: AtomId },
    NonFinitePosition { atom: AtomId },
    Unit(UnitError),
    Position(PositionError),
}

impl fmt::Display for StagedCoordinateError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CapacityOverflow => formatter.write_str("source coordinate capacity exceeded"),
            Self::MissingPosition { atom } => {
                write!(formatter, "missing source position for {atom}")
            }
            Self::NonFinitePosition { atom } => {
                write!(formatter, "source position for {atom} is not finite")
            }
            Self::Unit(error) => error.fmt(formatter),
            Self::Position(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for StagedCoordinateError {}

impl From<UnitError> for StagedCoordinateError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

impl From<PositionError> for StagedCoordinateError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}
