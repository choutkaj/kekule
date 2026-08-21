use super::*;
use crate::geometry::Point3;
use crate::units::{Quantity, Unit, UnitError, ANGSTROM};
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Conformer {
    pub(crate) positions: Vec<Option<Point3>>,
    unit: Unit,
    pub(crate) props: PropMap,
}

impl Conformer {
    /// Creates an empty conformer whose stored coordinates use `unit`.
    pub fn new(unit: Unit) -> std::result::Result<Self, ConformerError> {
        Self::with_atom_capacity(0, unit)
    }

    /// Creates an empty conformer with space for `atom_capacity` positions.
    pub fn with_atom_capacity(
        atom_capacity: usize,
        unit: Unit,
    ) -> std::result::Result<Self, ConformerError> {
        checked_fixed_id_collection_len(0, atom_capacity)
            .map_err(|_| ConformerError::PositionCapacityExceeded)?;
        unit.conversion_factor_to(ANGSTROM)?;
        Ok(Self {
            positions: vec![None; atom_capacity],
            unit,
            props: PropMap::new(),
        })
    }

    /// Returns the unit used by the stored coordinate array.
    pub const fn unit(&self) -> Unit {
        self.unit
    }

    /// Stores a position, converting it to the conformer's coordinate unit.
    pub fn set_position(
        &mut self,
        atom: AtomId,
        point: Quantity<Point3>,
    ) -> std::result::Result<(), ConformerError> {
        let required_len = atom
            .index()
            .checked_add(1)
            .ok_or(ConformerError::PositionCapacityExceeded)?;
        checked_fixed_id_collection_len(0, required_len)
            .map_err(|_| ConformerError::PositionCapacityExceeded)?;
        let point = point.to_unit(self.unit)?.to_value();
        if self.positions.len() <= atom.index() {
            self.positions.resize(required_len, None);
        }
        self.positions[atom.index()] = Some(point);
        Ok(())
    }

    pub fn clear_position(&mut self, atom: AtomId) {
        if let Some(position) = self.positions.get_mut(atom.index()) {
            *position = None;
        }
    }

    pub fn position(&self, atom: AtomId) -> Option<Quantity<Point3>> {
        self.position_value(atom)
            .map(|point| Quantity::new(point, self.unit))
    }

    pub fn positions(&self) -> impl Iterator<Item = (AtomId, Quantity<Point3>)> + '_ {
        self.positions_values()
            .map(|(atom, point)| (atom, Quantity::new(point, self.unit)))
    }

    pub(crate) fn position_value(&self, atom: AtomId) -> Option<Point3> {
        self.positions.get(atom.index()).copied().flatten()
    }

    pub(crate) fn positions_values(&self) -> impl Iterator<Item = (AtomId, Point3)> + '_ {
        (0..=u32::MAX)
            .zip(self.positions.iter())
            .filter_map(|(raw, point)| point.map(|point| (AtomId::new(raw), point)))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }
}

/// Errors from constructing or updating molecule-local conformer coordinates.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConformerError {
    Unit(UnitError),
    PositionCapacityExceeded,
}

impl fmt::Display for ConformerError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unit(error) => error.fmt(formatter),
            Self::PositionCapacityExceeded => {
                formatter.write_str("conformer position identifier capacity exceeded")
            }
        }
    }
}

impl std::error::Error for ConformerError {}

impl From<UnitError> for ConformerError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[cfg(all(test, target_pointer_width = "64"))]
mod capacity_tests {
    use super::*;

    #[test]
    fn conformer_capacity_is_rejected_before_allocation() {
        let first_unsupported_len = usize::try_from(u64::from(u32::MAX) + 2).expect("64-bit usize");
        assert_eq!(
            Conformer::with_atom_capacity(first_unsupported_len, ANGSTROM),
            Err(ConformerError::PositionCapacityExceeded)
        );
    }
}
