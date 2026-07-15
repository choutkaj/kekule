use super::*;
use crate::units::{Quantity, ScaleValue, Unit, UnitError, ANGSTROM};

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

impl ScaleValue for Point3 {
    fn scaled(self, factor: f64) -> Self {
        Self::new(self.x * factor, self.y * factor, self.z * factor)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Conformer {
    pub(crate) positions: Vec<Option<Point3>>,
    unit: Unit,
    pub(crate) props: PropMap,
}

impl Conformer {
    /// Creates an empty conformer whose stored coordinates use `unit`.
    pub fn new(unit: Unit) -> std::result::Result<Self, UnitError> {
        Self::with_atom_capacity(0, unit)
    }

    /// Creates an empty conformer with space for `atom_capacity` positions.
    pub fn with_atom_capacity(
        atom_capacity: usize,
        unit: Unit,
    ) -> std::result::Result<Self, UnitError> {
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
    ) -> std::result::Result<(), UnitError> {
        let point = point.into_unit(self.unit)?.into_value();
        if self.positions.len() <= atom.index() {
            self.positions.resize(atom.index() + 1, None);
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
        self.positions
            .iter()
            .enumerate()
            .filter_map(|(index, point)| point.map(|point| (AtomId::new(index as u32), point)))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }
}
