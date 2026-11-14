use crate::math::{Lattice, PeriodicCellError, Point3, PointGrid};

/// Points in a periodic grid of the lattice, measured by exact Cartesian
/// nearest-image distance, which skewed cells need.
pub(super) struct Neighbors<'a> {
    geometry: &'a Lattice,
    grid: PointGrid<usize>,
}
impl<'a> Neighbors<'a> {
    pub(super) fn new(geometry: &'a Lattice, cutoff: f64) -> Result<Self, PeriodicCellError> {
        Ok(Self {
            geometry,
            grid: PointGrid::periodic(geometry, cutoff)?,
        })
    }
    pub(super) fn insert(&mut self, point: Point3) -> Result<(), PeriodicCellError> {
        self.grid.insert(point, self.grid.len())
    }
    /// Predicate gets insertion index and exact periodic distance in nm.
    pub(super) fn any(
        &self,
        point: Point3,
        mut predicate: impl FnMut(usize, f64) -> bool,
    ) -> Result<bool, PeriodicCellError> {
        self.grid.try_any_near(&point, |other, index| {
            let distance = self
                .geometry
                .minimum_image_canonical(point - *other)?
                .norm();
            if !distance.is_finite() {
                return Err(PeriodicCellError::NumericalFailure);
            }
            Ok(predicate(index, distance))
        })
    }
}

#[cfg(test)]
mod tests;
