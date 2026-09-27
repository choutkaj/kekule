use crate::geometry::{PeriodicGeometry, PeriodicGeometryError, Point3};
use std::collections::BTreeMap;

/// Fractional bins sized from reciprocal-vector norms. Every displacement within
/// cutoff changes each wrapped fractional coordinate by at most one bin. Exact
/// Cartesian nearest-image checks are still necessary for skewed cells.
pub(super) struct Neighbors<'a> {
    geometry: &'a PeriodicGeometry,
    bins: [usize; 3],
    cells: BTreeMap<[usize; 3], Vec<usize>>,
    points: Vec<Point3>,
}
impl<'a> Neighbors<'a> {
    pub(super) fn new(
        geometry: &'a PeriodicGeometry,
        cutoff: f64,
    ) -> Result<Self, PeriodicGeometryError> {
        let spans = geometry.reciprocal_norms().map(|r| r * cutoff);
        if cutoff <= 0.0 || spans.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(PeriodicGeometryError::NumericalFailure);
        }
        let bins = spans.map(|span| {
            (1.0 / (span * (1.0 + 128.0 * f64::EPSILON)))
                .floor()
                .clamp(1.0, 1_000_000.0) as usize
        });
        Ok(Self {
            geometry,
            bins,
            cells: BTreeMap::new(),
            points: Vec::new(),
        })
    }
    fn cell(&self, point: Point3) -> Result<[usize; 3], PeriodicGeometryError> {
        let fractional = self.geometry.fractional(point - Point3::origin())?;
        if fractional.iter().any(|v| v.abs() >= 2.0_f64.powi(48)) {
            return Err(PeriodicGeometryError::NumericalFailure);
        }
        Ok(std::array::from_fn(|i| {
            ((fractional[i].rem_euclid(1.0) * self.bins[i] as f64).floor() as usize)
                .min(self.bins[i] - 1)
        }))
    }
    pub(super) fn insert(&mut self, point: Point3) -> Result<(), PeriodicGeometryError> {
        let key = self.cell(point)?;
        self.cells.entry(key).or_default().push(self.points.len());
        self.points.push(point);
        Ok(())
    }
    /// Predicate gets insertion index and exact periodic distance in nm.
    pub(super) fn any(
        &self,
        point: Point3,
        mut predicate: impl FnMut(usize, f64) -> bool,
    ) -> Result<bool, PeriodicGeometryError> {
        let key = self.cell(point)?;
        let axes: [[usize; 3]; 3] = std::array::from_fn(|i| {
            [
                (key[i] + self.bins[i] - 1) % self.bins[i],
                key[i],
                (key[i] + 1) % self.bins[i],
            ]
        });
        for a in 0..3 {
            if axes[0][..a].contains(&axes[0][a]) {
                continue;
            }
            for b in 0..3 {
                if axes[1][..b].contains(&axes[1][b]) {
                    continue;
                }
                for c in 0..3 {
                    if axes[2][..c].contains(&axes[2][c]) {
                        continue;
                    }
                    if let Some(indices) = self.cells.get(&[axes[0][a], axes[1][b], axes[2][c]]) {
                        for &index in indices {
                            let distance = self
                                .geometry
                                .minimum_image(point - self.points[index])?
                                .norm();
                            if !distance.is_finite() {
                                return Err(PeriodicGeometryError::NumericalFailure);
                            }
                            if predicate(index, distance) {
                                return Ok(true);
                            }
                        }
                    }
                }
            }
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        geometry::{PeriodicCell, Vector3},
        units::{Quantity, NANOMETER},
    };
    #[test]
    fn periodic_bins_match_exhaustive_search() {
        for basis in [
            [
                Vector3::new(2.0, 0.0, 0.0),
                Vector3::new(0.0, 3.0, 0.0),
                Vector3::new(0.0, 0.0, 4.0),
            ],
            [
                Vector3::new(2.0, 0.2, 0.1),
                Vector3::new(1.7, 0.8, 0.2),
                Vector3::new(0.2, 0.3, 2.0),
            ],
            [
                Vector3::new(0.0, 2.0, 0.0),
                Vector3::new(3.0, 0.0, 0.0),
                Vector3::new(0.0, 0.0, 4.0),
            ],
        ] {
            let geometry = PeriodicGeometry::new(
                PeriodicCell::new(Quantity::new(basis, NANOMETER), [true; 3]).unwrap(),
            )
            .unwrap();
            for cutoff in [0.05, 0.4, 2.5] {
                let mut search = Neighbors::new(&geometry, cutoff).unwrap();
                let points: Vec<_> = (0..80)
                    .map(|i| {
                        Point3::origin()
                            + geometry.cartesian([
                                (i * 7 % 79) as f64 / 79.0,
                                (i * 31 % 79) as f64 / 79.0,
                                (i * 13 % 79) as f64 / 79.0,
                            ])
                    })
                    .collect();
                for &point in &points {
                    search.insert(point).unwrap();
                }
                for i in 0..95 {
                    let query = Point3::origin()
                        + geometry.cartesian([
                            i as f64 / 95.0 - 1.0,
                            (i * 17 % 95) as f64 / 95.0,
                            1.999999,
                        ]);
                    let mut observed = vec![];
                    search
                        .any(query, |j, d| {
                            if d < cutoff {
                                observed.push(j);
                            }
                            false
                        })
                        .unwrap();
                    observed.sort_unstable();
                    let expected: Vec<_> = points
                        .iter()
                        .enumerate()
                        .filter_map(|(j, p)| {
                            (geometry.minimum_image(query - *p).unwrap().norm() < cutoff)
                                .then_some(j)
                        })
                        .collect();
                    assert_eq!(observed, expected, "cutoff={cutoff}, query={query:?}");
                }
            }
        }
    }
}
