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
mod tests {
    use super::*;
    use crate::{
        math::{PeriodicCell, Vector3},
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
            let geometry = Lattice::new(
                PeriodicCell::new(Quantity::new(basis, NANOMETER), [true; 3]).unwrap(),
            )
            .unwrap();
            for cutoff in [0.05, 0.4, 2.5] {
                let mut search = Neighbors::new(&geometry, cutoff).unwrap();
                let points: Vec<_> = (0..80)
                    .map(|i| {
                        Point3::origin()
                            + geometry.cartesian_canonical([
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
                        + geometry.cartesian_canonical([
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
                            (geometry.minimum_image_canonical(query - *p).unwrap().norm() < cutoff)
                                .then_some(j)
                        })
                        .collect();
                    assert_eq!(observed, expected, "cutoff={cutoff}, query={query:?}");
                }
            }
        }
    }

    /// Packing's per-point cutoffs and incremental insertion, against an
    /// exhaustive search over seeded points in several cell images.
    #[test]
    fn incremental_packing_matches_exhaustive_search() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut unit = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let geometry = Lattice::new(
            PeriodicCell::new(
                Quantity::new(
                    [
                        Vector3::new(2.2, 0.0, 0.0),
                        Vector3::new(0.9, 2.0, 0.0),
                        Vector3::new(0.3, -0.4, 1.8),
                    ],
                    NANOMETER,
                ),
                [true; 3],
            )
            .unwrap(),
        )
        .unwrap();
        let cutoffs: Vec<f64> = (0..150).map(|_| 0.1 + unit() * 0.3).collect();
        let max_cutoff = cutoffs.iter().copied().fold(0.0, f64::max);
        let mut search = Neighbors::new(&geometry, max_cutoff).unwrap();
        let mut accepted: Vec<Point3> = Vec::new();
        for _ in 0..600 {
            let point = Point3::origin()
                + geometry.cartesian_canonical([
                    unit() * 3.0 - 1.0,
                    unit() * 2.0 - 0.5,
                    unit() * 3.0 - 2.0,
                ]);
            let observed = search
                .any(point, |index, distance| distance < cutoffs[index])
                .unwrap();
            let expected = accepted.iter().enumerate().any(|(index, other)| {
                geometry
                    .minimum_image_canonical(point - *other)
                    .unwrap()
                    .norm()
                    < cutoffs[index]
            });
            assert_eq!(observed, expected, "{point:?}");
            if !observed && accepted.len() < cutoffs.len() {
                search.insert(point).unwrap();
                accepted.push(point);
            }
        }
        assert!(accepted.len() > 20);
    }
}
