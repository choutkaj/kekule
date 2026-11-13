//! A uniform grid of points with payloads, for searches within a cutoff.

use super::Point3;

type Cell = (i64, i64, i64);

/// Points binned into cubic cells at least as wide as the cutoff, sorted by
/// cell, so the points within the cutoff of a point lie in the 27 cells
/// around its own. Each point carries a payload, such as its index.
///
/// Contact occupancy of atom groups and non-covalent interaction search use
/// it; a grid is rebuilt in place per structure, reusing its storage.
#[derive(Debug, Clone)]
pub(crate) struct PointGrid<T> {
    low: Point3,
    scale: f64,
    /// The highest cell on each axis.
    high: Cell,
    entries: Vec<(Cell, Point3, T)>,
}

impl<T> Default for PointGrid<T> {
    fn default() -> Self {
        Self {
            low: Point3::origin(),
            scale: 1.0,
            high: (0, 0, 0),
            entries: Vec::new(),
        }
    }
}

impl<T: Copy> PointGrid<T> {
    /// Rebuilds the grid over `points`, whose bounding box is `bounds`.
    pub(crate) fn build(
        &mut self,
        points: impl IntoIterator<Item = (Point3, T)>,
        (low, high): (Point3, Point3),
        cutoff: f64,
    ) {
        // Cells slightly wider than the cutoff keep every contact within the
        // adjacent cells despite rounding, and at least a millionth of the
        // box, so cell coordinates stay small enough to round alike. Any
        // positive edge finds coincident points.
        let extent = (high.x - low.x).max(high.y - low.y).max(high.z - low.z);
        let edge = (cutoff * (1.0 + 1e-6)).max(extent * 1e-6);
        let edge = if edge > 0.0 { edge } else { 1.0 };
        self.low = low;
        self.scale = edge.recip();
        let mut entries = std::mem::take(&mut self.entries);
        entries.clear();
        entries.extend(
            points
                .into_iter()
                .map(|(point, payload)| (self.cell(&point), point, payload)),
        );
        self.entries = entries;
        // A stable sort keeps the input order within a cell.
        self.entries.sort_by_key(|(cell, _, _)| *cell);
        self.high = self.cell(&high);
    }

    /// The cell of `point`; one monotone computation places every point, so
    /// points nearer than an edge along an axis lie in adjacent cells.
    fn cell(&self, point: &Point3) -> Cell {
        let axis = |value: f64, origin: f64| ((value - origin) * self.scale).floor() as i64;
        (
            axis(point.x, self.low.x),
            axis(point.y, self.low.y),
            axis(point.z, self.low.z),
        )
    }

    /// Passes the entries of the 27 cells around `candidate`'s to `each`, a
    /// run of cells at a time, stopping at the first run for which `each`
    /// returns true; returns whether one did.
    fn around(
        &self,
        candidate: &Point3,
        mut each: impl FnMut(&[(Cell, Point3, T)]) -> bool,
    ) -> bool {
        let (x, y, z) = self.cell(candidate);
        if x < -1
            || y < -1
            || z < -1
            || x > self.high.0.saturating_add(1)
            || y > self.high.1.saturating_add(1)
            || z > self.high.2.saturating_add(1)
        {
            return false;
        }
        // Cells sort x-major, so the three cells along z around each (x, y)
        // column are one run.
        (x - 1..=x + 1).any(|x| {
            (y - 1..=y + 1).any(|y| {
                let first = self
                    .entries
                    .partition_point(|(cell, _, _)| *cell < (x, y, z - 1));
                let end = self
                    .entries
                    .partition_point(|(cell, _, _)| *cell <= (x, y, z + 1));
                each(&self.entries[first..end])
            })
        })
    }

    /// Whether any of `candidates` is within `cutoff` of a grid point, by
    /// `within`.
    pub(crate) fn any_within(
        &self,
        candidates: &[Point3],
        within: impl Fn(&Point3, &Point3) -> bool,
    ) -> bool {
        candidates.iter().any(|candidate| {
            self.around(candidate, |entries| {
                entries.iter().any(|(_, point, _)| within(candidate, point))
            })
        })
    }

    /// Passes every grid point within the grid's cell reach of `candidate`
    /// for which `within` holds, with its payload, to `each`, in grid order.
    pub(crate) fn for_each_near(
        &self,
        candidate: &Point3,
        within: impl Fn(&Point3, &Point3) -> bool,
        mut each: impl FnMut(&Point3, T),
    ) {
        self.around(candidate, |entries| {
            for (_, point, payload) in entries {
                if within(candidate, point) {
                    each(point, *payload);
                }
            }
            false
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_exactly_the_points_within_the_cutoff() {
        let points = (0..200)
            .map(|i| {
                let i = f64::from(i);
                Point3::new((i * 0.37).sin() * 2.0, (i * 0.11).cos() * 2.0, i * 0.01)
            })
            .collect::<Vec<_>>();
        let low = Point3::new(-2.0, -2.0, 0.0);
        let high = Point3::new(2.0, 2.0, 2.0);
        let mut grid = PointGrid::default();
        grid.build(points.iter().copied().zip(0usize..), (low, high), 0.3);
        let within = |a: &Point3, b: &Point3| (*a - *b).norm() <= 0.3;
        for candidate in &points {
            let mut found = Vec::new();
            grid.for_each_near(candidate, within, |_, index| found.push(index));
            found.sort_unstable();
            let expected = (0..points.len())
                .filter(|&index| within(candidate, &points[index]))
                .collect::<Vec<_>>();
            assert_eq!(found, expected);
            assert!(grid.any_within(std::slice::from_ref(candidate), within));
        }
        assert!(!grid.any_within(&[Point3::new(9.0, 9.0, 9.0)], within));
    }
}
