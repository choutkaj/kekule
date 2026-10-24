//! Uniform-grid neighbor search for `within`.

use std::collections::HashMap;

use crate::geometry::Point3;

/// Marks every candidate within `cutoff` (inclusive, Cartesian) of any
/// reference point.
pub(crate) fn within(
    candidates: &[Point3],
    references: &[Point3],
    cutoff: f64,
    mark: impl FnMut(usize),
) {
    Grid::new(references, cutoff).mark(references, candidates.iter().copied().enumerate(), mark);
}

/// A grid stores every cell of its box when the box has at most this many
/// cells per reference point, plus [`DENSE_SPARE_CELLS`]; a larger box
/// stores only its occupied cells.
const DENSE_CELLS_PER_POINT: usize = 8;
const DENSE_SPARE_CELLS: usize = 4096;

/// Reference points binned into cubic cells of edge `cutoff`, so each
/// candidate inside the references' cutoff-expanded bounding box examines
/// only the 27 cells around its own; candidates outside the box are rejected
/// at once.
///
/// The grid holds indices into the references it was built from; queries
/// pass the same points back.
#[derive(Debug, Default)]
pub(crate) struct Grid {
    low: Point3,
    high: Point3,
    edge: f64,
    cutoff: f64,
    cells: Cells,
}

type Cell = (i64, i64, i64);

#[derive(Debug)]
enum Cells {
    /// Every cell of the box, in x-major order: the members of slot `s` are
    /// `members[starts[s]..starts[s + 1]]`.
    Dense {
        dims: Cell,
        starts: Vec<u32>,
        members: Vec<u32>,
    },
    /// Only the occupied cells, for a box much larger than its points.
    Sparse {
        cells: HashMap<Cell, Vec<u32>>,
        /// A filter of the occupied cells, probed before each lookup, so a
        /// candidate far from every point costs 27 probes rather than 27
        /// hash lookups; a false positive costs only the lookup.
        occupied: CellFilter,
    },
}

impl Default for Cells {
    fn default() -> Self {
        Self::Sparse {
            cells: HashMap::new(),
            occupied: CellFilter::new(0),
        }
    }
}

/// A one-probe Bloom filter of cells, with at least 32 bits for each cell
/// it holds.
#[derive(Debug)]
struct CellFilter {
    words: Vec<u64>,
    shift: u32,
}

impl CellFilter {
    fn new(cells: usize) -> Self {
        let bits = cells.saturating_mul(32).max(64).next_power_of_two();
        Self {
            words: vec![0; bits / 64],
            shift: 64 - bits.trailing_zeros(),
        }
    }

    fn bit(&self, (x, y, z): Cell) -> usize {
        let hash = (x as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15)
            ^ (y as u64).wrapping_mul(0xc2b2_ae3d_27d4_eb4f)
            ^ (z as u64).wrapping_mul(0x1656_67b1_9e37_79f9);
        (hash.wrapping_mul(0x9e37_79b9_7f4a_7c15) >> self.shift) as usize
    }

    fn insert(&mut self, cell: Cell) {
        let bit = self.bit(cell);
        self.words[bit / 64] |= 1 << (bit % 64);
    }

    fn may_contain(&self, cell: Cell) -> bool {
        let bit = self.bit(cell);
        self.words[bit / 64] & (1 << (bit % 64)) != 0
    }
}

/// The 27 cells within one of `cell`, itself included.
fn around((x, y, z): Cell) -> impl Iterator<Item = Cell> {
    (-1..=1).flat_map(move |dx| {
        (-1..=1).flat_map(move |dy| {
            (-1..=1).map(move |dz| {
                (
                    x.saturating_add(dx),
                    y.saturating_add(dy),
                    z.saturating_add(dz),
                )
            })
        })
    })
}

impl Grid {
    pub(crate) fn new(references: &[Point3], cutoff: f64) -> Self {
        let Some(&first) = references.first() else {
            return Self::default();
        };
        let (mut low, mut high) = (first, first);
        for point in references {
            low = Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z));
            high = Point3::new(
                high.x.max(point.x),
                high.y.max(point.y),
                high.z.max(point.z),
            );
        }
        // Cells slightly wider than the cutoff keep every neighbor within the
        // adjacent cells despite rounding. A zero cutoff still needs a
        // positive cell; any edge finds coincident points in the candidate's
        // own cell.
        let edge = if cutoff > 0.0 {
            cutoff * (1.0 + 1e-6)
        } else {
            1.0
        };
        let mut grid = Self {
            low: Point3::new(low.x - cutoff, low.y - cutoff, low.z - cutoff),
            high: Point3::new(high.x + cutoff, high.y + cutoff, high.z + cutoff),
            edge,
            cutoff,
            cells: Cells::default(),
        };
        // Cells count from the low corner, so every point of the box lies
        // between the origin cell and the high corner's cell.
        let (x, y, z) = grid.cell(&grid.high);
        let dims = (
            x.saturating_add(1),
            y.saturating_add(1),
            z.saturating_add(1),
        );
        let limit = references
            .len()
            .saturating_mul(DENSE_CELLS_PER_POINT)
            .saturating_add(DENSE_SPARE_CELLS);
        let count = [dims.0, dims.1, dims.2]
            .into_iter()
            .try_fold(1usize, |count, side| {
                count.checked_mul(usize::try_from(side).ok()?)
            })
            .filter(|&count| count <= limit);
        grid.cells = match count {
            Some(count) => {
                // A counting sort keeps each cell's members in index order.
                let slots = references
                    .iter()
                    .map(|point| slot(dims, grid.cell(point)))
                    .collect::<Vec<_>>();
                let mut starts = vec![0u32; count + 1];
                for &slot in &slots {
                    starts[slot + 1] += 1;
                }
                for index in 1..starts.len() {
                    starts[index] += starts[index - 1];
                }
                let mut next = starts.clone();
                let mut members = vec![0u32; references.len()];
                for (index, &slot) in slots.iter().enumerate() {
                    members[next[slot] as usize] = index as u32;
                    next[slot] += 1;
                }
                Cells::Dense {
                    dims,
                    starts,
                    members,
                }
            }
            None => {
                let mut cells: HashMap<Cell, Vec<u32>> = HashMap::new();
                for (index, point) in references.iter().enumerate() {
                    cells
                        .entry(grid.cell(point))
                        .or_default()
                        .push(index as u32);
                }
                let mut occupied = CellFilter::new(cells.len());
                for &cell in cells.keys() {
                    occupied.insert(cell);
                }
                Cells::Sparse { cells, occupied }
            }
        };
        grid
    }

    fn cell(&self, point: &Point3) -> Cell {
        let axis = |value: f64, origin: f64| ((value - origin) / self.edge).floor() as i64;
        (
            axis(point.x, self.low.x),
            axis(point.y, self.low.y),
            axis(point.z, self.low.z),
        )
    }

    /// The references in one cell.
    fn members(&self, cell: Cell) -> &[u32] {
        match &self.cells {
            Cells::Dense {
                dims,
                starts,
                members,
            } => {
                let (x, y, z) = cell;
                if x < 0 || y < 0 || z < 0 || x >= dims.0 || y >= dims.1 || z >= dims.2 {
                    return &[];
                }
                let slot = slot(*dims, cell);
                &members[starts[slot] as usize..starts[slot + 1] as usize]
            }
            Cells::Sparse { cells, occupied } => {
                if !occupied.may_contain(cell) {
                    return &[];
                }
                cells.get(&cell).map_or(&[], Vec::as_slice)
            }
        }
    }

    fn is_empty(&self) -> bool {
        match &self.cells {
            Cells::Dense { members, .. } => members.is_empty(),
            Cells::Sparse { cells, .. } => cells.is_empty(),
        }
    }

    /// Marks the index of every candidate within the grid's cutoff of any of
    /// `references`, the points the grid was built from.
    pub(crate) fn mark(
        &self,
        references: &[Point3],
        candidates: impl IntoIterator<Item = (usize, Point3)>,
        mut mark: impl FnMut(usize),
    ) {
        if self.is_empty() {
            return;
        }
        let (low, high) = (self.low, self.high);
        let squared = self.cutoff * self.cutoff;
        for (index, point) in candidates {
            if point.x < low.x
                || point.y < low.y
                || point.z < low.z
                || point.x > high.x
                || point.y > high.y
                || point.z > high.z
            {
                continue;
            }
            let near = around(self.cell(&point)).any(|cell| {
                self.members(cell).iter().any(|&member| {
                    let other = references[member as usize];
                    let (ex, ey, ez) = (point.x - other.x, point.y - other.y, point.z - other.z);
                    ex * ex + ey * ey + ez * ez <= squared
                })
            });
            if near {
                mark(index);
            }
        }
    }
}

/// The dense slot of a cell inside `dims`.
fn slot(dims: Cell, (x, y, z): Cell) -> usize {
    ((x * dims.1 + y) * dims.2 + z) as usize
}

#[cfg(test)]
mod tests {
    use super::{within, Cells, Grid};
    use crate::geometry::Point3;

    fn brute(points: &[Point3], references: &[Point3], cutoff: f64) -> Vec<usize> {
        points
            .iter()
            .enumerate()
            .filter(|(_, point)| {
                references.iter().any(|other| {
                    let (x, y, z) = (point.x - other.x, point.y - other.y, point.z - other.z);
                    x * x + y * y + z * z <= cutoff * cutoff
                })
            })
            .map(|(index, _)| index)
            .collect()
    }

    #[test]
    fn grid_agrees_with_brute_force() {
        let mut points = (0..400)
            .map(|index| {
                let index = f64::from(index);
                Point3::new(
                    (index * 0.37).sin() * 3.0,
                    (index * 0.11).cos() * 3.0,
                    index * 0.01,
                )
            })
            .collect::<Vec<_>>();
        // A far point makes small cutoffs' boxes too large to store densely.
        points.push(Point3::new(1.0e4, -2.0e3, 5.0e2));
        let mut references = points[..17].to_vec();
        for far in [false, true] {
            if far {
                references.push(points[400]);
            }
            for cutoff in [0.0, 0.05, 0.4, 1.3, 10.0] {
                let mut grid = Vec::new();
                within(&points, &references, cutoff, |index| grid.push(index));
                assert_eq!(
                    grid,
                    brute(&points, &references, cutoff),
                    "cutoff {cutoff}, far {far}"
                );
            }
        }
        let dense = |references: &[Point3], cutoff| {
            matches!(Grid::new(references, cutoff).cells, Cells::Dense { .. })
        };
        assert!(dense(&references[..17], 1.3));
        assert!(!dense(&references[..17], 0.05));
        assert!(!dense(&references, 1.3));
    }
}
