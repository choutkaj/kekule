//! Uniform-grid neighbor search for `within`.

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::{Arc, Mutex, PoisonError};

use crate::geometry::Point3;

/// A grid maps every cell of its box, at under two bits a cell, when the box
/// has at most this many cells per reference point, plus
/// [`DENSE_CELLS_PER_CANDIDATE`] per candidate it is built for and
/// [`DENSE_SPARE_CELLS`]; a larger box keeps only its occupied cells.
const DENSE_CELLS_PER_POINT: usize = 256;
const DENSE_CELLS_PER_CANDIDATE: usize = 16;
const DENSE_SPARE_CELLS: usize = 1 << 16;

/// A grid measured once searches around candidates, rather than gathering
/// lists for them, when it has at most this many candidates per distinct
/// point.
const SEARCH_CANDIDATES_PER_POINT: usize = 4;

/// Reference points binned into cubic cells of edge `cutoff`, so each
/// candidate inside the references' cutoff-expanded bounding box examines
/// only the 27 cells around its own; candidates outside the box are rejected
/// at once.
///
/// The grid keeps a copy of each distinct point, cell by cell: repeated
/// points, as of copies stacked in one place, are measured once. A grid
/// shared by many measurements also keeps the near lists it gathers, so
/// measurements after the first gather only the cells no earlier one visited.
#[derive(Debug, Default)]
pub(crate) struct Grid {
    low: Point3,
    high: Point3,
    edge: f64,
    /// The reciprocal of `edge`.
    scale: f64,
    cutoff: f64,
    /// The distinct points, contiguous within each cell.
    points: Vec<Point3>,
    cells: Cells,
    /// The near lists kept across measurements, for a shared grid.
    shared: Option<SharedNear>,
    /// Whether a candidate that misses the point last matched searches the
    /// cells around its own rather than gathering their list: a grid
    /// measured once by few candidates for its points uses each list little.
    search: bool,
}

type Cell = (i64, i64, i64);

#[derive(Debug)]
enum Cells {
    /// Every cell of the box, in x-major order, as one bit each: set for a
    /// cell with points. The points of the cells before slot `s` are
    /// `points[..starts[rank(s)]]`, where `rank(s)` counts the occupied
    /// cells before `s`: `ranks` holds that count at each word.
    Dense {
        dims: Cell,
        occupied: Vec<u64>,
        ranks: Vec<u32>,
        starts: Vec<u32>,
    },
    /// Only the occupied cells and the ranges of their points, for a box much
    /// larger than its points.
    Sparse {
        cells: HashMap<Cell, (u32, u32)>,
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

/// The exact bits of a point, so equal coordinates sort and compare equal.
fn bits(point: &Point3) -> (u64, u64, u64) {
    (point.x.to_bits(), point.y.to_bits(), point.z.to_bits())
}

impl Grid {
    /// Grids `references` for measuring `candidates` points once; a grid
    /// shared by many measurements passes zero.
    pub(crate) fn new(references: &[Point3], cutoff: f64, candidates: usize) -> Self {
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
            scale: edge.recip(),
            cutoff,
            points: Vec::new(),
            cells: Cells::default(),
            shared: (candidates == 0).then(SharedNear::default),
            search: false,
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
            .saturating_add(candidates.saturating_mul(DENSE_CELLS_PER_CANDIDATE))
            .saturating_add(DENSE_SPARE_CELLS);
        let count = [dims.0, dims.1, dims.2]
            .into_iter()
            .try_fold(1usize, |count, side| {
                count.checked_mul(usize::try_from(side).ok()?)
            })
            .filter(|&count| count <= limit);
        // Sorting by cell, x-major as slots are, then by exact coordinates
        // makes each cell's points contiguous and repeats adjacent.
        let mut order = references
            .iter()
            .enumerate()
            .map(|(index, point)| (grid.cell(point), index as u32))
            .collect::<Vec<_>>();
        order.sort_unstable_by(|(cell, index), (other, next)| {
            cell.cmp(other).then_with(|| {
                bits(&references[*index as usize]).cmp(&bits(&references[*next as usize]))
            })
        });
        order.dedup_by(|(cell, index), (other, next)| {
            cell == other && bits(&references[*index as usize]) == bits(&references[*next as usize])
        });
        let points = order
            .iter()
            .map(|&(_, index)| references[index as usize])
            .collect::<Vec<_>>();
        let groups = order.chunk_by(|a, b| a.0 == b.0);
        grid.cells = match count {
            Some(count) => {
                let mut occupied = vec![0u64; count.div_ceil(64)];
                let mut starts = Vec::new();
                let mut start = 0;
                for group in groups {
                    let slot = slot(dims, group[0].0);
                    occupied[slot / 64] |= 1 << (slot % 64);
                    starts.push(start as u32);
                    start += group.len();
                }
                starts.push(start as u32);
                let mut ranks = Vec::with_capacity(occupied.len());
                let mut rank = 0;
                for word in &occupied {
                    ranks.push(rank);
                    rank += word.count_ones();
                }
                Cells::Dense {
                    dims,
                    occupied,
                    ranks,
                    starts,
                }
            }
            None => {
                let mut cells = HashMap::new();
                let mut start = 0;
                for group in groups {
                    let end = start + group.len() as u32;
                    cells.insert(group[0].0, (start, end));
                    start = end;
                }
                let mut occupied = CellFilter::new(cells.len());
                for &cell in cells.keys() {
                    occupied.insert(cell);
                }
                Cells::Sparse { cells, occupied }
            }
        };
        grid.search = candidates > 0
            && candidates <= points.len().saturating_mul(SEARCH_CANDIDATES_PER_POINT);
        grid.points = points;
        grid
    }

    fn cell(&self, point: &Point3) -> Cell {
        // Points of the box lie at or above its low corner, so truncation
        // rounds down; every point is placed by this one computation.
        let axis = |value: f64, origin: f64| ((value - origin) * self.scale) as i64;
        (
            axis(point.x, self.low.x),
            axis(point.y, self.low.y),
            axis(point.z, self.low.z),
        )
    }

    /// Passes the points of the 27 cells around `cell` to `each`, a slice
    /// at a time: a dense grid's three cells along z are one slice. Stops,
    /// returning true, at the first slice for which `each` returns true.
    fn around(&self, cell: Cell, mut each: impl FnMut(&[Point3]) -> bool) -> bool {
        let range = |start: u32, end: u32| &self.points[start as usize..end as usize];
        match &self.cells {
            Cells::Dense {
                dims,
                occupied,
                ranks,
                starts,
            } => {
                // The occupied cells before `slot`.
                let rank = |slot: usize| {
                    let (word, bit) = (slot / 64, slot % 64);
                    let below = occupied.get(word).map_or(0, |word| word & ((1 << bit) - 1));
                    ranks
                        .get(word)
                        .map_or(starts.len() - 1, |&rank| rank as usize)
                        + below.count_ones() as usize
                };
                let (x, y, z) = cell;
                let (low, high) = (
                    z.saturating_sub(1).max(0),
                    z.saturating_add(1).min(dims.2 - 1),
                );
                if low > high {
                    return false;
                }
                for x in x.saturating_sub(1).max(0)..=x.saturating_add(1).min(dims.0 - 1) {
                    for y in y.saturating_sub(1).max(0)..=y.saturating_add(1).min(dims.1 - 1) {
                        let first = rank(slot(*dims, (x, y, low)));
                        let end = rank(slot(*dims, (x, y, high)) + 1);
                        if first < end && each(range(starts[first], starts[end])) {
                            return true;
                        }
                    }
                }
                false
            }
            Cells::Sparse { cells, occupied } => {
                for cell in around(cell) {
                    if occupied.may_contain(cell) {
                        if let Some(&(start, end)) = cells.get(&cell) {
                            if each(range(start, end)) {
                                return true;
                            }
                        }
                    }
                }
                false
            }
        }
    }

    /// Marks the index of every candidate within the grid's cutoff of any of
    /// the points the grid was built from.
    ///
    /// Candidates in one cell measure against one list: the points of the
    /// cells around it that can lie within the cutoff of the cell. Atoms in
    /// dense order run through few cells in turn and come back to them, so
    /// the lists of recent cells are kept, and a shared grid keeps them for
    /// every measurement. A candidate tries the point that last matched
    /// first, then searches the cells around its own, which gathers their
    /// list when the search finds no point or the cell is searched often.
    pub(crate) fn mark(
        &self,
        candidates: impl IntoIterator<Item = (usize, Point3)>,
        mark: impl FnMut(usize),
    ) {
        if self.points.is_empty() {
            return;
        }
        match &self.shared {
            Some(shared) => self.mark_with(&mut SharedLists::new(shared), candidates, mark),
            None => self.mark_with(&mut NearLists::new(), candidates, mark),
        }
    }

    fn mark_with(
        &self,
        lists: &mut impl Near,
        candidates: impl IntoIterator<Item = (usize, Point3)>,
        mut mark: impl FnMut(usize),
    ) {
        let (low, high) = (self.low, self.high);
        let squared = self.cutoff * self.cutoff;
        let mut current = None;
        let mut near = None;
        let mut previous: Option<Point3> = None;
        let mut probes = None;
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
            let within = |other: &Point3| {
                let (ex, ey, ez) = (point.x - other.x, point.y - other.y, point.z - other.z);
                ex * ex + ey * ey + ez * ez <= squared
            };
            // Consecutive atoms mostly lie near one point.
            if previous.as_ref().is_some_and(within) {
                mark(index);
                continue;
            }
            let cell = self.cell(&point);
            if current != Some(cell) {
                current = Some(cell);
                near = lists.kept(cell);
            }
            let found = match &near {
                Some(list) => lists
                    .points(list)
                    .iter()
                    .find(|other| within(other))
                    .copied(),
                None if !self.search => {
                    let list = lists.near(self, cell);
                    let found = lists
                        .points(&list)
                        .iter()
                        .find(|other| within(other))
                        .copied();
                    near = Some(list);
                    found
                }
                // A search of the cells around stops at the first point
                // within the cutoff; one that finds none has gathered the
                // cell's list. A cell searched often gathers its list.
                None => match lists.search(self, cell, within) {
                    Ok(found) => {
                        if probes.get_or_insert_with(Probes::new).count(cell) >= PROBES_PER_LIST {
                            near = Some(lists.near(self, cell));
                        }
                        Some(found)
                    }
                    Err(list) => {
                        near = Some(list);
                        None
                    }
                },
            };
            if let Some(found) = found {
                previous = Some(found);
                mark(index);
            }
        }
    }

    /// Appends to `near` the points of the cells around `cell` that can lie
    /// within the cutoff of a point the grid places in `cell`, stopping at
    /// the first of them for which `stop` holds and returning it.
    fn gather(
        &self,
        cell: Cell,
        near: &mut Vec<Point3>,
        stop: impl Fn(&Point3) -> bool,
    ) -> Option<Point3> {
        let (low, high) = (self.low, self.high);
        // A point within the cutoff of a candidate lies within this reach of
        // the candidate's cell as computed, whatever the rounding of the
        // cell's bounds and of the distance to them.
        let magnitude = [low.x, low.y, low.z, high.x, high.y, high.z]
            .into_iter()
            .fold(0.0f64, |largest, value| largest.max(value.abs()));
        let reach = self.cutoff + self.edge * 1e-6 + magnitude * 1e-9;
        let reach = reach * reach;
        let corner = |origin: f64, step: i64| origin + step as f64 * self.edge;
        let from = Point3::new(
            corner(low.x, cell.0),
            corner(low.y, cell.1),
            corner(low.z, cell.2),
        );
        let to = Point3::new(from.x + self.edge, from.y + self.edge, from.z + self.edge);
        let gap = |value: f64, from: f64, to: f64| (from - value).max(value - to).max(0.0);
        let mut found = None;
        self.around(cell, |members| {
            for other in members {
                let (x, y, z) = (
                    gap(other.x, from.x, to.x),
                    gap(other.y, from.y, to.y),
                    gap(other.z, from.z, to.z),
                );
                if x * x + y * y + z * z <= reach {
                    if stop(other) {
                        found = Some(*other);
                        return true;
                    }
                    near.push(*other);
                }
            }
            false
        });
        found
    }
}

/// The searches of one cell that find a point, after which it gathers its
/// list.
const PROBES_PER_LIST: u8 = 16;

/// Approximate counts of the searches of recent cells, in a table of fixed
/// size: a cell sharing a slot with another restarts its count.
struct Probes {
    slots: Vec<(Cell, u8)>,
}

impl Probes {
    const SLOTS: usize = 1 << 10;

    fn new() -> Self {
        Self {
            slots: vec![((i64::MIN, i64::MIN, i64::MIN), 0); Self::SLOTS],
        }
    }

    /// Counts a search of `cell`, returning the searches counted so far.
    fn count(&mut self, cell: Cell) -> u8 {
        let mut hasher = CellHasher::default();
        hasher.write_i64(cell.0);
        hasher.write_i64(cell.1);
        hasher.write_i64(cell.2);
        let slot = &mut self.slots[hasher.finish() as usize % Self::SLOTS];
        if slot.0 != cell {
            *slot = (cell, 0);
        }
        slot.1 = slot.1.saturating_add(1);
        slot.1
    }
}

/// The near lists of cells.
trait Near {
    /// A list, held while it is measured against.
    type List;

    /// The near list of `cell`, when it was gathered and is still kept.
    fn kept(&mut self, cell: Cell) -> Option<Self::List>;

    /// The first point of the near list of `cell` in `grid` for which
    /// `within` holds, gathering the list as far as that point, or the
    /// whole list, kept, when there is none.
    fn search(
        &mut self,
        grid: &Grid,
        cell: Cell,
        within: impl Fn(&Point3) -> bool,
    ) -> Result<Point3, Self::List>;

    /// The near list of `cell` in `grid`, gathered and kept.
    fn near(&mut self, grid: &Grid, cell: Cell) -> Self::List {
        match self.search(grid, cell, |_| false) {
            Ok(_) => unreachable!("no point stops gathering a whole list"),
            Err(list) => list,
        }
    }

    /// The points of `list`.
    fn points<'a>(&'a self, list: &'a Self::List) -> &'a [Point3];
}

/// The near lists of the cells visited so far.
struct NearLists {
    /// The range of each visited cell's list in `points`.
    kept: HashMap<Cell, std::ops::Range<usize>, BuildHasherDefault<CellHasher>>,
    points: Vec<Point3>,
}

/// The points kept across near lists before they are dropped together,
/// which bounds the memory of a measurement whatever its points. Unit tests
/// drop them often.
#[cfg(not(test))]
const NEAR_POINTS: usize = 1 << 16;
#[cfg(test)]
const NEAR_POINTS: usize = 16;

impl NearLists {
    fn new() -> Self {
        Self {
            kept: HashMap::default(),
            points: Vec::new(),
        }
    }
}

impl Near for NearLists {
    /// The range of the list in [`Self::points`].
    type List = std::ops::Range<usize>;

    fn kept(&mut self, cell: Cell) -> Option<Self::List> {
        self.kept.get(&cell).cloned()
    }

    fn search(
        &mut self,
        grid: &Grid,
        cell: Cell,
        within: impl Fn(&Point3) -> bool,
    ) -> Result<Point3, Self::List> {
        if self.points.len() > NEAR_POINTS {
            self.points.clear();
            self.kept.clear();
        }
        let start = self.points.len();
        if let Some(found) = grid.gather(cell, &mut self.points, within) {
            self.points.truncate(start);
            return Ok(found);
        }
        let range = start..self.points.len();
        self.kept.insert(cell, range.clone());
        Err(range)
    }

    fn points<'a>(&'a self, list: &'a Self::List) -> &'a [Point3] {
        &self.points[list.clone()]
    }
}

/// The points a shared grid keeps across its near lists, which bounds its
/// memory; lists gathered after that serve their own measurement alone.
/// Unit tests reach it often.
#[cfg(not(test))]
const SHARED_NEAR_POINTS: usize = 1 << 20;
#[cfg(test)]
const SHARED_NEAR_POINTS: usize = 64;

/// The near lists a shared grid keeps for every measurement, each gathered
/// by the first measurement that needs it.
#[derive(Debug, Default)]
struct SharedNear {
    /// A poisoned lock holds only complete lists.
    kept: Mutex<KeptNear>,
}

#[derive(Debug, Default)]
struct KeptNear {
    lists: HashMap<Cell, Arc<[Point3]>, BuildHasherDefault<CellHasher>>,
    /// The points of `lists`.
    points: usize,
}

impl SharedNear {
    fn kept(&self, cell: Cell) -> Option<Arc<[Point3]>> {
        let kept = self.kept.lock().unwrap_or_else(PoisonError::into_inner);
        kept.lists.get(&cell).map(Arc::clone)
    }

    /// Keeps `list`, gathered outside the lock, as the near list of `cell`,
    /// unless another measurement kept one first, which it returns, or the
    /// grid keeps enough points.
    fn keep(&self, cell: Cell, list: Arc<[Point3]>) -> Arc<[Point3]> {
        let mut kept = self.kept.lock().unwrap_or_else(PoisonError::into_inner);
        let full = kept.points + list.len() > SHARED_NEAR_POINTS;
        match kept.lists.entry(cell) {
            Entry::Occupied(earlier) => Arc::clone(earlier.get()),
            Entry::Vacant(_) if full => list,
            Entry::Vacant(entry) => {
                entry.insert(Arc::clone(&list));
                kept.points += list.len();
                list
            }
        }
    }
}

/// The near lists of the cells one measurement of a shared grid visited so
/// far, read without the grid's lock.
struct SharedLists<'a> {
    shared: &'a SharedNear,
    kept: HashMap<Cell, Arc<[Point3]>, BuildHasherDefault<CellHasher>>,
    /// The points of `kept`, which bounds those of lists the grid does not
    /// keep.
    points: usize,
    /// The list being gathered.
    gathering: Vec<Point3>,
}

impl<'a> SharedLists<'a> {
    fn new(shared: &'a SharedNear) -> Self {
        Self {
            shared,
            kept: HashMap::default(),
            points: 0,
            gathering: Vec::new(),
        }
    }

    fn keep(&mut self, cell: Cell, list: &Arc<[Point3]>) {
        if self.points > NEAR_POINTS {
            self.kept.clear();
            self.points = 0;
        }
        self.points += list.len();
        self.kept.insert(cell, Arc::clone(list));
    }
}

impl Near for SharedLists<'_> {
    type List = Arc<[Point3]>;

    fn kept(&mut self, cell: Cell) -> Option<Self::List> {
        if let Some(list) = self.kept.get(&cell) {
            return Some(Arc::clone(list));
        }
        let list = self.shared.kept(cell)?;
        self.keep(cell, &list);
        Some(list)
    }

    fn search(
        &mut self,
        grid: &Grid,
        cell: Cell,
        within: impl Fn(&Point3) -> bool,
    ) -> Result<Point3, Self::List> {
        self.gathering.clear();
        if let Some(found) = grid.gather(cell, &mut self.gathering, within) {
            return Ok(found);
        }
        let list = self.shared.keep(cell, Arc::from(&self.gathering[..]));
        self.keep(cell, &list);
        Err(list)
    }

    fn points<'a>(&'a self, list: &'a Self::List) -> &'a [Point3] {
        list
    }
}

/// A fast hash of cells, whose coordinates are small integers.
#[derive(Default)]
struct CellHasher(u64);

impl Hasher for CellHasher {
    fn finish(&self) -> u64 {
        self.0 ^ (self.0 >> 29)
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(u64::from(byte));
        }
    }

    fn write_u64(&mut self, value: u64) {
        self.0 = (self.0.rotate_left(21) ^ value).wrapping_mul(0x9e37_79b9_7f4a_7c15);
    }

    fn write_i64(&mut self, value: i64) {
        self.write_u64(value as u64);
    }
}

/// The dense slot of a cell inside `dims`.
fn slot(dims: Cell, (x, y, z): Cell) -> usize {
    ((x * dims.1 + y) * dims.2 + z) as usize
}

#[cfg(test)]
mod tests {
    use super::{Cells, Grid};
    use crate::geometry::Point3;

    /// Marks every candidate within `cutoff` (inclusive, Cartesian) of any
    /// reference point.
    fn within(candidates: &[Point3], references: &[Point3], cutoff: f64, mark: impl FnMut(usize)) {
        Grid::new(references, cutoff, candidates.len())
            .mark(candidates.iter().copied().enumerate(), mark);
    }

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
            for cutoff in [0.0, 0.01, 0.05, 0.4, 1.3, 10.0] {
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
            matches!(Grid::new(references, cutoff, 0).cells, Cells::Dense { .. })
        };
        assert!(dense(&references[..17], 1.3));
        assert!(!dense(&references[..17], 0.01));
        assert!(!dense(&references, 1.3));
        // Grids measured once are dense for boxes up to their candidates'.
        assert!(matches!(
            Grid::new(&references[..17], 0.01, 1_000_000).cells,
            Cells::Dense { .. }
        ));
    }

    #[test]
    fn grids_measured_by_few_candidates_search_alike() {
        // Few candidates per point search around them; repeated candidates
        // search some cells often enough to gather their lists, and
        // candidates far from every point gather empty ones.
        let references = (0..300)
            .map(|index| {
                let index = f64::from(index);
                Point3::new(
                    (index * 0.53).sin() * 7.0,
                    (index * 0.19).cos() * 7.0,
                    index * 0.03,
                )
            })
            .collect::<Vec<_>>();
        let mut candidates = (0..600)
            .map(|index| {
                let index = f64::from(index);
                Point3::new(
                    (index * 0.31).cos() * 12.0,
                    (index * 0.23).sin() * 12.0,
                    (index * 0.07).sin() * 9.0 + 4.0,
                )
            })
            .collect::<Vec<_>>();
        candidates.extend(candidates[..30].repeat(20));
        for cutoff in [0.0, 0.3, 1.2, 3.0] {
            let grid = Grid::new(&references, cutoff, candidates.len());
            assert!(grid.search, "{cutoff}");
            let mut marked = Vec::new();
            grid.mark(candidates.iter().copied().enumerate(), |index| {
                marked.push(index)
            });
            let expected = brute(&candidates, &references, cutoff);
            assert!(!expected.is_empty() || cutoff == 0.0);
            assert_eq!(marked, expected, "{cutoff}");
        }
        assert!(!Grid::new(&references, 1.2, 0).search);
        assert!(!Grid::new(&references, 1.2, 1_201).search);
    }

    #[test]
    fn repeated_points_and_scattered_candidates_change_nothing() {
        // Copies stacked in one place repeat every point; candidates come in
        // an order that jumps between cells and returns to them.
        let base = (0..60)
            .map(|index| {
                let index = f64::from(index);
                Point3::new(
                    (index * 1.7).sin() * 6.0,
                    (index * 0.3).cos() * 6.0,
                    index * 0.1,
                )
            })
            .collect::<Vec<_>>();
        let stacked = base
            .iter()
            .cycle()
            .take(base.len() * 5)
            .copied()
            .collect::<Vec<_>>();
        let candidates = (0..300)
            .map(|index| base[(index * 37) % base.len()])
            .map(|point| Point3::new(point.x + 0.9, point.y - 0.4, point.z))
            .collect::<Vec<_>>();
        let references = &stacked[..stacked.len() / 3];
        for cutoff in [0.0, 0.5, 1.5, 4.0] {
            let grid = Grid::new(references, cutoff, 0);
            assert!(grid.points.len() <= base.len());
            let mut marked = Vec::new();
            grid.mark(candidates.iter().copied().enumerate(), |index| {
                marked.push(index)
            });
            assert_eq!(marked, brute(&candidates, references, cutoff), "{cutoff}");
        }
    }

    #[test]
    fn shared_near_lists_serve_every_measurement() {
        // Many small structures measured against one shared grid, one after
        // another and at once; under `cfg(test)` the grid keeps few points,
        // so later lists serve their measurement alone.
        let references = (0..500)
            .map(|index| {
                let index = f64::from(index);
                Point3::new(
                    (index * 0.73).sin() * 8.0,
                    (index * 0.29).cos() * 8.0,
                    (index * 0.05).sin() * 4.0,
                )
            })
            .collect::<Vec<_>>();
        let structures = (0..12)
            .map(|structure| {
                (0..80)
                    .map(|index| {
                        let index = f64::from(index * 13 + structure);
                        Point3::new(
                            (index * 0.41).cos() * 10.0,
                            (index * 0.17).sin() * 10.0,
                            (index * 0.07).cos() * 6.0,
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        for cutoff in [0.0, 0.8, 2.5] {
            let grid = Grid::new(&references, cutoff, 0);
            let measure = |candidates: &[Point3]| {
                let mut marked = Vec::new();
                grid.mark(candidates.iter().copied().enumerate(), |index| {
                    marked.push(index)
                });
                marked
            };
            let expected = structures
                .iter()
                .map(|structure| brute(structure, &references, cutoff))
                .collect::<Vec<_>>();
            for _ in 0..2 {
                for (structure, expected) in structures.iter().zip(&expected) {
                    assert_eq!(&measure(structure), expected, "{cutoff}");
                }
            }
            std::thread::scope(|scope| {
                let handles = structures
                    .iter()
                    .map(|structure| scope.spawn(|| measure(structure)))
                    .collect::<Vec<_>>();
                for (handle, expected) in handles.into_iter().zip(&expected) {
                    assert_eq!(&handle.join().unwrap(), expected, "{cutoff}");
                }
            });
            let kept = grid.shared.as_ref().unwrap().kept.lock().unwrap();
            assert!(!kept.lists.is_empty() || cutoff == 0.0);
            assert!(kept.points <= super::SHARED_NEAR_POINTS);
        }
    }

    #[test]
    fn far_coordinates_keep_points_on_cell_edges() {
        // Far from the origin, the rounding of cell bounds exceeds the
        // relative slack of the cell edge; points exactly at the cutoff stay.
        for offset in [0.0, 1.0e6, -3.0e9, 7.0e12] {
            let references = [Point3::new(offset, offset, offset)];
            let candidates = (0..64)
                .map(|index| {
                    let step = f64::from(index) * 0.125;
                    Point3::new(offset + step, offset - step * 0.5, offset)
                })
                .collect::<Vec<_>>();
            for cutoff in [1.0, 2.5, 5.0] {
                let mut marked = Vec::new();
                within(&candidates, &references, cutoff, |index| marked.push(index));
                assert_eq!(
                    marked,
                    brute(&candidates, &references, cutoff),
                    "{offset} {cutoff}"
                );
            }
        }
    }
}
