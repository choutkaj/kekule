//! A uniform grid of points with payloads, for searches within a cutoff.

use std::collections::HashMap;
use std::convert::Infallible;
use std::hash::{BuildHasherDefault, Hasher};

use super::{Lattice, PeriodicCellError, Point3};

/// A grid cell, by its integer coordinates along each axis.
pub(crate) type Cell = (i64, i64, i64);

/// A map keyed by cells.
type CellMap<V> = HashMap<Cell, V, BuildHasherDefault<CellHasher>>;

/// An open grid indexes every cell of its points' box, at about two bits a
/// cell, when the box has at most this many cells per point, plus
/// [`DENSE_CELLS_PER_QUERY`] per expected search and [`DENSE_SPARE_CELLS`];
/// a larger box indexes only its occupied cells.
const DENSE_CELLS_PER_POINT: usize = 256;
const DENSE_CELLS_PER_QUERY: usize = 16;
const DENSE_SPARE_CELLS: usize = 4096;

/// A box of at most this many cells per point, plus [`DENSE_SPARE_CELLS`],
/// keeps the first entry of each cell, so a search reads a run's bounds
/// directly rather than counting occupied cells.
const DIRECT_CELLS_PER_POINT: usize = 8;

/// Points with payloads, binned into cells at least as wide as a cutoff, so
/// the points within the cutoff of a point lie in the 27 cells around its
/// own.
///
/// An open grid bins Cartesian points into cubic cells and is built from
/// all its points at once, rebuilt in place to reuse its storage. A
/// periodic grid bins a lattice's wrapped fractional coordinates and grows
/// a point at a time. A search visits the cells around a candidate's in a
/// fixed order, x-major for an open grid, and each cell's points in the
/// order they were added. The grid only finds candidates; its callers
/// measure them.
#[derive(Debug, Clone)]
pub(crate) struct PointGrid<T> {
    space: Space,
    index: Index<T>,
    /// The points of an open grid, in cell order.
    entries: Vec<(Point3, T)>,
    /// A build's points with their cells, kept for its storage.
    keyed: Vec<(Cell, Point3, T)>,
    /// A counting sort's next free entry of each cell, kept for its storage.
    next: Vec<usize>,
    len: usize,
    /// How many cells along each axis a search reaches: 1, or 2 where
    /// rounding of cell coordinates far from the box's corner can move a
    /// point by up to a cell.
    reach: i64,
}

#[derive(Debug, Clone)]
enum Space {
    /// Cubic cells of edge `edge` counted from `low`, for points in the box
    /// from `low` to `high`.
    Open {
        low: Point3,
        high: Point3,
        edge: f64,
        /// The reciprocal of `edge`.
        scale: f64,
        cutoff: f64,
    },
    /// `bins` fractional bins along each lattice vector.
    Periodic { lattice: Lattice, bins: [i64; 3] },
}

#[derive(Debug, Clone)]
enum Index<T> {
    /// Every cell of the box of `dims` cells from `origin`, in x-major
    /// order: the points of the cells before slot `s` are
    /// `entries[..starts[s]]`.
    Direct {
        origin: Cell,
        dims: Cell,
        starts: Vec<usize>,
    },
    /// Every cell of the box of `dims` cells from `origin`, in x-major
    /// order, as one bit each: set for a cell with points. The points of
    /// the cells before slot `s` are `entries[..starts[rank(s)]]`, where
    /// `rank(s)` counts the occupied cells before `s`: `ranks` holds that
    /// count at each word.
    Dense {
        origin: Cell,
        dims: Cell,
        occupied: Vec<u64>,
        ranks: Vec<usize>,
        starts: Vec<usize>,
    },
    /// Only the occupied cells, from `origin` to `last`, with their ranges
    /// of `entries`.
    Sparse {
        origin: Cell,
        last: Cell,
        ranges: CellMap<(usize, usize)>,
        /// Probed before each lookup, so a candidate far from every point
        /// costs 27 probes rather than 27 hash lookups.
        filter: CellFilter,
    },
    /// The points of each occupied cell of a periodic grid.
    Lists(CellMap<Vec<(Point3, T)>>),
}

impl<T> Default for PointGrid<T> {
    fn default() -> Self {
        Self {
            space: Space::Open {
                low: Point3::origin(),
                high: Point3::origin(),
                edge: 1.0,
                scale: 1.0,
                cutoff: 0.0,
            },
            index: Index::Sparse {
                origin: (0, 0, 0),
                last: (0, 0, 0),
                ranges: CellMap::default(),
                filter: CellFilter::new(0),
            },
            entries: Vec::new(),
            keyed: Vec::new(),
            next: Vec::new(),
            len: 0,
            reach: 1,
        }
    }
}

/// The smallest box holding `points`; the origin for none.
pub(crate) fn bounding_box(points: impl IntoIterator<Item = Point3>) -> (Point3, Point3) {
    let mut points = points.into_iter();
    let Some(first) = points.next() else {
        return (Point3::origin(), Point3::origin());
    };
    points.fold((first, first), |(low, high), point| {
        (
            Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z)),
            Point3::new(
                high.x.max(point.x),
                high.y.max(point.y),
                high.z.max(point.z),
            ),
        )
    })
}

impl<T: Copy> PointGrid<T> {
    /// Rebuilds the grid as an open grid over `points`, whose bounding box
    /// is `bounds`, for searches within `cutoff`; `queries` estimates the
    /// searches it will serve, which buys it a denser index.
    pub(crate) fn build(
        &mut self,
        points: impl IntoIterator<Item = (Point3, T)>,
        (low, high): (Point3, Point3),
        cutoff: f64,
        queries: usize,
    ) {
        // Cells are slightly wider than the cutoff, so a contact lies in an
        // adjacent cell, and never narrower than a 2^50th of the box or
        // the smallest normal width, so cell coordinates and the reciprocal
        // stay finite. Placing a point rounds its cell coordinate by at
        // most `keys * 2^-51` cells; where that can carry a contact past
        // an adjacent cell, searches reach two cells. Any positive edge
        // finds coincident points.
        let extent = (high.x - low.x).max(high.y - low.y).max(high.z - low.z);
        let edge = (cutoff * (1.0 + 1e-6))
            .max(extent * 2.0_f64.powi(-50))
            .max(4.0 * f64::MIN_POSITIVE);
        let scale = edge.recip();
        let keys = extent * scale;
        self.reach = if keys * 2.0_f64.powi(-51) < 4e-7 {
            1
        } else {
            2
        };
        self.space = Space::Open {
            low,
            high,
            edge,
            scale,
            cutoff,
        };
        let mut keyed = std::mem::take(&mut self.keyed);
        keyed.clear();
        keyed.extend(
            points
                .into_iter()
                .map(|(point, payload)| (open_cell(low, scale, &point), point, payload)),
        );
        self.len = keyed.len();
        let (origin, last) = keyed.iter().fold(
            (
                (i64::MAX, i64::MAX, i64::MAX),
                (i64::MIN, i64::MIN, i64::MIN),
            ),
            |(low, high), &(cell, _, _)| {
                (
                    (low.0.min(cell.0), low.1.min(cell.1), low.2.min(cell.2)),
                    (high.0.max(cell.0), high.1.max(cell.1), high.2.max(cell.2)),
                )
            },
        );
        let side = |low: i64, high: i64| high.checked_sub(low)?.checked_add(1);
        let budget = self
            .len
            .saturating_mul(DENSE_CELLS_PER_POINT)
            .saturating_add(queries.saturating_mul(DENSE_CELLS_PER_QUERY))
            .saturating_add(DENSE_SPARE_CELLS);
        // The box's cells, when it has points and few enough to index all.
        let dense = (|| {
            let dims = (
                side(origin.0, last.0)?,
                side(origin.1, last.1)?,
                side(origin.2, last.2)?,
            );
            let count = [dims.0, dims.1, dims.2]
                .into_iter()
                .try_fold(1usize, |count, side| {
                    count.checked_mul(usize::try_from(side).ok()?)
                })?;
            (self.len > 0 && count <= budget).then_some((dims, count))
        })();
        let direct = self
            .len
            .saturating_mul(DIRECT_CELLS_PER_POINT)
            .saturating_add(DENSE_SPARE_CELLS);
        let previous = std::mem::replace(&mut self.index, Index::Lists(CellMap::default()));
        self.index = match (dense, previous) {
            (Some((dims, count)), previous) if count <= direct => {
                let mut starts = match previous {
                    Index::Direct { starts, .. } | Index::Dense { starts, .. } => starts,
                    _ => Vec::new(),
                };
                self.sort_by_bucket(&keyed, count, |cell| slot(origin, dims, cell), &mut starts);
                Index::Direct {
                    origin,
                    dims,
                    starts,
                }
            }
            (Some((dims, count)), previous) => {
                let (mut occupied, mut ranks, mut starts) = match previous {
                    Index::Dense {
                        occupied,
                        ranks,
                        starts,
                        ..
                    } => (occupied, ranks, starts),
                    Index::Direct { starts, .. } => (Vec::new(), Vec::new(), starts),
                    _ => Default::default(),
                };
                occupied.clear();
                occupied.resize(count.div_ceil(64), 0);
                for &(cell, _, _) in &keyed {
                    let slot = slot(origin, dims, cell);
                    occupied[slot / 64] |= 1 << (slot % 64);
                }
                ranks.clear();
                let mut rank = 0;
                for word in &occupied {
                    ranks.push(rank);
                    rank += word.count_ones() as usize;
                }
                self.sort_by_bucket(
                    &keyed,
                    rank,
                    |cell| dense_rank(&occupied, &ranks, slot(origin, dims, cell)),
                    &mut starts,
                );
                Index::Dense {
                    origin,
                    dims,
                    occupied,
                    ranks,
                    starts,
                }
            }
            (None, previous) => {
                let mut ranges = match previous {
                    Index::Sparse { ranges, .. } => ranges,
                    _ => CellMap::default(),
                };
                // A stable sort keeps the input order within a cell.
                keyed.sort_by_key(|(cell, _, _)| *cell);
                self.entries.clear();
                self.entries
                    .extend(keyed.iter().map(|&(_, point, payload)| (point, payload)));
                ranges.clear();
                let mut start = 0;
                for group in keyed.chunk_by(|a, b| a.0 == b.0) {
                    ranges.insert(group[0].0, (start, start + group.len()));
                    start += group.len();
                }
                let mut filter = CellFilter::new(ranges.len());
                for &cell in ranges.keys() {
                    filter.insert(cell);
                }
                Index::Sparse {
                    origin,
                    last,
                    ranges,
                    filter,
                }
            }
        };
        keyed.clear();
        self.keyed = keyed;
    }

    /// Places the points of `keyed`, which is not empty, in `entries` by a
    /// counting sort into `buckets` buckets, which keeps the input order
    /// within a bucket. Sets `starts[b]` to the first entry of bucket `b`,
    /// and `starts[buckets]` to the number of points.
    fn sort_by_bucket(
        &mut self,
        keyed: &[(Cell, Point3, T)],
        buckets: usize,
        bucket: impl Fn(Cell) -> usize,
        starts: &mut Vec<usize>,
    ) {
        starts.clear();
        starts.resize(buckets + 1, 0);
        for &(cell, _, _) in keyed {
            starts[bucket(cell) + 1] += 1;
        }
        for index in 1..starts.len() {
            starts[index] += starts[index - 1];
        }
        let (_, point, payload) = keyed[0];
        self.entries.clear();
        self.entries.resize(keyed.len(), (point, payload));
        let mut next = std::mem::take(&mut self.next);
        next.clear();
        next.extend_from_slice(&starts[..buckets]);
        for &(cell, point, payload) in keyed {
            let free = &mut next[bucket(cell)];
            self.entries[*free] = (point, payload);
            *free += 1;
        }
        self.next = next;
    }

    /// An empty periodic grid over `lattice` for searches within `cutoff`,
    /// whose fractional bins are sized from the reciprocal-vector norms, so
    /// that a displacement within the cutoff changes each wrapped fractional
    /// coordinate by at most one bin. A skewed cell still needs exact
    /// nearest-image measurement of the candidates.
    ///
    /// Fails with [`PeriodicCellError::NumericalFailure`] for a cutoff that is not
    /// positive or bins that are not finite.
    pub(crate) fn periodic(lattice: &Lattice, cutoff: f64) -> Result<Self, PeriodicCellError> {
        let spans = lattice.reciprocal_norms().map(|norm| norm * cutoff);
        if cutoff <= 0.0 || spans.iter().any(|span| !span.is_finite() || *span <= 0.0) {
            return Err(PeriodicCellError::NumericalFailure);
        }
        let bins = spans.map(|span| {
            (1.0 / (span * (1.0 + 128.0 * f64::EPSILON)))
                .floor()
                .clamp(1.0, 1_000_000.0) as i64
        });
        Ok(Self {
            space: Space::Periodic {
                lattice: lattice.clone(),
                bins,
            },
            index: Index::Lists(CellMap::default()),
            ..Self::default()
        })
    }

    /// Adds `point` to a periodic grid, after the points already there.
    ///
    /// Fails with [`PeriodicCellError::NumericalFailure`] when its fractional
    /// coordinates are not finite or too large to wrap exactly.
    pub(crate) fn insert(&mut self, point: Point3, payload: T) -> Result<(), PeriodicCellError> {
        let cell = self.locate(&point)?;
        let Index::Lists(lists) = &mut self.index else {
            unreachable!("only a periodic grid grows a point at a time");
        };
        lists.entry(cell).or_default().push((point, payload));
        self.len += 1;
        Ok(())
    }

    /// The number of points.
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// Whether the grid holds no points.
    pub(crate) fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Whether the grid indexes every cell of its box.
    #[cfg(test)]
    pub(crate) fn is_dense(&self) -> bool {
        matches!(self.index, Index::Direct { .. } | Index::Dense { .. })
    }

    /// Whether the grid keeps the first entry of every cell of its box.
    #[cfg(test)]
    pub(crate) fn is_direct(&self) -> bool {
        matches!(self.index, Index::Direct { .. })
    }

    /// The cell of `point`; `None` when a periodic grid cannot wrap it.
    pub(crate) fn cell(&self, point: &Point3) -> Option<Cell> {
        self.locate(point).ok()
    }

    fn locate(&self, point: &Point3) -> Result<Cell, PeriodicCellError> {
        match &self.space {
            Space::Open { low, scale, .. } => Ok(open_cell(*low, *scale, point)),
            Space::Periodic { lattice, bins } => {
                let fractional = lattice.fractional_canonical(*point - Point3::origin())?;
                if fractional
                    .iter()
                    .any(|value| value.abs() >= 2.0_f64.powi(48))
                {
                    return Err(PeriodicCellError::NumericalFailure);
                }
                let axis = |axis: usize| {
                    ((fractional[axis].rem_euclid(1.0) * bins[axis] as f64).floor() as i64)
                        .min(bins[axis] - 1)
                };
                Ok((axis(0), axis(1), axis(2)))
            }
        }
    }

    /// Passes the points of the cells around `cell`, within [`Self::reach`]
    /// along each axis, to `each`, a run of cells at a time, stopping at the first run for which `each` returns
    /// true; returns whether one did.
    fn around(&self, cell: Cell, mut each: impl FnMut(&[(Point3, T)]) -> bool) -> bool {
        if self.len == 0 {
            return false;
        }
        let (x, y, z) = cell;
        let reach = self.reach;
        match &self.index {
            Index::Direct {
                origin,
                dims,
                starts,
            } => {
                let Some((x, y, z)) = relative(cell, *origin, *dims, reach) else {
                    return false;
                };
                let (low, high) = ((z - reach).max(0), (z + reach).min(dims.2 - 1));
                if low > high {
                    return false;
                }
                // Cells sort x-major, so the cells along z around each (x, y)
                // column are one run.
                for x in (x - reach).max(0)..=(x + reach).min(dims.0 - 1) {
                    for y in (y - reach).max(0)..=(y + reach).min(dims.1 - 1) {
                        let first = starts[dense_slot(*dims, (x, y, low))];
                        let end = starts[dense_slot(*dims, (x, y, high)) + 1];
                        if first < end && each(&self.entries[first..end]) {
                            return true;
                        }
                    }
                }
                false
            }
            Index::Dense {
                origin,
                dims,
                occupied,
                ranks,
                starts,
            } => {
                let Some((x, y, z)) = relative(cell, *origin, *dims, reach) else {
                    return false;
                };
                // The occupied cells before `slot`.
                let rank = |slot: usize| {
                    let (word, bit) = (slot / 64, slot % 64);
                    let below = occupied.get(word).map_or(0, |word| word & ((1 << bit) - 1));
                    ranks.get(word).map_or(starts.len() - 1, |&rank| rank)
                        + below.count_ones() as usize
                };
                let (low, high) = ((z - reach).max(0), (z + reach).min(dims.2 - 1));
                if low > high {
                    return false;
                }
                // Cells sort x-major, so the cells along z around each (x, y)
                // column are one run.
                for x in (x - reach).max(0)..=(x + reach).min(dims.0 - 1) {
                    for y in (y - reach).max(0)..=(y + reach).min(dims.1 - 1) {
                        let first = rank(dense_slot(*dims, (x, y, low)));
                        let end = rank(dense_slot(*dims, (x, y, high)) + 1);
                        if first < end && each(&self.entries[starts[first]..starts[end]]) {
                            return true;
                        }
                    }
                }
                false
            }
            Index::Sparse {
                origin,
                last,
                ranges,
                filter,
            } => {
                let reaches = |value: i64, origin: i64, last: i64| {
                    value >= origin.saturating_sub(reach) && value <= last.saturating_add(reach)
                };
                if !(reaches(x, origin.0, last.0)
                    && reaches(y, origin.1, last.1)
                    && reaches(z, origin.2, last.2))
                {
                    return false;
                }
                let around =
                    |value: i64| (-reach..=reach).filter_map(move |step| value.checked_add(step));
                for x in around(x) {
                    for y in around(y) {
                        for z in around(z) {
                            if !filter.may_contain((x, y, z)) {
                                continue;
                            }
                            if let Some(&(start, end)) = ranges.get(&(x, y, z)) {
                                if each(&self.entries[start..end]) {
                                    return true;
                                }
                            }
                        }
                    }
                }
                false
            }
            Index::Lists(lists) => {
                let Space::Periodic { bins, .. } = &self.space else {
                    unreachable!("only a periodic grid keeps lists");
                };
                // A lattice of fewer than three bins along an axis wraps
                // onto a cell already visited, which is visited once.
                let axes: [[i64; 3]; 3] = std::array::from_fn(|axis| {
                    let (key, bins) = ([x, y, z][axis], bins[axis]);
                    [(key + bins - 1) % bins, key, (key + 1) % bins]
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
                            if let Some(list) = lists.get(&(axes[0][a], axes[1][b], axes[2][c])) {
                                if each(list) {
                                    return true;
                                }
                            }
                        }
                    }
                }
                false
            }
        }
    }

    /// Passes the points around `cell` to `test` in grid order until it
    /// returns true or fails.
    fn search<E>(
        &self,
        cell: Cell,
        mut test: impl FnMut(&Point3, T) -> Result<bool, E>,
    ) -> Result<bool, E> {
        let mut outcome = Ok(false);
        self.around(cell, |run| {
            for (point, payload) in run {
                match test(point, *payload) {
                    Ok(false) => {}
                    stop => {
                        outcome = stop;
                        return true;
                    }
                }
            }
            false
        });
        outcome
    }

    /// Whether `test` holds for a point around `candidate`, passing the
    /// points to it in grid order until it does or fails.
    ///
    /// Fails as `test` does, or with [`PeriodicCellError::NumericalFailure`]
    /// when a periodic grid cannot wrap `candidate`.
    pub(crate) fn try_any_near<E: From<PeriodicCellError>>(
        &self,
        candidate: &Point3,
        test: impl FnMut(&Point3, T) -> Result<bool, E>,
    ) -> Result<bool, E> {
        let cell = self.locate(candidate)?;
        self.search(cell, test)
    }

    /// Whether `test` holds for a point around `candidate`, passing the
    /// points to it in grid order until it does. A candidate a periodic
    /// grid cannot wrap has none.
    pub(crate) fn any_near(
        &self,
        candidate: &Point3,
        mut test: impl FnMut(&Point3, T) -> bool,
    ) -> bool {
        let Ok(cell) = self.locate(candidate) else {
            return false;
        };
        match self.search::<Infallible>(cell, |point, payload| Ok(test(point, payload))) {
            Ok(found) => found,
            Err(never) => match never {},
        }
    }

    /// Passes every point around `candidate`, with its payload, to `each`,
    /// in grid order. A candidate a periodic grid cannot wrap has none.
    pub(crate) fn for_each_near(&self, candidate: &Point3, mut each: impl FnMut(&Point3, T)) {
        self.any_near(candidate, |point, payload| {
            each(point, payload);
            false
        });
    }

    /// Passes to `each`, in grid order, the points of the cells around
    /// `cell` that can lie within the cutoff of a point an open grid places
    /// in `cell`, stopping at the first for which `each` returns true;
    /// returns whether one did. A periodic grid has none.
    pub(crate) fn near_cell(&self, cell: Cell, mut each: impl FnMut(&Point3) -> bool) -> bool {
        let Space::Open {
            low,
            high,
            edge,
            cutoff,
            ..
        } = self.space
        else {
            return false;
        };
        // A point within the cutoff of a candidate lies within this reach of
        // the candidate's cell as computed, whatever the rounding of the
        // cell's bounds and of the distance to them.
        let magnitude = [low.x, low.y, low.z, high.x, high.y, high.z]
            .into_iter()
            .fold(0.0f64, |largest, value| largest.max(value.abs()));
        let reach = cutoff + edge * 1e-6 + magnitude * 1e-9;
        let reach = reach * reach;
        let corner = |origin: f64, step: i64| origin + step as f64 * edge;
        let from = Point3::new(
            corner(low.x, cell.0),
            corner(low.y, cell.1),
            corner(low.z, cell.2),
        );
        let to = Point3::new(from.x + edge, from.y + edge, from.z + edge);
        let gap = |value: f64, from: f64, to: f64| (from - value).max(value - to).max(0.0);
        self.around(cell, |run| {
            run.iter().any(|(point, _)| {
                let (x, y, z) = (
                    gap(point.x, from.x, to.x),
                    gap(point.y, from.y, to.y),
                    gap(point.z, from.z, to.z),
                );
                x * x + y * y + z * z <= reach && each(point)
            })
        })
    }
}

/// The open-grid cell of `point`; one monotone computation places every
/// point, so points nearer than an edge along an axis lie in adjacent cells.
fn open_cell(low: Point3, scale: f64, point: &Point3) -> Cell {
    let axis = |value: f64, origin: f64| ((value - origin) * scale).floor() as i64;
    (
        axis(point.x, low.x),
        axis(point.y, low.y),
        axis(point.z, low.z),
    )
}

/// `cell` relative to the box of `dims` cells from `origin`, when the cells
/// within `reach` of it can meet the box.
fn relative(cell: Cell, origin: Cell, dims: Cell, reach: i64) -> Option<Cell> {
    let axis = |value: i64, origin: i64, side: i64| {
        let offset = value.checked_sub(origin)?;
        (-reach..side + reach).contains(&offset).then_some(offset)
    };
    Some((
        axis(cell.0, origin.0, dims.0)?,
        axis(cell.1, origin.1, dims.1)?,
        axis(cell.2, origin.2, dims.2)?,
    ))
}

/// The dense slot of `cell`, inside the box of `dims` cells from `origin`.
fn slot(origin: Cell, dims: Cell, cell: Cell) -> usize {
    dense_slot(
        dims,
        (cell.0 - origin.0, cell.1 - origin.1, cell.2 - origin.2),
    )
}

/// The occupied cells before `slot` of a dense index, a slot inside its box.
fn dense_rank(occupied: &[u64], ranks: &[usize], slot: usize) -> usize {
    let (word, bit) = (slot / 64, slot % 64);
    ranks[word] + (occupied[word] & ((1 << bit) - 1)).count_ones() as usize
}

/// The dense slot of a cell relative to its box of `dims` cells.
fn dense_slot(dims: Cell, (x, y, z): Cell) -> usize {
    ((x * dims.1 + y) * dims.2 + z) as usize
}

/// A one-probe Bloom filter of cells, with at least 32 bits for each cell
/// it holds.
#[derive(Debug, Clone)]
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

/// A fast hash of cells, whose coordinates are small integers.
#[derive(Default)]
pub(crate) struct CellHasher(u64);

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

#[cfg(test)]
mod tests;
