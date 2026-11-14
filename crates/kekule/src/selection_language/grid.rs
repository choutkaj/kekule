//! Neighbor search for `within`, over the shared point grid.

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};
use std::sync::{Arc, Mutex, PoisonError};

use crate::math::{bounding_box, Cell, CellHasher, Point3, PointGrid};

/// A grid measured once searches around candidates, rather than gathering
/// lists for them, when it has at most this many candidates per distinct
/// point.
const SEARCH_CANDIDATES_PER_POINT: usize = 4;

/// Reference points binned into cubic cells at least as wide as `cutoff`, so
/// each candidate inside the references' cutoff-expanded bounding box
/// examines only the 27 cells around its own; candidates outside the box are
/// rejected at once.
///
/// The grid keeps a copy of each distinct point: repeated points, as of
/// copies stacked in one place, are measured once. A grid shared by many
/// measurements also keeps the near lists it gathers, so measurements after
/// the first gather only the cells no earlier one visited.
#[derive(Debug, Default)]
pub(crate) struct Grid {
    low: Point3,
    high: Point3,
    cutoff: f64,
    /// The distinct points.
    points: PointGrid<()>,
    /// The near lists kept across measurements, for a shared grid.
    shared: Option<SharedNear>,
    /// Whether a candidate that misses the point last matched searches the
    /// cells around its own rather than gathering their list: a grid
    /// measured once by few candidates for its points uses each list little.
    search: bool,
}

/// The exact bits of a point, so equal coordinates compare equal.
fn bits(point: &Point3) -> (u64, u64, u64) {
    (point.x.to_bits(), point.y.to_bits(), point.z.to_bits())
}

impl Grid {
    /// Grids `references` for measuring `candidates` points once; a grid
    /// shared by many measurements passes zero.
    pub(crate) fn new(references: &[Point3], cutoff: f64, candidates: usize) -> Self {
        if references.is_empty() {
            return Self::default();
        }
        let (low, high) = bounding_box(references.iter().copied());
        let low = Point3::new(low.x - cutoff, low.y - cutoff, low.z - cutoff);
        let high = Point3::new(high.x + cutoff, high.y + cutoff, high.z + cutoff);
        // Points, as of many structures, often repeat: the grid keeps each
        // once.
        let mut seen = HashSet::with_hasher(BuildHasherDefault::<CellHasher>::default());
        let mut points = PointGrid::default();
        points.build(
            references
                .iter()
                .filter(|point| seen.insert(bits(point)))
                .map(|&point| (point, ())),
            (low, high),
            cutoff,
            candidates,
        );
        let search = candidates > 0
            && candidates <= points.len().saturating_mul(SEARCH_CANDIDATES_PER_POINT);
        Self {
            low,
            high,
            cutoff,
            points,
            shared: (candidates == 0).then(SharedNear::default),
            search,
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
            let Some(cell) = self.points.cell(&point) else {
                continue;
            };
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
        let mut found = None;
        self.points.near_cell(cell, |other| {
            if stop(other) {
                found = Some(*other);
                return true;
            }
            near.push(*other);
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

/// The points a shared grid keeps across its near lists, each list counting
/// one more so that empty lists count too, which bounds its memory; lists
/// gathered after that serve their own measurement alone. Unit tests reach
/// it often.
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
    /// The points of `lists`, each list counting one more, so empty lists
    /// count too.
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
        let full = kept.points + list.len() + 1 > SHARED_NEAR_POINTS;
        match kept.lists.entry(cell) {
            Entry::Occupied(earlier) => Arc::clone(earlier.get()),
            Entry::Vacant(_) if full => list,
            Entry::Vacant(entry) => {
                entry.insert(Arc::clone(&list));
                kept.points += list.len() + 1;
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
    /// The points of `kept`, each list counting one more, which bounds those
    /// of lists the grid does not keep.
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
        self.points += list.len() + 1;
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

#[cfg(test)]
mod tests {
    use super::Grid;
    use crate::math::Point3;

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
        let dense =
            |references: &[Point3], cutoff| Grid::new(references, cutoff, 0).points.is_dense();
        assert!(dense(&references[..17], 1.3));
        assert!(!dense(&references[..17], 0.01));
        assert!(!dense(&references, 1.3));
        // Grids measured once are dense for boxes up to their candidates'.
        assert!(Grid::new(&references[..17], 0.01, 1_000_000)
            .points
            .is_dense());
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
    fn shared_grids_bound_the_empty_lists_they_keep() {
        // Two far points make a box of empty cells; each measurement visits
        // cells no other did, whose lists are empty.
        let references = [Point3::new(0.0, 0.0, 0.0), Point3::new(40.0, 40.0, 40.0)];
        let grid = Grid::new(&references, 1.0, 0);
        for structure in 0..100 {
            let candidates = (0..20)
                .map(|index| {
                    let cell = structure * 20 + index;
                    Point3::new(
                        5.5 + f64::from(cell % 30),
                        5.5 + f64::from(cell / 30 % 30),
                        5.5 + f64::from(cell / 900),
                    )
                })
                .collect::<Vec<_>>();
            let mut marked = Vec::new();
            grid.mark(candidates.iter().copied().enumerate(), |index| {
                marked.push(index)
            });
            assert!(marked.is_empty());
        }
        let kept = grid.shared.as_ref().unwrap().kept.lock().unwrap();
        assert!(!kept.lists.is_empty());
        assert!(kept.lists.len() <= super::SHARED_NEAR_POINTS);
        assert!(kept.points <= super::SHARED_NEAR_POINTS);
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
