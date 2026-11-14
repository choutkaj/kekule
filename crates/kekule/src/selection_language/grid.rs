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
mod tests;
