use super::*;
use crate::math::{PeriodicCell, Vector3};
use crate::units::{Quantity, NANOMETER};

/// A deterministic linear congruential generator.
struct Lcg(u64);

impl Lcg {
    /// A value in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn within(cutoff: f64) -> impl Fn(&Point3, &Point3) -> bool {
    move |a, b| (*a - *b).norm() <= cutoff
}

/// The indexes `for_each_near` passes for `candidate`, in its order, of
/// the points `within` accepts.
fn near(
    grid: &PointGrid<usize>,
    candidate: &Point3,
    within: &impl Fn(&Point3, &Point3) -> bool,
) -> Vec<usize> {
    let mut found = Vec::new();
    grid.for_each_near(candidate, |point, index| {
        if within(candidate, point) {
            found.push(index);
        }
    });
    found
}

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
    grid.build(points.iter().copied().zip(0usize..), (low, high), 0.3, 0);
    let within = within(0.3);
    for candidate in &points {
        let mut found = near(&grid, candidate, &within);
        found.sort_unstable();
        let expected = (0..points.len())
            .filter(|&index| within(candidate, &points[index]))
            .collect::<Vec<_>>();
        assert_eq!(found, expected);
        assert!(grid.any_near(candidate, |point, _| within(candidate, point)));
    }
    let far = Point3::new(9.0, 9.0, 9.0);
    assert!(!grid.any_near(&far, |point, _| within(&far, point)));
}

/// Dense and sparse grids find the same points in the same order: cells
/// x-major, then each cell's points in input order. Points sit on cell
/// boundaries and exactly at the cutoff.
#[test]
fn dense_and_sparse_grids_agree_with_brute_force_in_grid_order() {
    let mut random = Lcg(7);
    for (round, cutoff) in [0.0_f64, 0.25, 0.5, 1.0, 3.0].into_iter().enumerate() {
        let mut points = (0..300)
            .map(|_| {
                Point3::new(
                    random.unit() * 4.0,
                    random.unit() * 4.0,
                    random.unit() * 2.0,
                )
            })
            .collect::<Vec<_>>();
        // Grid-aligned points, their copies, and partners at the cutoff.
        for step in 0..20 {
            let at = f64::from(step) * cutoff.max(0.125);
            points.push(Point3::new(at, at * 0.5, 0.0));
            points.push(Point3::new(at, at * 0.5, 0.0));
            points.push(Point3::new(at + cutoff, at * 0.5, 0.0));
        }
        let bounds = bounding_box(points.iter().copied());
        let within = within(cutoff);
        let mut dense = PointGrid::default();
        dense.build(points.iter().copied().zip(0..), bounds, cutoff, 1 << 24);
        // A zero cutoff bins by a 2^50th of the box, too fine to index
        // densely.
        assert_eq!(dense.is_dense(), cutoff > 0.0, "{round}");
        // A far point makes the box too large to index densely.
        let mut sparse_points = points.clone();
        sparse_points.push(Point3::new(1.0e7, -1.0e7, 1.0e7));
        let mut sparse = PointGrid::default();
        sparse.build(
            sparse_points.iter().copied().zip(0..),
            bounding_box(sparse_points.iter().copied()),
            cutoff,
            0,
        );
        assert!(!sparse.is_dense(), "{round}");
        for candidate in points.iter().chain(&[Point3::new(-0.5, 5.0, 1.0)]) {
            let found = near(&dense, candidate, &within);
            let mut expected = (0..points.len())
                .filter(|&index| within(candidate, &points[index]))
                .collect::<Vec<_>>();
            let cell = |index: usize| dense.cell(&points[index]).unwrap();
            expected.sort_by_key(|&index| (cell(index), index));
            assert_eq!(found, expected, "{round} {candidate:?}");
            // The far point never lies within the cutoff, so the sparse
            // grid's cells order the same points alike.
            let mut found = near(&sparse, candidate, &within);
            let cell = |index: usize| sparse.cell(&sparse_points[index]).unwrap();
            assert!(found
                .windows(2)
                .all(|pair| (cell(pair[0]), pair[0]) < (cell(pair[1]), pair[1])));
            found.sort_unstable();
            expected.sort_unstable();
            assert_eq!(found, expected, "{round} {candidate:?}");
        }
    }
}

/// A distant point keeps cells cutoff-sized: a candidate meets only its
/// neighbours, not every point the outlier's box would share a cell
/// with.
#[test]
fn a_distant_point_keeps_the_search_local() {
    let mut points = (0..10_000)
        .map(|i| Point3::new(f64::from(i) * 10.0, 0.0, 0.0))
        .collect::<Vec<_>>();
    points.push(Point3::new(1.0e15, 0.0, 0.0));
    let mut grid = PointGrid::default();
    grid.build(
        points.iter().copied().zip(0usize..),
        bounding_box(points.iter().copied()),
        5.0,
        points.len(),
    );
    for candidate in &points[..100] {
        let mut visited = 0;
        grid.for_each_near(candidate, |_, _| visited += 1);
        assert!(visited <= 5, "{visited}");
    }
}

/// Far from the box's corner, rounding of cell coordinates can carry a
/// contact past an adjacent cell; searches reach far enough to find it.
#[test]
fn far_contacts_are_found_despite_rounding() {
    let mut random = Lcg(11);
    for offset in [1.0e9, 1.0e12, 1.0e14] {
        let mut points = vec![Point3::origin()];
        for _ in 0..200 {
            let x = offset + (random.unit() * 64.0).floor() * 1.000_001;
            points.push(Point3::new(x, random.unit(), 0.0));
            points.push(Point3::new(x + 1.0, random.unit(), 0.0));
        }
        let within = within(1.0);
        let mut grid = PointGrid::default();
        grid.build(
            points.iter().copied().zip(0usize..),
            bounding_box(points.iter().copied()),
            1.0,
            0,
        );
        for candidate in &points {
            let mut found = near(&grid, candidate, &within);
            found.sort_unstable();
            let expected = (0..points.len())
                .filter(|&index| within(candidate, &points[index]))
                .collect::<Vec<_>>();
            assert_eq!(found, expected, "{offset} {candidate:?}");
        }
    }
}

/// A subnormal cutoff still bins every point at a finite cell.
#[test]
fn subnormal_cutoffs_find_their_contacts() {
    let points = [Point3::origin(), Point3::new(5.0e-321, 0.0, 0.0)];
    let mut grid = PointGrid::default();
    grid.build(
        points.iter().copied().zip(0usize..),
        bounding_box(points.iter().copied()),
        1.0e-320,
        0,
    );
    let within = within(1.0e-320);
    assert_eq!(near(&grid, &points[0], &within), [0, 1]);
    assert_eq!(near(&grid, &points[1], &within), [0, 1]);
}

#[test]
fn rebuilding_reuses_the_grid() {
    let mut grid = PointGrid::default();
    for (round, count) in [50usize, 0, 400, 3].into_iter().enumerate() {
        let points = (0..count)
            .map(|index| {
                let index = index as f64;
                Point3::new(index.sin(), index.cos() * 3.0, index * 1e-3)
            })
            .collect::<Vec<_>>();
        grid.build(
            points.iter().copied().zip(0..),
            bounding_box(points.iter().copied()),
            0.4,
            0,
        );
        assert_eq!(grid.len(), count);
        let within = within(0.4);
        for candidate in &points {
            let mut found = near(&grid, candidate, &within);
            found.sort_unstable();
            let expected = (0..count)
                .filter(|&index| within(candidate, &points[index]))
                .collect::<Vec<_>>();
            assert_eq!(found, expected, "{round}");
        }
    }
}

#[test]
fn near_cells_hold_every_point_within_the_cutoff_of_the_cell() {
    let mut random = Lcg(11);
    let points = (0..500)
        .map(|_| Point3::new(random.unit() * 6.0, random.unit() * 6.0, random.unit()))
        .collect::<Vec<_>>();
    for cutoff in [0.1, 0.7, 2.0] {
        let mut grid = PointGrid::default();
        grid.build(
            points.iter().map(|&point| (point, ())),
            bounding_box(points.iter().copied()),
            cutoff,
            0,
        );
        for candidate in &points {
            let mut near = Vec::new();
            grid.near_cell(grid.cell(candidate).unwrap(), |point| {
                near.push(*point);
                false
            });
            for point in &points {
                if (*candidate - *point).norm() <= cutoff {
                    assert!(near.contains(point), "{cutoff}");
                }
            }
        }
    }
}

fn lattices() -> Vec<Lattice> {
    [
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
    ]
    .into_iter()
    .map(|basis| {
        Lattice::new(PeriodicCell::new(Quantity::new(basis, NANOMETER), [true; 3]).unwrap())
            .unwrap()
    })
    .collect()
}

/// Periodic searches find every minimum-image neighbor, in the order of
/// the wrapped cells around the candidate's, each visited once, then of
/// insertion. Points sit on bin boundaries and outside the home cell.
#[test]
fn periodic_grids_agree_with_brute_force() {
    let mut random = Lcg(3);
    for lattice in lattices() {
        // Cutoffs from many bins down to a single bin per axis.
        for cutoff in [0.05, 0.4, 0.9, 2.5] {
            let mut grid = PointGrid::periodic(&lattice, cutoff).unwrap();
            let Space::Periodic { bins, .. } = grid.space.clone() else {
                unreachable!()
            };
            let mut points = Vec::new();
            for index in 0..120 {
                let fractional = if index % 4 == 0 {
                    // On a bin boundary, in another image of the cell.
                    std::array::from_fn(|axis| {
                        (index as f64 % bins[axis] as f64) / bins[axis] as f64 + (index % 3) as f64
                            - 1.0
                    })
                } else {
                    [
                        random.unit() * 3.0 - 1.0,
                        random.unit(),
                        random.unit() * 2.0 - 1.0,
                    ]
                };
                let point = Point3::origin() + lattice.cartesian_canonical(fractional);
                grid.insert(point, points.len()).unwrap();
                points.push(point);
            }
            assert_eq!(grid.len(), points.len());
            let distance =
                |a: Point3, b: Point3| lattice.minimum_image_canonical(a - b).unwrap().norm();
            for (index, &candidate) in points.iter().enumerate() {
                let query = if index % 2 == 0 {
                    candidate
                } else {
                    Point3::origin()
                        + lattice.cartesian_canonical([random.unit() * 4.0 - 2.0, 0.5, 1.0])
                };
                let mut visited = Vec::new();
                let found = grid
                    .try_any_near(&query, |point, index| {
                        assert_eq!(*point, points[index]);
                        visited.push(index);
                        Ok::<_, PeriodicCellError>(false)
                    })
                    .unwrap();
                assert!(!found);
                // Each point is visited at most once, cells in order.
                let mut unique = visited.clone();
                unique.sort_unstable();
                unique.dedup();
                assert_eq!(unique.len(), visited.len());
                let mut near = visited
                    .iter()
                    .copied()
                    .filter(|&other| distance(query, points[other]) < cutoff)
                    .collect::<Vec<_>>();
                near.sort_unstable();
                let expected = (0..points.len())
                    .filter(|&other| distance(query, points[other]) < cutoff)
                    .collect::<Vec<_>>();
                assert_eq!(near, expected, "{cutoff} {query:?}");
                // An early exit stops at the first match in that order.
                let first = visited
                    .iter()
                    .copied()
                    .find(|&other| distance(query, points[other]) < cutoff);
                let mut stopped = None;
                let any = grid
                    .try_any_near(&query, |point, index| {
                        let hit = distance(query, *point) < cutoff;
                        if hit {
                            stopped = Some(index);
                        }
                        Ok::<_, PeriodicCellError>(hit)
                    })
                    .unwrap();
                assert_eq!(any, first.is_some());
                assert_eq!(stopped, first);
            }
        }
    }
}

#[test]
fn periodic_grids_reject_points_they_cannot_wrap() {
    let lattice = lattices().remove(0);
    assert_eq!(
        PointGrid::<()>::periodic(&lattice, 0.0).unwrap_err(),
        PeriodicCellError::NumericalFailure
    );
    let mut grid = PointGrid::periodic(&lattice, 0.5).unwrap();
    let far = Point3::new(1.0e20, 0.0, 0.0);
    assert_eq!(
        grid.insert(far, ()).unwrap_err(),
        PeriodicCellError::NumericalFailure
    );
    assert!(grid.is_empty());
    grid.insert(Point3::origin(), ()).unwrap();
    assert_eq!(
        grid.try_any_near(&far, |_, ()| Ok::<_, PeriodicCellError>(true)),
        Err(PeriodicCellError::NumericalFailure)
    );
    assert!(!grid.any_near(&far, |_, ()| true));
    assert!(grid.any_near(&Point3::new(0.1, 0.0, 0.0), |_, ()| true));
}
