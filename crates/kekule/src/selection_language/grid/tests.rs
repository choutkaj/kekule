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
    let dense = |references: &[Point3], cutoff| Grid::new(references, cutoff, 0).points.is_dense();
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
