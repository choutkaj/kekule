//! Uniform-grid neighbor search for `within`.

use std::collections::HashMap;

use crate::geometry::Point3;

/// Marks every candidate within `cutoff` (inclusive, Cartesian) of any
/// reference point.
///
/// References are hashed into cubic cells of edge `cutoff`, so each candidate
/// inside the references' cutoff-expanded bounding box examines only the 27
/// cells around its own; candidates outside the box are rejected at once.
pub(crate) fn within(
    candidates: &[Point3],
    references: &[Point3],
    cutoff: f64,
    mut mark: impl FnMut(usize),
) {
    if references.is_empty() {
        return;
    }
    let mut low = references[0];
    let mut high = references[0];
    for point in references {
        low = Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z));
        high = Point3::new(
            high.x.max(point.x),
            high.y.max(point.y),
            high.z.max(point.z),
        );
    }
    low = Point3::new(low.x - cutoff, low.y - cutoff, low.z - cutoff);
    high = Point3::new(high.x + cutoff, high.y + cutoff, high.z + cutoff);
    // Cells slightly wider than the cutoff keep every neighbor within the
    // adjacent cells despite rounding. A zero cutoff still needs a positive
    // cell; any edge finds coincident points in the candidate's own cell.
    let edge = if cutoff > 0.0 {
        cutoff * (1.0 + 1e-6)
    } else {
        1.0
    };
    let cell = |point: &Point3| {
        let axis = |value: f64, origin: f64| ((value - origin) / edge).floor() as i64;
        (
            axis(point.x, low.x),
            axis(point.y, low.y),
            axis(point.z, low.z),
        )
    };
    let mut cells: HashMap<(i64, i64, i64), Vec<u32>> = HashMap::new();
    for (index, point) in references.iter().enumerate() {
        cells.entry(cell(point)).or_default().push(index as u32);
    }
    let squared = cutoff * cutoff;
    for (index, point) in candidates.iter().enumerate() {
        if point.x < low.x
            || point.y < low.y
            || point.z < low.z
            || point.x > high.x
            || point.y > high.y
            || point.z > high.z
        {
            continue;
        }
        let (x, y, z) = cell(point);
        let near = (-1..=1).any(|dx| {
            (-1..=1).any(|dy| {
                (-1..=1).any(|dz| {
                    let key = (
                        x.saturating_add(dx),
                        y.saturating_add(dy),
                        z.saturating_add(dz),
                    );
                    cells.get(&key).is_some_and(|members| {
                        members.iter().any(|&member| {
                            let other = references[member as usize];
                            let (ex, ey, ez) =
                                (point.x - other.x, point.y - other.y, point.z - other.z);
                            ex * ex + ey * ey + ez * ez <= squared
                        })
                    })
                })
            })
        });
        if near {
            mark(index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::within;
    use crate::geometry::Point3;

    #[test]
    fn grid_agrees_with_brute_force() {
        let points = (0..400)
            .map(|index| {
                let index = f64::from(index);
                Point3::new(
                    (index * 0.37).sin() * 3.0,
                    (index * 0.11).cos() * 3.0,
                    index * 0.01,
                )
            })
            .collect::<Vec<_>>();
        let references = &points[..17];
        for cutoff in [0.0, 0.05, 0.4, 1.3, 10.0] {
            let mut grid = Vec::new();
            within(&points, references, cutoff, |index| grid.push(index));
            let brute = points
                .iter()
                .enumerate()
                .filter(|(_, point)| {
                    references.iter().any(|other| {
                        let (x, y, z) = (point.x - other.x, point.y - other.y, point.z - other.z);
                        x * x + y * y + z * z <= cutoff * cutoff
                    })
                })
                .map(|(index, _)| index)
                .collect::<Vec<_>>();
            assert_eq!(grid, brute, "cutoff {cutoff}");
        }
    }
}
