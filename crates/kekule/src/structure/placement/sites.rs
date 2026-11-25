//! Ideal substituent directions around one atom: electron-domain templates
//! oriented on the atom's existing neighbors, staggered against the
//! substituents of a single neighbor, and fitted to all its neighbors for
//! settling.
use super::AtomShape;
use crate::math::{Point3, Vector3};

/// The tetrahedral angle, `acos(-1/3)`.
const TETRAHEDRAL: f64 = 1.910_633_236_249_018_5;
/// Half the tetrahedral angle, between either of two new tetrahedral bonds
/// and their bisector: `acos(1/sqrt(3))`.
const HALF_TETRAHEDRAL: f64 = 0.955_316_618_124_509_3;
/// Below this length a direction is undefined.
const DEGENERATE: f64 = 1.0e-6;

/// What the site search needs to know about one atom.
#[derive(Debug, Clone, Copy)]
pub(crate) struct SiteQuery<'a> {
    pub(crate) center: Point3,
    /// Existing neighbors, in a deterministic order.
    pub(crate) neighbors: &'a [Point3],
    /// With exactly one neighbor: that neighbor and the substituent of it the
    /// first new site is placed anti to.
    pub(crate) reference: Option<Point3>,
    pub(crate) shape: AtomShape,
    /// Lone pairs that occupy sites of the shape.
    pub(crate) lone_pairs: usize,
}

/// The neighbors coincide with the atom or with each other's directions so
/// that no orientation follows from them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct DegenerateSites;

/// The free substituent directions of an atom, unit vectors in order of
/// preference: the site anti to the reference first. Sites taken by
/// neighbors or lone pairs are excluded.
pub(crate) fn free_sites(query: SiteQuery<'_>) -> Result<Vec<Vector3>, DegenerateSites> {
    let domains = query.shape.domains();
    let directions = query
        .neighbors
        .iter()
        .map(|neighbor| unit(*neighbor - query.center).ok_or(DegenerateSites))
        .collect::<Result<Vec<_>, _>>()?;
    let k = directions.len();
    if k + query.lone_pairs >= domains {
        return Ok(Vec::new());
    }
    let wanted = domains - k - query.lone_pairs;
    let template = template(query.shape);
    let (sites, occupied) = match (k, query.shape) {
        (0, _) => (template, Vec::new()),
        (1, _) => {
            let z = directions[0];
            let x = query
                .reference
                .and_then(|reference| {
                    let along = reference - query.neighbors[0];
                    unit(along - z * along.dot(z)).map(|p| -p)
                })
                .unwrap_or_else(|| perpendicular(z));
            let y = z.cross(x);
            let sites = template
                .iter()
                .map(|site| x * site.x + y * site.y + z * site.z)
                .collect();
            (sites, vec![0])
        }
        // Opposite neighbors, as of an atom bent from linear, leave the
        // plane of the new sites to a deterministic perpendicular.
        (2, AtomShape::TrigonalPlanar) => {
            let site = unit(-(directions[0] + directions[1]))
                .unwrap_or_else(|| perpendicular(directions[0]));
            return Ok(vec![site]);
        }
        (2, AtomShape::Tetrahedral) => {
            let bisector = unit(-(directions[0] + directions[1]))
                .unwrap_or_else(|| perpendicular(directions[0]));
            let normal = unit(directions[0].cross(directions[1])).unwrap_or_else(|| {
                unit(directions[0].cross(bisector)).unwrap_or_else(|| perpendicular(bisector))
            });
            let (sin, cos) = HALF_TETRAHEDRAL.sin_cos();
            let sites = [bisector * cos + normal * sin, bisector * cos - normal * sin];
            return Ok(sites[..wanted].to_vec());
        }
        (3, AtomShape::Tetrahedral) => {
            // Equiangular to the three neighbors: normal to the plane through
            // the tips of their directions, away from them.
            let away = -(directions[0] + directions[1] + directions[2]);
            let site =
                match unit((directions[1] - directions[0]).cross(directions[2] - directions[0])) {
                    Some(normal) if normal.dot(away) < 0.0 => -normal,
                    Some(normal) => normal,
                    None => unit(away).ok_or(DegenerateSites)?,
                };
            return Ok(vec![site]);
        }
        _ => fit_template(&template, &directions)?,
    };
    let mut free = (0..domains)
        .filter(|index| !occupied.contains(index))
        .collect::<Vec<_>>();
    let lone_pairs = lone_pair_sites(query.shape, &free, &occupied, query.lone_pairs);
    free.retain(|index| !lone_pairs.contains(index));
    Ok(free
        .into_iter()
        .take(wanted)
        .map(|index| sites[index])
        .collect())
}

/// The shape's sites in a frame whose first site is `+z`; with one neighbor
/// on `+z`, the first free site has a positive `x` component, which the
/// frame turns anti to the reference.
fn template(shape: AtomShape) -> Vec<Vector3> {
    let z = Vector3::new(0.0, 0.0, 1.0);
    match shape {
        AtomShape::Terminal => vec![z],
        AtomShape::Linear => vec![z, -z],
        AtomShape::TrigonalPlanar => {
            let (sin, cos) = (2.0 * std::f64::consts::FRAC_PI_3).sin_cos();
            vec![z, Vector3::new(sin, 0.0, cos), Vector3::new(-sin, 0.0, cos)]
        }
        AtomShape::Tetrahedral => {
            let (sin, cos) = TETRAHEDRAL.sin_cos();
            let mut sites = vec![z];
            for turn in 0..3 {
                let (sin_phi, cos_phi) =
                    (f64::from(turn) * 2.0 * std::f64::consts::FRAC_PI_3).sin_cos();
                sites.push(Vector3::new(sin * cos_phi, sin * sin_phi, cos));
            }
            sites
        }
        AtomShape::TrigonalBipyramidal => {
            let mut sites = vec![z, -z];
            for turn in 0..3 {
                let (sin_phi, cos_phi) =
                    (f64::from(turn) * 2.0 * std::f64::consts::FRAC_PI_3).sin_cos();
                sites.push(Vector3::new(cos_phi, sin_phi, 0.0));
            }
            sites
        }
        AtomShape::Octahedral => vec![
            z,
            -z,
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(-1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(0.0, -1.0, 0.0),
        ],
    }
}

/// The sites of `free` that lone pairs take: equatorial sites of a trigonal
/// bipyramid, sites trans to each other or to a neighbor of an octahedron,
/// and otherwise the last free sites.
fn lone_pair_sites(
    shape: AtomShape,
    free: &[usize],
    occupied: &[usize],
    count: usize,
) -> Vec<usize> {
    let mut chosen = Vec::with_capacity(count);
    match shape {
        AtomShape::TrigonalBipyramidal => {
            let order = free
                .iter()
                .rev()
                .filter(|index| **index >= 2)
                .chain(free.iter().rev().filter(|index| **index < 2));
            chosen.extend(order.take(count).copied());
        }
        AtomShape::Octahedral => {
            let trans = |index: usize| index ^ 1;
            while chosen.len() < count {
                let open = |index: &usize| !chosen.contains(index);
                let pick = free
                    .iter()
                    .rev()
                    .copied()
                    .filter(open)
                    .find(|index| chosen.contains(&trans(*index)))
                    .or_else(|| {
                        free.iter()
                            .rev()
                            .copied()
                            .filter(open)
                            .find(|index| occupied.contains(&trans(*index)))
                    })
                    .or_else(|| free.iter().rev().copied().find(open));
                match pick {
                    Some(index) => chosen.push(index),
                    None => break,
                }
            }
        }
        _ => chosen.extend(free.iter().rev().take(count).copied()),
    }
    chosen
}

/// Orients `template` on two or more neighbor directions: one neighbor on a
/// template site exactly, a second as close as the template allows, the
/// rest on their nearest remaining sites, choosing the orientation with the
/// smallest total deviation. Returns the oriented sites and the indices the
/// neighbors take.
fn fit_template(
    template: &[Vector3],
    directions: &[Vector3],
) -> Result<(Vec<Vector3>, Vec<usize>), DegenerateSites> {
    let first = directions[0];
    let second = directions[1..]
        .iter()
        .copied()
        .find(|direction| first.cross(*direction).norm() > 1.0e-3);
    let mut best: Option<(f64, Vec<Vector3>, Vec<usize>)> = None;
    for a in 0..template.len() {
        let orientations: Vec<[Vector3; 3]> = match second {
            Some(second) => (0..template.len())
                .filter(|b| *b != a)
                .filter_map(|b| frame_rotation(template[a], template[b], first, second))
                .collect(),
            None => {
                let target = perpendicular(first);
                let source = perpendicular(template[a]);
                frame_rotation(template[a], source, first, target)
                    .into_iter()
                    .collect()
            }
        };
        for rotation in orientations {
            let sites = template
                .iter()
                .map(|site| rotation[0] * site.x + rotation[1] * site.y + rotation[2] * site.z)
                .collect::<Vec<_>>();
            let mut taken = Vec::with_capacity(directions.len());
            let mut deviation = 0.0;
            for direction in directions {
                let Some((index, dot)) = sites
                    .iter()
                    .enumerate()
                    .filter(|(index, _)| !taken.contains(index))
                    .map(|(index, site)| (index, site.dot(*direction)))
                    .max_by(|left, right| left.1.total_cmp(&right.1))
                else {
                    break;
                };
                taken.push(index);
                deviation += 1.0 - dot;
            }
            if taken.len() == directions.len()
                && best
                    .as_ref()
                    .is_none_or(|(score, _, _)| deviation < *score - 1.0e-9)
            {
                best = Some((deviation, sites, taken));
            }
        }
    }
    best.map(|(_, sites, taken)| (sites, taken))
        .ok_or(DegenerateSites)
}

/// A held neighbor's deviation from its site weighs this much more than a
/// turning neighbor's in a shape fit.
const HELD_WEIGHT: f64 = 1.0e3;

/// A shape oriented on an atom's neighbors.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ShapeFit {
    /// The sites that bonds take, unit vectors: the shape's sites less those
    /// of its lone pairs.
    pub(crate) sites: Vec<Vector3>,
    /// The site of each neighbor direction, in their order.
    pub(crate) taken: Vec<usize>,
}

impl ShapeFit {
    /// The sites no neighbor takes, in order of preference.
    pub(crate) fn free(&self) -> Vec<Vector3> {
        (0..self.sites.len())
            .filter(|index| !self.taken.contains(index))
            .map(|index| self.sites[index])
            .collect()
    }
}

/// Orients the sites of `shape` that bonds take beside `lone_pairs` (lone
/// pairs equatorial on a trigonal bipyramid, trans on an octahedron) on the
/// unit neighbor `directions`, the first `held` of which do not move. The
/// first lies on a site exactly, and any second direction not along it
/// lies in the plane of its site; when both are held, the orientation may
/// instead split their deviations alike. The other directions take the
/// nearest free sites. The orientation with the smallest deviation wins,
/// held directions weighing most, the first found on a tie.
pub(crate) fn fit_shape(
    shape: AtomShape,
    lone_pairs: usize,
    directions: &[Vector3],
    held: usize,
) -> Result<ShapeFit, DegenerateSites> {
    let template = template(shape);
    let all = (0..template.len()).collect::<Vec<_>>();
    let lone = lone_pair_sites(shape, &all, &[], lone_pairs);
    let sites = all
        .into_iter()
        .filter(|index| !lone.contains(index))
        .map(|index| template[index])
        .collect::<Vec<_>>();
    let Some(&first) = directions.first() else {
        return Err(DegenerateSites);
    };
    if directions.len() > sites.len() {
        return Err(DegenerateSites);
    }
    let weight = |index: usize| if index < held { HELD_WEIGHT } else { 1.0 };
    let mut best: Option<(f64, ShapeFit)> = None;
    for a in 0..sites.len() {
        // A second direction not along the first lies in the plane of its
        // site; when both are held, the frame may instead split their
        // deviations alike. Every direction may be that second, so nearly
        // opposite neighbors leave the frame to another one.
        let mut frames = Vec::new();
        for (index, other) in directions.iter().enumerate().skip(1) {
            if first.cross(*other).norm() <= 1.0e-3 {
                continue;
            }
            for b in (0..sites.len()).filter(|b| *b != a) {
                frames.extend(
                    frame_rotation(sites[a], sites[b], first, *other)
                        .map(|frame| (frame, Some((index, b)))),
                );
                if index < held {
                    if let (Some(from), Some(to)) =
                        (unit(sites[a] + sites[b]), unit(first + *other))
                    {
                        frames.extend(
                            frame_rotation(from, sites[a], to, first)
                                .map(|frame| (frame, Some((index, b)))),
                        );
                    }
                }
            }
        }
        if frames.is_empty() {
            // All directions lie along the first, or their sites along its.
            frames.extend(
                frame_rotation(
                    sites[a],
                    perpendicular(sites[a]),
                    first,
                    perpendicular(first),
                )
                .map(|frame| (frame, None)),
            );
        }
        for (frame, second) in frames {
            let oriented = sites
                .iter()
                .map(|site| frame[0] * site.x + frame[1] * site.y + frame[2] * site.z)
                .collect::<Vec<_>>();
            let mut taken = vec![usize::MAX; directions.len()];
            let mut used = vec![false; sites.len()];
            taken[0] = a;
            used[a] = true;
            if let Some((index, b)) = second {
                taken[index] = b;
                used[b] = true;
            }
            // The rest take the free sites nearest them, best matches first.
            let mut pairs = Vec::new();
            for (index, direction) in directions.iter().enumerate() {
                if taken[index] != usize::MAX {
                    continue;
                }
                for (site_index, site) in oriented.iter().enumerate() {
                    if !used[site_index] {
                        pairs.push((site.dot(*direction), index, site_index));
                    }
                }
            }
            pairs.sort_by(|x, y| y.0.total_cmp(&x.0).then((x.1, x.2).cmp(&(y.1, y.2))));
            for (_, index, site_index) in pairs {
                if taken[index] == usize::MAX && !used[site_index] {
                    taken[index] = site_index;
                    used[site_index] = true;
                }
            }
            let deviation = directions
                .iter()
                .zip(&taken)
                .enumerate()
                .map(|(index, (direction, site))| {
                    weight(index) * (1.0 - oriented[*site].dot(*direction))
                })
                .sum::<f64>();
            if best
                .as_ref()
                .is_none_or(|(score, _)| deviation < *score - 1.0e-9)
            {
                best = Some((
                    deviation,
                    ShapeFit {
                        sites: oriented,
                        taken,
                    },
                ));
            }
        }
    }
    best.map(|(_, fit)| fit).ok_or(DegenerateSites)
}

/// The smallest angle between two domains of `shape`, in degrees.
pub(crate) fn smallest_angle(shape: AtomShape) -> f64 {
    match shape {
        AtomShape::Terminal | AtomShape::Linear => 180.0,
        AtomShape::TrigonalPlanar => 120.0,
        AtomShape::Tetrahedral => TETRAHEDRAL.to_degrees(),
        AtomShape::TrigonalBipyramidal | AtomShape::Octahedral => 90.0,
    }
}

/// The rotation, as the images of the unit axes, that takes `from` onto `to`
/// and turns `from_plane` toward `to_plane` about it; `None` when either
/// pair is parallel.
fn frame_rotation(
    from: Vector3,
    from_plane: Vector3,
    to: Vector3,
    to_plane: Vector3,
) -> Option<[Vector3; 3]> {
    let e1 = unit(from)?;
    let e2 = unit(from_plane - e1 * from_plane.dot(e1))?;
    let e3 = e1.cross(e2);
    let f1 = unit(to)?;
    let f2 = unit(to_plane - f1 * to_plane.dot(f1))?;
    let f3 = f1.cross(f2);
    let axis = |component: fn(Vector3) -> f64| {
        f1 * component(e1) + f2 * component(e2) + f3 * component(e3)
    };
    Some([axis(|v| v.x), axis(|v| v.y), axis(|v| v.z)])
}

/// A deterministic unit vector perpendicular to the unit vector `v`.
pub(crate) fn perpendicular(v: Vector3) -> Vector3 {
    let axis = if v.x.abs() <= v.y.abs() && v.x.abs() <= v.z.abs() {
        Vector3::new(1.0, 0.0, 0.0)
    } else if v.y.abs() <= v.z.abs() {
        Vector3::new(0.0, 1.0, 0.0)
    } else {
        Vector3::new(0.0, 0.0, 1.0)
    };
    unit(axis - v * axis.dot(v)).expect("an axis least aligned with a unit vector is not parallel")
}

/// A proper rotation, as the images of the unit axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Rotation([Vector3; 3]);

impl Rotation {
    pub(crate) const IDENTITY: Self = Self([
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    ]);

    /// The rotation by `angle` radians about the unit vector `axis`.
    pub(crate) fn about(axis: Vector3, angle: f64) -> Self {
        let (sin, cos) = angle.sin_cos();
        let turn = |v: Vector3| v * cos + axis.cross(v) * sin + axis * (axis.dot(v) * (1.0 - cos));
        Self([
            turn(Vector3::new(1.0, 0.0, 0.0)),
            turn(Vector3::new(0.0, 1.0, 0.0)),
            turn(Vector3::new(0.0, 0.0, 1.0)),
        ])
    }

    /// The smallest rotation taking the unit vector `from` onto the unit
    /// vector `to`.
    pub(crate) fn between(from: Vector3, to: Vector3) -> Self {
        let axis = from.cross(to);
        let (sin, cos) = (axis.norm(), from.dot(to));
        match unit(axis) {
            Some(axis) if sin > 1.0e-9 => Self::about(axis, sin.atan2(cos)),
            _ if cos > 0.0 => Self::IDENTITY,
            _ => Self::about(perpendicular(from), std::f64::consts::PI),
        }
    }

    pub(crate) fn apply(&self, v: Vector3) -> Vector3 {
        self.0[0] * v.x + self.0[1] * v.y + self.0[2] * v.z
    }

    /// This rotation applied after `first`.
    pub(crate) fn after(&self, first: &Self) -> Self {
        Self(first.0.map(|axis| self.apply(axis)))
    }
}

pub(crate) fn unit(v: Vector3) -> Option<Vector3> {
    let norm = v.norm();
    (norm.is_finite() && norm > DEGENERATE).then(|| v / norm)
}

#[cfg(test)]
mod tests;
