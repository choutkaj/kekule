//! Ideal substituent directions around one atom: electron-domain templates
//! oriented on the atom's existing neighbors, staggered against the
//! substituents of a single neighbor.
use super::AtomShape;
use crate::geometry::{Point3, Vector3};

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
        (2, AtomShape::TrigonalPlanar) => {
            let site = unit(-(directions[0] + directions[1])).ok_or(DegenerateSites)?;
            return Ok(vec![site]);
        }
        (2, AtomShape::Tetrahedral) => {
            let bisector = unit(-(directions[0] + directions[1])).ok_or(DegenerateSites)?;
            let normal = unit(directions[0].cross(directions[1])).ok_or(DegenerateSites)?;
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
fn perpendicular(v: Vector3) -> Vector3 {
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
mod tests {
    use super::*;

    fn angle(a: Vector3, b: Vector3) -> f64 {
        a.dot(b).clamp(-1.0, 1.0).acos().to_degrees()
    }

    fn query(shape: AtomShape, neighbors: &[Point3], lone_pairs: usize) -> SiteQuery<'_> {
        SiteQuery {
            center: Point3::origin(),
            neighbors,
            reference: None,
            shape,
            lone_pairs,
        }
    }

    #[test]
    fn isolated_shapes_have_ideal_angles() {
        for (shape, ideal) in [
            (AtomShape::Linear, vec![180.0]),
            (AtomShape::TrigonalPlanar, vec![120.0]),
            (AtomShape::Tetrahedral, vec![109.471_220_634]),
            (AtomShape::TrigonalBipyramidal, vec![90.0, 120.0, 180.0]),
            (AtomShape::Octahedral, vec![90.0, 180.0]),
        ] {
            let sites = free_sites(query(shape, &[], 0)).unwrap();
            assert_eq!(sites.len(), shape.domains());
            for (i, a) in sites.iter().enumerate() {
                assert!((a.norm() - 1.0).abs() < 1e-12);
                for b in &sites[i + 1..] {
                    let found = angle(*a, *b);
                    assert!(
                        ideal.iter().any(|ideal| (found - ideal).abs() < 1e-6),
                        "{shape:?}: {found}"
                    );
                }
            }
        }
    }

    #[test]
    fn one_neighbor_stagger_puts_the_first_site_anti() {
        let neighbor = Point3::new(0.0, 0.0, 0.153);
        let reference = Point3::new(0.144, 0.0, 0.204);
        let mut q = query(AtomShape::Tetrahedral, std::slice::from_ref(&neighbor), 0);
        q.reference = Some(reference);
        let sites = free_sites(q).unwrap();
        assert_eq!(sites.len(), 3);
        let bond = unit(neighbor - Point3::origin()).unwrap();
        for site in &sites {
            assert!((angle(*site, bond) - 109.471_220_634).abs() < 1e-6);
        }
        // Anti: the first site and the reference lie on opposite sides.
        assert!(sites[0].x < 0.0 && sites[0].y.abs() < 1e-12);
    }

    #[test]
    fn two_and_three_neighbors_fill_the_remaining_sites() {
        let sites = free_sites(query(AtomShape::Tetrahedral, &[], 0)).unwrap();
        let points = sites
            .iter()
            .map(|s| Point3::origin() + *s * 0.15)
            .collect::<Vec<_>>();
        let two = free_sites(query(AtomShape::Tetrahedral, &points[..2], 0)).unwrap();
        assert_eq!(two.len(), 2);
        for site in &two {
            assert!(sites[2..].iter().any(|ideal| angle(*ideal, *site) < 1e-6));
        }
        let three = free_sites(query(AtomShape::Tetrahedral, &points[..3], 0)).unwrap();
        assert!(angle(three[0], sites[3]) < 1e-6);
        let lone = free_sites(query(AtomShape::Tetrahedral, &points[..2], 1)).unwrap();
        assert_eq!(lone.len(), 1);
    }

    #[test]
    fn hypervalent_lone_pairs_take_equatorial_and_trans_sites() {
        // Sulfur tetrafluoride: a seesaw, lone pair equatorial.
        let sites = free_sites(query(AtomShape::TrigonalBipyramidal, &[], 1)).unwrap();
        assert_eq!(sites.len(), 4);
        let axial = sites.iter().filter(|s| s.z.abs() > 0.99).count();
        assert_eq!(axial, 2);
        // Xenon tetrafluoride: square planar, lone pairs trans.
        let sites = free_sites(query(AtomShape::Octahedral, &[], 2)).unwrap();
        assert_eq!(sites.len(), 4);
        let sum = sites.iter().fold(Vector3::zero(), |sum, s| sum + *s);
        assert!(sum.norm() < 1e-12);
        // Square pyramid around a fitted pair of neighbors.
        let neighbors = [Point3::new(0.17, 0.0, 0.0), Point3::new(0.0, 0.17, 0.0)];
        let sites = free_sites(query(AtomShape::Octahedral, &neighbors, 1)).unwrap();
        assert_eq!(sites.len(), 3);
        for site in &sites {
            let to_first = angle(*site, Vector3::new(1.0, 0.0, 0.0));
            assert!((to_first - 90.0).abs() < 1e-6 || (to_first - 180.0).abs() < 1e-6);
        }
    }

    #[test]
    fn rotations_align_and_turn() {
        let x = Vector3::new(1.0, 0.0, 0.0);
        let y = Vector3::new(0.0, 1.0, 0.0);
        let z = Vector3::new(0.0, 0.0, 1.0);
        assert!((Rotation::about(z, std::f64::consts::FRAC_PI_2).apply(x) - y).norm() < 1e-12);
        for (from, to) in [(x, y), (x, x), (x, -x), (unit(x + y + z).unwrap(), -z)] {
            assert!((Rotation::between(from, to).apply(from) - to).norm() < 1e-12);
        }
        let composed = Rotation::about(z, 1.0).after(&Rotation::about(z, 0.5));
        assert!((composed.apply(x) - Rotation::about(z, 1.5).apply(x)).norm() < 1e-12);
    }

    #[test]
    fn coincident_neighbors_are_degenerate() {
        let neighbors = [Point3::origin()];
        assert_eq!(
            free_sites(query(AtomShape::Tetrahedral, &neighbors, 0)),
            Err(DegenerateSites)
        );
    }
}
