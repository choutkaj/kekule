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

#[test]
fn a_fitted_seesaw_keeps_its_lone_pair_equatorial() {
    // Three tetrahedral directions on a trigonal bipyramid with one lone
    // pair: the first stays, and the four bond sites hold one axial pair.
    let tetrahedron = free_sites(query(AtomShape::Tetrahedral, &[], 0)).unwrap();
    let fit = fit_shape(AtomShape::TrigonalBipyramidal, 1, &tetrahedron[..3], 1).unwrap();
    assert_eq!(fit.sites.len(), 4);
    assert!(angle(fit.sites[fit.taken[0]], tetrahedron[0]) < 1e-9);
    let free = fit.free();
    assert_eq!(free.len(), 1);
    let straight = fit
        .sites
        .iter()
        .enumerate()
        .flat_map(|(i, a)| fit.sites[i + 1..].iter().map(move |b| angle(*a, *b)))
        .filter(|found| (found - 180.0).abs() < 1e-6)
        .count();
    assert_eq!(straight, 1);
}

#[test]
fn two_held_directions_share_their_deviation() {
    // Two ring bonds at 110° on a seesaw take its equatorial pair, each 5°
    // off.
    let (sin, cos) = 55.0_f64.to_radians().sin_cos();
    let directions = [Vector3::new(sin, 0.0, cos), Vector3::new(-sin, 0.0, cos)];
    let fit = fit_shape(AtomShape::TrigonalBipyramidal, 1, &directions, 2).unwrap();
    for (direction, site) in directions.iter().zip(&fit.taken) {
        assert!((angle(*direction, fit.sites[*site]) - 5.0).abs() < 1e-6);
    }
    assert!((angle(fit.sites[fit.taken[0]], fit.sites[fit.taken[1]]) - 120.0).abs() < 1e-6);
    assert_eq!(
        fit_shape(
            AtomShape::Tetrahedral,
            1,
            &free_sites(query(AtomShape::Tetrahedral, &[], 0)).unwrap(),
            1
        ),
        Err(DegenerateSites)
    );
}

#[test]
fn nearly_opposite_neighbors_leave_the_frame_to_another() {
    // A trigonal bipyramid whose axial pair is bent to 179°: every
    // neighbor stays within a degree of its site.
    let (sin, cos) = 0.5_f64.to_radians().sin_cos();
    let mut directions = vec![Vector3::new(sin, 0.0, cos), Vector3::new(sin, 0.0, -cos)];
    for turn in 0..3 {
        let (sin, cos) = (f64::from(turn) * 120.0_f64.to_radians()).sin_cos();
        directions.push(Vector3::new(cos, sin, 0.0));
    }
    let fit = fit_shape(AtomShape::TrigonalBipyramidal, 0, &directions, 1).unwrap();
    for (direction, site) in directions.iter().zip(&fit.taken) {
        assert!(angle(*direction, fit.sites[*site]) < 1.0 + 1e-6);
    }
    // With three lone pairs, the same axial pair alone still fits.
    let fit = fit_shape(AtomShape::TrigonalBipyramidal, 3, &directions[..2], 1).unwrap();
    for (direction, site) in directions.iter().zip(&fit.taken) {
        assert!(angle(*direction, fit.sites[*site]) < 1.0 + 1e-6);
    }
}

#[test]
fn clearest_sites_leave_planar_neighbors_alike() {
    // Three neighbors in a plane: one site goes normal to it.
    let plane = free_sites(query(AtomShape::TrigonalPlanar, &[], 0)).unwrap();
    let start = [unit(plane[0] + Vector3::new(0.0, 1.0, 0.0) * 0.7).unwrap()];
    assert!(clearance(&plane, &start) < 60.0);
    let clearest = clearest_sites(&plane, &start);
    assert!((clearance(&plane, &clearest) - 90.0).abs() < 1e-3);
    // Two at 120°: two sites straddle their plane alike.
    let start = [Vector3::new(0.3, 1.0, -1.0), Vector3::new(0.1, -1.0, -0.6)]
        .map(|site| unit(site).unwrap());
    let clearest = clearest_sites(&plane[1..], &start);
    assert!((clearest[0].y + clearest[1].y).abs() < 1e-3);
    assert!(clearance(&plane[1..], &clearest) > 106.0);
}
