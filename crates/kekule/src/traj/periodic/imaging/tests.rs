use super::*;
use crate::units::NANOMETER;

#[test]
fn nearest_images_match_exhaustive_cartesian_search_for_all_axis_flags() {
    let base = [
        Vector3::new(1.0, 0.2, -0.1),
        Vector3::new(0.4, 1.2, 0.3),
        Vector3::new(0.2, -0.1, 0.8),
    ];
    for basis in [
        base,
        base.map(|v| Vector3::new(-v.y, v.x, v.z)),
        [base[1], base[0], base[2]],
    ] {
        for mask in 1..8 {
            let axes = std::array::from_fn(|axis| mask & (1 << axis) != 0);
            let cell = PeriodicCell::new(Quantity::new(basis, NANOMETER), axes).unwrap();
            let lattice = frame_lattice(Some(cell), 7).unwrap();
            for geometry in 0..19 {
                let delta = Vector3::new(
                    f64::from((geometry * 17) % 31 - 15) / 11.0,
                    f64::from((geometry * 11) % 29 - 14) / 9.0,
                    f64::from((geometry * 7) % 23 - 11) / 8.0,
                );
                let (nearest, images) = lattice.nearest_image_canonical(delta).unwrap();
                let mut exhaustive = f64::INFINITY;
                for a in -4..=4 {
                    for b in -4..=4 {
                        for c in -4..=4 {
                            if (!axes[0] && a != 0) || (!axes[1] && b != 0) || (!axes[2] && c != 0)
                            {
                                continue;
                            }
                            let residual = delta
                                - basis[0] * f64::from(a)
                                - basis[1] * f64::from(b)
                                - basis[2] * f64::from(c);
                            exhaustive = exhaustive.min(residual.norm_squared());
                        }
                    }
                }
                assert!(
                    (nearest.norm_squared() - exhaustive).abs() < 1.0e-12,
                    "basis={basis:?}, axes={axes:?}, delta={delta:?}, images={images:?}"
                );
                for axis in 0..3 {
                    assert!(axes[axis] || images[axis] == 0);
                }
                assert!(
                    (nearest - (delta - lattice.cartesian_canonical(images.map(|n| n as f64))))
                        .norm()
                        < 1.0e-12
                );
            }
        }
    }
}

#[test]
fn rectangular_cells_use_independent_images_even_with_extreme_aspect_ratios() {
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0e-6, 1.0e6), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let lattice = frame_lattice(Some(cell), 0).unwrap();
    let (nearest, images) = lattice
        .nearest_image_canonical(Vector3::new(1.25, 1.25e-6, 1.25e6))
        .unwrap();
    assert_eq!(images, [1; 3]);
    assert!((nearest.x - 0.25).abs() < 1.0e-14);
    assert!((nearest.y - 0.25e-6).abs() < 1.0e-20);
    assert!((nearest.z - 0.25e6).abs() < 1.0e-8);
}

#[test]
fn nearest_image_fails_explicitly_on_excessive_search_or_lost_integer_precision() {
    let cell = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.99999, 1.0e-6, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            NANOMETER,
        ),
        [true; 3],
    )
    .unwrap();
    let lattice = frame_lattice(Some(cell), 7).unwrap();
    let nearest_image = |delta| {
        lattice
            .nearest_image_canonical(delta)
            .map_err(|error| cell_error(error, 7))
    };
    assert!(matches!(
        nearest_image(Vector3::new(0.49, 0.49, 0.0)),
        Err(PeriodicError::LimitExceeded { frame: 7, limit })
            if limit.resource() == crate::diagnostics::Resource::LatticeImages
    ));
    assert!(matches!(
        nearest_image(Vector3::new(1.0e100, 0.0, 0.0)),
        Err(PeriodicError::NumericalFailure { frame: 7 })
    ));
}
