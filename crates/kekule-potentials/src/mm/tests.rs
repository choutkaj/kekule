use super::*;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Central differences of `energy` along every coordinate, Richardson-extrapolated.
fn numerical_gradient(
    points: &[Point3],
    h: f64,
    energy: impl Fn(&[Point3]) -> f64,
) -> Vec<[f64; 3]> {
    let estimate = |atom: usize, axis: usize, h: f64| {
        let shifted = |sign: f64| {
            let mut moved = points.to_vec();
            match axis {
                0 => moved[atom].x += sign * h,
                1 => moved[atom].y += sign * h,
                _ => moved[atom].z += sign * h,
            }
            energy(&moved)
        };
        (shifted(1.0) - shifted(-1.0)) / (2.0 * h)
    };
    (0..points.len())
        .map(|atom| {
            std::array::from_fn(|axis| {
                (4.0 * estimate(atom, axis, h / 2.0) - estimate(atom, axis, h)) / 3.0
            })
        })
        .collect()
}

fn assert_gradient(analytic: &[Vector3], numerical: &[[f64; 3]], tolerance: f64) {
    for (atom, (a, n)) in analytic.iter().zip(numerical).enumerate() {
        for (axis, (a, n)) in [a.x, a.y, a.z].into_iter().zip(n).enumerate() {
            assert!(
                (a - n).abs() <= tolerance * (1.0 + n.abs()),
                "atom {atom} axis {axis}: analytic {a} numerical {n}"
            );
        }
    }
}

fn total(gradients: Gradients) -> Vec<Vector3> {
    match gradients {
        Gradients::Total(g) => g,
        _ => unreachable!(),
    }
}

#[test]
fn signed_dihedral_and_phase_are_preserved() {
    let phi = dihedral([p(1., 0., 0.), p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.)]).unwrap();
    assert!((phi - std::f64::consts::FRAC_PI_2).abs() < 1e-14);
    let torsion = PeriodicTorsion {
        atoms: [0, 1, 2, 3],
        terms: Box::new([FourierTerm {
            k: 2.0,
            periodicity: 1.0,
            phase: std::f64::consts::FRAC_PI_2,
        }]),
    };
    let points = [p(1., 0., 0.), p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.)];
    let energy = torsions(
        0,
        TermKind::ProperTorsion,
        &[torsion],
        &points,
        &mut Gradients::None,
    )
    .unwrap();
    assert!((energy - 4.0).abs() < 1e-14);
    assert!(dihedral([p(0., 0., 0.); 4]).is_none());
    // An outer atom on the axis leaves the dihedral undefined.
    assert!(dihedral([p(0., 0., -1.), p(0., 0., 0.), p(0., 0., 1.), p(0., 1., 1.)]).is_none());
}

#[test]
fn analytic_gradients_differentiate_every_kernel() {
    let x = [
        p(0.01, -0.02, 0.03),
        p(0.15, 0.01, -0.01),
        p(0.21, 0.13, 0.02),
        p(0.35, 0.16, 0.11),
        p(0.05, 0.31, -0.12),
    ];
    let bond = [HarmonicBond {
        atoms: [0, 1],
        k: 3.0e5,
        length: 0.11,
    }];
    let angle = [HarmonicAngle {
        atoms: [0, 1, 2],
        k: 500.0,
        angle: 1.9,
    }];
    let torsion = [PeriodicTorsion {
        atoms: [0, 1, 2, 3],
        terms: Box::new([
            FourierTerm {
                k: 1.3,
                periodicity: 3.0,
                phase: 0.0,
            },
            FourierTerm {
                k: 0.7,
                periodicity: 2.0,
                phase: std::f64::consts::PI,
            },
        ]),
    }];
    let nonbonded = Nonbonded {
        sigma: vec![0.3, 0.25, 0.32, 0.28, 0.31],
        epsilon: vec![0.4, 0.1, 0.5, 0.2, 0.3],
        charge: vec![-0.4, 0.2, 0.1, -0.3, 0.4],
        exceptions: vec![
            vec![(1, 0.0, 0.0), (3, 0.5, 0.8333)],
            vec![],
            vec![],
            vec![],
            vec![],
        ],
    };
    type Kernel<'a> = Box<dyn Fn(&[Point3], &mut Gradients) -> f64 + 'a>;
    let kernels: [Kernel; 4] = [
        Box::new(|x, g| bonds(&bond, x, g).unwrap()),
        Box::new(|x, g| angles(&angle, x, g).unwrap()),
        Box::new(|x, g| torsions(0, TermKind::ProperTorsion, &torsion, x, g).unwrap()),
        Box::new(|x, g| {
            let (v, q) = nonbonded.evaluate(x, g).unwrap();
            v + q
        }),
    ];
    for kernel in kernels {
        let mut gradients = Gradients::total(x.len());
        kernel(&x, &mut gradients);
        let numerical = numerical_gradient(&x, 1e-5, |x| kernel(x, &mut Gradients::None));
        assert_gradient(&total(gradients), &numerical, 1e-6);
    }
}

#[test]
fn near_linear_angle_gradient_differentiates_its_energy() {
    // Mapped atoms 16-17-18 of PubChem CID 443915, first geometry of the
    // externally supplied OpenFF validation panel, with the reference angle
    // parameter: sin(theta) is about 5.2e-5 and theta0 about 112 degrees.
    let x = [
        p(
            -0.18678127896333585,
            0.046354342980343965,
            -0.038190366434239835,
        ),
        p(
            -0.06602942179699119,
            0.004475868142689026,
            -0.04519745169972588,
        ),
        p(
            0.05763445072849957,
            -0.038419437413943415,
            -0.052371549722386224,
        ),
    ];
    let angle = [HarmonicAngle {
        atoms: [0, 1, 2],
        k: 985.271616332088,
        angle: 1.953427004874672,
    }];
    let mut gradients = Gradients::total(3);
    angles(&angle, &x, &mut gradients).unwrap();
    let numerical = numerical_gradient(&x, 1e-7, |x| {
        angles(&angle, x, &mut Gradients::None).unwrap()
    });
    assert_gradient(&total(gradients), &numerical, 1e-6);
}

#[test]
fn undefined_quantities_are_rejected_only_when_requested() {
    let linear = [p(0., 0., 0.), p(0.1, 0., 0.), p(0.2, 0., 0.)];
    let angle = [HarmonicAngle {
        atoms: [0, 1, 2],
        k: 100.0,
        angle: 2.0,
    }];
    let energy = angles(&angle, &linear, &mut Gradients::None).unwrap();
    let expected = 0.5 * 100.0 * (std::f64::consts::PI - 2.0).powi(2);
    assert!((energy - expected).abs() < 1e-12);
    let error = angles(&angle, &linear, &mut Gradients::total(3)).unwrap_err();
    assert_eq!(error.kind, SingularGeometry::LinearAngle);
    // Short but nonzero arms still define the angle and its gradient.
    let short = [p(1e-13, 0., 0.), p(0., 0., 0.), p(0., 2e-13, 0.)];
    let mut gradients = Gradients::total(3);
    let energy = angles(&angle, &short, &mut gradients).unwrap();
    let expected = 0.5 * 100.0 * (std::f64::consts::FRAC_PI_2 - 2.0).powi(2);
    assert!((energy - expected).abs() < 1e-12);
    let numerical = numerical_gradient(&short, 1e-15, |x| {
        angles(&angle, x, &mut Gradients::None).unwrap()
    });
    assert_gradient(&total(gradients), &numerical, 1e-6);
    let zero_arm = [p(0., 0., 0.), p(0., 0., 0.), p(0.1, 0., 0.)];
    for mut gradients in [Gradients::None, Gradients::total(3)] {
        let error = angles(&angle, &zero_arm, &mut gradients).unwrap_err();
        assert_eq!(error.kind, SingularGeometry::DegenerateAngle);
    }
    // A torsion without amplitude is zero at any geometry, including a
    // collinear one where the dihedral is undefined; any amplitude needs it.
    let collinear = [
        p(0., 0., 0.),
        p(0.1, 0., 0.),
        p(0.22, 0., 0.),
        p(0.32, 0., 0.),
    ];
    let torsion = |k| {
        [PeriodicTorsion {
            atoms: [0, 1, 2, 3],
            terms: Box::new([
                FourierTerm {
                    k: 0.0,
                    periodicity: 1.0,
                    phase: 0.0,
                },
                FourierTerm {
                    k,
                    periodicity: 3.0,
                    phase: std::f64::consts::PI,
                },
            ]),
        }]
    };
    let energy = torsions(
        0,
        TermKind::ProperTorsion,
        &torsion(0.0),
        &collinear,
        &mut Gradients::None,
    )
    .unwrap();
    assert_eq!(energy, 0.0);
    let mut gradients = Gradients::total(4);
    let energy = torsions(
        0,
        TermKind::ProperTorsion,
        &torsion(0.0),
        &collinear,
        &mut gradients,
    )
    .unwrap();
    assert_eq!(energy, 0.0);
    assert_eq!(total(gradients), vec![Vector3::zero(); 4]);
    for mut gradients in [Gradients::None, Gradients::total(4)] {
        let error = torsions(
            0,
            TermKind::ProperTorsion,
            &torsion(1.0),
            &collinear,
            &mut gradients,
        )
        .unwrap_err();
        assert_eq!(
            (error.interaction, error.atoms, error.kind),
            (
                TermKind::ProperTorsion,
                vec![0, 1, 2, 3],
                SingularGeometry::DegenerateDihedral
            )
        );
    }
    // A linear equilibrium is stationary, so its zero gradient is defined.
    let straight = [HarmonicAngle {
        atoms: [0, 1, 2],
        k: 100.0,
        angle: std::f64::consts::PI,
    }];
    angles(&straight, &linear, &mut Gradients::total(3)).unwrap();

    let coincident = [p(0., 0., 0.), p(0., 0., 0.)];
    let bond = [HarmonicBond {
        atoms: [0, 1],
        k: 2.0,
        length: 0.1,
    }];
    assert!((bonds(&bond, &coincident, &mut Gradients::None).unwrap() - 0.01).abs() < 1e-15);
    assert_eq!(
        bonds(&bond, &coincident, &mut Gradients::total(2))
            .unwrap_err()
            .kind,
        SingularGeometry::CoincidentAtoms
    );
    let nonbonded = Nonbonded {
        sigma: vec![0.3; 2],
        epsilon: vec![0.1; 2],
        charge: vec![0.0; 2],
        exceptions: vec![vec![], vec![]],
    };
    let error = nonbonded
        .evaluate(&coincident, &mut Gradients::None)
        .unwrap_err();
    assert_eq!(
        (error.atoms, error.kind),
        (vec![0, 1], SingularGeometry::CoincidentAtoms)
    );
    // Fully excluded pairs are never evaluated.
    let excluded = Nonbonded {
        exceptions: vec![vec![(1, 0.0, 0.0)], vec![]],
        ..nonbonded
    };
    assert_eq!(
        excluded
            .evaluate(&coincident, &mut Gradients::None)
            .unwrap(),
        (0.0, 0.0)
    );
}

#[test]
fn lennard_jones_and_coulomb_match_closed_forms() {
    let x = [p(0., 0., 0.), p(0.4, 0., 0.)];
    let nonbonded = Nonbonded {
        sigma: vec![0.3, 0.34],
        epsilon: vec![0.5, 0.2],
        charge: vec![0.5, -0.25],
        exceptions: vec![vec![(1, 0.5, 0.25)], vec![]],
    };
    let (vdw, electrostatics) = nonbonded.evaluate(&x, &mut Gradients::None).unwrap();
    let sigma: f64 = 0.32;
    let epsilon = 0.1_f64.sqrt();
    let sr6 = (sigma / 0.4).powi(6);
    assert!((vdw - 0.5 * 4.0 * epsilon * (sr6 * sr6 - sr6)).abs() < 1e-15);
    assert!((electrostatics - 0.25 * COULOMB * -0.125 / 0.4).abs() < 1e-13);
}
