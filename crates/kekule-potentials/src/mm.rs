//! Functional-form kernels shared by classical molecular-mechanics backends.
//!
//! Kernels operate on dense coordinates in nm and parameters already converted
//! to canonical units (kJ/mol, nm, rad, e). They know nothing about force-field
//! rules; a backend lowers its parameters into these terms once at preparation.
//!
//! Every kernel follows the crate-level singular-geometry policy: gradients
//! are exact, uncapped derivatives, and coordinates are rejected only where a
//! requested quantity is undefined.

use kekule::geometry::{Point3, Vector3};

use crate::potential::SingularGeometry;

/// Coulomb conversion factor, kJ mol^-1 nm e^-2, matching OpenMM.
pub(crate) const COULOMB: f64 = 138.935_457_644_381_98;
const MIN_DISTANCE_SQUARED: f64 = 1e-24;

pub(crate) const BONDS: usize = 0;
pub(crate) const ANGLES: usize = 1;
pub(crate) const PROPER_TORSIONS: usize = 2;
pub(crate) const IMPROPER_TORSIONS: usize = 3;
pub(crate) const VAN_DER_WAALS: usize = 4;
pub(crate) const ELECTROSTATICS: usize = 5;
pub(crate) const COMPONENTS: usize = 6;

/// A singular interaction, reported with dense atom indices.
#[derive(Debug)]
pub(crate) struct Singular {
    pub interaction: &'static str,
    pub atoms: Vec<usize>,
    pub kind: SingularGeometry,
}

fn singular(interaction: &'static str, atoms: &[usize], kind: SingularGeometry) -> Singular {
    Singular {
        interaction,
        atoms: atoms.to_vec(),
        kind,
    }
}

/// Requested gradient storage.
pub(crate) enum Gradients {
    None,
    Total(Vec<Vector3>),
    Components(Vec<Vec<Vector3>>),
}

impl Gradients {
    pub(crate) fn total(atoms: usize) -> Self {
        Self::Total(vec![Vector3::zero(); atoms])
    }

    pub(crate) fn components(atoms: usize) -> Self {
        Self::Components(vec![vec![Vector3::zero(); atoms]; COMPONENTS])
    }

    fn enabled(&self) -> bool {
        !matches!(self, Self::None)
    }

    #[inline]
    fn add(&mut self, component: usize, atom: usize, vector: Vector3, scale: f64) {
        let target = match self {
            Self::None => return,
            Self::Total(gradient) => &mut gradient[atom],
            Self::Components(gradients) => &mut gradients[component][atom],
        };
        target.x += scale * vector.x;
        target.y += scale * vector.y;
        target.z += scale * vector.z;
    }
}

fn delta(a: Point3, b: Point3) -> Vector3 {
    Vector3::new(a.x - b.x, a.y - b.y, a.z - b.z)
}

fn scaled(v: Vector3, s: f64) -> Vector3 {
    Vector3::new(v.x * s, v.y * s, v.z * s)
}

/// `k / 2 * (r - length)^2`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct HarmonicBond {
    pub atoms: [usize; 2],
    pub k: f64,
    pub length: f64,
}

/// `k / 2 * (theta - angle)^2`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct HarmonicAngle {
    pub atoms: [usize; 3],
    pub k: f64,
    pub angle: f64,
}

/// `k * (1 + cos(periodicity * phi - phase))`, with any divisor already applied to `k`.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FourierTerm {
    pub k: f64,
    pub periodicity: f64,
    pub phase: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct PeriodicTorsion {
    pub atoms: [usize; 4],
    pub terms: Box<[FourierTerm]>,
}

pub(crate) fn bonds(
    terms: &[HarmonicBond],
    x: &[Point3],
    gradients: &mut Gradients,
) -> Result<f64, Singular> {
    let mut energy = 0.0;
    for term in terms {
        let [a, b] = term.atoms;
        let r = delta(x[a], x[b]);
        let r2 = r.norm_squared();
        let distance = r2.sqrt();
        let stretch = distance - term.length;
        energy += 0.5 * term.k * stretch * stretch;
        if gradients.enabled() {
            if r2 < MIN_DISTANCE_SQUARED {
                return Err(singular(
                    "harmonic bond",
                    &term.atoms,
                    SingularGeometry::CoincidentAtoms,
                ));
            }
            let scale = term.k * stretch / distance;
            gradients.add(BONDS, a, r, scale);
            gradients.add(BONDS, b, r, -scale);
        }
    }
    Ok(energy)
}

pub(crate) fn angles(
    terms: &[HarmonicAngle],
    x: &[Point3],
    gradients: &mut Gradients,
) -> Result<f64, Singular> {
    let mut energy = 0.0;
    for term in terms {
        let [a, b, c] = term.atoms;
        let u = delta(x[a], x[b]);
        let v = delta(x[c], x[b]);
        let r1 = u.norm_squared();
        let r2 = v.norm_squared();
        let normal = u.cross(v);
        let cosine = u.dot(v);
        // A zero arm has no direction. Arms so short that both the sine and
        // cosine terms underflow are indistinguishable from zero.
        if r1 == 0.0 || r2 == 0.0 || (normal == Vector3::zero() && cosine == 0.0) {
            return Err(singular(
                "harmonic angle",
                &term.atoms,
                SingularGeometry::DegenerateAngle,
            ));
        }
        let norm = normal.norm_squared().sqrt();
        // atan2 avoids acos cancellation near linearity. Neither the normal nor
        // the derivative is capped, so the gradient matches the energy exactly.
        let offset = norm.atan2(cosine) - term.angle;
        energy += 0.5 * term.k * offset * offset;
        if !gradients.enabled() {
            continue;
        }
        if norm == 0.0 {
            if offset.abs() > 1e-12 && term.k != 0.0 {
                return Err(singular(
                    "harmonic angle",
                    &term.atoms,
                    SingularGeometry::LinearAngle,
                ));
            }
            continue;
        }
        let normal = scaled(normal, 1.0 / norm);
        let da = scaled(u.cross(normal), term.k * offset / r1);
        let dc = scaled(normal.cross(v), term.k * offset / r2);
        gradients.add(ANGLES, a, da, 1.0);
        gradients.add(ANGLES, c, dc, 1.0);
        gradients.add(ANGLES, b, da, -1.0);
        gradients.add(ANGLES, b, dc, -1.0);
    }
    Ok(energy)
}

pub(crate) fn torsions(
    component: usize,
    interaction: &'static str,
    terms: &[PeriodicTorsion],
    x: &[Point3],
    gradients: &mut Gradients,
) -> Result<f64, Singular> {
    let mut energy = 0.0;
    for term in terms {
        let points = term.atoms.map(|atom| x[atom]);
        let Some(phi) = dihedral(points) else {
            return Err(singular(
                interaction,
                &term.atoms,
                SingularGeometry::DegenerateDihedral,
            ));
        };
        let mut derivative = 0.0;
        for t in term.terms.iter() {
            let argument = t.periodicity * phi - t.phase;
            energy += t.k * (1.0 + argument.cos());
            derivative -= t.k * t.periodicity * argument.sin();
        }
        if gradients.enabled() {
            for (atom, vector) in term.atoms.into_iter().zip(dihedral_gradient(points)) {
                gradients.add(component, atom, vector, derivative);
            }
        }
    }
    Ok(energy)
}

/// Signed IUPAC dihedral angle in `(-pi, pi]`, or `None` when undefined.
fn dihedral([a, b, c, d]: [Point3; 4]) -> Option<f64> {
    let axis = delta(c, b);
    let norm = axis.norm_squared().sqrt();
    if norm < 1e-12 {
        return None;
    }
    let v = axis.cross(delta(a, b));
    let w = axis.cross(delta(d, c));
    let scale = (v.norm_squared() * w.norm_squared()).sqrt();
    if scale < 1e-18 * norm * norm {
        return None;
    }
    let cosine = (v.dot(w) / scale).clamp(-1.0, 1.0);
    let sine = (axis.dot(v.cross(w)) / (norm * scale)).clamp(-1.0, 1.0);
    Some(sine.atan2(cosine))
}

/// Cartesian derivative of [`dihedral`] for each of the four atoms.
fn dihedral_gradient([a, b, c, d]: [Point3; 4]) -> [Vector3; 4] {
    let left = delta(a, b);
    let axis = delta(c, b);
    let right = delta(d, c);
    let axis_sq = axis.norm_squared();
    let axis_norm = axis_sq.sqrt();
    let n0 = left.cross(axis);
    let n3 = axis.cross(right);
    let g0 = scaled(n0, axis_norm / n0.norm_squared());
    let g3 = scaled(n3, axis_norm / n3.norm_squared());
    let alpha = left.dot(axis) / axis_sq;
    let beta = right.dot(axis) / axis_sq;
    let combine = |p: f64, q: f64| {
        Vector3::new(
            p * g0.x + q * g3.x,
            p * g0.y + q * g3.y,
            p * g0.z + q * g3.z,
        )
    };
    [
        g0,
        combine(alpha - 1.0, beta),
        combine(-alpha, -(beta + 1.0)),
        g3,
    ]
}

/// Lennard-Jones 12-6 with Lorentz-Berthelot mixing plus Coulomb, over every
/// atom pair without cutoff. Listed exceptions scale individual pairs.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Nonbonded {
    pub sigma: Vec<f64>,
    pub epsilon: Vec<f64>,
    pub charge: Vec<f64>,
    /// For each atom `i`, exceptions `(j, vdw scale, electrostatics scale)`
    /// with `j > i`, sorted by `j`.
    pub exceptions: Vec<Vec<(usize, f64, f64)>>,
}

impl Nonbonded {
    pub(crate) fn evaluate(
        &self,
        x: &[Point3],
        gradients: &mut Gradients,
    ) -> Result<(f64, f64), Singular> {
        let (mut vdw, mut electrostatics) = (0.0, 0.0);
        let n = x.len();
        for i in 0..n {
            let exceptions = &self.exceptions[i];
            let mut next = 0;
            for j in i + 1..n {
                let (vs, qs) = match exceptions.get(next) {
                    Some(&(atom, vs, qs)) if atom == j => {
                        next += 1;
                        (vs, qs)
                    }
                    _ => (1.0, 1.0),
                };
                let lj = vs != 0.0 && self.epsilon[i] != 0.0 && self.epsilon[j] != 0.0;
                let coulomb = qs != 0.0 && self.charge[i] != 0.0 && self.charge[j] != 0.0;
                if !lj && !coulomb {
                    continue;
                }
                let r = delta(x[i], x[j]);
                let r2 = r.norm_squared();
                if r2 < MIN_DISTANCE_SQUARED {
                    return Err(singular(
                        "nonbonded pair",
                        &[i, j],
                        SingularGeometry::CoincidentAtoms,
                    ));
                }
                if lj {
                    let sigma = (self.sigma[i] + self.sigma[j]) / 2.0;
                    let epsilon = (self.epsilon[i] * self.epsilon[j]).sqrt();
                    let s2 = sigma * sigma / r2;
                    let s6 = s2 * s2 * s2;
                    let s12 = s6 * s6;
                    vdw += vs * 4.0 * epsilon * (s12 - s6);
                    // (dE/dr) / r
                    let scale = vs * 24.0 * epsilon * (s6 - 2.0 * s12) / r2;
                    gradients.add(VAN_DER_WAALS, i, r, scale);
                    gradients.add(VAN_DER_WAALS, j, r, -scale);
                }
                if coulomb {
                    let energy = qs * COULOMB * self.charge[i] * self.charge[j] / r2.sqrt();
                    electrostatics += energy;
                    let scale = -energy / r2;
                    gradients.add(ELECTROSTATICS, i, r, scale);
                    gradients.add(ELECTROSTATICS, j, r, -scale);
                }
            }
        }
        Ok((vdw, electrostatics))
    }
}

#[cfg(test)]
mod tests {
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
        let energy = torsions(0, "t", &[torsion], &points, &mut Gradients::None).unwrap();
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
            Box::new(|x, g| torsions(0, "t", &torsion, x, g).unwrap()),
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
}
