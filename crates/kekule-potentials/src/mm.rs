//! Functional-form kernels shared by classical molecular-mechanics backends.
//!
//! Kernels operate on dense coordinates in nm and parameters already converted
//! to canonical units (kJ/mol, nm, rad, e). They know nothing about force-field
//! rules; a backend lowers its parameters into these terms once at preparation.
//!
//! Every kernel follows the crate-level singular-geometry policy: gradients
//! are exact, uncapped derivatives, and coordinates are rejected only where a
//! requested quantity is undefined.

use kekule::math::{Point3, Vector3};

use crate::potential::{SingularGeometry, TermKind};

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
    pub interaction: TermKind,
    pub atoms: Vec<usize>,
    pub kind: SingularGeometry,
}

fn singular(interaction: TermKind, atoms: &[usize], kind: SingularGeometry) -> Singular {
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
        let r = x[a] - x[b];
        let r2 = r.norm_squared();
        let distance = r2.sqrt();
        let stretch = distance - term.length;
        energy += 0.5 * term.k * stretch * stretch;
        if gradients.enabled() {
            if r2 < MIN_DISTANCE_SQUARED {
                return Err(singular(
                    TermKind::Bond,
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
        let u = x[a] - x[b];
        let v = x[c] - x[b];
        let r1 = u.norm_squared();
        let r2 = v.norm_squared();
        let normal = u.cross(v);
        let cosine = u.dot(v);
        // A zero arm has no direction. Arms so short that both the sine and
        // cosine terms underflow are indistinguishable from zero.
        if r1 == 0.0 || r2 == 0.0 || (normal == Vector3::zero() && cosine == 0.0) {
            return Err(singular(
                TermKind::Angle,
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
                    TermKind::Angle,
                    &term.atoms,
                    SingularGeometry::LinearAngle,
                ));
            }
            continue;
        }
        let normal = normal * (1.0 / norm);
        let da = u.cross(normal) * (term.k * offset / r1);
        let dc = normal.cross(v) * (term.k * offset / r2);
        gradients.add(ANGLES, a, da, 1.0);
        gradients.add(ANGLES, c, dc, 1.0);
        gradients.add(ANGLES, b, da, -1.0);
        gradients.add(ANGLES, b, dc, -1.0);
    }
    Ok(energy)
}

pub(crate) fn torsions(
    component: usize,
    interaction: TermKind,
    terms: &[PeriodicTorsion],
    x: &[Point3],
    gradients: &mut Gradients,
) -> Result<f64, Singular> {
    let mut energy = 0.0;
    for term in terms {
        // Without amplitude the energy and its gradient vanish identically, so
        // the dihedral need not be defined; linear groups rely on this.
        if term.terms.iter().all(|t| t.k == 0.0) {
            continue;
        }
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
///
/// Not `kekule::structure::measure::dihedral`: that measurement normalizes
/// each bond first and fails on other quartets, so it can differ in the last
/// bits and in which geometries are singular; energies keep this kernel.
fn dihedral([a, b, c, d]: [Point3; 4]) -> Option<f64> {
    let axis = c - b;
    let norm = axis.norm_squared().sqrt();
    if norm < 1e-12 {
        return None;
    }
    let v = axis.cross(a - b);
    let w = axis.cross(d - c);
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
    let left = a - b;
    let axis = c - b;
    let right = d - c;
    let axis_sq = axis.norm_squared();
    let axis_norm = axis_sq.sqrt();
    let n0 = left.cross(axis);
    let n3 = axis.cross(right);
    let g0 = n0 * (axis_norm / n0.norm_squared());
    let g3 = n3 * (axis_norm / n3.norm_squared());
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
                let r = x[i] - x[j];
                let r2 = r.norm_squared();
                if r2 < MIN_DISTANCE_SQUARED {
                    return Err(singular(
                        TermKind::NonbondedPair,
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
mod tests;
