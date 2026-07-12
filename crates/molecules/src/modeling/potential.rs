use std::collections::BTreeSet;
use std::fmt;

use crate::core::BondOrder;

use super::{InstanceAtomId, InstanceBondId, MolecularModel};

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// Three-dimensional Cartesian vector used for energy gradients.
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub const fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    pub(crate) fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    pub(crate) fn add_scaled(&mut self, other: Self, scale: f64) {
        self.x += other.x * scale;
        self.y += other.y * scale;
        self.z += other.z * scale;
    }
}

#[derive(Debug, Clone, PartialEq)]
/// Validated energy and Cartesian gradient from a [`Potential`].
///
/// Energy is expressed in kJ/mol and gradient components in kJ/mol/angstrom.
pub struct PotentialEvaluation {
    energy: f64,
    gradient: Vec<Vector3>,
}

impl PotentialEvaluation {
    pub fn new(
        model: &MolecularModel,
        energy: f64,
        gradient: Vec<Vector3>,
    ) -> Result<Self, PotentialError> {
        if !energy.is_finite() {
            return Err(PotentialError::NonFiniteEnergy);
        }
        if gradient.len() != model.atom_count() {
            return Err(PotentialError::GradientLengthMismatch {
                expected: model.atom_count(),
                actual: gradient.len(),
            });
        }
        for (index, vector) in gradient.iter().copied().enumerate() {
            if !vector.is_finite() {
                return Err(PotentialError::NonFiniteGradient {
                    atom: model.topology().atom_ids()[index],
                });
            }
        }
        Ok(Self { energy, gradient })
    }

    pub fn energy(&self) -> f64 {
        self.energy
    }

    pub fn gradient(&self) -> &[Vector3] {
        &self.gradient
    }

    pub fn gradient_for(&self, model: &MolecularModel, atom: InstanceAtomId) -> Option<Vector3> {
        let index = model.topology().atom_index(atom)?;
        self.gradient.get(index.index()).copied()
    }
}

/// Energy-and-gradient evaluator for a fixed-topology molecular model.
///
/// Implementations may retain mutable caches between calls. Every returned
/// evaluation must contain one finite gradient vector per model atom.
pub trait Potential {
    fn evaluate(&mut self, model: &MolecularModel) -> Result<PotentialEvaluation, PotentialError>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Explicit parameter for one harmonic bond term.
pub struct HarmonicBondParameter {
    /// Bond in the model topology.
    pub bond: InstanceBondId,
    /// Equilibrium bond length in angstroms.
    pub equilibrium_length: f64,
    /// Harmonic force constant in kJ/mol/angstrom squared.
    pub force_constant: f64,
}

impl HarmonicBondParameter {
    pub const fn new(bond: InstanceBondId, equilibrium_length: f64, force_constant: f64) -> Self {
        Self {
            bond,
            equilibrium_length,
            force_constant,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
/// Caller-parameterized harmonic bond potential.
///
/// Each term contributes `0.5 * k * (r - r0)^2`. No parameters are inferred,
/// and angle, torsion, and nonbonded interactions are intentionally absent.
pub struct HarmonicBondPotential {
    terms: Vec<HarmonicBondTerm>,
}

#[derive(Debug, Clone, PartialEq)]
struct HarmonicBondTerm {
    bond: InstanceBondId,
    a: InstanceAtomId,
    b: InstanceAtomId,
    order: BondOrder,
    equilibrium_length: f64,
    force_constant: f64,
}

impl HarmonicBondPotential {
    pub fn new(
        model: &MolecularModel,
        parameters: impl IntoIterator<Item = HarmonicBondParameter>,
    ) -> Result<Self, PotentialError> {
        let mut seen = BTreeSet::new();
        let mut terms = Vec::new();
        for parameter in parameters {
            if !seen.insert(parameter.bond) {
                return Err(PotentialError::DuplicateBondParameter(parameter.bond));
            }
            if !parameter.equilibrium_length.is_finite() || parameter.equilibrium_length <= 0.0 {
                return Err(PotentialError::InvalidBondParameter {
                    bond: parameter.bond,
                    parameter: "equilibrium length must be finite and positive",
                });
            }
            if !parameter.force_constant.is_finite() || parameter.force_constant <= 0.0 {
                return Err(PotentialError::InvalidBondParameter {
                    bond: parameter.bond,
                    parameter: "force constant must be finite and positive",
                });
            }
            let bond = model
                .topology()
                .bond(parameter.bond)
                .map_err(|_| PotentialError::InvalidBondId(parameter.bond))?;
            let (a, b) = bond.endpoints();
            let a = InstanceAtomId::new(parameter.bond.molecule(), a);
            let b = InstanceAtomId::new(parameter.bond.molecule(), b);
            terms.push(HarmonicBondTerm {
                bond: parameter.bond,
                a,
                b,
                order: bond.order,
                equilibrium_length: parameter.equilibrium_length,
                force_constant: parameter.force_constant,
            });
        }
        Ok(Self { terms })
    }
}

impl Potential for HarmonicBondPotential {
    fn evaluate(&mut self, model: &MolecularModel) -> Result<PotentialEvaluation, PotentialError> {
        let mut energy = 0.0;
        let mut gradient = vec![Vector3::zero(); model.atom_count()];
        for term in &self.terms {
            let bond = model
                .topology()
                .bond(term.bond)
                .map_err(|_| PotentialError::ModelTopologyMismatch(term.bond))?;
            if bond.endpoints() != (term.a.atom(), term.b.atom()) || bond.order != term.order {
                return Err(PotentialError::ModelTopologyMismatch(term.bond));
            }
            let a = model
                .position(term.a)
                .map_err(|_| PotentialError::ModelTopologyMismatch(term.bond))?;
            let b = model
                .position(term.b)
                .map_err(|_| PotentialError::ModelTopologyMismatch(term.bond))?;
            let displacement = Vector3::new(a.x - b.x, a.y - b.y, a.z - b.z);
            let distance = displacement.norm();
            if distance == 0.0 {
                return Err(PotentialError::CoincidentBondAtoms(term.bond));
            }
            let extension = distance - term.equilibrium_length;
            energy += 0.5 * term.force_constant * extension * extension;
            let scale = term.force_constant * extension / distance;
            let a_index = model
                .topology()
                .atom_index(term.a)
                .expect("validated harmonic atom");
            let b_index = model
                .topology()
                .atom_index(term.b)
                .expect("validated harmonic atom");
            gradient[a_index.index()].add_scaled(displacement, scale);
            gradient[b_index.index()].add_scaled(displacement, -scale);
        }
        PotentialEvaluation::new(model, energy, gradient)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PotentialError {
    InvalidBondId(InstanceBondId),
    DuplicateBondParameter(InstanceBondId),
    InvalidBondParameter {
        bond: InstanceBondId,
        parameter: &'static str,
    },
    ModelTopologyMismatch(InstanceBondId),
    CoincidentBondAtoms(InstanceBondId),
    NonFiniteEnergy,
    GradientLengthMismatch {
        expected: usize,
        actual: usize,
    },
    NonFiniteGradient {
        atom: InstanceAtomId,
    },
    Custom(String),
}

impl PotentialError {
    pub fn custom(message: impl Into<String>) -> Self {
        Self::Custom(message.into())
    }
}

impl fmt::Display for PotentialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBondId(bond) => write!(f, "invalid harmonic bond id: {bond}"),
            Self::DuplicateBondParameter(bond) => {
                write!(f, "duplicate harmonic parameter for bond {bond}")
            }
            Self::InvalidBondParameter { bond, parameter } => {
                write!(f, "invalid harmonic parameter for bond {bond}: {parameter}")
            }
            Self::ModelTopologyMismatch(bond) => {
                write!(f, "model topology does not match harmonic bond {bond}")
            }
            Self::CoincidentBondAtoms(bond) => {
                write!(
                    f,
                    "bond {bond} has coincident atoms and an undefined gradient"
                )
            }
            Self::NonFiniteEnergy => write!(f, "potential returned a non-finite energy"),
            Self::GradientLengthMismatch { expected, actual } => write!(
                f,
                "potential returned {actual} gradients for a model with {expected} atoms"
            ),
            Self::NonFiniteGradient { atom } => {
                write!(
                    f,
                    "potential returned a non-finite gradient for atom {atom}"
                )
            }
            Self::Custom(message) => write!(f, "potential evaluation failed: {message}"),
        }
    }
}

impl std::error::Error for PotentialError {}
