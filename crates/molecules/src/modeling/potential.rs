use std::collections::BTreeSet;
use std::fmt;

use crate::core::{AtomId, BondId};

use super::{ModelDefinitionKey, MolecularModel};

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
                    atom: AtomId::new(index as u32),
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

    pub fn gradient_for(&self, atom: AtomId) -> Option<Vector3> {
        self.gradient.get(atom.index()).copied()
    }
}

/// Energy-and-gradient evaluator for a fixed-topology molecular model.
///
/// Implementations may retain mutable caches between calls. Every returned
/// evaluation must contain one finite gradient vector per model atom. Prepared
/// implementations should bind to [`MolecularModel::definition_key`] and return
/// [`PotentialError::IncompatibleModel`] for a different definition.
pub trait Potential {
    fn evaluate(&mut self, model: &MolecularModel) -> Result<PotentialEvaluation, PotentialError>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Explicit parameter for one harmonic bond term.
pub struct HarmonicBondParameter {
    /// Bond in the model topology.
    pub bond: BondId,
    /// Equilibrium bond length in angstroms.
    pub equilibrium_length: f64,
    /// Harmonic force constant in kJ/mol/angstrom squared.
    pub force_constant: f64,
}

impl HarmonicBondParameter {
    pub const fn new(bond: BondId, equilibrium_length: f64, force_constant: f64) -> Self {
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
    definition: ModelDefinitionKey,
    terms: Vec<HarmonicBondTerm>,
}

#[derive(Debug, Clone, PartialEq)]
struct HarmonicBondTerm {
    a: AtomId,
    b: AtomId,
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
            terms.push(HarmonicBondTerm {
                a,
                b,
                equilibrium_length: parameter.equilibrium_length,
                force_constant: parameter.force_constant,
            });
        }
        Ok(Self {
            definition: model.definition_key().clone(),
            terms,
        })
    }
}

impl Potential for HarmonicBondPotential {
    fn evaluate(&mut self, model: &MolecularModel) -> Result<PotentialEvaluation, PotentialError> {
        if &self.definition != model.definition_key() {
            return Err(PotentialError::IncompatibleModel);
        }
        let mut energy = 0.0;
        let mut gradient = vec![Vector3::zero(); model.atom_count()];
        for term in &self.terms {
            let a = model
                .position(term.a)
                .map_err(|_| PotentialError::IncompatibleModel)?;
            let b = model
                .position(term.b)
                .map_err(|_| PotentialError::IncompatibleModel)?;
            let displacement = Vector3::new(a.x - b.x, a.y - b.y, a.z - b.z);
            let distance = displacement.norm();
            if distance == 0.0 {
                return Err(PotentialError::invalid_geometry(
                    "harmonic bond",
                    [term.a, term.b],
                    PotentialGeometryError::CoincidentAtoms,
                ));
            }
            let extension = distance - term.equilibrium_length;
            energy += 0.5 * term.force_constant * extension * extension;
            let scale = term.force_constant * extension / distance;
            gradient[term.a.index()].add_scaled(displacement, scale);
            gradient[term.b.index()].add_scaled(displacement, -scale);
        }
        PotentialEvaluation::new(model, energy, gradient)
    }
}

/// Coordinate singularity reported by a potential evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PotentialGeometryError {
    CoincidentAtoms,
    DegenerateAngle,
    DegenerateDihedral,
    DegenerateInversion,
}

impl fmt::Display for PotentialGeometryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CoincidentAtoms => f.write_str("coincident atoms"),
            Self::DegenerateAngle => f.write_str("a degenerate angle"),
            Self::DegenerateDihedral => f.write_str("a degenerate dihedral"),
            Self::DegenerateInversion => f.write_str("a degenerate inversion"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum PotentialError {
    InvalidBondId(BondId),
    DuplicateBondParameter(BondId),
    InvalidBondParameter {
        bond: BondId,
        parameter: &'static str,
    },
    IncompatibleModel,
    InvalidGeometry {
        interaction: &'static str,
        atoms: Vec<AtomId>,
        kind: PotentialGeometryError,
    },
    NonFiniteEnergy,
    GradientLengthMismatch {
        expected: usize,
        actual: usize,
    },
    NonFiniteGradient {
        atom: AtomId,
    },
    Backend {
        backend: &'static str,
        message: String,
    },
}

impl PotentialError {
    pub fn invalid_geometry(
        interaction: &'static str,
        atoms: impl IntoIterator<Item = AtomId>,
        kind: PotentialGeometryError,
    ) -> Self {
        Self::InvalidGeometry {
            interaction,
            atoms: atoms.into_iter().collect(),
            kind,
        }
    }

    pub fn backend(backend: &'static str, message: impl Into<String>) -> Self {
        Self::Backend {
            backend,
            message: message.into(),
        }
    }

    /// Returns whether the failure is caused only by the evaluated coordinates.
    pub const fn is_invalid_geometry(&self) -> bool {
        matches!(self, Self::InvalidGeometry { .. })
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
            Self::IncompatibleModel => write!(
                f,
                "model definition differs from the definition bound to the potential"
            ),
            Self::InvalidGeometry {
                interaction,
                atoms,
                kind,
            } => {
                write!(f, "{interaction} has {kind} for atoms [")?;
                for (index, atom) in atoms.iter().enumerate() {
                    if index != 0 {
                        f.write_str(", ")?;
                    }
                    write!(f, "{atom}")?;
                }
                f.write_str("]")
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
            Self::Backend { backend, message } => {
                write!(f, "{backend} potential evaluation failed: {message}")
            }
        }
    }
}

impl std::error::Error for PotentialError {}
