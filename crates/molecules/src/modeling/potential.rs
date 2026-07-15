use std::collections::BTreeSet;
use std::fmt;

use super::{InstanceAtomId, InstanceBondId, Model, ModelDefinitionKey};
use crate::units::{
    Quantity, ScaleValue, UnitError, MODEL_ENERGY_UNIT, MODEL_FORCE_CONSTANT_UNIT,
    MODEL_GRADIENT_UNIT, MODEL_LENGTH_UNIT,
};

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

impl ScaleValue for Vector3 {
    fn scaled(self, factor: f64) -> Self {
        Self::new(self.x * factor, self.y * factor, self.z * factor)
    }
}

#[derive(Debug, Clone, PartialEq)]
/// Validated energy and Cartesian gradient from a [`Potential`].
///
/// Values are converted once to the modelling kernel's explicit canonical
/// energy and gradient units.
pub struct PotentialEvaluation {
    energy: Quantity<f64>,
    gradient: Quantity<Vec<Vector3>>,
}

impl PotentialEvaluation {
    pub fn new(
        model: &Model,
        energy: Quantity<f64>,
        gradient: Quantity<Vec<Vector3>>,
    ) -> Result<Self, PotentialError> {
        let energy = energy.into_unit(MODEL_ENERGY_UNIT)?;
        let gradient = gradient.into_unit(MODEL_GRADIENT_UNIT)?;
        if !energy.value().is_finite() {
            return Err(PotentialError::NonFiniteEnergy);
        }
        if gradient.value().len() != model.atom_count() {
            return Err(PotentialError::GradientLengthMismatch {
                expected: model.atom_count(),
                actual: gradient.value().len(),
            });
        }
        for (index, vector) in gradient.value().iter().copied().enumerate() {
            if !vector.is_finite() {
                return Err(PotentialError::NonFiniteGradient {
                    atom: model.topology().atom_ids()[index],
                });
            }
        }
        Ok(Self { energy, gradient })
    }

    pub fn energy(&self) -> Quantity<f64> {
        self.energy
    }

    pub fn gradient(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.gradient.value().as_slice(), self.gradient.unit())
    }

    pub fn gradient_for(&self, model: &Model, atom: InstanceAtomId) -> Option<Quantity<Vector3>> {
        let index = model.topology().atom_index(atom)?;
        self.gradient
            .value()
            .get(index.index())
            .copied()
            .map(|vector| Quantity::new(vector, self.gradient.unit()))
    }
}

/// Energy-and-gradient evaluator for a fixed-topology molecular model.
///
/// Implementations may retain mutable caches between calls. Every returned
/// evaluation must contain one finite gradient vector per model atom. Prepared
/// implementations should bind to [`Model::definition_key`] and return
/// [`PotentialError::IncompatibleModel`] for a different definition.
pub trait Potential {
    fn evaluate(&mut self, model: &Model) -> Result<PotentialEvaluation, PotentialError>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Explicit parameter for one harmonic bond term.
pub struct HarmonicBondParameter {
    /// Bond in the model topology.
    pub bond: InstanceBondId,
    /// Equilibrium bond length.
    pub equilibrium_length: Quantity<f64>,
    /// Harmonic force constant (energy per squared length).
    pub force_constant: Quantity<f64>,
}

impl HarmonicBondParameter {
    pub const fn new(
        bond: InstanceBondId,
        equilibrium_length: Quantity<f64>,
        force_constant: Quantity<f64>,
    ) -> Self {
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
    a: InstanceAtomId,
    b: InstanceAtomId,
    equilibrium_length: f64,
    force_constant: f64,
}

impl HarmonicBondPotential {
    pub fn new(
        model: &Model,
        parameters: impl IntoIterator<Item = HarmonicBondParameter>,
    ) -> Result<Self, PotentialError> {
        let mut seen = BTreeSet::new();
        let mut terms = Vec::new();
        for parameter in parameters {
            if !seen.insert(parameter.bond) {
                return Err(PotentialError::DuplicateBondParameter(parameter.bond));
            }
            let equilibrium_length = parameter
                .equilibrium_length
                .into_unit(MODEL_LENGTH_UNIT)?
                .into_value();
            let force_constant = parameter
                .force_constant
                .into_unit(MODEL_FORCE_CONSTANT_UNIT)?
                .into_value();
            if !equilibrium_length.is_finite() || equilibrium_length <= 0.0 {
                return Err(PotentialError::InvalidBondParameter {
                    bond: parameter.bond,
                    parameter: "equilibrium length must be finite and positive",
                });
            }
            if !force_constant.is_finite() || force_constant <= 0.0 {
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
                a,
                b,
                equilibrium_length,
                force_constant,
            });
        }
        Ok(Self {
            definition: model.definition_key().clone(),
            terms,
        })
    }
}

impl Potential for HarmonicBondPotential {
    fn evaluate(&mut self, model: &Model) -> Result<PotentialEvaluation, PotentialError> {
        if &self.definition != model.definition_key() {
            return Err(PotentialError::IncompatibleModel);
        }
        let mut energy = 0.0;
        let mut gradient = vec![Vector3::zero(); model.atom_count()];
        for term in &self.terms {
            let a = model
                .position(term.a)
                .map_err(|_| PotentialError::IncompatibleModel)?
                .into_value();
            let b = model
                .position(term.b)
                .map_err(|_| PotentialError::IncompatibleModel)?
                .into_value();
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
        PotentialEvaluation::new(
            model,
            Quantity::new(energy, MODEL_ENERGY_UNIT),
            Quantity::new(gradient, MODEL_GRADIENT_UNIT),
        )
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

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PotentialError {
    InvalidBondId(InstanceBondId),
    DuplicateBondParameter(InstanceBondId),
    InvalidBondParameter {
        bond: InstanceBondId,
        parameter: &'static str,
    },
    IncompatibleModel,
    InvalidGeometry {
        interaction: &'static str,
        atoms: Vec<InstanceAtomId>,
        kind: PotentialGeometryError,
    },
    NonFiniteEnergy,
    GradientLengthMismatch {
        expected: usize,
        actual: usize,
    },
    NonFiniteGradient {
        atom: InstanceAtomId,
    },
    Unit(UnitError),
    Backend {
        backend: &'static str,
        message: String,
    },
}

impl PotentialError {
    pub fn invalid_geometry(
        interaction: &'static str,
        atoms: impl IntoIterator<Item = InstanceAtomId>,
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
            Self::Unit(error) => write!(f, "invalid potential quantity unit: {error}"),
            Self::Backend { backend, message } => {
                write!(f, "{backend} potential evaluation failed: {message}")
            }
        }
    }
}

impl std::error::Error for PotentialError {}

impl From<UnitError> for PotentialError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
