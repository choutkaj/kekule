use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::geometry::Vector3;
use crate::structure::ModelView;
use crate::topology::{InstanceAtomId, InstanceBondId, Topology};
use crate::units::{
    Quantity, UnitError, MODEL_ENERGY_UNIT, MODEL_FORCE_CONSTANT_UNIT, MODEL_GRADIENT_UNIT,
    MODEL_LENGTH_UNIT,
};

#[derive(Debug, Clone)]
/// Validated energy and Cartesian gradient from a [`Potential`].
///
/// Values are converted once to the modelling kernel's explicit canonical
/// energy and gradient units.
pub struct PotentialEvaluation {
    topology: Arc<Topology>,
    energy: Quantity<f64>,
    gradient: Quantity<Vec<Vector3>>,
}

impl PartialEq for PotentialEvaluation {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.energy == other.energy
            && self.gradient == other.gradient
    }
}

impl PotentialEvaluation {
    pub fn new(
        model: ModelView<'_>,
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
        Ok(Self {
            topology: model.shared_topology(),
            energy,
            gradient,
        })
    }

    pub fn energy(&self) -> Quantity<f64> {
        self.energy
    }

    pub fn gradient(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.gradient.value().as_slice(), self.gradient.unit())
    }

    pub fn gradient_for(
        &self,
        model: ModelView<'_>,
        atom: InstanceAtomId,
    ) -> Option<Quantity<Vector3>> {
        if !Arc::ptr_eq(&self.topology, model.topology_arc()) {
            return None;
        }
        let index = model.topology().atom_index(atom)?;
        self.gradient
            .value()
            .get(index.index())
            .copied()
            .map(|vector| Quantity::new(vector, self.gradient.unit()))
    }
}

/// Energy-and-gradient evaluator for a borrowed structural view.
///
/// Implementations may retain mutable caches between calls. Every returned
/// evaluation must contain one finite gradient vector per topology atom.
/// Prepared implementations bind to one shared `Arc<Topology>` allocation. Accepting a
/// [`ModelView`] does not imply support for every configuration field;
/// implementations must document capabilities such as periodic-cell support
/// and return a structured error for unsupported state.
pub trait Potential {
    fn evaluate(&mut self, model: ModelView<'_>) -> Result<PotentialEvaluation, PotentialError>;
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Explicit parameter for one harmonic bond term.
pub struct HarmonicBondParameter {
    /// Bond in the exact topology used to prepare the potential.
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

#[derive(Debug, Clone)]
/// Caller-parameterized harmonic bond potential.
///
/// Each term contributes `0.5 * k * (r - r0)^2`. No parameters are inferred,
/// and angle, torsion, and nonbonded interactions are intentionally absent.
/// This potential is nonperiodic and rejects any evaluated configuration with
/// a periodic cell.
pub struct HarmonicBondPotential {
    topology: Arc<Topology>,
    terms: Vec<HarmonicBondTerm>,
}

impl PartialEq for HarmonicBondPotential {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.terms == other.terms
    }
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
        topology: &Arc<Topology>,
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
            let bond = topology
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
            topology: Arc::clone(topology),
            terms,
        })
    }
}

impl Potential for HarmonicBondPotential {
    fn evaluate(&mut self, model: ModelView<'_>) -> Result<PotentialEvaluation, PotentialError> {
        if !Arc::ptr_eq(&self.topology, model.topology_arc()) {
            return Err(PotentialError::IncompatibleTopology);
        }
        if model.cell().is_some() {
            return Err(PotentialError::UnsupportedPeriodicCell);
        }
        let mut energy = 0.0;
        let mut gradient = vec![Vector3::zero(); model.atom_count()];
        for term in &self.terms {
            let a = model
                .position(term.a)
                .map_err(|_| PotentialError::IncompatibleTopology)?
                .into_value();
            let b = model
                .position(term.b)
                .map_err(|_| PotentialError::IncompatibleTopology)?
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
    /// A supplied bond does not exist in the topology being parameterized.
    InvalidBondId(InstanceBondId),
    /// More than one parameter was supplied for the same topology bond.
    DuplicateBondParameter(InstanceBondId),
    /// A harmonic bond parameter is numerically or dimensionally invalid.
    InvalidBondParameter {
        bond: InstanceBondId,
        parameter: &'static str,
    },
    /// The evaluated view does not share the potential's exact topology.
    IncompatibleTopology,
    /// The potential has no complete periodic-coordinate evaluation policy.
    UnsupportedPeriodicCell,
    /// The evaluated coordinates are singular for a required interaction.
    InvalidGeometry {
        interaction: &'static str,
        atoms: Vec<InstanceAtomId>,
        kind: PotentialGeometryError,
    },
    /// The implementation returned a non-finite energy.
    NonFiniteEnergy,
    /// The implementation returned the wrong number of gradient vectors.
    GradientLengthMismatch { expected: usize, actual: usize },
    /// The implementation returned a non-finite gradient vector.
    NonFiniteGradient { atom: InstanceAtomId },
    /// A public potential quantity used incompatible units.
    Unit(UnitError),
    /// A potential backend reported a non-geometric evaluation failure.
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
            Self::IncompatibleTopology => write!(
                f,
                "model view belongs to a different exact topology than the potential"
            ),
            Self::UnsupportedPeriodicCell => {
                f.write_str("potential does not support periodic-cell configurations")
            }
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
