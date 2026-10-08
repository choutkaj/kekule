//! The evaluation contract shared by every potential and optimizer.
//!
//! [`Potential`] makes no assumption about how an energy is built, so pairwise
//! force fields, machine-learned models, and quantum-chemical backends can all
//! implement it. A prepared potential binds one exact shared `Arc<Topology>` and
//! evaluates any [`ModelView`] of that snapshot: a [`kekule::structure::Model`],
//! an ensemble member, or a trajectory frame. Results are validated at this
//! boundary and stored in Kekule's canonical units (kJ/mol and kJ/mol/nm).

use std::fmt;
use std::sync::Arc;

use kekule::geometry::Vector3;
use kekule::structure::ModelView;
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::{Quantity, UnitError, CANONICAL_ENERGY_UNIT, CANONICAL_GRADIENT_UNIT};

/// Energy-and-gradient evaluator bound to one exact topology snapshot.
///
/// A prepared potential is immutable: evaluation never changes it, and it can
/// be shared between threads, for example to evaluate trajectory frames in
/// parallel. Coordinate-dependent caches belong to the caller driving the
/// evaluations, not to the potential.
///
/// Accepting a [`ModelView`] does not imply support for every realization
/// field. Implementations document capabilities such as periodic-cell support
/// and return a structured [`EvaluationError`] for unsupported state. Every
/// evaluation of a view from a different topology snapshot, including an
/// independently equal one, returns [`EvaluationError::IncompatibleTopology`].
pub trait Potential: Send + Sync {
    /// The exact topology snapshot this potential was prepared for.
    fn topology(&self) -> &Arc<Topology>;

    /// Evaluates the energy and its Cartesian gradient.
    fn evaluate(&self, model: ModelView<'_>) -> Result<Evaluation, EvaluationError>;

    /// Evaluates only the energy.
    ///
    /// Implementations should override this when an energy-only path is
    /// cheaper. The default evaluates and discards the gradient. Some
    /// coordinates, such as an exactly linear strained angle, have a defined
    /// energy but no defined gradient; an energy-only override may accept them.
    fn energy(&self, model: ModelView<'_>) -> Result<Energy, EvaluationError> {
        self.evaluate(model).map(Evaluation::into_energy)
    }
}

/// Named contribution to a decomposable energy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum ComponentKind {
    Bonds,
    Angles,
    ProperTorsions,
    ImproperTorsions,
    VanDerWaals,
    Electrostatics,
}

impl ComponentKind {
    pub const fn name(self) -> &'static str {
        match self {
            Self::Bonds => "bonds",
            Self::Angles => "angles",
            Self::ProperTorsions => "proper torsions",
            Self::ImproperTorsions => "improper torsions",
            Self::VanDerWaals => "van der Waals",
            Self::Electrostatics => "electrostatics",
        }
    }
}

impl fmt::Display for ComponentKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.pad(self.name())
    }
}

/// One component of an [`Energy`], in the canonical energy unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnergyComponent {
    pub kind: ComponentKind,
    pub energy: Quantity<f64>,
}

/// A validated finite energy, optionally decomposed into components.
///
/// When components are present, the total is their sum in the order given.
/// Potentials without a meaningful decomposition report no components.
#[derive(Debug, Clone, PartialEq)]
pub struct Energy {
    total: Quantity<f64>,
    components: Vec<EnergyComponent>,
}

impl Energy {
    /// An energy without a decomposition.
    pub fn new(total: Quantity<f64>) -> Result<Self, EvaluationError> {
        let total = total.into_unit(CANONICAL_ENERGY_UNIT)?;
        if !total.value().is_finite() {
            return Err(EvaluationError::NonFiniteEnergy);
        }
        Ok(Self {
            total,
            components: Vec::new(),
        })
    }

    /// An energy whose total is the ordered sum of distinct components.
    pub fn from_components(
        components: impl IntoIterator<Item = EnergyComponent>,
    ) -> Result<Self, EvaluationError> {
        let mut total = 0.0;
        let mut checked: Vec<EnergyComponent> = Vec::new();
        for component in components {
            let energy = component.energy.into_unit(CANONICAL_ENERGY_UNIT)?;
            if !energy.value().is_finite() {
                return Err(EvaluationError::NonFiniteEnergy);
            }
            if checked.iter().any(|c| c.kind == component.kind) {
                return Err(EvaluationError::DuplicateComponent(component.kind));
            }
            total += energy.value();
            checked.push(EnergyComponent {
                kind: component.kind,
                energy,
            });
        }
        if !total.is_finite() {
            return Err(EvaluationError::NonFiniteEnergy);
        }
        Ok(Self {
            total: Quantity::new(total, CANONICAL_ENERGY_UNIT),
            components: checked,
        })
    }

    pub fn total(&self) -> Quantity<f64> {
        self.total
    }

    pub fn components(&self) -> &[EnergyComponent] {
        &self.components
    }

    pub fn component(&self, kind: ComponentKind) -> Option<Quantity<f64>> {
        self.components
            .iter()
            .find(|c| c.kind == kind)
            .map(|c| c.energy)
    }
}

/// A validated energy and Cartesian gradient `dE/dx` for one realization.
///
/// The gradient has one finite vector per topology atom in dense order. It is
/// the derivative of the potential energy; forces are its negation.
#[derive(Debug, Clone)]
pub struct Evaluation {
    topology: Arc<Topology>,
    energy: Energy,
    gradient: Vec<Vector3>,
}

impl PartialEq for Evaluation {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.energy == other.energy
            && self.gradient == other.gradient
    }
}

impl Evaluation {
    pub fn new(
        model: ModelView<'_>,
        energy: Energy,
        gradient: Quantity<Vec<Vector3>>,
    ) -> Result<Self, EvaluationError> {
        let gradient = gradient.into_unit(CANONICAL_GRADIENT_UNIT)?.into_value();
        if gradient.len() != model.atom_count() {
            return Err(EvaluationError::GradientLengthMismatch {
                expected: model.atom_count(),
                actual: gradient.len(),
            });
        }
        if let Some(index) = gradient.iter().position(|g| !g.is_finite()) {
            return Err(EvaluationError::NonFiniteGradient {
                atom: model.topology().atom_ids()[index],
            });
        }
        Ok(Self {
            topology: model.shared_topology(),
            energy,
            gradient,
        })
    }

    /// The exact topology snapshot of the evaluated realization.
    pub fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    pub fn energy(&self) -> &Energy {
        &self.energy
    }

    pub fn into_energy(self) -> Energy {
        self.energy
    }

    pub fn gradient(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.gradient.as_slice(), CANONICAL_GRADIENT_UNIT)
    }

    /// The gradient of one atom, if `model` shares this evaluation's snapshot.
    pub fn gradient_for(
        &self,
        model: ModelView<'_>,
        atom: InstanceAtomId,
    ) -> Option<Quantity<Vector3>> {
        if !std::ptr::eq(self.topology.as_ref(), model.topology()) {
            return None;
        }
        let index = model.topology().atom_index(atom)?;
        self.gradient
            .get(index.index())
            .map(|g| Quantity::new(*g, CANONICAL_GRADIENT_UNIT))
    }

    /// The largest per-atom gradient norm, the usual convergence measure.
    pub fn max_gradient_norm(&self) -> Quantity<f64> {
        Quantity::new(max_norm(&self.gradient), CANONICAL_GRADIENT_UNIT)
    }
}

pub(crate) fn max_norm(vectors: &[Vector3]) -> f64 {
    vectors.iter().map(|v| v.norm()).fold(0.0, f64::max)
}

/// Coordinates at which a required interaction is mathematically undefined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum SingularGeometry {
    /// Two interacting atoms share a position.
    CoincidentAtoms,
    /// An angle arm has zero (or underflowing) length.
    DegenerateAngle,
    /// A strained angle is exactly linear, so its gradient direction is undefined.
    LinearAngle,
    /// A torsion axis has zero length or an outer atom lies on the axis.
    DegenerateDihedral,
}

impl fmt::Display for SingularGeometry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::CoincidentAtoms => "coincident atoms",
            Self::DegenerateAngle => "a degenerate angle",
            Self::LinearAngle => "a strained linear angle",
            Self::DegenerateDihedral => "a degenerate dihedral",
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum EvaluationError {
    /// The evaluated view does not share the potential's exact topology.
    IncompatibleTopology,
    /// The potential has no evaluation policy for periodic cells.
    UnsupportedPeriodicCell,
    /// The coordinates are singular for a required interaction.
    InvalidGeometry {
        interaction: &'static str,
        atoms: Vec<InstanceAtomId>,
        kind: SingularGeometry,
    },
    NonFiniteEnergy,
    /// An energy decomposition named one component twice.
    DuplicateComponent(ComponentKind),
    GradientLengthMismatch {
        expected: usize,
        actual: usize,
    },
    NonFiniteGradient {
        atom: InstanceAtomId,
    },
    Unit(UnitError),
    /// A backend reported a non-geometric failure.
    Backend {
        backend: &'static str,
        message: String,
    },
}

impl EvaluationError {
    /// Whether the failure is caused only by the evaluated coordinates.
    pub const fn is_invalid_geometry(&self) -> bool {
        matches!(self, Self::InvalidGeometry { .. })
    }
}

impl fmt::Display for EvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IncompatibleTopology => {
                f.write_str("model view belongs to a different exact topology than the potential")
            }
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
            Self::NonFiniteEnergy => f.write_str("potential produced a non-finite energy"),
            Self::DuplicateComponent(kind) => {
                write!(f, "energy decomposition repeats the {kind} component")
            }
            Self::GradientLengthMismatch { expected, actual } => write!(
                f,
                "potential produced {actual} gradients for a model with {expected} atoms"
            ),
            Self::NonFiniteGradient { atom } => {
                write!(
                    f,
                    "potential produced a non-finite gradient for atom {atom}"
                )
            }
            Self::Unit(error) => write!(f, "invalid potential quantity unit: {error}"),
            Self::Backend { backend, message } => {
                write!(f, "{backend} potential evaluation failed: {message}")
            }
        }
    }
}

impl std::error::Error for EvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for EvaluationError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
