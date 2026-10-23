//! The evaluation contract shared by every potential and optimizer.
//!
//! [`Potential`] makes no assumption about how an energy is built, so pairwise
//! force fields, machine-learned models, and quantum-chemical backends can all
//! implement it. A prepared potential binds one topology layout and evaluates
//! any [`ModelView`] sharing it, including perceived snapshots: a
//! [`kekule::structure::Model`], an ensemble sample, or a trajectory frame.
//! Results are validated at this boundary and stored in Kekule's canonical
//! units (kJ/mol and kJ/mol/nm).

use std::fmt;
use std::sync::Arc;

use kekule::geometry::Vector3;
use kekule::structure::ModelView;
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::{Quantity, UnitError, CANONICAL_ENERGY_UNIT, CANONICAL_GRADIENT_UNIT};

/// Energy-and-gradient evaluator bound to one topology layout.
///
/// A prepared potential is immutable: evaluation never changes it, and it can
/// be shared between threads, for example to evaluate trajectory frames in
/// parallel. Coordinate-dependent caches belong to the caller driving the
/// evaluations, not to the potential.
///
/// Accepting a [`ModelView`] does not imply support for every structure
/// field. Implementations document capabilities such as periodic-cell support
/// and return a structured [`EvaluationError`] for unsupported state.
///
/// # Implementing this trait
///
/// This trait is the extension point for energy backends. An implementation
/// must keep this contract:
///
/// - [`Self::topology`] returns the snapshot the potential was prepared for,
///   unchanged for the potential's whole life.
/// - [`Self::evaluate`] and [`Self::energy`] accept every view whose topology
///   shares that layout, including perceived snapshots, and return
///   [`EvaluationError::IncompatibleTopology`] for any other view, including
///   one of an independently equal topology.
/// - State the potential cannot evaluate, such as an unsupported periodic
///   cell, is rejected with an error rather than ignored.
/// - Results are built with [`Evaluation::new`] and [`Energy`], so they are
///   finite and in canonical units; gradients are exact derivatives.
/// - Evaluation does not mutate the potential, which is `Send + Sync`.
pub trait Potential: Send + Sync {
    /// The topology snapshot this potential was prepared for. Evaluation
    /// accepts any snapshot sharing its layout.
    fn topology(&self) -> &Topology;

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
    /// Bond-stretch energy.
    Bonds,
    /// Angle-bend energy.
    Angles,
    /// Proper torsion energy.
    ProperTorsions,
    /// Improper torsion energy.
    ImproperTorsions,
    /// Van der Waals (dispersion and repulsion) energy.
    VanDerWaals,
    /// Electrostatic energy.
    Electrostatics,
}

impl ComponentKind {
    /// Lowercase human-readable name, also used by [`fmt::Display`].
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
    pub(crate) kind: ComponentKind,
    pub(crate) energy: Quantity<f64>,
}

impl EnergyComponent {
    /// Which contribution this is.
    pub const fn kind(&self) -> ComponentKind {
        self.kind
    }

    /// The contribution's energy, finite and in kJ/mol when part of an
    /// [`Energy`].
    pub const fn energy(&self) -> Quantity<f64> {
        self.energy
    }
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
    ///
    /// # Errors
    ///
    /// [`EvaluationError::Unit`] if `total` is not an energy, and
    /// [`EvaluationError::NonFiniteEnergy`] if it is not finite.
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
    ///
    /// # Errors
    ///
    /// [`EvaluationError::Unit`] if a component is not an energy,
    /// [`EvaluationError::NonFiniteEnergy`] if a component or the total is not
    /// finite, and [`EvaluationError::DuplicateComponent`] if a kind repeats.
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

    /// The finite total energy in kJ/mol.
    pub fn total(&self) -> Quantity<f64> {
        self.total
    }

    /// The decomposition in the order given, each in kJ/mol; empty when the
    /// potential reports none.
    pub fn components(&self) -> &[EnergyComponent] {
        &self.components
    }

    /// The energy of one component in kJ/mol, or `None` if the decomposition
    /// does not contain `kind`.
    pub fn component(&self, kind: ComponentKind) -> Option<Quantity<f64>> {
        self.components
            .iter()
            .find(|c| c.kind == kind)
            .map(|c| c.energy)
    }
}

/// A validated energy and Cartesian gradient `dE/dx` for one structure.
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
        self.topology.shares_layout(&other.topology)
            && self.energy == other.energy
            && self.gradient == other.gradient
    }
}

impl Evaluation {
    /// Validates a backend's result for `model`, which it binds by snapshot.
    ///
    /// `gradient` holds one vector per topology atom in dense order.
    ///
    /// # Errors
    ///
    /// - [`EvaluationError::Unit`] if `gradient` is not an energy per length.
    /// - [`EvaluationError::GradientLengthMismatch`] if it does not have one
    ///   vector per atom of `model`.
    /// - [`EvaluationError::NonFiniteGradient`] naming the first atom whose
    ///   gradient is not finite.
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

    /// The exact topology snapshot of the evaluated structure.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Clones the handle to the evaluated snapshot, without copying topology
    /// data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// The validated energy.
    pub fn energy(&self) -> &Energy {
        &self.energy
    }

    /// The energy, discarding the gradient.
    pub fn into_energy(self) -> Energy {
        self.energy
    }

    /// The finite gradient `dE/dx` in kJ/(mol nm), one vector per topology atom
    /// in dense order.
    pub fn gradient(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.gradient.as_slice(), CANONICAL_GRADIENT_UNIT)
    }

    /// The gradient of one atom in kJ/(mol nm), or `None` if `model` does not
    /// share this evaluation's topology layout or has no atom `atom`.
    pub fn gradient_for(
        &self,
        model: ModelView<'_>,
        atom: InstanceAtomId,
    ) -> Option<Quantity<Vector3>> {
        if !self.topology.shares_layout(model.topology()) {
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

/// The kind of interaction whose coordinates are singular.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum InteractionKind {
    /// A harmonic bond stretch.
    Bond,
    /// A harmonic angle bend.
    Angle,
    /// A proper torsion.
    ProperTorsion,
    /// An improper torsion.
    ImproperTorsion,
    /// A nonbonded (van der Waals or electrostatic) atom pair.
    NonbondedPair,
    /// A positional or distance restraint, such as a confining wall.
    Restraint,
}

impl fmt::Display for InteractionKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Bond => "harmonic bond",
            Self::Angle => "harmonic angle",
            Self::ProperTorsion => "proper torsion",
            Self::ImproperTorsion => "improper torsion",
            Self::NonbondedPair => "nonbonded pair",
            Self::Restraint => "restraint",
        })
    }
}

/// A failed energy or gradient evaluation.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum EvaluationError {
    /// The evaluated view does not share the potential's topology layout.
    IncompatibleTopology,
    /// The potential has no evaluation policy for periodic cells.
    UnsupportedPeriodicCell,
    /// The coordinates are singular for a required interaction.
    InvalidGeometry {
        /// The interaction that is undefined.
        interaction: InteractionKind,
        /// The atoms of the interaction.
        atoms: Vec<InstanceAtomId>,
        /// Why the interaction is undefined.
        kind: SingularGeometry,
    },
    /// The potential produced a non-finite energy.
    NonFiniteEnergy,
    /// An energy decomposition named one component twice.
    DuplicateComponent(ComponentKind),
    /// The potential produced a gradient array of the wrong length.
    GradientLengthMismatch {
        /// The model's atom count.
        expected: usize,
        /// The number of gradients produced.
        actual: usize,
    },
    /// The potential produced a non-finite gradient.
    NonFiniteGradient {
        /// The atom whose gradient is not finite.
        atom: InstanceAtomId,
    },
    /// A quantity could not be converted to the canonical unit; see the
    /// source.
    Unit(UnitError),
    /// A backend reported a non-geometric failure; see the source.
    Backend(Arc<dyn std::error::Error + Send + Sync>),
}

impl PartialEq for EvaluationError {
    /// Backend failures compare equal when their messages are equal.
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::IncompatibleTopology, Self::IncompatibleTopology)
            | (Self::UnsupportedPeriodicCell, Self::UnsupportedPeriodicCell)
            | (Self::NonFiniteEnergy, Self::NonFiniteEnergy) => true,
            (
                Self::InvalidGeometry {
                    interaction,
                    atoms,
                    kind,
                },
                Self::InvalidGeometry {
                    interaction: other_interaction,
                    atoms: other_atoms,
                    kind: other_kind,
                },
            ) => interaction == other_interaction && atoms == other_atoms && kind == other_kind,
            (Self::DuplicateComponent(kind), Self::DuplicateComponent(other)) => kind == other,
            (
                Self::GradientLengthMismatch { expected, actual },
                Self::GradientLengthMismatch {
                    expected: other_expected,
                    actual: other_actual,
                },
            ) => expected == other_expected && actual == other_actual,
            (Self::NonFiniteGradient { atom }, Self::NonFiniteGradient { atom: other }) => {
                atom == other
            }
            (Self::Unit(error), Self::Unit(other)) => error == other,
            (Self::Backend(error), Self::Backend(other)) => {
                Arc::ptr_eq(error, other) || error.to_string() == other.to_string()
            }
            _ => false,
        }
    }
}

impl EvaluationError {
    /// A backend failure caused by `error`.
    pub fn backend(error: impl std::error::Error + Send + Sync + 'static) -> Self {
        Self::Backend(Arc::new(error))
    }

    /// Whether the failure is caused only by the evaluated coordinates.
    pub const fn is_invalid_geometry(&self) -> bool {
        matches!(self, Self::InvalidGeometry { .. })
    }
}

impl fmt::Display for EvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IncompatibleTopology => {
                f.write_str("model view belongs to a different topology layout than the potential")
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
                write!(f, "potential produced a non-finite gradient for {atom}")
            }
            Self::Unit(_) => f.write_str("invalid potential quantity unit"),
            Self::Backend(_) => f.write_str("potential backend evaluation failed"),
        }
    }
}

impl std::error::Error for EvaluationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            Self::Backend(error) => Some(error.as_ref()),
            _ => None,
        }
    }
}

impl From<UnitError> for EvaluationError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn backend_failures_keep_their_cause_as_the_source() {
        let cause = std::io::Error::other("model server unavailable");
        let error = EvaluationError::backend(cause);
        assert_eq!(error.to_string(), "potential backend evaluation failed");
        let source = std::error::Error::source(&error).unwrap();
        assert_eq!(source.to_string(), "model server unavailable");
        assert!(source.downcast_ref::<std::io::Error>().is_some());
        assert_eq!(
            error,
            EvaluationError::backend(std::io::Error::other("model server unavailable"))
        );
        assert_ne!(
            error,
            EvaluationError::backend(std::io::Error::other("timeout"))
        );
    }
}
