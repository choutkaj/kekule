//! OpenFF (SMIRNOFF) energies and gradients from a [`ParameterizedTopology`].
//!
//! [`OpenFfPotential::new`] lowers an assignment produced by `kekule-openff`
//! into dense harmonic, Fourier, Lennard-Jones, and Coulomb terms once. The
//! potential then evaluates any [`ModelView`] sharing that topology layout
//! snapshot through [`Potential`].
//!
//! # Hamiltonian and capabilities
//!
//! - Vacuum, nonperiodic evaluation: SMIRNOFF's nonperiodic `no-cutoff` vdW and
//!   plain Coulomb methods. Every non-excluded atom pair is evaluated, which
//!   costs `O(N^2)` per call. A force field whose nonperiodic vdW method is a
//!   cutoff is rejected at preparation. Views with a periodic cell are rejected
//!   with [`EvaluationError::UnsupportedPeriodicCell`].
//! - Lennard-Jones 12-6 with Lorentz-Berthelot mixing, Coulomb with fixed
//!   assigned charges, and the assigned 1-2/1-3/1-4 pair exceptions.
//! - Distance constraints are not applied: every bond term, including those of
//!   constrained bonds, contributes to the energy and gradient. Constraints
//!   remain in the [`ParameterizedTopology`] for dynamics consumers.
//! - Improper torsions are the three assigned trefoil terms, each contributing
//!   `k / idivf * (1 + cos(n * phi - phase))`.
//!
//! See the crate-level singular-geometry policy for coordinates at which an
//! energy or gradient is undefined.

use std::fmt;
use std::sync::Arc;

use kekule::geometry::Vector3;
use kekule::structure::ModelView;
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::{
    Quantity, UnitError, CANONICAL_ANGLE_UNIT, CANONICAL_CHARGE_UNIT, CANONICAL_ENERGY_UNIT,
    CANONICAL_FORCE_CONSTANT_UNIT, CANONICAL_GRADIENT_UNIT, CANONICAL_LENGTH_UNIT,
};
use kekule_openff::{
    ElectrostaticsMethod, Interaction, ParameterizedTopology, TorsionParameter, VdwMethod,
};

use crate::mm::{
    self, FourierTerm, Gradients, HarmonicAngle, HarmonicBond, Nonbonded, PeriodicTorsion, Singular,
};
use crate::potential::{
    ComponentKind, Energy, EnergyComponent, Evaluation, EvaluationError, InteractionKind, Potential,
};

const KINDS: [ComponentKind; mm::COMPONENTS] = [
    ComponentKind::Bonds,
    ComponentKind::Angles,
    ComponentKind::ProperTorsions,
    ComponentKind::ImproperTorsions,
    ComponentKind::VanDerWaals,
    ComponentKind::Electrostatics,
];

/// Prepared OpenFF potential bound to one topology layout.
#[derive(Debug, Clone)]
pub struct OpenFfPotential {
    topology: Arc<Topology>,
    bonds: Vec<HarmonicBond>,
    angles: Vec<HarmonicAngle>,
    propers: Vec<PeriodicTorsion>,
    impropers: Vec<PeriodicTorsion>,
    nonbonded: Nonbonded,
}

/// The energy and gradient of one OpenFF component.
#[derive(Debug, Clone, PartialEq)]
pub struct ComponentEvaluation {
    pub(crate) kind: ComponentKind,
    pub(crate) energy: Quantity<f64>,
    pub(crate) gradient: Quantity<Vec<Vector3>>,
}

impl ComponentEvaluation {
    pub const fn kind(&self) -> ComponentKind {
        self.kind
    }

    pub const fn energy(&self) -> Quantity<f64> {
        self.energy
    }

    pub fn gradient(&self) -> &Quantity<Vec<Vector3>> {
        &self.gradient
    }
}

impl OpenFfPotential {
    /// Lowers every assigned interaction into dense evaluation terms.
    ///
    /// The potential binds `parameters.topology()`; it does not copy rule
    /// identities, constraints, or charge provenance.
    pub fn new(parameters: &ParameterizedTopology) -> Result<Self, OpenFfPotentialError> {
        let topology = parameters.shared_topology();
        let settings = parameters.nonbonded_settings();
        if settings.vdw_nonperiodic_method() != VdwMethod::NoCutoff {
            return Err(OpenFfPotentialError::UnsupportedVdwMethod(
                settings.vdw_nonperiodic_method(),
            ));
        }
        if settings.electrostatics_nonperiodic_method() != ElectrostaticsMethod::Coulomb {
            return Err(OpenFfPotentialError::UnsupportedElectrostaticsMethod(
                settings.electrostatics_nonperiodic_method(),
            ));
        }
        let mut bonds = Vec::with_capacity(parameters.bonds().len());
        for term in parameters.bonds() {
            bonds.push(HarmonicBond {
                atoms: dense(&topology, term.atoms())?,
                k: finite(
                    OpenFfParameter::BondK,
                    term.parameter()
                        .k()
                        .value_in(CANONICAL_FORCE_CONSTANT_UNIT)?,
                )?,
                length: finite(
                    OpenFfParameter::BondLength,
                    term.parameter().length().value_in(CANONICAL_LENGTH_UNIT)?,
                )?,
            });
        }
        let mut angles = Vec::with_capacity(parameters.angles().len());
        let angle_k = CANONICAL_ENERGY_UNIT.try_div(CANONICAL_ANGLE_UNIT.try_powi(2)?)?;
        for term in parameters.angles() {
            angles.push(HarmonicAngle {
                atoms: dense(&topology, term.atoms())?,
                k: finite(
                    OpenFfParameter::AngleK,
                    term.parameter().k().value_in(angle_k)?,
                )?,
                angle: finite(
                    OpenFfParameter::Angle,
                    term.parameter().angle().value_in(CANONICAL_ANGLE_UNIT)?,
                )?,
            });
        }
        let propers = torsions(&topology, parameters.proper_torsions())?;
        let impropers = torsions(&topology, parameters.improper_torsions())?;

        let n = topology.atom_count();
        if parameters.vdw().len() != n {
            return Err(OpenFfPotentialError::AtomCountMismatch {
                array: PerAtomArray::VdwParameters,
                expected: n,
                actual: parameters.vdw().len(),
            });
        }
        let mut sigma = Vec::with_capacity(n);
        let mut epsilon = Vec::with_capacity(n);
        for vdw in parameters.vdw() {
            sigma.push(non_negative(
                OpenFfParameter::Sigma,
                vdw.sigma().value_in(CANONICAL_LENGTH_UNIT)?,
            )?);
            epsilon.push(non_negative(
                OpenFfParameter::Epsilon,
                vdw.epsilon().value_in(CANONICAL_ENERGY_UNIT)?,
            )?);
        }
        let mut exceptions = vec![Vec::new(); n];
        for exception in parameters.pair_exceptions() {
            let mut pair = dense(&topology, exception.atoms())?;
            pair.sort_unstable();
            if pair[0] == pair[1] {
                return Err(OpenFfPotentialError::SelfPairException);
            }
            exceptions[pair[0]].push((
                pair[1],
                finite(OpenFfParameter::VdwExceptionScale, exception.vdw_scale())?,
                finite(
                    OpenFfParameter::ElectrostaticsExceptionScale,
                    exception.electrostatics_scale(),
                )?,
            ));
        }
        for list in &mut exceptions {
            list.sort_unstable_by_key(|&(j, _, _)| j);
            if list.windows(2).any(|w| w[0].0 == w[1].0) {
                return Err(OpenFfPotentialError::DuplicatePairException);
            }
        }
        let mut potential = Self {
            topology,
            bonds,
            angles,
            propers,
            impropers,
            nonbonded: Nonbonded {
                sigma,
                epsilon,
                charge: Vec::new(),
                exceptions,
            },
        };
        potential.set_charges(parameters.charges().clone())?;
        Ok(potential)
    }

    /// Replaces the assigned partial charges, in dense topology atom order.
    ///
    /// Use this to evaluate the same valence and vdW assignment with charges
    /// from another source. Charges are not checked against formal charges.
    pub fn with_charges(
        mut self,
        charges: Quantity<Vec<f64>>,
    ) -> Result<Self, OpenFfPotentialError> {
        self.set_charges(charges)?;
        Ok(self)
    }

    fn set_charges(&mut self, charges: Quantity<Vec<f64>>) -> Result<(), OpenFfPotentialError> {
        let charges = charges.into_unit(CANONICAL_CHARGE_UNIT)?.into_value();
        if charges.len() != self.topology.atom_count() {
            return Err(OpenFfPotentialError::AtomCountMismatch {
                array: PerAtomArray::Charges,
                expected: self.topology.atom_count(),
                actual: charges.len(),
            });
        }
        for &charge in &charges {
            finite(OpenFfParameter::Charge, charge)?;
        }
        self.nonbonded.charge = charges;
        Ok(())
    }

    /// Partial charges used for electrostatics, in dense topology atom order.
    pub fn charges(&self) -> Quantity<&[f64]> {
        Quantity::new(self.nonbonded.charge.as_slice(), CANONICAL_CHARGE_UNIT)
    }

    /// Evaluates every component's energy and gradient separately.
    ///
    /// Components are returned in [`Energy::components`] order. This is a
    /// diagnostic path: it stores one gradient array per component.
    pub fn evaluate_components(
        &self,
        model: ModelView<'_>,
    ) -> Result<Vec<ComponentEvaluation>, EvaluationError> {
        let (energies, gradients) = self.run(model, Gradients::components(model.atom_count()))?;
        let Gradients::Components(gradients) = gradients else {
            unreachable!("component gradients were requested");
        };
        let mut result = Vec::with_capacity(mm::COMPONENTS);
        for ((kind, energy), gradient) in KINDS.into_iter().zip(energies).zip(gradients) {
            if !energy.is_finite() {
                return Err(EvaluationError::NonFiniteEnergy);
            }
            if let Some(index) = gradient.iter().position(|g| !g.is_finite()) {
                return Err(EvaluationError::NonFiniteGradient {
                    atom: model.topology().atom_ids()[index],
                });
            }
            result.push(ComponentEvaluation {
                kind,
                energy: Quantity::new(energy, CANONICAL_ENERGY_UNIT),
                gradient: Quantity::new(gradient, CANONICAL_GRADIENT_UNIT),
            });
        }
        Ok(result)
    }

    fn run(
        &self,
        model: ModelView<'_>,
        mut gradients: Gradients,
    ) -> Result<([f64; mm::COMPONENTS], Gradients), EvaluationError> {
        if !self.topology.shares_layout(model.topology()) {
            return Err(EvaluationError::IncompatibleTopology);
        }
        if model.cell().is_some() {
            return Err(EvaluationError::UnsupportedPeriodicCell);
        }
        let x = model.positions().values().into_value();
        let located = |error: Singular| EvaluationError::InvalidGeometry {
            interaction: error.interaction,
            atoms: error
                .atoms
                .iter()
                .map(|&i| self.topology.atom_ids()[i])
                .collect(),
            kind: error.kind,
        };
        let mut energies = [0.0; mm::COMPONENTS];
        energies[mm::BONDS] = mm::bonds(&self.bonds, x, &mut gradients).map_err(located)?;
        energies[mm::ANGLES] = mm::angles(&self.angles, x, &mut gradients).map_err(located)?;
        energies[mm::PROPER_TORSIONS] = mm::torsions(
            mm::PROPER_TORSIONS,
            InteractionKind::ProperTorsion,
            &self.propers,
            x,
            &mut gradients,
        )
        .map_err(located)?;
        energies[mm::IMPROPER_TORSIONS] = mm::torsions(
            mm::IMPROPER_TORSIONS,
            InteractionKind::ImproperTorsion,
            &self.impropers,
            x,
            &mut gradients,
        )
        .map_err(located)?;
        let (vdw, electrostatics) = self
            .nonbonded
            .evaluate(x, &mut gradients)
            .map_err(located)?;
        energies[mm::VAN_DER_WAALS] = vdw;
        energies[mm::ELECTROSTATICS] = electrostatics;
        Ok((energies, gradients))
    }

    fn energy_from(energies: [f64; mm::COMPONENTS]) -> Result<Energy, EvaluationError> {
        Energy::from_components(KINDS.into_iter().zip(energies).map(|(kind, energy)| {
            EnergyComponent {
                kind,
                energy: Quantity::new(energy, CANONICAL_ENERGY_UNIT),
            }
        }))
    }
}

impl OpenFfPotential {
    /// Clones the handle to the bound snapshot, without copying topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }
}

impl Potential for OpenFfPotential {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn evaluate(&self, model: ModelView<'_>) -> Result<Evaluation, EvaluationError> {
        let (energies, gradients) = self.run(model, Gradients::total(model.atom_count()))?;
        let Gradients::Total(gradient) = gradients else {
            unreachable!("a total gradient was requested");
        };
        Evaluation::new(
            model,
            Self::energy_from(energies)?,
            Quantity::new(gradient, CANONICAL_GRADIENT_UNIT),
        )
    }

    /// Evaluates without a gradient. Bonds of zero length and exactly linear
    /// strained angles have a defined energy and are accepted here.
    fn energy(&self, model: ModelView<'_>) -> Result<Energy, EvaluationError> {
        let (energies, _) = self.run(model, Gradients::None)?;
        Self::energy_from(energies)
    }
}

/// Dense topology indices of qualified atoms.
fn dense<const N: usize>(
    topology: &Topology,
    atoms: [InstanceAtomId; N],
) -> Result<[usize; N], OpenFfPotentialError> {
    let mut indices = [0; N];
    for (index, atom) in indices.iter_mut().zip(atoms) {
        *index = topology
            .atom_index(atom)
            .ok_or(OpenFfPotentialError::InvalidAtom(atom))?
            .index();
    }
    Ok(indices)
}

fn torsions(
    topology: &Topology,
    interactions: &[Interaction<4, TorsionParameter>],
) -> Result<Vec<PeriodicTorsion>, OpenFfPotentialError> {
    let mut result = Vec::with_capacity(interactions.len());
    for term in interactions {
        let mut terms = Vec::with_capacity(term.parameter().terms().len());
        for t in term.parameter().terms() {
            let idivf = finite(OpenFfParameter::TorsionIdivf, t.idivf())?;
            if idivf == 0.0 {
                return Err(OpenFfPotentialError::ZeroTorsionDivisor);
            }
            terms.push(FourierTerm {
                k: finite(
                    OpenFfParameter::TorsionK,
                    t.k().value_in(CANONICAL_ENERGY_UNIT)? / idivf,
                )?,
                periodicity: f64::from(t.periodicity()),
                phase: finite(
                    OpenFfParameter::TorsionPhase,
                    t.phase().value_in(CANONICAL_ANGLE_UNIT)?,
                )?,
            });
        }
        result.push(PeriodicTorsion {
            atoms: dense(topology, term.atoms())?,
            terms: terms.into_boxed_slice(),
        });
    }
    Ok(result)
}

fn finite(parameter: OpenFfParameter, value: f64) -> Result<f64, OpenFfPotentialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(OpenFfPotentialError::InvalidParameter { parameter, value })
    }
}

fn non_negative(parameter: OpenFfParameter, value: f64) -> Result<f64, OpenFfPotentialError> {
    if finite(parameter, value)? < 0.0 {
        return Err(OpenFfPotentialError::InvalidParameter { parameter, value });
    }
    Ok(value)
}

/// A lowered OpenFF parameter that must be finite (and, for Lennard-Jones
/// parameters, non-negative).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum OpenFfParameter {
    /// A bond force constant.
    BondK,
    /// A bond equilibrium length.
    BondLength,
    /// An angle force constant.
    AngleK,
    /// An angle equilibrium value.
    Angle,
    /// A torsion barrier height.
    TorsionK,
    /// A torsion phase.
    TorsionPhase,
    /// A torsion `idivf` divisor.
    TorsionIdivf,
    /// A Lennard-Jones sigma.
    Sigma,
    /// A Lennard-Jones epsilon.
    Epsilon,
    /// A pair exception's van der Waals scale.
    VdwExceptionScale,
    /// A pair exception's electrostatics scale.
    ElectrostaticsExceptionScale,
    /// An atomic partial charge.
    Charge,
}

impl fmt::Display for OpenFfParameter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::BondK => "bond k",
            Self::BondLength => "bond length",
            Self::AngleK => "angle k",
            Self::Angle => "angle",
            Self::TorsionK => "torsion k",
            Self::TorsionPhase => "torsion phase",
            Self::TorsionIdivf => "torsion idivf",
            Self::Sigma => "sigma",
            Self::Epsilon => "epsilon",
            Self::VdwExceptionScale => "vdW exception scale",
            Self::ElectrostaticsExceptionScale => "electrostatics exception scale",
            Self::Charge => "charge",
        })
    }
}

/// A per-atom parameter array that must cover the bound topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum PerAtomArray {
    /// The assigned van der Waals parameters.
    VdwParameters,
    /// The partial charges.
    Charges,
}

impl fmt::Display for PerAtomArray {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::VdwParameters => "vdW parameters",
            Self::Charges => "charges",
        })
    }
}

/// A parameterization that cannot be lowered into an evaluable potential.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum OpenFfPotentialError {
    /// The nonperiodic vdW method has no implementation in this potential;
    /// only [`VdwMethod::NoCutoff`] is implemented.
    UnsupportedVdwMethod(VdwMethod),
    /// The nonperiodic electrostatics method has no implementation in this
    /// potential; only [`ElectrostaticsMethod::Coulomb`] is implemented.
    UnsupportedElectrostaticsMethod(ElectrostaticsMethod),
    /// An interaction refers to an atom outside the bound topology.
    InvalidAtom(InstanceAtomId),
    /// A per-atom array does not cover the topology.
    AtomCountMismatch {
        /// The array that does not cover the topology.
        array: PerAtomArray,
        /// The topology's atom count.
        expected: usize,
        /// The array's length.
        actual: usize,
    },
    /// A parameter is not finite, or a Lennard-Jones parameter is negative.
    InvalidParameter {
        /// The parameter at fault.
        parameter: OpenFfParameter,
        /// Its value in the canonical unit.
        value: f64,
    },
    /// A torsion term's `idivf` divisor is zero.
    ZeroTorsionDivisor,
    /// A pair exception names the same atom twice.
    SelfPairException,
    /// An atom pair has more than one pair exception.
    DuplicatePairException,
    /// A parameter could not be converted to the canonical unit; see the
    /// source.
    Unit(UnitError),
}

impl fmt::Display for OpenFfPotentialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedVdwMethod(method) => write!(
                f,
                "unsupported nonperiodic vdW method {:?}; only no-cutoff vdW and Coulomb electrostatics are implemented",
                method.smirnoff_name()
            ),
            Self::UnsupportedElectrostaticsMethod(method) => write!(
                f,
                "unsupported nonperiodic Electrostatics method {:?}; only no-cutoff vdW and Coulomb electrostatics are implemented",
                method.smirnoff_name()
            ),
            Self::InvalidAtom(atom) => write!(f, "interaction refers to unknown {atom}"),
            Self::AtomCountMismatch {
                array,
                expected,
                actual,
            } => write!(
                f,
                "{array} cover {actual} atoms, but the topology has {expected}"
            ),
            Self::InvalidParameter { parameter, value } if value.is_finite() => {
                write!(f, "invalid {parameter}: {value} is negative")
            }
            Self::InvalidParameter { parameter, value } => {
                write!(f, "invalid {parameter}: {value} is not finite")
            }
            Self::ZeroTorsionDivisor => {
                f.write_str("invalid torsion idivf: the divisor must be nonzero")
            }
            Self::SelfPairException => {
                f.write_str("invalid pair exception: an exception must name two distinct atoms")
            }
            Self::DuplicatePairException => {
                f.write_str("invalid pair exception: an atom pair has more than one exception")
            }
            Self::Unit(_) => f.write_str("invalid parameter unit"),
        }
    }
}

impl std::error::Error for OpenFfPotentialError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for OpenFfPotentialError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
