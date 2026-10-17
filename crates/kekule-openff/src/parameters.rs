use kekule::{
    molecule::AtomId,
    topology::{InstanceAtomId, Topology},
    units::Quantity,
};
use std::{collections::BTreeMap, sync::Arc};

/// Source rule retained for reproducibility and inspection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParameterIdentity {
    pub(crate) id: String,
    pub(crate) smirks: String,
}

impl ParameterIdentity {
    /// OFFXML rule ID, or an empty string when the optional `id` was omitted.
    /// IDs are metadata; ordered SMIRKS rules determine assignment precedence.
    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn smirks(&self) -> &str {
        &self.smirks
    }
}

/// Harmonic bond: `k / 2 * (r - length)^2`.
#[derive(Debug, Clone, PartialEq)]
pub struct BondParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) length: Quantity<f64>,
    pub(crate) k: Quantity<f64>,
}

impl BondParameter {
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    pub const fn length(&self) -> Quantity<f64> {
        self.length
    }

    pub const fn k(&self) -> Quantity<f64> {
        self.k
    }
}
/// Harmonic angle: `k / 2 * (theta - angle)^2`.
#[derive(Debug, Clone, PartialEq)]
pub struct AngleParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) angle: Quantity<f64>,
    pub(crate) k: Quantity<f64>,
}

impl AngleParameter {
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    pub const fn angle(&self) -> Quantity<f64> {
        self.angle
    }

    pub const fn k(&self) -> Quantity<f64> {
        self.k
    }
}
/// One Fourier term: `k / idivf * (1 + cos(periodicity * theta - phase))`.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionTerm {
    pub(crate) periodicity: u32,
    pub(crate) phase: Quantity<f64>,
    pub(crate) k: Quantity<f64>,
    pub(crate) idivf: f64,
}

impl TorsionTerm {
    pub const fn periodicity(&self) -> u32 {
        self.periodicity
    }

    pub const fn phase(&self) -> Quantity<f64> {
        self.phase
    }

    pub const fn k(&self) -> Quantity<f64> {
        self.k
    }

    pub const fn idivf(&self) -> f64 {
        self.idivf
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) terms: Vec<TorsionTerm>,
}

impl TorsionParameter {
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    pub fn terms(&self) -> &[TorsionTerm] {
        &self.terms
    }
}
/// Lennard-Jones 12-6 parameters. Combine sigma arithmetically and epsilon geometrically.
#[derive(Debug, Clone, PartialEq)]
pub struct VdwParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) sigma: Quantity<f64>,
    pub(crate) epsilon: Quantity<f64>,
}

impl VdwParameter {
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    pub const fn sigma(&self) -> Quantity<f64> {
        self.sigma
    }

    pub const fn epsilon(&self) -> Quantity<f64> {
        self.epsilon
    }
}
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) distance: Quantity<f64>,
}

impl ConstraintParameter {
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    pub const fn distance(&self) -> Quantity<f64> {
        self.distance
    }
}
/// Parameter assignment on atoms qualified by an instance in the owning topology.
///
/// Interactions assigned by the same rule share one parameter allocation, so
/// rule identities are not copied for every interaction of a large system.
#[derive(Debug, Clone, PartialEq)]
pub struct Interaction<const N: usize, P> {
    pub(crate) atoms: [InstanceAtomId; N],
    pub(crate) parameter: Arc<P>,
}

impl<const N: usize, P> Interaction<N, P> {
    /// The interacting atoms, in the order the parameter's functional form uses.
    pub const fn atoms(&self) -> [InstanceAtomId; N] {
        self.atoms
    }

    /// The assigned parameter, shared by every interaction of the same rule.
    pub fn parameter(&self) -> &Arc<P> {
        &self.parameter
    }
}
/// Explicit pair exceptions for graph distances 1, 2, and 3 (1-2/1-3/1-4).
#[derive(Debug, Clone, PartialEq)]
pub struct PairException {
    pub(crate) atoms: [InstanceAtomId; 2],
    pub(crate) vdw_scale: f64,
    pub(crate) electrostatics_scale: f64,
}

impl PairException {
    pub const fn atoms(&self) -> [InstanceAtomId; 2] {
        self.atoms
    }

    pub const fn vdw_scale(&self) -> f64 {
        self.vdw_scale
    }

    pub const fn electrostatics_scale(&self) -> f64 {
        self.electrostatics_scale
    }
}
/// SMIRNOFF vdW treatment for one boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum VdwMethod {
    /// Truncate at the handler cutoff, with its optional switching width.
    Cutoff,
    /// Evaluate every non-excluded pair without truncation.
    NoCutoff,
}

impl VdwMethod {
    /// The method's spelling in OFFXML.
    pub const fn smirnoff_name(self) -> &'static str {
        match self {
            Self::Cutoff => "cutoff",
            Self::NoCutoff => "no-cutoff",
        }
    }
}

/// SMIRNOFF electrostatics treatment for one boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum ElectrostaticsMethod {
    /// Plain Coulomb interactions between every non-excluded pair.
    Coulomb,
    /// Ewald summation with conducting (tin-foil) boundary conditions; PME is
    /// one accepted implementation.
    Ewald3DConductingBoundary,
}

impl ElectrostaticsMethod {
    /// The method's spelling in OFFXML.
    pub const fn smirnoff_name(self) -> &'static str {
        match self {
            Self::Coulomb => "Coulomb",
            Self::Ewald3DConductingBoundary => "Ewald3D-ConductingBoundary",
        }
    }
}

/// Nonbonded policy, preserved independently of any simulation backend.
///
/// Scales are indexed by graph distance: 1-2, 1-3, 1-4, and 1-5.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedSettings {
    pub(crate) vdw_cutoff: Quantity<f64>,
    pub(crate) vdw_switch_width: Quantity<f64>,
    pub(crate) electrostatics_cutoff: Quantity<f64>,
    pub(crate) electrostatics_switch_width: Quantity<f64>,
    pub(crate) vdw_scales: [f64; 4],
    pub(crate) electrostatics_scales: [f64; 4],
    pub(crate) vdw_periodic_method: VdwMethod,
    pub(crate) vdw_nonperiodic_method: VdwMethod,
    pub(crate) electrostatics_periodic_method: ElectrostaticsMethod,
    pub(crate) electrostatics_nonperiodic_method: ElectrostaticsMethod,
}

impl NonbondedSettings {
    pub const fn vdw_cutoff(&self) -> Quantity<f64> {
        self.vdw_cutoff
    }

    pub const fn vdw_switch_width(&self) -> Quantity<f64> {
        self.vdw_switch_width
    }

    pub const fn electrostatics_cutoff(&self) -> Quantity<f64> {
        self.electrostatics_cutoff
    }

    pub const fn electrostatics_switch_width(&self) -> Quantity<f64> {
        self.electrostatics_switch_width
    }

    pub const fn vdw_scales(&self) -> [f64; 4] {
        self.vdw_scales
    }

    pub const fn electrostatics_scales(&self) -> [f64; 4] {
        self.electrostatics_scales
    }

    pub const fn vdw_periodic_method(&self) -> VdwMethod {
        self.vdw_periodic_method
    }

    pub const fn vdw_nonperiodic_method(&self) -> VdwMethod {
        self.vdw_nonperiodic_method
    }

    pub const fn electrostatics_periodic_method(&self) -> ElectrostaticsMethod {
        self.electrostatics_periodic_method
    }

    pub const fn electrostatics_nonperiodic_method(&self) -> ElectrostaticsMethod {
        self.electrostatics_nonperiodic_method
    }
}
/// Complete parameterization tied to one exact topology snapshot.
///
/// Atom parameters and charges follow `topology().atom_ids()` dense order.
/// Bond terms remain present for constrained bonds; an evaluation backend must
/// decide how to handle constrained degrees of freedom. Atom vdW parameters
/// assigned by the same rule share one allocation. Improper interactions
/// are the three SMIRNOFF trefoil terms, with the central atom FIRST, matching
/// Interchange's exported four-atom dihedral convention. This differs from the
/// central atom's second position in an improper SMIRKS pattern or label key.
#[derive(Debug, Clone)]
pub struct ParameterizedTopology {
    pub(crate) topology: Arc<Topology>,
    pub(crate) bonds: Vec<Interaction<2, BondParameter>>,
    pub(crate) angles: Vec<Interaction<3, AngleParameter>>,
    pub(crate) propers: Vec<Interaction<4, TorsionParameter>>,
    pub(crate) impropers: Vec<Interaction<4, TorsionParameter>>,
    pub(crate) constraints: Vec<Interaction<2, ConstraintParameter>>,
    pub(crate) vdw: Vec<Arc<VdwParameter>>,
    pub(crate) charges: Quantity<Vec<f64>>,
    pub(crate) charge_sources: Vec<crate::ChargeSource>,
    pub(crate) exceptions: Vec<PairException>,
    pub(crate) settings: NonbondedSettings,
}
impl ParameterizedTopology {
    /// The exact topology snapshot that was parameterized.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Clones the handle to the parameterized snapshot, without copying
    /// topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }
    pub fn bonds(&self) -> &[Interaction<2, BondParameter>] {
        &self.bonds
    }
    pub fn angles(&self) -> &[Interaction<3, AngleParameter>] {
        &self.angles
    }
    pub fn proper_torsions(&self) -> &[Interaction<4, TorsionParameter>] {
        &self.propers
    }
    pub fn improper_torsions(&self) -> &[Interaction<4, TorsionParameter>] {
        &self.impropers
    }
    pub fn constraints(&self) -> &[Interaction<2, ConstraintParameter>] {
        &self.constraints
    }
    pub fn vdw(&self) -> &[Arc<VdwParameter>] {
        &self.vdw
    }
    pub fn charges(&self) -> &Quantity<Vec<f64>> {
        &self.charges
    }
    /// One source per molecule instance, in topology instance order.
    pub fn charge_sources(&self) -> &[crate::ChargeSource] {
        &self.charge_sources
    }
    pub fn pair_exceptions(&self) -> &[PairException] {
        &self.exceptions
    }
    pub fn nonbonded_settings(&self) -> &NonbondedSettings {
        &self.settings
    }
}

/// Rule labels before numerical assignment; keys use molecule-local identities.
/// Improper keys retain the center in position 1 and sort the three outer atoms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MoleculeLabels {
    pub(crate) bonds: BTreeMap<[AtomId; 2], ParameterIdentity>,
    pub(crate) angles: BTreeMap<[AtomId; 3], ParameterIdentity>,
    pub(crate) proper_torsions: BTreeMap<[AtomId; 4], ParameterIdentity>,
    pub(crate) improper_torsions: BTreeMap<[AtomId; 4], ParameterIdentity>,
    pub(crate) constraints: BTreeMap<[AtomId; 2], ParameterIdentity>,
    pub(crate) vdw: BTreeMap<[AtomId; 1], ParameterIdentity>,
}

impl MoleculeLabels {
    pub fn bonds(&self) -> &BTreeMap<[AtomId; 2], ParameterIdentity> {
        &self.bonds
    }

    pub fn angles(&self) -> &BTreeMap<[AtomId; 3], ParameterIdentity> {
        &self.angles
    }

    pub fn proper_torsions(&self) -> &BTreeMap<[AtomId; 4], ParameterIdentity> {
        &self.proper_torsions
    }

    pub fn improper_torsions(&self) -> &BTreeMap<[AtomId; 4], ParameterIdentity> {
        &self.improper_torsions
    }

    pub fn constraints(&self) -> &BTreeMap<[AtomId; 2], ParameterIdentity> {
        &self.constraints
    }

    pub fn vdw(&self) -> &BTreeMap<[AtomId; 1], ParameterIdentity> {
        &self.vdw
    }
}
