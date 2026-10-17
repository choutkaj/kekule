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

    /// The rule's SMIRKS pattern, exactly as written in the OFFXML source.
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
    /// The rule that assigned this parameter.
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    /// Equilibrium bond length, positive, in nm.
    pub const fn length(&self) -> Quantity<f64> {
        self.length
    }

    /// Non-negative force constant in kJ/(mol nm^2).
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
    /// The rule that assigned this parameter.
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    /// Equilibrium angle in rad, within `[0, pi]`.
    pub const fn angle(&self) -> Quantity<f64> {
        self.angle
    }

    /// Non-negative force constant in kJ/(mol rad^2).
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
    /// Positive integer multiplicity of the cosine term.
    pub const fn periodicity(&self) -> u32 {
        self.periodicity
    }

    /// Phase offset in rad.
    pub const fn phase(&self) -> Quantity<f64> {
        self.phase
    }

    /// Barrier height in kJ/mol, before division by [`Self::idivf`].
    pub const fn k(&self) -> Quantity<f64> {
        self.k
    }

    /// Positive divisor applied to [`Self::k`].
    ///
    /// An `auto` proper-torsion divisor is resolved at assignment to the number
    /// of proper torsions around the central bond; `auto` impropers use 3.
    pub const fn idivf(&self) -> f64 {
        self.idivf
    }
}
/// Periodic torsion: the sum of its [`TorsionTerm`]s.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) terms: Vec<TorsionTerm>,
}

impl TorsionParameter {
    /// The rule that assigned this parameter.
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    /// The Fourier terms in OFFXML index order; never empty.
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
    /// The rule that assigned this parameter.
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    /// Positive Lennard-Jones sigma in nm; an OFFXML `rmin_half` is converted to sigma.
    pub const fn sigma(&self) -> Quantity<f64> {
        self.sigma
    }

    /// Non-negative well depth epsilon in kJ/mol.
    pub const fn epsilon(&self) -> Quantity<f64> {
        self.epsilon
    }
}
/// Fixed-distance constraint between two atoms.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintParameter {
    pub(crate) source: ParameterIdentity,
    pub(crate) distance: Quantity<f64>,
}

impl ConstraintParameter {
    /// The constraint rule that assigned this parameter.
    pub fn source(&self) -> &ParameterIdentity {
        &self.source
    }

    /// Positive constrained distance in nm. A rule without an explicit
    /// distance uses the equilibrium length of the bond assigned to the pair.
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
    /// The two atoms of one molecule instance, lower molecule-local atom ID first.
    pub const fn atoms(&self) -> [InstanceAtomId; 2] {
        self.atoms
    }

    /// Factor applied to this pair's Lennard-Jones energy.
    pub const fn vdw_scale(&self) -> f64 {
        self.vdw_scale
    }

    /// Factor applied to this pair's Coulomb energy.
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
    /// Lennard-Jones cutoff distance in nm (OFFXML default 9 Å).
    pub const fn vdw_cutoff(&self) -> Quantity<f64> {
        self.vdw_cutoff
    }

    /// Width in nm of the Lennard-Jones switching region ending at the cutoff
    /// (OFFXML default 1 Å); zero disables switching.
    pub const fn vdw_switch_width(&self) -> Quantity<f64> {
        self.vdw_switch_width
    }

    /// Electrostatics cutoff distance in nm (OFFXML default 9 Å).
    pub const fn electrostatics_cutoff(&self) -> Quantity<f64> {
        self.electrostatics_cutoff
    }

    /// Width in nm of the electrostatics switching region (OFFXML default 0 Å).
    pub const fn electrostatics_switch_width(&self) -> Quantity<f64> {
        self.electrostatics_switch_width
    }

    /// Lennard-Jones scale factors for 1-2, 1-3, 1-4 and 1-5 pairs; the 1-5
    /// factor is always 1.
    pub const fn vdw_scales(&self) -> [f64; 4] {
        self.vdw_scales
    }

    /// Coulomb scale factors for 1-2, 1-3, 1-4 and 1-5 pairs; the 1-5 factor
    /// is always 1.
    pub const fn electrostatics_scales(&self) -> [f64; 4] {
        self.electrostatics_scales
    }

    /// Lennard-Jones treatment for periodic systems.
    pub const fn vdw_periodic_method(&self) -> VdwMethod {
        self.vdw_periodic_method
    }

    /// Lennard-Jones treatment for non-periodic systems.
    pub const fn vdw_nonperiodic_method(&self) -> VdwMethod {
        self.vdw_nonperiodic_method
    }

    /// Electrostatics treatment for periodic systems.
    pub const fn electrostatics_periodic_method(&self) -> ElectrostaticsMethod {
        self.electrostatics_periodic_method
    }

    /// Electrostatics treatment for non-periodic systems; always
    /// [`ElectrostaticsMethod::Coulomb`].
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

    /// Harmonic bond terms, including bonds that are also constrained.
    pub fn bonds(&self) -> &[Interaction<2, BondParameter>] {
        &self.bonds
    }

    /// Harmonic angle terms.
    pub fn angles(&self) -> &[Interaction<3, AngleParameter>] {
        &self.angles
    }

    /// Proper torsion terms, with `auto` divisors resolved.
    pub fn proper_torsions(&self) -> &[Interaction<4, TorsionParameter>] {
        &self.propers
    }

    /// Improper torsion terms: three trefoil interactions per matched improper,
    /// each with the central atom first.
    pub fn improper_torsions(&self) -> &[Interaction<4, TorsionParameter>] {
        &self.impropers
    }

    /// Distance constraints.
    pub fn constraints(&self) -> &[Interaction<2, ConstraintParameter>] {
        &self.constraints
    }

    /// One Lennard-Jones parameter per topology atom, in dense atom order.
    pub fn vdw(&self) -> &[Arc<VdwParameter>] {
        &self.vdw
    }

    /// One partial charge per topology atom in e, in dense atom order.
    pub fn charges(&self) -> &Quantity<Vec<f64>> {
        &self.charges
    }

    /// One source per molecule instance, in topology instance order.
    pub fn charge_sources(&self) -> &[crate::ChargeSource] {
        &self.charge_sources
    }

    /// Scaled nonbonded pairs at graph distance 1, 2 or 3 within each instance;
    /// every other pair uses the full 1-5 interaction.
    pub fn pair_exceptions(&self) -> &[PairException] {
        &self.exceptions
    }

    /// The force field's nonbonded policy.
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
    /// Bond rule per bonded atom pair.
    pub fn bonds(&self) -> &BTreeMap<[AtomId; 2], ParameterIdentity> {
        &self.bonds
    }

    /// Angle rule per angle, keyed with the central atom in position 1.
    pub fn angles(&self) -> &BTreeMap<[AtomId; 3], ParameterIdentity> {
        &self.angles
    }

    /// Proper torsion rule per bonded atom quadruple.
    pub fn proper_torsions(&self) -> &BTreeMap<[AtomId; 4], ParameterIdentity> {
        &self.proper_torsions
    }

    /// Improper torsion rule per center and sorted outer atoms.
    pub fn improper_torsions(&self) -> &BTreeMap<[AtomId; 4], ParameterIdentity> {
        &self.improper_torsions
    }

    /// Constraint rule per constrained atom pair; empty without a
    /// Constraints handler.
    pub fn constraints(&self) -> &BTreeMap<[AtomId; 2], ParameterIdentity> {
        &self.constraints
    }

    /// Lennard-Jones rule per atom.
    pub fn vdw(&self) -> &BTreeMap<[AtomId; 1], ParameterIdentity> {
        &self.vdw
    }
}
