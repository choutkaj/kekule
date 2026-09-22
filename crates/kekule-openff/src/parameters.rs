use kekule::{
    core::AtomId,
    topology::{InstanceAtomId, Topology},
    units::Quantity,
};
use std::{collections::BTreeMap, sync::Arc};

/// Source rule retained for reproducibility and inspection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParameterIdentity {
    /// OFFXML rule ID, or an empty string when the optional `id` was omitted.
    /// IDs are metadata; ordered SMIRKS rules determine assignment precedence.
    pub id: String,
    pub smirks: String,
}

/// Harmonic bond: `k / 2 * (r - length)^2`.
#[derive(Debug, Clone, PartialEq)]
pub struct BondParameter {
    pub source: ParameterIdentity,
    pub length: Quantity<f64>,
    pub k: Quantity<f64>,
}
/// Harmonic angle: `k / 2 * (theta - angle)^2`.
#[derive(Debug, Clone, PartialEq)]
pub struct AngleParameter {
    pub source: ParameterIdentity,
    pub angle: Quantity<f64>,
    pub k: Quantity<f64>,
}
/// One Fourier term: `k / idivf * (1 + cos(periodicity * theta - phase))`.
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionTerm {
    pub periodicity: u32,
    pub phase: Quantity<f64>,
    pub k: Quantity<f64>,
    pub idivf: f64,
}
#[derive(Debug, Clone, PartialEq)]
pub struct TorsionParameter {
    pub source: ParameterIdentity,
    pub terms: Vec<TorsionTerm>,
}
/// Lennard-Jones 12-6 parameters. Combine sigma arithmetically and epsilon geometrically.
#[derive(Debug, Clone, PartialEq)]
pub struct VdwParameter {
    pub source: ParameterIdentity,
    pub sigma: Quantity<f64>,
    pub epsilon: Quantity<f64>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct ConstraintParameter {
    pub source: ParameterIdentity,
    pub distance: Quantity<f64>,
}
/// Parameter assignment on atoms qualified by an instance in the owning topology.
#[derive(Debug, Clone, PartialEq)]
pub struct Interaction<const N: usize, P> {
    pub atoms: [InstanceAtomId; N],
    pub parameter: P,
}
/// Explicit pair exceptions for graph distances 1, 2, and 3 (1-2/1-3/1-4).
#[derive(Debug, Clone, PartialEq)]
pub struct PairException {
    pub atoms: [InstanceAtomId; 2],
    pub vdw_scale: f64,
    pub electrostatics_scale: f64,
}
/// Nonbonded policy, preserved independently of any simulation backend.
#[derive(Debug, Clone, PartialEq)]
pub struct NonbondedSettings {
    pub vdw_cutoff: Quantity<f64>,
    pub vdw_switch_width: Quantity<f64>,
    pub electrostatics_cutoff: Quantity<f64>,
    pub electrostatics_switch_width: Quantity<f64>,
    pub vdw_scales: [f64; 4],
    pub electrostatics_scales: [f64; 4],
    pub vdw_periodic_method: String,
    pub vdw_nonperiodic_method: String,
    pub electrostatics_periodic_method: String,
    pub electrostatics_nonperiodic_method: String,
}
/// Complete parameterization tied to one exact topology snapshot.
///
/// Atom parameters and charges follow `topology().atom_ids()` dense order.
/// Bond terms remain present for constrained bonds; an evaluation backend must
/// decide how to handle constrained degrees of freedom. Improper interactions
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
    pub(crate) vdw: Vec<VdwParameter>,
    pub(crate) charges: Quantity<Vec<f64>>,
    pub(crate) charge_sources: Vec<crate::ChargeSource>,
    pub(crate) exceptions: Vec<PairException>,
    pub(crate) settings: NonbondedSettings,
}
impl ParameterizedTopology {
    pub fn topology(&self) -> &Arc<Topology> {
        &self.topology
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
    pub fn vdw(&self) -> &[VdwParameter] {
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
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MoleculeLabels {
    pub bonds: BTreeMap<[AtomId; 2], ParameterIdentity>,
    pub angles: BTreeMap<[AtomId; 3], ParameterIdentity>,
    pub proper_torsions: BTreeMap<[AtomId; 4], ParameterIdentity>,
    pub improper_torsions: BTreeMap<[AtomId; 4], ParameterIdentity>,
    pub constraints: BTreeMap<[AtomId; 2], ParameterIdentity>,
    pub vdw: BTreeMap<[AtomId; 1], ParameterIdentity>,
}
