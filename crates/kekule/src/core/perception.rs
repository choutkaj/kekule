use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use super::{
    checked_fixed_id_collection_len, AtomId, BondId, Molecule, StereoDescriptor, StereoElementId,
};

static EMPTY_CIP_DESCRIPTORS: BTreeMap<StereoElementId, StereoDescriptor> = BTreeMap::new();

/// The valence model used to produce installed valence state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValenceModel {
    RdkitLike,
}

/// The aromaticity model used to produce installed aromaticity state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    RdkitLike,
}

/// The algorithm used to select an installed ring basis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingBasisModel {
    /// Kekule's deterministic Figueras/SSSR-like basis selection.
    FiguerasSssrLike,
}

/// Cycle membership over the stable atom and bond slots of a molecule.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingMembership {
    pub(crate) atom_flags: Vec<bool>,
    pub(crate) bond_flags: Vec<bool>,
}

impl RingMembership {
    /// Constructs detached ring membership over complete stable atom and bond slots.
    ///
    /// Slot lengths and live references are checked when the containing
    /// [`PerceptionState`] is installed on a molecule.
    pub fn from_slot_flags(atom_flags: Vec<bool>, bond_flags: Vec<bool>) -> Self {
        Self {
            atom_flags,
            bond_flags,
        }
    }

    /// Returns the complete stable atom-slot flags, including tombstones.
    pub fn atom_slot_flags(&self) -> &[bool] {
        &self.atom_flags
    }

    /// Returns the complete stable bond-slot flags, including tombstones.
    pub fn bond_slot_flags(&self) -> &[bool] {
        &self.bond_flags
    }

    pub fn atom_in_ring(&self, atom: AtomId) -> bool {
        self.atom_flags.get(atom.index()).copied().unwrap_or(false)
    }

    pub fn bond_in_ring(&self, bond: BondId) -> bool {
        self.bond_flags.get(bond.index()).copied().unwrap_or(false)
    }

    pub fn ring_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        (0..=u32::MAX)
            .zip(self.atom_flags.iter())
            .filter_map(|(raw, in_ring)| in_ring.then_some(AtomId::new(raw)))
    }

    pub fn ring_bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        (0..=u32::MAX)
            .zip(self.bond_flags.iter())
            .filter_map(|(raw, in_ring)| in_ring.then_some(BondId::new(raw)))
    }
}

/// One ring in an installed deterministic ring basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

/// A deterministic ring basis installed in molecule perception state.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingSet {
    pub(crate) rings: Vec<Ring>,
}

impl RingSet {
    /// Constructs a detached deterministic ring basis.
    ///
    /// Ring references and graph coherence are checked when the containing
    /// [`PerceptionState`] is installed on a molecule.
    pub fn from_rings(rings: Vec<Ring>) -> Self {
        Self { rings }
    }

    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    pub fn len(&self) -> usize {
        self.rings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }
}

/// One installed ring basis and its optional algorithm provenance.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingBasisState {
    pub(super) model: Option<RingBasisModel>,
    pub(super) rings: RingSet,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValencePerceptionState {
    pub(super) model: Option<ValenceModel>,
    pub(super) implicit_hydrogens: BTreeMap<AtomId, u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingPerceptionState {
    pub(super) membership: RingMembership,
    pub(super) basis: Option<RingBasisState>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AromaticityPerceptionState {
    pub(super) model: AromaticityModel,
    pub(super) atoms: BTreeSet<AtomId>,
    pub(super) bonds: BTreeSet<BondId>,
}

/// Installed stereo-derived perception, currently CIP descriptor assignments.
///
/// Presence of this section records that CIP derivation completed successfully,
/// including when no descriptors were assigned.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StereoPerceptionState {
    pub(super) cip_descriptors: BTreeMap<StereoElementId, StereoDescriptor>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PerceptionState {
    pub(super) valence: Option<ValencePerceptionState>,
    pub(super) rings: Option<RingPerceptionState>,
    pub(super) aromaticity: Option<AromaticityPerceptionState>,
    pub(super) stereo: Option<StereoPerceptionState>,
}

impl ValencePerceptionState {
    /// Returns the named model, or `None` for an installed model-neutral section.
    pub const fn model(&self) -> Option<ValenceModel> {
        self.model
    }

    /// Iterates every installed implicit-hydrogen assignment.
    pub fn implicit_hydrogens(
        &self,
    ) -> impl ExactSizeIterator<Item = (AtomId, u8)> + DoubleEndedIterator + '_ {
        self.implicit_hydrogens
            .iter()
            .map(|(atom, count)| (*atom, *count))
    }
}

impl RingPerceptionState {
    /// Returns complete ring membership over stable atom and bond slots.
    pub const fn membership(&self) -> &RingMembership {
        &self.membership
    }

    /// Returns the selected ring basis and its provenance, when installed.
    pub const fn basis(&self) -> Option<&RingBasisState> {
        self.basis.as_ref()
    }

    /// Returns the installed deterministic ring basis, when present.
    pub const fn ring_set(&self) -> Option<&RingSet> {
        match &self.basis {
            Some(basis) => Some(&basis.rings),
            None => None,
        }
    }
}

impl RingBasisState {
    /// Constructs a detached ring basis with optional algorithm provenance.
    pub const fn new(model: Option<RingBasisModel>, rings: RingSet) -> Self {
        Self { model, rings }
    }

    /// Returns the named basis-selection model, if provenance is known.
    pub const fn model(&self) -> Option<RingBasisModel> {
        self.model
    }

    /// Returns the selected deterministic ring set.
    pub const fn ring_set(&self) -> &RingSet {
        &self.rings
    }
}

impl AromaticityPerceptionState {
    /// Returns the model that perceived this semantic aromaticity membership.
    pub const fn model(&self) -> AromaticityModel {
        self.model
    }

    /// Iterates every installed aromatic atom.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomId> + DoubleEndedIterator + '_ {
        self.atoms.iter().copied()
    }

    /// Iterates every installed aromatic bond.
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = BondId> + DoubleEndedIterator + '_ {
        self.bonds.iter().copied()
    }
}

impl StereoPerceptionState {
    /// Iterates every installed CIP descriptor assignment.
    pub fn cip_descriptors(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoElementId, StereoDescriptor)> + DoubleEndedIterator + '_
    {
        self.cip_descriptors
            .iter()
            .map(|(element, descriptor)| (*element, *descriptor))
    }

    /// Returns the installed descriptor for one stereo element, when assigned.
    pub fn cip_descriptor(&self, element: StereoElementId) -> Option<StereoDescriptor> {
        self.cip_descriptors.get(&element).copied()
    }
}

impl PerceptionState {
    /// Starts construction of a detached complete perception state.
    pub fn builder() -> PerceptionStateBuilder {
        PerceptionStateBuilder::default()
    }

    /// Returns the exact installed valence section, including model-neutral state.
    pub const fn valence_state(&self) -> Option<&ValencePerceptionState> {
        self.valence.as_ref()
    }

    /// Returns the exact installed ring section.
    pub const fn ring_state(&self) -> Option<&RingPerceptionState> {
        self.rings.as_ref()
    }

    /// Returns the exact installed aromaticity section.
    pub const fn aromaticity_state(&self) -> Option<&AromaticityPerceptionState> {
        self.aromaticity.as_ref()
    }

    /// Returns the exact installed stereo-perception section.
    pub const fn stereo_state(&self) -> Option<&StereoPerceptionState> {
        self.stereo.as_ref()
    }

    /// Iterates every installed CIP descriptor assignment.
    pub fn cip_descriptors(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoElementId, StereoDescriptor)> + DoubleEndedIterator + '_
    {
        self.stereo
            .as_ref()
            .map(|state| &state.cip_descriptors)
            .unwrap_or(&EMPTY_CIP_DESCRIPTORS)
            .iter()
            .map(|(element, descriptor)| (*element, *descriptor))
    }

    pub fn has_valence(&self) -> bool {
        self.valence.is_some()
    }

    pub fn has_rings(&self) -> bool {
        self.rings.is_some()
    }

    pub fn has_aromaticity(&self) -> bool {
        self.aromaticity.is_some()
    }

    /// Returns whether stereo perception, including an empty CIP result, is installed.
    pub fn has_stereo(&self) -> bool {
        self.stereo.is_some()
    }

    /// Returns whether the installed stereo section contains any CIP assignments.
    ///
    /// This is a payload query, not a section-presence query. Use
    /// [`Self::has_stereo`] to distinguish an absent section from a successful
    /// CIP derivation with no assignments.
    pub fn has_cip_descriptors(&self) -> bool {
        self.stereo
            .as_ref()
            .is_some_and(|state| !state.cip_descriptors.is_empty())
    }

    /// Returns the named valence model, if the installed section has provenance.
    ///
    /// Use [`Self::has_valence`] or [`Self::valence_state`] to distinguish an
    /// absent section from installed model-neutral valence.
    pub fn valence_model(&self) -> Option<ValenceModel> {
        self.valence.as_ref().and_then(|state| state.model)
    }

    pub fn implicit_hydrogens(&self, atom: AtomId) -> Option<u8> {
        self.valence
            .as_ref()
            .and_then(|state| state.implicit_hydrogens.get(&atom).copied())
    }

    pub fn ring_membership(&self) -> Option<&RingMembership> {
        self.rings.as_ref().map(|state| &state.membership)
    }

    pub fn ring_set(&self) -> Option<&RingSet> {
        self.rings.as_ref().and_then(RingPerceptionState::ring_set)
    }

    /// Returns the named ring-basis model, if a named basis is installed.
    ///
    /// Use [`Self::ring_state`] and [`RingPerceptionState::basis`] to
    /// distinguish absent ring perception, membership-only state, and an
    /// installed model-neutral basis.
    pub fn ring_basis_model(&self) -> Option<RingBasisModel> {
        self.rings
            .as_ref()
            .and_then(RingPerceptionState::basis)
            .and_then(RingBasisState::model)
    }

    pub fn aromaticity_model(&self) -> Option<AromaticityModel> {
        self.aromaticity.as_ref().map(|state| state.model)
    }

    pub fn atom_is_aromatic(&self, atom: AtomId) -> Option<bool> {
        self.aromaticity
            .as_ref()
            .map(|state| state.atoms.contains(&atom))
    }

    pub fn bond_is_aromatic(&self, bond: BondId) -> Option<bool> {
        self.aromaticity
            .as_ref()
            .map(|state| state.bonds.contains(&bond))
    }

    pub fn cip_descriptor(&self, element: StereoElementId) -> Option<StereoDescriptor> {
        self.stereo
            .as_ref()
            .and_then(|state| state.cip_descriptor(element))
    }
}

/// Constructs one detached perception state without mutating a molecule.
#[derive(Debug, Clone, Default)]
pub struct PerceptionStateBuilder {
    state: PerceptionState,
}

impl PerceptionStateBuilder {
    /// Installs an exact valence section on the detached state.
    pub fn with_valence(
        mut self,
        model: Option<ValenceModel>,
        assignments: Vec<(AtomId, u8)>,
    ) -> std::result::Result<Self, PerceptionStateBuildError> {
        check_perception_component_capacity(
            assignments.len(),
            PerceptionStateComponent::ImplicitHydrogens,
        )?;
        let mut implicit_hydrogens = BTreeMap::new();
        for (atom, count) in assignments {
            if implicit_hydrogens.insert(atom, count).is_some() {
                return Err(PerceptionStateBuildError::DuplicateImplicitHydrogen(atom));
            }
        }
        self.state.valence = Some(ValencePerceptionState {
            model,
            implicit_hydrogens,
        });
        Ok(self)
    }

    /// Installs exact ring membership and an optional deterministic ring basis.
    pub fn with_rings(mut self, membership: RingMembership, basis: Option<RingBasisState>) -> Self {
        self.state.rings = Some(RingPerceptionState { membership, basis });
        self
    }

    /// Installs exact model-perceived aromaticity membership.
    pub fn with_aromaticity(
        mut self,
        model: AromaticityModel,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
    ) -> std::result::Result<Self, PerceptionStateBuildError> {
        check_perception_component_capacity(atoms.len(), PerceptionStateComponent::AromaticAtoms)?;
        check_perception_component_capacity(bonds.len(), PerceptionStateComponent::AromaticBonds)?;
        let mut atom_set = BTreeSet::new();
        for atom in atoms {
            if !atom_set.insert(atom) {
                return Err(PerceptionStateBuildError::DuplicateAromaticAtom(atom));
            }
        }
        let mut bond_set = BTreeSet::new();
        for bond in bonds {
            if !bond_set.insert(bond) {
                return Err(PerceptionStateBuildError::DuplicateAromaticBond(bond));
            }
        }
        self.state.aromaticity = Some(AromaticityPerceptionState {
            model,
            atoms: atom_set,
            bonds: bond_set,
        });
        Ok(self)
    }

    /// Installs the complete stereo-perception section.
    ///
    /// An empty assignment list intentionally installs a present, empty
    /// section, recording that CIP derivation completed without assignments.
    pub fn with_cip_descriptors(
        mut self,
        assignments: Vec<(StereoElementId, StereoDescriptor)>,
    ) -> std::result::Result<Self, PerceptionStateBuildError> {
        check_perception_component_capacity(
            assignments.len(),
            PerceptionStateComponent::CipDescriptors,
        )?;
        let mut descriptors = BTreeMap::new();
        for (element, descriptor) in assignments {
            if descriptors.insert(element, descriptor).is_some() {
                return Err(PerceptionStateBuildError::DuplicateCipDescriptor(element));
            }
        }
        self.state.stereo = Some(StereoPerceptionState {
            cip_descriptors: descriptors,
        });
        Ok(self)
    }

    /// Finishes detached construction. Molecule-specific validation occurs on install.
    pub fn build(self) -> PerceptionState {
        self.state
    }
}

/// Capacity-bounded components accepted by [`PerceptionStateBuilder`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerceptionStateComponent {
    /// Complete atom-slot ring-membership flags.
    RingAtomSlots,
    /// Complete bond-slot ring-membership flags.
    RingBondSlots,
    /// Rings in one installed deterministic basis.
    Rings,
    /// Atom references in one installed ring.
    RingAtoms,
    /// Bond references in one installed ring.
    RingBonds,
    /// Atom-wise implicit-hydrogen assignments.
    ImplicitHydrogens,
    /// Aromatic atom references.
    AromaticAtoms,
    /// Aromatic bond references.
    AromaticBonds,
    /// Stereo-element CIP descriptor assignments.
    CipDescriptors,
}

impl fmt::Display for PerceptionStateComponent {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::RingAtomSlots => "ring atom-slot flags",
            Self::RingBondSlots => "ring bond-slot flags",
            Self::Rings => "installed rings",
            Self::RingAtoms => "installed ring atom references",
            Self::RingBonds => "installed ring bond references",
            Self::ImplicitHydrogens => "implicit-hydrogen assignments",
            Self::AromaticAtoms => "aromatic atom references",
            Self::AromaticBonds => "aromatic bond references",
            Self::CipDescriptors => "CIP descriptor assignments",
        })
    }
}

/// Errors constructing a detached perception state.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionStateBuildError {
    /// A component contains more entries than fixed-width stable IDs can address.
    ComponentCapacityExceeded(PerceptionStateComponent),
    /// One atom has more than one implicit-hydrogen assignment.
    DuplicateImplicitHydrogen(AtomId),
    /// One aromatic atom is listed more than once.
    DuplicateAromaticAtom(AtomId),
    /// One aromatic bond is listed more than once.
    DuplicateAromaticBond(BondId),
    /// One stereo element has more than one CIP descriptor assignment.
    DuplicateCipDescriptor(StereoElementId),
}

impl fmt::Display for PerceptionStateBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ComponentCapacityExceeded(component) => {
                write!(formatter, "{component} capacity exceeded")
            }
            Self::DuplicateImplicitHydrogen(atom) => {
                write!(
                    formatter,
                    "duplicate implicit-hydrogen assignment for {atom}"
                )
            }
            Self::DuplicateAromaticAtom(atom) => {
                write!(formatter, "duplicate aromatic atom reference {atom}")
            }
            Self::DuplicateAromaticBond(bond) => {
                write!(formatter, "duplicate aromatic bond reference {bond}")
            }
            Self::DuplicateCipDescriptor(element) => {
                write!(
                    formatter,
                    "duplicate CIP descriptor assignment for {element}"
                )
            }
        }
    }
}

impl std::error::Error for PerceptionStateBuildError {}

/// Structural reasons an installed ring is not a simple graph cycle.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MalformedRingReason {
    /// Fewer than three atoms were supplied.
    TooFewAtoms,
    /// The atom and bond reference counts differ.
    AtomBondCountMismatch,
    /// One atom occurs more than once.
    DuplicateAtom,
    /// One bond occurs more than once.
    DuplicateBond,
    /// A ring bond endpoint is absent from the ring atom set.
    BondEndpointOutsideRing,
    /// Ring-local degrees do not describe one simple cycle.
    NotSimpleCycle,
    /// Ring-local bonds form more than one connected component.
    Disconnected,
}

impl fmt::Display for MalformedRingReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::TooFewAtoms => "fewer than three atoms",
            Self::AtomBondCountMismatch => "atom and bond counts differ",
            Self::DuplicateAtom => "duplicate atom reference",
            Self::DuplicateBond => "duplicate bond reference",
            Self::BondEndpointOutsideRing => "bond endpoint outside the ring atom set",
            Self::NotSimpleCycle => "ring-local degrees do not form a simple cycle",
            Self::Disconnected => "ring-local bonds are disconnected",
        })
    }
}

/// Errors validating a detached perception state against one molecule.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionStateInstallError {
    /// One component exceeds the fixed-width stable-ID capacity.
    ComponentCapacityExceeded(PerceptionStateComponent),
    /// An atom reference is not live.
    InvalidAtomId(AtomId),
    /// A bond reference is not live.
    InvalidBondId(BondId),
    /// A stereo-element reference is not live.
    InvalidStereoElementId(StereoElementId),
    /// Ring atom flags do not cover the molecule's exact stable atom slots.
    RingAtomSlotCountMismatch {
        /// Required molecule atom-slot count.
        expected: usize,
        /// Supplied flag count.
        actual: usize,
    },
    /// Ring bond flags do not cover the molecule's exact stable bond slots.
    RingBondSlotCountMismatch {
        /// Required molecule bond-slot count.
        expected: usize,
        /// Supplied flag count.
        actual: usize,
    },
    /// An installed ring is not a coherent simple graph cycle.
    MalformedRing {
        /// Zero-based ring position in the detached basis.
        ring: usize,
        /// Structural failure.
        reason: MalformedRingReason,
    },
    /// A ring atom is absent from installed membership.
    InconsistentRingAtomMembership(AtomId),
    /// A ring bond is absent from installed membership.
    InconsistentRingBondMembership(BondId),
    /// Membership marks an atom not covered by the installed basis.
    RingAtomMissingFromBasis(AtomId),
    /// Membership marks a bond not covered by the installed basis.
    RingBondMissingFromBasis(BondId),
}

impl fmt::Display for PerceptionStateInstallError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ComponentCapacityExceeded(component) => {
                write!(formatter, "{component} capacity exceeded")
            }
            Self::InvalidAtomId(atom) => write!(formatter, "invalid perception atom id: {atom}"),
            Self::InvalidBondId(bond) => write!(formatter, "invalid perception bond id: {bond}"),
            Self::InvalidStereoElementId(element) => {
                write!(formatter, "invalid perception stereo-element id: {element}")
            }
            Self::RingAtomSlotCountMismatch { expected, actual } => write!(
                formatter,
                "ring atom-slot count mismatch: expected {expected}, got {actual}"
            ),
            Self::RingBondSlotCountMismatch { expected, actual } => write!(
                formatter,
                "ring bond-slot count mismatch: expected {expected}, got {actual}"
            ),
            Self::MalformedRing { ring, reason } => {
                write!(formatter, "malformed installed ring {ring}: {reason}")
            }
            Self::InconsistentRingAtomMembership(atom) => {
                write!(
                    formatter,
                    "installed ring atom {atom} is not marked in membership"
                )
            }
            Self::InconsistentRingBondMembership(bond) => {
                write!(
                    formatter,
                    "installed ring bond {bond} is not marked in membership"
                )
            }
            Self::RingAtomMissingFromBasis(atom) => {
                write!(
                    formatter,
                    "ring membership atom {atom} is absent from the ring basis"
                )
            }
            Self::RingBondMissingFromBasis(bond) => {
                write!(
                    formatter,
                    "ring membership bond {bond} is absent from the ring basis"
                )
            }
        }
    }
}

impl std::error::Error for PerceptionStateInstallError {}

fn check_perception_component_capacity(
    length: usize,
    component: PerceptionStateComponent,
) -> std::result::Result<(), PerceptionStateBuildError> {
    checked_fixed_id_collection_len(0, length)
        .map_err(|_| PerceptionStateBuildError::ComponentCapacityExceeded(component))
}

fn check_install_component_capacity(
    length: usize,
    component: PerceptionStateComponent,
) -> std::result::Result<(), PerceptionStateInstallError> {
    checked_fixed_id_collection_len(0, length)
        .map_err(|_| PerceptionStateInstallError::ComponentCapacityExceeded(component))
}

pub(super) fn validate_perception_state(
    molecule: &Molecule,
    state: &PerceptionState,
) -> std::result::Result<(), PerceptionStateInstallError> {
    if let Some(valence) = &state.valence {
        for atom in valence.implicit_hydrogens.keys().copied() {
            if molecule
                .atoms
                .get(atom.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionStateInstallError::InvalidAtomId(atom));
            }
        }
    }

    if let Some(ring_state) = &state.rings {
        validate_ring_state(molecule, ring_state)?;
    }

    if let Some(aromaticity) = &state.aromaticity {
        for atom in aromaticity.atoms.iter().copied() {
            if molecule
                .atoms
                .get(atom.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionStateInstallError::InvalidAtomId(atom));
            }
        }
        for bond in aromaticity.bonds.iter().copied() {
            if molecule
                .bonds
                .get(bond.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionStateInstallError::InvalidBondId(bond));
            }
        }
    }

    if let Some(stereo) = &state.stereo {
        for element in stereo.cip_descriptors.keys().copied() {
            if molecule
                .stereo_elements
                .get(element.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionStateInstallError::InvalidStereoElementId(element));
            }
        }
    }
    Ok(())
}

fn validate_ring_state(
    molecule: &Molecule,
    state: &RingPerceptionState,
) -> std::result::Result<(), PerceptionStateInstallError> {
    let membership = &state.membership;
    let atom_slots = membership.atom_slot_flags().len();
    let bond_slots = membership.bond_slot_flags().len();
    check_install_component_capacity(atom_slots, PerceptionStateComponent::RingAtomSlots)?;
    check_install_component_capacity(bond_slots, PerceptionStateComponent::RingBondSlots)?;
    if atom_slots != molecule.atoms.len() {
        return Err(PerceptionStateInstallError::RingAtomSlotCountMismatch {
            expected: molecule.atoms.len(),
            actual: atom_slots,
        });
    }
    if bond_slots != molecule.bonds.len() {
        return Err(PerceptionStateInstallError::RingBondSlotCountMismatch {
            expected: molecule.bonds.len(),
            actual: bond_slots,
        });
    }
    for (raw, in_ring) in (0..=u32::MAX).zip(membership.atom_slot_flags()) {
        if *in_ring
            && molecule
                .atoms
                .get(AtomId::new(raw).index())
                .and_then(Option::as_ref)
                .is_none()
        {
            return Err(PerceptionStateInstallError::InvalidAtomId(AtomId::new(raw)));
        }
    }
    for (raw, in_ring) in (0..=u32::MAX).zip(membership.bond_slot_flags()) {
        if *in_ring
            && molecule
                .bonds
                .get(BondId::new(raw).index())
                .and_then(Option::as_ref)
                .is_none()
        {
            return Err(PerceptionStateInstallError::InvalidBondId(BondId::new(raw)));
        }
    }

    let Some(basis) = &state.basis else {
        return Ok(());
    };
    let ring_set = &basis.rings;
    check_install_component_capacity(ring_set.rings().len(), PerceptionStateComponent::Rings)?;

    let mut covered_atoms = BTreeSet::new();
    let mut covered_bonds = BTreeSet::new();
    for (ring_index, ring) in ring_set.rings().iter().enumerate() {
        validate_installed_ring(molecule, ring_index, ring)?;
        for atom in ring.atoms.iter().copied() {
            if !membership.atom_in_ring(atom) {
                return Err(PerceptionStateInstallError::InconsistentRingAtomMembership(
                    atom,
                ));
            }
            covered_atoms.insert(atom);
        }
        for bond in ring.bonds.iter().copied() {
            if !membership.bond_in_ring(bond) {
                return Err(PerceptionStateInstallError::InconsistentRingBondMembership(
                    bond,
                ));
            }
            covered_bonds.insert(bond);
        }
    }
    for atom in membership.ring_atom_ids() {
        if !covered_atoms.contains(&atom) {
            return Err(PerceptionStateInstallError::RingAtomMissingFromBasis(atom));
        }
    }
    for bond in membership.ring_bond_ids() {
        if !covered_bonds.contains(&bond) {
            return Err(PerceptionStateInstallError::RingBondMissingFromBasis(bond));
        }
    }
    Ok(())
}

fn validate_installed_ring(
    molecule: &Molecule,
    ring_index: usize,
    ring: &Ring,
) -> std::result::Result<(), PerceptionStateInstallError> {
    check_install_component_capacity(ring.atoms.len(), PerceptionStateComponent::RingAtoms)?;
    check_install_component_capacity(ring.bonds.len(), PerceptionStateComponent::RingBonds)?;
    if ring.atoms.len() < 3 {
        return Err(PerceptionStateInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::TooFewAtoms,
        });
    }
    if ring.atoms.len() != ring.bonds.len() {
        return Err(PerceptionStateInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::AtomBondCountMismatch,
        });
    }
    let atom_set = ring.atoms.iter().copied().collect::<BTreeSet<_>>();
    if atom_set.len() != ring.atoms.len() {
        return Err(PerceptionStateInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::DuplicateAtom,
        });
    }
    let bond_set = ring.bonds.iter().copied().collect::<BTreeSet<_>>();
    if bond_set.len() != ring.bonds.len() {
        return Err(PerceptionStateInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::DuplicateBond,
        });
    }

    let mut degrees = BTreeMap::<AtomId, usize>::new();
    let mut adjacency = BTreeMap::<AtomId, Vec<AtomId>>::new();
    for atom in atom_set.iter().copied() {
        if molecule
            .atoms
            .get(atom.index())
            .and_then(Option::as_ref)
            .is_none()
        {
            return Err(PerceptionStateInstallError::InvalidAtomId(atom));
        }
        degrees.insert(atom, 0);
        adjacency.insert(atom, Vec::new());
    }
    for bond_id in bond_set {
        let Some(bond) = molecule.bonds.get(bond_id.index()).and_then(Option::as_ref) else {
            return Err(PerceptionStateInstallError::InvalidBondId(bond_id));
        };
        if !atom_set.contains(&bond.a) || !atom_set.contains(&bond.b) {
            return Err(PerceptionStateInstallError::MalformedRing {
                ring: ring_index,
                reason: MalformedRingReason::BondEndpointOutsideRing,
            });
        }
        *degrees
            .get_mut(&bond.a)
            .expect("validated ring endpoint should have a degree") += 1;
        *degrees
            .get_mut(&bond.b)
            .expect("validated ring endpoint should have a degree") += 1;
        adjacency
            .get_mut(&bond.a)
            .expect("validated ring endpoint should have adjacency")
            .push(bond.b);
        adjacency
            .get_mut(&bond.b)
            .expect("validated ring endpoint should have adjacency")
            .push(bond.a);
    }
    if degrees.values().any(|degree| *degree != 2) {
        return Err(PerceptionStateInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::NotSimpleCycle,
        });
    }
    let start = ring.atoms[0];
    let mut seen = BTreeSet::from([start]);
    let mut stack = vec![start];
    while let Some(atom) = stack.pop() {
        for neighbor in &adjacency[&atom] {
            if seen.insert(*neighbor) {
                stack.push(*neighbor);
            }
        }
    }
    if seen.len() != atom_set.len() {
        return Err(PerceptionStateInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::Disconnected,
        });
    }
    Ok(())
}
