use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use super::{
    checked_fixed_id_collection_len, AtomId, BondId, Molecule, StereoDescriptor, StereoElementId,
    StereoElementKind,
};
use super::{ConjugationModel, ConjugationPerception, ResonanceGroup, ResonancePerception};

/// The valence model used to produce installed valence state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ValenceModel {
    /// RDKit-like allowed valences and hydrogen inference on localized chemistry.
    ///
    /// Zero-order and dative bonds contribute zero at both endpoints. Dative
    /// acceptor valence therefore differs from RDKit's directional bond model.
    RdkitLike,
}

/// The aromaticity model used to produce installed aromaticity state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum AromaticityModel {
    /// RDKit's default aromaticity model (RDKit 2026.03.3), in which fused
    /// rings are searched only as single rings and pairs.
    RdkitLike,
    /// MDL model as implemented by RDKit (the OpenFF OEAroModel_MDL profile).
    Mdl,
}

/// How valence perception treats states RDKit's strict property cache rejects.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum ValenceValidation {
    /// Reject invalid valence and implicit-hydrogen states.
    #[default]
    Strict,
    /// Follow RDKit's permissive rules: excess valence implies zero
    /// additional hydrogens instead of an error.
    Permissive,
}

/// The requested algorithm for selecting an installed ring set.
///
/// Whether the selection fell back to depth-first cycles is recorded
/// separately, by [`RingBasisState::used_depth_first_fallback`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum RingBasisModel {
    /// RDKit-like Figueras selection and symmetrization. When the candidate
    /// search is incomplete, a deterministic depth-first cycle fallback
    /// chooses the rings instead.
    ///
    /// This follows the selected-ring model rather than promising a minimum
    /// cycle basis or that the selected rings span the complete cycle space.
    FiguerasSssrLike,
}

impl AromaticityModel {
    /// Whether the model applies the MDL candidate and electron-count rules.
    pub(crate) const fn is_mdl(self) -> bool {
        match self {
            Self::RdkitLike => false,
            Self::Mdl => true,
        }
    }
}

/// Cycle membership indexed by a molecule's atom and bond IDs.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingMembership {
    pub(crate) atom_flags: Vec<bool>,
    pub(crate) bond_flags: Vec<bool>,
}

impl RingMembership {
    /// Constructs detached ring membership with one flag per atom and bond ID.
    ///
    /// Lengths and references are checked when the containing [`Perception`]
    /// is installed on a molecule.
    pub fn from_flags(atom_flags: Vec<bool>, bond_flags: Vec<bool>) -> Self {
        Self {
            atom_flags,
            bond_flags,
        }
    }

    /// Returns one ring-membership flag per atom, indexed by [`AtomId::index`].
    pub fn atom_flags(&self) -> &[bool] {
        &self.atom_flags
    }

    /// Returns one ring-membership flag per bond, indexed by [`BondId::index`].
    pub fn bond_flags(&self) -> &[bool] {
        &self.bond_flags
    }

    /// Whether an atom is in a ring; `false` for an ID beyond the flags.
    pub fn atom_in_ring(&self, atom: AtomId) -> bool {
        self.atom_flags.get(atom.index()).copied().unwrap_or(false)
    }

    /// Whether a bond is in a ring; `false` for an ID beyond the flags.
    pub fn bond_in_ring(&self, bond: BondId) -> bool {
        self.bond_flags.get(bond.index()).copied().unwrap_or(false)
    }

    /// Ring atoms in ascending ID order.
    pub fn ring_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        (0..=u32::MAX)
            .zip(self.atom_flags.iter())
            .filter_map(|(raw, in_ring)| in_ring.then_some(AtomId::new(raw)))
    }

    /// Ring bonds in ascending ID order.
    pub fn ring_bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        (0..=u32::MAX)
            .zip(self.bond_flags.iter())
            .filter_map(|(raw, in_ring)| in_ring.then_some(BondId::new(raw)))
    }
}

/// One ring in an installed deterministic selected ring set.
///
/// An installed ring is one simple cycle of at least three atoms, with as many
/// bonds as atoms.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ring {
    /// The ring atoms; perception lists them in order around the cycle.
    pub atoms: Vec<AtomId>,
    /// The ring bonds; perception lists them in ascending ID order.
    pub bonds: Vec<BondId>,
}

/// A deterministic selected ring set installed in molecule perception state.
///
/// The RDKit-like model adds symmetric alternatives to its Figueras selection
/// and may fall back to depth-first cycles. The result is not guaranteed to be
/// a minimum cycle basis or to span the complete cycle space. Cycle membership
/// is available separately through [`RingMembership`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingSet {
    pub(crate) rings: Vec<Ring>,
}

impl RingSet {
    /// Constructs a detached selected ring set.
    ///
    /// Ring references and graph coherence are checked when the containing
    /// [`Perception`] is installed on a molecule.
    pub fn from_rings(rings: Vec<Ring>) -> Self {
        Self { rings }
    }

    /// The selected rings.
    pub fn rings(&self) -> &[Ring] {
        &self.rings
    }

    /// The number of selected rings.
    pub fn len(&self) -> usize {
        self.rings.len()
    }

    /// Whether no rings are selected, as for an acyclic molecule.
    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }
}

/// One installed ring basis and its provenance: the requested selection model
/// and whether the depth-first fallback replaced an incomplete candidate
/// search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingBasisState {
    pub(super) model: RingBasisModel,
    pub(super) depth_first_fallback: bool,
    pub(super) rings: RingSet,
}

/// Installed valence: the inferred implicit-hydrogen counts and the model and
/// validation policy that produced them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValencePerception {
    pub(super) model: ValenceModel,
    pub(super) validation: ValenceValidation,
    pub(super) inferred_hydrogens: BTreeMap<AtomId, u8>,
}

/// Installed ring membership and, when requested, a selected ring basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingPerception {
    pub(super) membership: RingMembership,
    pub(super) basis: Option<RingBasisState>,
}

/// Installed aromatic atoms and bonds and the model that perceived them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AromaticityPerception {
    pub(super) model: AromaticityModel,
    pub(super) atoms: BTreeSet<AtomId>,
    pub(super) bonds: BTreeSet<BondId>,
}

/// Installed stereo-derived perception, currently CIP descriptor assignments.
///
/// Presence of this section records that CIP derivation completed successfully,
/// including when no descriptors were assigned.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CipPerception {
    pub(super) cip_descriptors: BTreeMap<StereoElementId, StereoDescriptor>,
}

/// Reconstructible chemistry derived from a molecule's represented graph.
///
/// Each perception step installs one section, and each section owns its
/// queries and records the model or policy that produced it:
/// [`Self::valence`], [`Self::rings`], [`Self::aromaticity`],
/// [`Self::conjugation`], [`Self::resonance`], and the CIP state
/// [`Self::cip`]. An absent section has not been perceived, or a chemical
/// edit invalidated it.
///
/// ```
/// use kekule::smiles;
/// let mut benzene = smiles::to_molecules("c1ccccc1")?.remove(0);
/// benzene.perceive()?;
/// let perception = benzene.perception();
/// let aromaticity = perception.aromaticity().expect("perceived");
/// assert_eq!(aromaticity.atoms().len(), 6);
/// assert_eq!(perception.rings().and_then(|rings| rings.ring_set()).map(|set| set.len()), Some(1));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Perception {
    pub(super) valence: Option<ValencePerception>,
    pub(super) rings: Option<RingPerception>,
    pub(super) aromaticity: Option<AromaticityPerception>,
    pub(crate) conjugation: Option<ConjugationPerception>,
    pub(crate) resonance: Option<ResonancePerception>,
    pub(super) stereo: Option<CipPerception>,
}

impl ValencePerception {
    /// The valence model that inferred these counts.
    pub const fn model(&self) -> ValenceModel {
        self.model
    }

    /// The validation policy the counts were inferred under: strict
    /// perception rejected invalid valence states, permissive perception
    /// gave them zero additional hydrogens.
    pub const fn validation(&self) -> ValenceValidation {
        self.validation
    }

    /// The model's inferred count for one atom. It is the implicit count of
    /// an atom whose count is inferred and zero for a fixed count; chemical
    /// callers use [`Molecule::implicit_hydrogens`], which resolves both.
    pub fn inferred_hydrogen_count(&self, atom: AtomId) -> Option<u8> {
        self.inferred_hydrogens.get(&atom).copied()
    }

    /// Iterates every installed inferred-hydrogen assignment.
    pub fn inferred_hydrogens(
        &self,
    ) -> impl ExactSizeIterator<Item = (AtomId, u8)> + DoubleEndedIterator + '_ {
        self.inferred_hydrogens
            .iter()
            .map(|(atom, count)| (*atom, *count))
    }
}

impl RingPerception {
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
    /// A ring basis selected by `model` from a complete candidate search.
    pub const fn new(model: RingBasisModel, rings: RingSet) -> Self {
        Self {
            model,
            depth_first_fallback: false,
            rings,
        }
    }

    /// A ring basis requested from `model` whose incomplete candidate search
    /// fell back to deterministic depth-first cycles.
    pub const fn with_depth_first_fallback(model: RingBasisModel, rings: RingSet) -> Self {
        Self {
            model,
            depth_first_fallback: true,
            rings,
        }
    }

    /// The requested basis-selection model.
    pub const fn model(&self) -> RingBasisModel {
        self.model
    }

    /// Whether the candidate search was incomplete and the depth-first
    /// fallback chose the rings. A fallback basis is unsuitable for MDL
    /// aromaticity and exact SMARTS ring counts.
    pub const fn used_depth_first_fallback(&self) -> bool {
        self.depth_first_fallback
    }

    /// Returns the selected deterministic ring set.
    pub const fn ring_set(&self) -> &RingSet {
        &self.rings
    }

    /// Whether the selected rings follow the SSSR-like selection that MDL
    /// aromaticity and SMARTS ring counts require.
    pub(crate) const fn is_sssr_like(&self) -> bool {
        match self.model {
            RingBasisModel::FiguerasSssrLike => !self.depth_first_fallback,
        }
    }
}

impl AromaticityPerception {
    /// Returns the model that perceived this semantic aromaticity membership.
    pub const fn model(&self) -> AromaticityModel {
        self.model
    }

    /// Whether one atom is aromatic.
    pub fn atom_is_aromatic(&self, atom: AtomId) -> bool {
        self.atoms.contains(&atom)
    }

    /// Whether one bond is aromatic.
    pub fn bond_is_aromatic(&self, bond: BondId) -> bool {
        self.bonds.contains(&bond)
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

impl CipPerception {
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

impl Perception {
    /// Starts construction of a detached complete perception state.
    pub fn builder() -> PerceptionBuilder {
        PerceptionBuilder::default()
    }

    /// Installed valence, or `None` before valence perception.
    pub const fn valence(&self) -> Option<&ValencePerception> {
        self.valence.as_ref()
    }

    /// Installed ring membership and basis, or `None` before ring perception.
    pub const fn rings(&self) -> Option<&RingPerception> {
        self.rings.as_ref()
    }

    /// Installed aromaticity, or `None` before aromaticity perception.
    pub const fn aromaticity(&self) -> Option<&AromaticityPerception> {
        self.aromaticity.as_ref()
    }

    /// Installed conjugation, including a successfully computed empty set.
    pub const fn conjugation(&self) -> Option<&ConjugationPerception> {
        self.conjugation.as_ref()
    }

    /// The explicitly prepared conjugated-group partition.
    pub const fn resonance(&self) -> Option<&ResonancePerception> {
        self.resonance.as_ref()
    }

    /// Installed CIP state, including a successful derivation that assigned
    /// no descriptors.
    pub const fn cip(&self) -> Option<&CipPerception> {
        self.stereo.as_ref()
    }
}

/// Constructs one detached perception state without mutating a molecule.
#[derive(Debug, Clone, Default)]
pub struct PerceptionBuilder {
    state: Perception,
}

impl PerceptionBuilder {
    /// Adds exact conjugation membership. Graph coherence is checked at installation.
    ///
    /// # Errors
    ///
    /// [`PerceptionBuildError::DuplicateConjugationMember`] when an atom or bond
    /// repeats, or [`PerceptionBuildError::ComponentCapacityExceeded`] when a
    /// list exceeds ID capacity.
    pub fn with_conjugation(
        mut self,
        model: ConjugationModel,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
    ) -> Result<Self, PerceptionBuildError> {
        check_perception_component_capacity(atoms.len(), PerceptionComponent::Conjugation)?;
        check_perception_component_capacity(bonds.len(), PerceptionComponent::Conjugation)?;
        let atom_set: BTreeSet<_> = atoms.iter().copied().collect();
        let bond_set: BTreeSet<_> = bonds.iter().copied().collect();
        if atom_set.len() != atoms.len() || bond_set.len() != bonds.len() {
            return Err(PerceptionBuildError::DuplicateConjugationMember);
        }
        self.state.conjugation = Some(ConjugationPerception {
            model,
            atoms: atom_set,
            bonds: bond_set,
        });
        Ok(self)
    }
    /// Adds connected groups, ordered by their smallest bond ID.
    /// Installation checks the complete partition, including empty or duplicate members.
    ///
    /// # Errors
    ///
    /// [`PerceptionBuildError::ComponentCapacityExceeded`] when the groups or a
    /// group's members exceed ID capacity.
    pub fn with_resonance_groups(
        mut self,
        mut groups: Vec<ResonanceGroup>,
    ) -> Result<Self, PerceptionBuildError> {
        check_perception_component_capacity(groups.len(), PerceptionComponent::Conjugation)?;
        for group in &groups {
            check_perception_component_capacity(
                group.atoms.len(),
                PerceptionComponent::Conjugation,
            )?;
            check_perception_component_capacity(
                group.bonds.len(),
                PerceptionComponent::Conjugation,
            )?;
        }
        groups.sort_by_key(|group| group.bonds.iter().min().copied());
        self.state.resonance = Some(ResonancePerception { groups });
        Ok(self)
    }
    /// Installs exact inferred counts with the model and validation policy
    /// that produced them. Assignments exclude graph hydrogen neighbors.
    /// Installation rejects a nonzero count on an atom whose implicit count is
    /// fixed.
    pub fn with_valence(
        mut self,
        model: ValenceModel,
        validation: ValenceValidation,
        assignments: Vec<(AtomId, u8)>,
    ) -> std::result::Result<Self, PerceptionBuildError> {
        check_perception_component_capacity(
            assignments.len(),
            PerceptionComponent::InferredHydrogens,
        )?;
        let mut inferred_hydrogens = BTreeMap::new();
        for (atom, count) in assignments {
            if inferred_hydrogens.insert(atom, count).is_some() {
                return Err(PerceptionBuildError::DuplicateInferredHydrogen(atom));
            }
        }
        self.state.valence = Some(ValencePerception {
            model,
            validation,
            inferred_hydrogens,
        });
        Ok(self)
    }

    /// Installs exact ring membership and an optional deterministic ring basis.
    pub fn with_rings(mut self, membership: RingMembership, basis: Option<RingBasisState>) -> Self {
        self.state.rings = Some(RingPerception { membership, basis });
        self
    }

    /// Installs exact model-perceived aromaticity membership.
    ///
    /// # Errors
    ///
    /// [`PerceptionBuildError::DuplicateAromaticAtom`] or
    /// [`PerceptionBuildError::DuplicateAromaticBond`] when a member repeats, or
    /// [`PerceptionBuildError::ComponentCapacityExceeded`] when a list exceeds
    /// ID capacity.
    pub fn with_aromaticity(
        mut self,
        model: AromaticityModel,
        atoms: Vec<AtomId>,
        bonds: Vec<BondId>,
    ) -> std::result::Result<Self, PerceptionBuildError> {
        check_perception_component_capacity(atoms.len(), PerceptionComponent::AromaticAtoms)?;
        check_perception_component_capacity(bonds.len(), PerceptionComponent::AromaticBonds)?;
        let mut atom_set = BTreeSet::new();
        for atom in atoms {
            if !atom_set.insert(atom) {
                return Err(PerceptionBuildError::DuplicateAromaticAtom(atom));
            }
        }
        let mut bond_set = BTreeSet::new();
        for bond in bonds {
            if !bond_set.insert(bond) {
                return Err(PerceptionBuildError::DuplicateAromaticBond(bond));
            }
        }
        self.state.aromaticity = Some(AromaticityPerception {
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
    ///
    /// # Errors
    ///
    /// [`PerceptionBuildError::DuplicateCipDescriptor`] when an element repeats,
    /// or [`PerceptionBuildError::ComponentCapacityExceeded`] when the list
    /// exceeds ID capacity.
    pub fn with_cip_descriptors(
        mut self,
        assignments: Vec<(StereoElementId, StereoDescriptor)>,
    ) -> std::result::Result<Self, PerceptionBuildError> {
        check_perception_component_capacity(
            assignments.len(),
            PerceptionComponent::CipDescriptors,
        )?;
        let mut descriptors = BTreeMap::new();
        for (element, descriptor) in assignments {
            if descriptors.insert(element, descriptor).is_some() {
                return Err(PerceptionBuildError::DuplicateCipDescriptor(element));
            }
        }
        self.state.stereo = Some(CipPerception {
            cip_descriptors: descriptors,
        });
        Ok(self)
    }

    /// Finishes detached construction. Molecule-specific validation occurs on install.
    pub fn build(self) -> Perception {
        self.state
    }
}

/// Capacity-bounded components accepted by [`PerceptionBuilder`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerceptionComponent {
    /// Conjugated atoms, bonds, or connected groups.
    Conjugation,
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
    /// Atom-wise inferred-hydrogen assignments.
    InferredHydrogens,
    /// Aromatic atom references.
    AromaticAtoms,
    /// Aromatic bond references.
    AromaticBonds,
    /// Stereo-element CIP descriptor assignments.
    CipDescriptors,
}

impl fmt::Display for PerceptionComponent {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Conjugation => "conjugation membership",
            Self::RingAtomSlots => "ring atom-slot flags",
            Self::RingBondSlots => "ring bond-slot flags",
            Self::Rings => "installed rings",
            Self::RingAtoms => "installed ring atom references",
            Self::RingBonds => "installed ring bond references",
            Self::InferredHydrogens => "inferred-hydrogen assignments",
            Self::AromaticAtoms => "aromatic atom references",
            Self::AromaticBonds => "aromatic bond references",
            Self::CipDescriptors => "CIP descriptor assignments",
        })
    }
}

/// Errors constructing a detached perception state.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionBuildError {
    /// An atom or bond was repeated in conjugation membership.
    DuplicateConjugationMember,
    /// A component contains more entries than fixed-width stable IDs can address.
    ComponentCapacityExceeded(PerceptionComponent),
    /// One atom has more than one inferred-hydrogen assignment.
    DuplicateInferredHydrogen(AtomId),
    /// One aromatic atom is listed more than once.
    DuplicateAromaticAtom(AtomId),
    /// One aromatic bond is listed more than once.
    DuplicateAromaticBond(BondId),
    /// One stereo element has more than one CIP descriptor assignment.
    DuplicateCipDescriptor(StereoElementId),
}

impl fmt::Display for PerceptionBuildError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DuplicateConjugationMember => formatter.write_str("duplicate conjugation member"),
            Self::ComponentCapacityExceeded(component) => {
                write!(formatter, "{component} capacity exceeded")
            }
            Self::DuplicateInferredHydrogen(atom) => {
                write!(
                    formatter,
                    "duplicate inferred-hydrogen assignment for {atom}"
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

impl std::error::Error for PerceptionBuildError {}

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

/// Why installed conjugation or resonance groups are inconsistent.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DelocalizationIssue {
    /// Resonance groups are installed without conjugation.
    GroupsWithoutConjugation,
    /// Conjugation is installed without aromaticity.
    ConjugationWithoutAromaticity,
    /// The conjugated atoms are not exactly the endpoints of the conjugated
    /// bonds.
    ConjugatedAtomsNotBondEndpoints,
    /// An aromatic bond is not conjugated.
    AromaticBondNotConjugated,
    /// The zero-based resonance group `group` has no bonds.
    EmptyGroup {
        /// The group position.
        group: usize,
    },
    /// A bond of group `group` is not conjugated or belongs to an earlier
    /// group.
    GroupBondsNotPartition {
        /// The group position.
        group: usize,
    },
    /// An atom of group `group` is repeated or belongs to an earlier group.
    OverlappingGroupAtom {
        /// The group position.
        group: usize,
    },
    /// The atoms of group `group` are not exactly the endpoints of its bonds.
    GroupAtomsNotBondEndpoints {
        /// The group position.
        group: usize,
    },
    /// The bonds of group `group` are not connected.
    DisconnectedGroup {
        /// The group position.
        group: usize,
    },
    /// The groups do not cover every conjugated atom and bond.
    IncompleteGroupPartition,
}

impl fmt::Display for DelocalizationIssue {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::GroupsWithoutConjugation => formatter.write_str("groups require conjugation"),
            Self::ConjugationWithoutAromaticity => {
                formatter.write_str("conjugation requires aromaticity")
            }
            Self::ConjugatedAtomsNotBondEndpoints => {
                formatter.write_str("conjugated atoms must be exactly bond endpoints")
            }
            Self::AromaticBondNotConjugated => {
                formatter.write_str("aromatic bonds must be conjugated")
            }
            Self::EmptyGroup { group } => write!(formatter, "group {group} is empty"),
            Self::GroupBondsNotPartition { group } => write!(
                formatter,
                "group {group} bonds do not partition conjugation"
            ),
            Self::OverlappingGroupAtom { group } => {
                write!(
                    formatter,
                    "group {group} has a duplicate or overlapping atom"
                )
            }
            Self::GroupAtomsNotBondEndpoints { group } => write!(
                formatter,
                "group {group} atoms are not exactly its bond endpoints"
            ),
            Self::DisconnectedGroup { group } => write!(formatter, "group {group} is disconnected"),
            Self::IncompleteGroupPartition => {
                formatter.write_str("incomplete conjugated-group partition")
            }
        }
    }
}

/// Errors validating a detached perception state against one molecule.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerceptionInstallError {
    /// Conjugation or its group partition is structurally inconsistent.
    InvalidDelocalization(DelocalizationIssue),
    /// One component exceeds the fixed-width stable-ID capacity.
    ComponentCapacityExceeded(PerceptionComponent),
    /// An atom reference is not live.
    InvalidAtomId(AtomId),
    /// Valence assigns inferred hydrogens to an atom whose implicit count is fixed.
    InferredHydrogensOnFixedAtom(AtomId),
    /// A bond reference is not live.
    InvalidBondId(BondId),
    /// A stereo-element reference is not live.
    InvalidStereoElementId(StereoElementId),
    /// A descriptor uses a different stereo geometry than its represented focus.
    IncompatibleStereoDescriptor {
        /// The stereo element being described.
        element: StereoElementId,
        /// The incompatible descriptor.
        descriptor: StereoDescriptor,
    },
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

impl fmt::Display for PerceptionInstallError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDelocalization(issue) => {
                write!(formatter, "invalid delocalization state: {issue}")
            }
            Self::ComponentCapacityExceeded(component) => {
                write!(formatter, "{component} capacity exceeded")
            }
            Self::InvalidAtomId(atom) => write!(formatter, "invalid perception atom id: {atom}"),
            Self::InferredHydrogensOnFixedAtom(atom) => write!(
                formatter,
                "perception infers hydrogens on {atom}, whose implicit count is fixed"
            ),
            Self::InvalidBondId(bond) => write!(formatter, "invalid perception bond id: {bond}"),
            Self::InvalidStereoElementId(element) => {
                write!(formatter, "invalid perception stereo-element id: {element}")
            }
            Self::IncompatibleStereoDescriptor {
                element,
                descriptor,
            } => write!(
                formatter,
                "CIP descriptor {descriptor:?} is incompatible with {element}"
            ),
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
                    "installed ring {atom} is not marked in membership"
                )
            }
            Self::InconsistentRingBondMembership(bond) => {
                write!(
                    formatter,
                    "installed ring {bond} is not marked in membership"
                )
            }
            Self::RingAtomMissingFromBasis(atom) => {
                write!(
                    formatter,
                    "ring membership {atom} is absent from the ring basis"
                )
            }
            Self::RingBondMissingFromBasis(bond) => {
                write!(
                    formatter,
                    "ring membership {bond} is absent from the ring basis"
                )
            }
        }
    }
}

impl std::error::Error for PerceptionInstallError {}

fn check_perception_component_capacity(
    length: usize,
    component: PerceptionComponent,
) -> std::result::Result<(), PerceptionBuildError> {
    checked_fixed_id_collection_len(0, length)
        .map_err(|_| PerceptionBuildError::ComponentCapacityExceeded(component))
}

fn check_install_component_capacity(
    length: usize,
    component: PerceptionComponent,
) -> std::result::Result<(), PerceptionInstallError> {
    checked_fixed_id_collection_len(0, length)
        .map_err(|_| PerceptionInstallError::ComponentCapacityExceeded(component))
}

pub(super) fn validate_perception(
    molecule: &Molecule,
    state: &Perception,
) -> std::result::Result<(), PerceptionInstallError> {
    validate_delocalization(molecule, state)?;
    if let Some(valence) = &state.valence {
        for (&atom, &count) in &valence.inferred_hydrogens {
            let Some(payload) = molecule
                .graph
                .atoms
                .get(atom.index())
                .and_then(Option::as_ref)
            else {
                return Err(PerceptionInstallError::InvalidAtomId(atom));
            };
            if count != 0 && !payload.hydrogens.is_inferred() {
                return Err(PerceptionInstallError::InferredHydrogensOnFixedAtom(atom));
            }
        }
    }

    if let Some(ring_state) = &state.rings {
        validate_ring_state(molecule, ring_state)?;
    }

    if let Some(aromaticity) = &state.aromaticity {
        for atom in aromaticity.atoms.iter().copied() {
            if molecule
                .graph
                .atoms
                .get(atom.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionInstallError::InvalidAtomId(atom));
            }
        }
        for bond in aromaticity.bonds.iter().copied() {
            if molecule
                .graph
                .bonds
                .get(bond.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionInstallError::InvalidBondId(bond));
            }
        }
    }

    if let Some(stereo) = &state.stereo {
        for (&element, &descriptor) in &stereo.cip_descriptors {
            let represented = molecule
                .graph
                .stereo_elements
                .get(element.index())
                .and_then(Option::as_ref)
                .ok_or(PerceptionInstallError::InvalidStereoElementId(element))?;
            let compatible = matches!(
                (&represented.kind, descriptor),
                (
                    StereoElementKind::Tetrahedral(_),
                    StereoDescriptor::R
                        | StereoDescriptor::S
                        | StereoDescriptor::LowerR
                        | StereoDescriptor::LowerS
                ) | (
                    StereoElementKind::DoubleBond(_),
                    StereoDescriptor::E
                        | StereoDescriptor::Z
                        | StereoDescriptor::SeqCis
                        | StereoDescriptor::SeqTrans
                ) | (
                    StereoElementKind::Axis(_),
                    StereoDescriptor::M
                        | StereoDescriptor::P
                        | StereoDescriptor::LowerM
                        | StereoDescriptor::LowerP
                )
            );
            if !compatible {
                return Err(PerceptionInstallError::IncompatibleStereoDescriptor {
                    element,
                    descriptor,
                });
            }
        }
    }
    Ok(())
}

fn validate_delocalization(
    mol: &Molecule,
    state: &Perception,
) -> Result<(), PerceptionInstallError> {
    let bad = PerceptionInstallError::InvalidDelocalization;
    let Some(conjugation) = &state.conjugation else {
        return if state.resonance.is_some() {
            Err(bad(DelocalizationIssue::GroupsWithoutConjugation))
        } else {
            Ok(())
        };
    };
    let aromaticity = state
        .aromaticity
        .as_ref()
        .ok_or(bad(DelocalizationIssue::ConjugationWithoutAromaticity))?;
    let mut endpoints = BTreeSet::new();
    for &atom in &conjugation.atoms {
        mol.atom(atom)
            .ok_or(PerceptionInstallError::InvalidAtomId(atom))?;
    }
    for &id in &conjugation.bonds {
        let bond = mol
            .bond(id)
            .ok_or(PerceptionInstallError::InvalidBondId(id))?;
        endpoints.extend([bond.a(), bond.b()]);
    }
    if endpoints != conjugation.atoms {
        return Err(bad(DelocalizationIssue::ConjugatedAtomsNotBondEndpoints));
    }
    if !aromaticity.bonds.is_subset(&conjugation.bonds) {
        return Err(bad(DelocalizationIssue::AromaticBondNotConjugated));
    }
    if let Some(resonance) = &state.resonance {
        let mut atoms = BTreeSet::new();
        let mut bonds = BTreeSet::new();
        for (index, group) in resonance.groups.iter().enumerate() {
            if group.bonds.is_empty() {
                return Err(bad(DelocalizationIssue::EmptyGroup { group: index }));
            }
            let mut group_ends = BTreeSet::new();
            for &id in &group.bonds {
                if !conjugation.bonds.contains(&id) || !bonds.insert(id) {
                    return Err(bad(DelocalizationIssue::GroupBondsNotPartition {
                        group: index,
                    }));
                }
                let bond = mol
                    .bond(id)
                    .ok_or(PerceptionInstallError::InvalidBondId(id))?;
                group_ends.extend([bond.a(), bond.b()]);
            }
            for &id in &group.atoms {
                if !atoms.insert(id) {
                    return Err(bad(DelocalizationIssue::OverlappingGroupAtom {
                        group: index,
                    }));
                }
            }
            if group_ends != group.atoms.iter().copied().collect() {
                return Err(bad(DelocalizationIssue::GroupAtomsNotBondEndpoints {
                    group: index,
                }));
            }
            let member_bonds: BTreeSet<_> = group.bonds.iter().copied().collect();
            let mut reached = BTreeSet::new();
            let mut todo = vec![group.atoms[0]];
            while let Some(atom) = todo.pop() {
                if reached.insert(atom) {
                    for (id, bond) in mol.incident_bonds(atom).expect("validated atom") {
                        if member_bonds.contains(&id) {
                            todo.push(bond.other_atom(atom));
                        }
                    }
                }
            }
            if reached != group_ends {
                return Err(bad(DelocalizationIssue::DisconnectedGroup { group: index }));
            }
        }
        if atoms != conjugation.atoms || bonds != conjugation.bonds {
            return Err(bad(DelocalizationIssue::IncompleteGroupPartition));
        }
    }
    Ok(())
}

fn validate_ring_state(
    molecule: &Molecule,
    state: &RingPerception,
) -> std::result::Result<(), PerceptionInstallError> {
    let membership = &state.membership;
    let atom_slots = membership.atom_flags().len();
    let bond_slots = membership.bond_flags().len();
    check_install_component_capacity(atom_slots, PerceptionComponent::RingAtomSlots)?;
    check_install_component_capacity(bond_slots, PerceptionComponent::RingBondSlots)?;
    if atom_slots != molecule.graph.atoms.len() {
        return Err(PerceptionInstallError::RingAtomSlotCountMismatch {
            expected: molecule.graph.atoms.len(),
            actual: atom_slots,
        });
    }
    if bond_slots != molecule.graph.bonds.len() {
        return Err(PerceptionInstallError::RingBondSlotCountMismatch {
            expected: molecule.graph.bonds.len(),
            actual: bond_slots,
        });
    }
    for (raw, in_ring) in (0..=u32::MAX).zip(membership.atom_flags()) {
        if *in_ring
            && molecule
                .graph
                .atoms
                .get(AtomId::new(raw).index())
                .and_then(Option::as_ref)
                .is_none()
        {
            return Err(PerceptionInstallError::InvalidAtomId(AtomId::new(raw)));
        }
    }
    for (raw, in_ring) in (0..=u32::MAX).zip(membership.bond_flags()) {
        if *in_ring
            && molecule
                .graph
                .bonds
                .get(BondId::new(raw).index())
                .and_then(Option::as_ref)
                .is_none()
        {
            return Err(PerceptionInstallError::InvalidBondId(BondId::new(raw)));
        }
    }

    let Some(basis) = &state.basis else {
        return Ok(());
    };
    let ring_set = &basis.rings;
    check_install_component_capacity(ring_set.rings().len(), PerceptionComponent::Rings)?;

    let mut covered_atoms = BTreeSet::new();
    let mut covered_bonds = BTreeSet::new();
    for (ring_index, ring) in ring_set.rings().iter().enumerate() {
        validate_installed_ring(molecule, ring_index, ring)?;
        for atom in ring.atoms.iter().copied() {
            if !membership.atom_in_ring(atom) {
                return Err(PerceptionInstallError::InconsistentRingAtomMembership(atom));
            }
            covered_atoms.insert(atom);
        }
        for bond in ring.bonds.iter().copied() {
            if !membership.bond_in_ring(bond) {
                return Err(PerceptionInstallError::InconsistentRingBondMembership(bond));
            }
            covered_bonds.insert(bond);
        }
    }
    for atom in membership.ring_atom_ids() {
        if !covered_atoms.contains(&atom) {
            return Err(PerceptionInstallError::RingAtomMissingFromBasis(atom));
        }
    }
    for bond in membership.ring_bond_ids() {
        if !covered_bonds.contains(&bond) {
            return Err(PerceptionInstallError::RingBondMissingFromBasis(bond));
        }
    }
    Ok(())
}

fn validate_installed_ring(
    molecule: &Molecule,
    ring_index: usize,
    ring: &Ring,
) -> std::result::Result<(), PerceptionInstallError> {
    check_install_component_capacity(ring.atoms.len(), PerceptionComponent::RingAtoms)?;
    check_install_component_capacity(ring.bonds.len(), PerceptionComponent::RingBonds)?;
    if ring.atoms.len() < 3 {
        return Err(PerceptionInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::TooFewAtoms,
        });
    }
    if ring.atoms.len() != ring.bonds.len() {
        return Err(PerceptionInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::AtomBondCountMismatch,
        });
    }
    let atom_set = ring.atoms.iter().copied().collect::<BTreeSet<_>>();
    if atom_set.len() != ring.atoms.len() {
        return Err(PerceptionInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::DuplicateAtom,
        });
    }
    let bond_set = ring.bonds.iter().copied().collect::<BTreeSet<_>>();
    if bond_set.len() != ring.bonds.len() {
        return Err(PerceptionInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::DuplicateBond,
        });
    }

    let mut degrees = BTreeMap::<AtomId, usize>::new();
    let mut adjacency = BTreeMap::<AtomId, Vec<AtomId>>::new();
    for atom in atom_set.iter().copied() {
        if molecule
            .graph
            .atoms
            .get(atom.index())
            .and_then(Option::as_ref)
            .is_none()
        {
            return Err(PerceptionInstallError::InvalidAtomId(atom));
        }
        degrees.insert(atom, 0);
        adjacency.insert(atom, Vec::new());
    }
    for bond_id in bond_set {
        let Some(bond) = molecule
            .graph
            .bonds
            .get(bond_id.index())
            .and_then(Option::as_ref)
        else {
            return Err(PerceptionInstallError::InvalidBondId(bond_id));
        };
        if !atom_set.contains(&bond.a) || !atom_set.contains(&bond.b) {
            return Err(PerceptionInstallError::MalformedRing {
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
        return Err(PerceptionInstallError::MalformedRing {
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
        return Err(PerceptionInstallError::MalformedRing {
            ring: ring_index,
            reason: MalformedRingReason::Disconnected,
        });
    }
    Ok(())
}
