use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::{Deref, DerefMut};

use super::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityProvenance {
    Imported,
    Perceived(AromaticityModel),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValencePerceptionState {
    model: Option<ValenceModel>,
    implicit_hydrogens: BTreeMap<AtomId, u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingPerceptionState {
    membership: RingMembership,
    rings: Option<RingSet>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AromaticityPerceptionState {
    provenance: AromaticityProvenance,
    atoms: BTreeSet<AtomId>,
    bonds: BTreeSet<BondId>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PerceptionState {
    valence: Option<ValencePerceptionState>,
    rings: Option<RingPerceptionState>,
    aromaticity: Option<AromaticityPerceptionState>,
    cip_descriptors: BTreeMap<StereoElementId, StereoDescriptor>,
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

    /// Returns the installed deterministic ring basis, when present.
    pub const fn ring_set(&self) -> Option<&RingSet> {
        self.rings.as_ref()
    }
}

impl AromaticityPerceptionState {
    /// Returns whether aromaticity was imported or perceived with a named model.
    pub const fn provenance(&self) -> AromaticityProvenance {
        self.provenance
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

    /// Iterates every installed CIP descriptor assignment.
    pub fn cip_descriptors(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoElementId, StereoDescriptor)> + DoubleEndedIterator + '_
    {
        self.cip_descriptors
            .iter()
            .map(|(element, descriptor)| (*element, *descriptor))
    }

    pub fn has_valence(&self) -> bool {
        self.valence
            .as_ref()
            .is_some_and(|state| state.model.is_some())
    }

    pub fn has_rings(&self) -> bool {
        self.rings.is_some()
    }

    pub fn has_aromaticity(&self) -> bool {
        self.aromaticity.is_some()
    }

    pub fn has_cip_descriptors(&self) -> bool {
        !self.cip_descriptors.is_empty()
    }

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
        self.rings.as_ref().and_then(|state| state.rings.as_ref())
    }

    pub fn aromaticity_provenance(&self) -> Option<AromaticityProvenance> {
        self.aromaticity.as_ref().map(|state| state.provenance)
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
        self.cip_descriptors.get(&element).copied()
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
    pub fn with_rings(mut self, membership: RingMembership, rings: Option<RingSet>) -> Self {
        self.state.rings = Some(RingPerceptionState { membership, rings });
        self
    }

    /// Installs exact aromaticity membership and provenance.
    pub fn with_aromaticity(
        mut self,
        provenance: AromaticityProvenance,
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
            provenance,
            atoms: atom_set,
            bonds: bond_set,
        });
        Ok(self)
    }

    /// Installs the complete CIP descriptor assignment map.
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
        self.state.cip_descriptors = descriptors;
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

/// One completed connected molecular graph.
///
/// Empty and single-atom values are valid boundary cases. Build nontrivial
/// graphs with [`Molecule::builder`] and change topology transactionally with
/// [`Molecule::edit`], both of which enforce connectedness before publication.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    pub(crate) atoms: Vec<Option<Atom>>,
    pub(crate) bonds: Vec<Option<Bond>>,
    pub(crate) adjacency: Vec<Vec<BondId>>,
    pub(crate) conformers: Vec<Option<Conformer>>,
    pub(crate) stereo_elements: Vec<Option<StereoElement>>,
    pub(crate) stereo_groups: Vec<Option<StereoGroup>>,
    pub(crate) stereo_bond_marks: Vec<StereoBondMark>,
    pub(crate) props: PropMap,
    pub(crate) perception: PerceptionState,
}

pub struct AtomMut<'a> {
    molecule: &'a mut Molecule,
    id: AtomId,
    original: AtomChemistry,
}

impl Deref for AtomMut<'_> {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        self.molecule.atoms[self.id.index()]
            .as_ref()
            .expect("validated atom must remain live while borrowed")
    }
}

impl DerefMut for AtomMut<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.molecule.atoms[self.id.index()]
            .as_mut()
            .expect("validated atom must remain live while borrowed")
    }
}

impl Drop for AtomMut<'_> {
    fn drop(&mut self) {
        if AtomChemistry::from(&**self) != self.original {
            self.molecule.invalidate_topology();
        }
    }
}

pub struct BondMut<'a> {
    molecule: &'a mut Molecule,
    id: BondId,
    original: BondChemistry,
}

impl Deref for BondMut<'_> {
    type Target = Bond;

    fn deref(&self) -> &Self::Target {
        self.molecule.bonds[self.id.index()]
            .as_ref()
            .expect("validated bond must remain live while borrowed")
    }
}

impl DerefMut for BondMut<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.molecule.bonds[self.id.index()]
            .as_mut()
            .expect("validated bond must remain live while borrowed")
    }
}

impl Drop for BondMut<'_> {
    fn drop(&mut self) {
        if BondChemistry::from(&**self) != self.original {
            self.molecule.invalidate_topology();
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct AtomChemistry {
    element: Element,
    isotope: Option<u16>,
    formal_charge: i8,
    radical: Option<AtomRadical>,
    explicit_hydrogens: u8,
    no_implicit_hydrogens: bool,
}

impl From<&Atom> for AtomChemistry {
    fn from(atom: &Atom) -> Self {
        Self {
            element: atom.element,
            isotope: atom.isotope,
            formal_charge: atom.formal_charge,
            radical: atom.radical,
            explicit_hydrogens: atom.explicit_hydrogens,
            no_implicit_hydrogens: atom.no_implicit_hydrogens,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
struct BondChemistry {
    order: BondOrder,
}

impl From<&Bond> for BondChemistry {
    fn from(bond: &Bond) -> Self {
        Self { order: bond.order }
    }
}

impl Molecule {
    /// Creates the valid empty molecule boundary case.
    pub fn new() -> Self {
        Self::default()
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.iter().flatten().count()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.iter().flatten().count()
    }

    /// Returns the sum of the asserted formal charges on all live atoms.
    ///
    /// This aggregate does not require sanitization or perception.
    pub fn formal_charge(&self) -> i64 {
        self.atoms()
            .map(|(_, atom)| i64::from(atom.formal_charge))
            .sum()
    }

    /// Inserts an atom into crate-private construction/edit state.
    ///
    /// Public molecule construction goes through [`Molecule::builder`], because
    /// adding an atom to an already populated graph would temporarily violate
    /// the connected-molecule invariant.
    pub(crate) fn add_atom(&mut self, atom: Atom) -> Result<AtomId> {
        self.add_atom_at_slot(atom, self.atoms.len())
    }

    fn add_atom_at_slot(&mut self, atom: Atom, slot: usize) -> Result<AtomId> {
        let id = checked_molecule_id(slot, MoleculeIdKind::Atom, AtomId::new)?;
        debug_assert_eq!(slot, self.atoms.len());
        self.atoms.push(Some(atom));
        self.adjacency.push(Vec::new());
        self.invalidate_topology();
        Ok(id)
    }

    /// Removes an atom only inside crate-private construction/edit state.
    pub(crate) fn delete_atom(&mut self, id: AtomId) -> Result<Atom> {
        self.atom(id)?;
        let incident = self.adjacency[id.index()].clone();
        for bond_id in incident {
            if self
                .bonds
                .get(bond_id.index())
                .and_then(Option::as_ref)
                .is_some()
            {
                self.delete_bond(bond_id)?;
            }
        }
        self.adjacency[id.index()].clear();
        let atom = self.atoms[id.index()]
            .take()
            .ok_or(MoleculeError::InvalidAtomId(id))?;
        for conformer in self.conformers.iter_mut().flatten() {
            conformer.clear_position(id);
        }
        self.prune_stereo_for_atom(id);
        self.invalidate_topology();
        Ok(atom)
    }

    pub fn atom(&self, id: AtomId) -> Result<&Atom> {
        self.atoms
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidAtomId(id))
    }

    pub fn atom_mut(&mut self, id: AtomId) -> Result<AtomMut<'_>> {
        let original = AtomChemistry::from(self.atom(id)?);
        Ok(AtomMut {
            molecule: self,
            id,
            original,
        })
    }

    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)> {
        (0..=u32::MAX)
            .zip(self.atoms.iter())
            .filter_map(|(raw, atom)| atom.as_ref().map(|atom| (AtomId::new(raw), atom)))
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms().map(|(id, _)| id)
    }

    /// Inserts a bond after validating its endpoints and identifier capacity.
    ///
    /// Adding a bond cannot disconnect a valid molecule, so it remains a safe
    /// direct public topology mutation. Every failure leaves the molecule unchanged.
    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.atom(a)?;
        self.atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.bond_between(a, b)?.is_some() {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        let id = checked_molecule_id(self.bonds.len(), MoleculeIdKind::Bond, BondId::new)?;
        self.bonds.push(Some(Bond::new(a, b, order)));
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        self.invalidate_topology();
        Ok(id)
    }

    /// Removes a bond only inside crate-private construction/edit state.
    pub(crate) fn delete_bond(&mut self, id: BondId) -> Result<Bond> {
        let bond = self
            .bonds
            .get_mut(id.index())
            .and_then(Option::take)
            .ok_or(MoleculeError::InvalidBondId(id))?;
        self.remove_incident_bond(bond.a, id);
        self.remove_incident_bond(bond.b, id);
        self.prune_stereo_for_bond(id);
        self.invalidate_topology();
        Ok(bond)
    }

    pub fn bond(&self, id: BondId) -> Result<&Bond> {
        self.bonds
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidBondId(id))
    }

    pub fn bond_mut(&mut self, id: BondId) -> Result<BondMut<'_>> {
        let original = BondChemistry::from(self.bond(id)?);
        Ok(BondMut {
            molecule: self,
            id,
            original,
        })
    }

    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)> {
        (0..=u32::MAX)
            .zip(self.bonds.iter())
            .filter_map(|(raw, bond)| bond.as_ref().map(|bond| (BondId::new(raw), bond)))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds().map(|(id, _)| id)
    }

    pub fn neighbors(&self, id: AtomId) -> Result<impl Iterator<Item = AtomId> + '_> {
        self.atom(id)?;
        Ok(self.adjacency[id.index()]
            .iter()
            .filter_map(|bond_id| self.bond(*bond_id).ok())
            .map(move |bond| bond.other_atom(id)))
    }

    /// Returns graph components for validation and graph algorithms.
    ///
    /// A completed nonempty public molecule has exactly one component. The
    /// general result shape also supports empty values and private builder,
    /// editor, and format-interpretation staging.
    pub fn connected_components(&self) -> Vec<Vec<AtomId>> {
        let mut seen = vec![false; self.atoms.len()];
        let mut components = Vec::new();
        for start in self.atom_ids() {
            if seen[start.index()] {
                continue;
            }
            seen[start.index()] = true;
            let mut stack = vec![start];
            let mut component = Vec::new();
            while let Some(atom) = stack.pop() {
                component.push(atom);
                let mut neighbors = self
                    .neighbors(atom)
                    .expect("live atom must have valid adjacency")
                    .filter(|neighbor| !seen[neighbor.index()])
                    .collect::<Vec<_>>();
                neighbors.sort_unstable_by(|left, right| right.cmp(left));
                for neighbor in neighbors {
                    if !seen[neighbor.index()] {
                        seen[neighbor.index()] = true;
                        stack.push(neighbor);
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        components
    }

    pub fn incident_bonds(&self, id: AtomId) -> Result<impl Iterator<Item = (BondId, &Bond)> + '_> {
        self.atom(id)?;
        Ok(self.adjacency[id.index()]
            .iter()
            .filter_map(|bond_id| self.bond(*bond_id).ok().map(|bond| (*bond_id, bond))))
    }

    pub fn bond_between(&self, a: AtomId, b: AtomId) -> Result<Option<BondId>> {
        self.atom(a)?;
        self.atom(b)?;
        Ok(self.adjacency[a.index()].iter().copied().find(|bond_id| {
            self.bond(*bond_id)
                .map(|bond| bond.connects(a, b))
                .unwrap_or(false)
        }))
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn perception(&self) -> &PerceptionState {
        &self.perception
    }

    /// Validates and atomically replaces the complete installed perception state.
    ///
    /// The detached state is checked against all stable graph and stereo slots
    /// before mutation. Failure leaves the previous state exactly unchanged.
    pub fn install_perception_state(
        &mut self,
        state: PerceptionState,
    ) -> std::result::Result<(), PerceptionStateInstallError> {
        self.validate_perception_state(&state)?;
        self.perception = state;
        #[cfg(test)]
        self.sync_legacy_perception_fields();
        Ok(())
    }

    pub fn implicit_hydrogens(&self, atom: AtomId) -> Result<Option<u8>> {
        self.atom(atom)?;
        let perceived = self.perception.implicit_hydrogens(atom);
        #[cfg(test)]
        return Ok(perceived.or(self.atom(atom)?.implicit_hydrogens));
        #[cfg(not(test))]
        Ok(perceived)
    }

    pub fn atom_is_aromatic(&self, atom: AtomId) -> Result<Option<bool>> {
        self.atom(atom)?;
        let perceived = self.perception.atom_is_aromatic(atom);
        #[cfg(test)]
        return Ok(perceived.or(self.atom(atom)?.aromatic.then_some(true)));
        #[cfg(not(test))]
        Ok(perceived)
    }

    pub fn bond_is_aromatic(&self, bond: BondId) -> Result<Option<bool>> {
        self.bond(bond)?;
        let perceived = self.perception.bond_is_aromatic(bond);
        #[cfg(test)]
        return Ok(perceived.or(self.bond(bond)?.aromatic.then_some(true)));
        #[cfg(not(test))]
        Ok(perceived)
    }

    pub fn cip_descriptor(&self, element: StereoElementId) -> Result<Option<StereoDescriptor>> {
        self.stereo_element(element)?;
        let perceived = self.perception.cip_descriptor(element);
        #[cfg(test)]
        return Ok(perceived.or(self.stereo_element(element)?.descriptor));
        #[cfg(not(test))]
        Ok(perceived)
    }

    pub fn ring_membership(&self) -> Option<&RingMembership> {
        self.perception.ring_membership()
    }

    pub fn ring_set(&self) -> Option<&RingSet> {
        self.perception.ring_set()
    }

    /// Attaches a conformer after validating its atom slots and identifier capacity.
    pub fn add_conformer(&mut self, mut conformer: Conformer) -> Result<ConformerId> {
        checked_fixed_id_collection_len(0, conformer.positions.len())
            .map_err(|_| MoleculeError::IdentifierCapacityExceeded(MoleculeIdKind::Atom))?;
        let id = checked_molecule_id(
            self.conformers.len(),
            MoleculeIdKind::Conformer,
            ConformerId::new,
        )?;
        for (raw, position) in (0..=u32::MAX).zip(conformer.positions.iter()) {
            if position.is_some()
                && self
                    .atoms
                    .get(AtomId::new(raw).index())
                    .and_then(Option::as_ref)
                    .is_none()
            {
                return Err(MoleculeError::InvalidAtomId(AtomId::new(raw)));
            }
        }
        if conformer.positions.len() < self.atoms.len() {
            conformer.positions.resize(self.atoms.len(), None);
        }
        self.conformers.push(Some(conformer));
        Ok(id)
    }

    pub fn conformer(&self, id: ConformerId) -> Result<&Conformer> {
        self.conformers
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidConformerId(id))
    }

    pub fn conformer_mut(&mut self, id: ConformerId) -> Result<&mut Conformer> {
        self.conformers
            .get_mut(id.index())
            .and_then(Option::as_mut)
            .ok_or(MoleculeError::InvalidConformerId(id))
    }

    pub fn conformers(&self) -> impl Iterator<Item = (ConformerId, &Conformer)> {
        (0..=u32::MAX)
            .zip(self.conformers.iter())
            .filter_map(|(raw, conformer)| {
                conformer
                    .as_ref()
                    .map(|conformer| (ConformerId::new(raw), conformer))
            })
    }

    pub fn first_conformer(&self) -> Option<(ConformerId, &Conformer)> {
        self.conformers().next()
    }

    /// Inserts an ungrouped validated stereo element without narrowing its collection slot.
    ///
    /// Group membership must be established separately through
    /// [`Self::add_stereo_group`]. A pre-grouped element is rejected before
    /// slot allocation or perception invalidation.
    pub fn add_stereo_element(&mut self, element: StereoElement) -> Result<StereoElementId> {
        if element.group.is_some() {
            return Err(MoleculeError::InvalidStereoReference(
                "stereo element group membership must be established through add_stereo_group",
            ));
        }
        self.validate_stereo_element_refs(&element)?;
        let id = checked_molecule_id(
            self.stereo_elements.len(),
            MoleculeIdKind::StereoElement,
            StereoElementId::new,
        )?;
        self.stereo_elements.push(Some(element));
        self.invalidate_stereo();
        Ok(id)
    }

    pub fn stereo_element(&self, id: StereoElementId) -> Result<&StereoElement> {
        self.stereo_elements
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidStereoElementId(id))
    }

    pub fn replace_stereo_element(
        &mut self,
        id: StereoElementId,
        replacement: StereoElement,
    ) -> Result<StereoElement> {
        let Some(current) = self
            .stereo_elements
            .get(id.index())
            .and_then(Option::as_ref)
        else {
            return Err(MoleculeError::InvalidStereoElementId(id));
        };
        if replacement.group != current.group {
            return Err(MoleculeError::InvalidStereoReference(
                "stereo group membership must be changed through stereo-group operations",
            ));
        }
        self.validate_stereo_element_refs(&replacement)?;
        let previous = std::mem::replace(
            self.stereo_elements[id.index()]
                .as_mut()
                .expect("validated stereo element should remain live"),
            replacement,
        );
        self.invalidate_stereo();
        Ok(previous)
    }

    /// Removes a stereo element and returns it detached from any relation group.
    pub fn remove_stereo_element(&mut self, id: StereoElementId) -> Result<StereoElement> {
        let mut element = self
            .stereo_elements
            .get_mut(id.index())
            .and_then(Option::take)
            .ok_or(MoleculeError::InvalidStereoElementId(id))?;
        self.remove_stereo_element_from_groups(id);
        self.invalidate_stereo();
        element.group = None;
        Ok(element)
    }

    pub fn stereo_elements(&self) -> impl Iterator<Item = (StereoElementId, &StereoElement)> {
        (0..=u32::MAX)
            .zip(self.stereo_elements.iter())
            .filter_map(|(raw, element)| {
                element
                    .as_ref()
                    .map(|element| (StereoElementId::new(raw), element))
            })
    }

    pub fn stereo_element_ids(&self) -> impl Iterator<Item = StereoElementId> + '_ {
        self.stereo_elements().map(|(id, _)| id)
    }

    /// Inserts a validated stereo group transactionally.
    pub fn add_stereo_group(&mut self, group: StereoGroup) -> Result<StereoGroupId> {
        if group.members.is_empty() {
            return Err(MoleculeError::InvalidStereoReference(
                "stereo group must contain at least one element",
            ));
        }
        if group.members.iter().copied().collect::<BTreeSet<_>>().len() != group.members.len() {
            return Err(MoleculeError::InvalidStereoReference(
                "stereo group members must be unique",
            ));
        }
        for member in &group.members {
            let element = self.stereo_element(*member)?;
            if element.group.is_some() {
                return Err(MoleculeError::InvalidStereoReference(
                    "stereo element already belongs to a group",
                ));
            }
        }
        let id = checked_molecule_id(
            self.stereo_groups.len(),
            MoleculeIdKind::StereoGroup,
            StereoGroupId::new,
        )?;
        for member in &group.members {
            self.stereo_elements[member.index()]
                .as_mut()
                .expect("validated stereo group member should remain live")
                .group = Some(id);
        }
        self.stereo_groups.push(Some(group));
        self.invalidate_stereo();
        Ok(id)
    }

    pub fn stereo_group(&self, id: StereoGroupId) -> Result<&StereoGroup> {
        self.stereo_groups
            .get(id.index())
            .and_then(Option::as_ref)
            .ok_or(MoleculeError::InvalidStereoGroupId(id))
    }

    pub fn remove_stereo_group(&mut self, id: StereoGroupId) -> Result<StereoGroup> {
        let group = self
            .stereo_groups
            .get_mut(id.index())
            .and_then(Option::take)
            .ok_or(MoleculeError::InvalidStereoGroupId(id))?;
        for member in &group.members {
            if let Some(element) = self
                .stereo_elements
                .get_mut(member.index())
                .and_then(Option::as_mut)
            {
                if element.group == Some(id) {
                    element.group = None;
                }
            }
        }
        self.invalidate_stereo();
        Ok(group)
    }

    pub fn stereo_groups(&self) -> impl Iterator<Item = (StereoGroupId, &StereoGroup)> {
        (0..=u32::MAX)
            .zip(self.stereo_groups.iter())
            .filter_map(|(raw, group)| group.as_ref().map(|group| (StereoGroupId::new(raw), group)))
    }

    /// Returns the complete stereo-group stable-slot count, including tombstones.
    pub fn stereo_group_slot_count(&self) -> usize {
        self.stereo_groups.len()
    }

    /// Iterates every stereo-group stable slot, including interior and trailing tombstones.
    pub fn stereo_group_slots(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoGroupId, Option<&StereoGroup>)> + DoubleEndedIterator + '_
    {
        self.stereo_groups.iter().enumerate().map(|(slot, group)| {
            let raw = u32::try_from(slot)
                .expect("stereo-group slot capacity is checked before insertion");
            (StereoGroupId::new(raw), group.as_ref())
        })
    }

    /// Appends one deleted stereo-group slot without changing live stereo or CIP state.
    pub fn append_stereo_group_tombstone(&mut self) -> Result<StereoGroupId> {
        let id = checked_molecule_id(
            self.stereo_groups.len(),
            MoleculeIdKind::StereoGroup,
            StereoGroupId::new,
        )?;
        self.stereo_groups.push(None);
        Ok(id)
    }

    pub fn set_stereo_bond_mark(&mut self, mark: StereoBondMark) -> Result<()> {
        self.bond(mark.bond)?;
        if let Some(existing) = self
            .stereo_bond_marks
            .iter_mut()
            .find(|existing| existing.bond == mark.bond)
        {
            *existing = mark;
        } else {
            self.stereo_bond_marks.push(mark);
        }
        self.invalidate_stereo();
        Ok(())
    }

    pub fn clear_stereo_bond_mark(&mut self, bond: BondId) -> Result<Option<StereoBondMark>> {
        self.bond(bond)?;
        let Some(index) = self
            .stereo_bond_marks
            .iter()
            .position(|mark| mark.bond == bond)
        else {
            return Ok(None);
        };
        self.invalidate_stereo();
        Ok(Some(self.stereo_bond_marks.remove(index)))
    }

    pub fn stereo_bond_mark(&self, bond: BondId) -> Option<&StereoBondMark> {
        self.stereo_bond_marks.iter().find(|mark| mark.bond == bond)
    }

    pub fn stereo_bond_marks(&self) -> impl Iterator<Item = &StereoBondMark> {
        self.stereo_bond_marks.iter()
    }

    pub fn invalidate_topology(&mut self) {
        self.perception = PerceptionState::default();
    }

    fn remove_incident_bond(&mut self, atom: AtomId, bond: BondId) {
        if let Some(incident) = self.adjacency.get_mut(atom.index()) {
            incident.retain(|id| *id != bond);
        }
    }

    pub(crate) fn invalidate_stereo(&mut self) {
        self.perception.cip_descriptors.clear();
    }

    pub(crate) fn install_valence(
        &mut self,
        model: ValenceModel,
        implicit_hydrogens: BTreeMap<AtomId, u8>,
    ) {
        #[cfg(test)]
        for (raw, atom) in (0..=u32::MAX).zip(self.atoms.iter_mut()) {
            if let Some(atom) = atom {
                atom.implicit_hydrogens = implicit_hydrogens.get(&AtomId::new(raw)).copied();
            }
        }
        self.perception.valence = Some(ValencePerceptionState {
            model: Some(model),
            implicit_hydrogens,
        });
        self.perception.aromaticity = None;
        self.perception.cip_descriptors.clear();
    }

    pub(crate) fn set_implicit_hydrogens(&mut self, atom: AtomId, count: u8) {
        self.perception
            .valence
            .get_or_insert_with(|| ValencePerceptionState {
                model: None,
                implicit_hydrogens: BTreeMap::new(),
            })
            .implicit_hydrogens
            .insert(atom, count);
        #[cfg(test)]
        if let Some(payload) = self.atoms.get_mut(atom.index()).and_then(Option::as_mut) {
            payload.implicit_hydrogens = Some(count);
        }
    }

    pub(crate) fn clear_valence(&mut self) {
        self.perception.valence = None;
        self.perception.aromaticity = None;
        self.perception.cip_descriptors.clear();
        #[cfg(test)]
        for atom in self.atoms.iter_mut().flatten() {
            atom.implicit_hydrogens = None;
        }
    }

    pub(crate) fn install_ring_membership(&mut self, membership: RingMembership) {
        self.perception.rings = Some(RingPerceptionState {
            membership,
            rings: None,
        });
    }

    pub(crate) fn install_rings(&mut self, membership: RingMembership, rings: RingSet) {
        self.perception.rings = Some(RingPerceptionState {
            membership,
            rings: Some(rings),
        });
    }

    pub(crate) fn clear_rings(&mut self) {
        self.perception.rings = None;
        self.perception.aromaticity = None;
        self.perception.cip_descriptors.clear();
    }

    pub(crate) fn discard_ring_results(&mut self) {
        self.perception.rings = None;
    }

    pub(crate) fn begin_aromaticity(&mut self, provenance: AromaticityProvenance) {
        self.perception.aromaticity = Some(AromaticityPerceptionState {
            provenance,
            atoms: BTreeSet::new(),
            bonds: BTreeSet::new(),
        });
        self.perception.cip_descriptors.clear();
        #[cfg(test)]
        {
            for atom in self.atoms.iter_mut().flatten() {
                atom.aromatic = false;
            }
            for bond in self.bonds.iter_mut().flatten() {
                bond.aromatic = false;
            }
        }
    }

    pub(crate) fn clear_aromaticity(&mut self) {
        self.perception.aromaticity = None;
        self.perception.cip_descriptors.clear();
        #[cfg(test)]
        {
            for atom in self.atoms.iter_mut().flatten() {
                atom.aromatic = false;
            }
            for bond in self.bonds.iter_mut().flatten() {
                bond.aromatic = false;
            }
        }
    }

    pub(crate) fn set_atom_aromatic(&mut self, atom: AtomId, aromatic: bool) {
        let Some(state) = self.perception.aromaticity.as_mut() else {
            return;
        };
        if aromatic {
            state.atoms.insert(atom);
        } else {
            state.atoms.remove(&atom);
        }
        #[cfg(test)]
        if let Some(payload) = self.atoms.get_mut(atom.index()).and_then(Option::as_mut) {
            payload.aromatic = aromatic;
        }
    }

    pub(crate) fn set_bond_aromatic(&mut self, bond: BondId, aromatic: bool) {
        let Some(state) = self.perception.aromaticity.as_mut() else {
            return;
        };
        if aromatic {
            state.bonds.insert(bond);
        } else {
            state.bonds.remove(&bond);
        }
        #[cfg(test)]
        if let Some(payload) = self.bonds.get_mut(bond.index()).and_then(Option::as_mut) {
            payload.aromatic = aromatic;
        }
    }

    pub(crate) fn install_cip_descriptor(
        &mut self,
        element: StereoElementId,
        descriptor: StereoDescriptor,
    ) {
        self.perception.cip_descriptors.insert(element, descriptor);
        #[cfg(test)]
        if let Some(payload) = self
            .stereo_elements
            .get_mut(element.index())
            .and_then(Option::as_mut)
        {
            payload.descriptor = Some(descriptor);
        }
    }

    pub(crate) fn clear_cip_descriptors(&mut self) {
        self.perception.cip_descriptors.clear();
        #[cfg(test)]
        for element in self.stereo_elements.iter_mut().flatten() {
            element.descriptor = None;
        }
    }

    fn validate_perception_state(
        &self,
        state: &PerceptionState,
    ) -> std::result::Result<(), PerceptionStateInstallError> {
        if let Some(valence) = &state.valence {
            for atom in valence.implicit_hydrogens.keys().copied() {
                if self
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
            self.validate_ring_state(ring_state)?;
        }

        if let Some(aromaticity) = &state.aromaticity {
            for atom in aromaticity.atoms.iter().copied() {
                if self
                    .atoms
                    .get(atom.index())
                    .and_then(Option::as_ref)
                    .is_none()
                {
                    return Err(PerceptionStateInstallError::InvalidAtomId(atom));
                }
            }
            for bond in aromaticity.bonds.iter().copied() {
                if self
                    .bonds
                    .get(bond.index())
                    .and_then(Option::as_ref)
                    .is_none()
                {
                    return Err(PerceptionStateInstallError::InvalidBondId(bond));
                }
            }
        }

        for element in state.cip_descriptors.keys().copied() {
            if self
                .stereo_elements
                .get(element.index())
                .and_then(Option::as_ref)
                .is_none()
            {
                return Err(PerceptionStateInstallError::InvalidStereoElementId(element));
            }
        }
        Ok(())
    }

    fn validate_ring_state(
        &self,
        state: &RingPerceptionState,
    ) -> std::result::Result<(), PerceptionStateInstallError> {
        let membership = &state.membership;
        let atom_slots = membership.atom_slot_flags().len();
        let bond_slots = membership.bond_slot_flags().len();
        check_install_component_capacity(atom_slots, PerceptionStateComponent::RingAtomSlots)?;
        check_install_component_capacity(bond_slots, PerceptionStateComponent::RingBondSlots)?;
        if atom_slots != self.atoms.len() {
            return Err(PerceptionStateInstallError::RingAtomSlotCountMismatch {
                expected: self.atoms.len(),
                actual: atom_slots,
            });
        }
        if bond_slots != self.bonds.len() {
            return Err(PerceptionStateInstallError::RingBondSlotCountMismatch {
                expected: self.bonds.len(),
                actual: bond_slots,
            });
        }
        for (raw, in_ring) in (0..=u32::MAX).zip(membership.atom_slot_flags()) {
            if *in_ring
                && self
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
                && self
                    .bonds
                    .get(BondId::new(raw).index())
                    .and_then(Option::as_ref)
                    .is_none()
            {
                return Err(PerceptionStateInstallError::InvalidBondId(BondId::new(raw)));
            }
        }

        let Some(ring_set) = &state.rings else {
            return Ok(());
        };
        check_install_component_capacity(ring_set.rings().len(), PerceptionStateComponent::Rings)?;

        let mut covered_atoms = BTreeSet::new();
        let mut covered_bonds = BTreeSet::new();
        for (ring_index, ring) in ring_set.rings().iter().enumerate() {
            self.validate_installed_ring(ring_index, ring)?;
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
        &self,
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
            if self
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
            let Some(bond) = self.bonds.get(bond_id.index()).and_then(Option::as_ref) else {
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

    #[cfg(test)]
    fn sync_legacy_perception_fields(&mut self) {
        for (raw, atom) in (0..=u32::MAX).zip(self.atoms.iter_mut()) {
            if let Some(atom) = atom {
                atom.implicit_hydrogens = self
                    .perception
                    .valence
                    .as_ref()
                    .and_then(|state| state.implicit_hydrogens.get(&AtomId::new(raw)).copied());
                atom.aromatic = self
                    .perception
                    .aromaticity
                    .as_ref()
                    .is_some_and(|state| state.atoms.contains(&AtomId::new(raw)));
            }
        }
        for (raw, bond) in (0..=u32::MAX).zip(self.bonds.iter_mut()) {
            if let Some(bond) = bond {
                bond.aromatic = self
                    .perception
                    .aromaticity
                    .as_ref()
                    .is_some_and(|state| state.bonds.contains(&BondId::new(raw)));
            }
        }
        for (raw, element) in (0..=u32::MAX).zip(self.stereo_elements.iter_mut()) {
            if let Some(element) = element {
                element.descriptor = self
                    .perception
                    .cip_descriptors
                    .get(&StereoElementId::new(raw))
                    .copied();
            }
        }
    }

    pub(crate) fn without_conformers(mut self) -> Self {
        self.conformers.clear();
        self
    }

    /// Clones coordinate-independent molecular state without cloning conformers.
    pub(crate) fn clone_without_conformers(&self) -> Self {
        Self {
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            adjacency: self.adjacency.clone(),
            conformers: Vec::new(),
            stereo_elements: self.stereo_elements.clone(),
            stereo_groups: self.stereo_groups.clone(),
            stereo_bond_marks: self.stereo_bond_marks.clone(),
            props: self.props.clone(),
            perception: self.perception.clone(),
        }
    }

    fn validate_stereo_element_refs(&self, element: &StereoElement) -> Result<()> {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                self.atom(stereo.center)?;
                self.validate_stereo_carriers(&stereo.carriers)?;
            }
            StereoElementKind::DoubleBond(stereo) => {
                let bond = self.bond(stereo.bond)?;
                if !bond.connects(stereo.left, stereo.right) {
                    return Err(MoleculeError::InvalidStereoReference(
                        "double-bond stereo focus does not match bond endpoints",
                    ));
                }
                self.validate_stereo_carriers(&[stereo.left_carrier, stereo.right_carrier])?;
            }
            StereoElementKind::Axis(stereo) => {
                self.bond(stereo.axis)?;
                self.validate_stereo_carriers(&stereo.carriers)?;
            }
        }
        Ok(())
    }

    fn validate_stereo_carriers(&self, carriers: &[StereoCarrier]) -> Result<()> {
        for carrier in carriers {
            if let StereoCarrier::Atom(atom) = carrier {
                self.atom(*atom)?;
            }
        }
        Ok(())
    }

    fn prune_stereo_for_atom(&mut self, atom: AtomId) {
        let removed = self
            .stereo_elements()
            .filter_map(|(id, element)| element.references_atom(atom).then_some(id))
            .collect::<Vec<_>>();
        for id in removed {
            self.stereo_elements[id.index()] = None;
            self.remove_stereo_element_from_groups(id);
        }
        self.invalidate_stereo();
    }

    fn prune_stereo_for_bond(&mut self, bond: BondId) {
        let removed = self
            .stereo_elements()
            .filter_map(|(id, element)| element.references_bond(bond).then_some(id))
            .collect::<Vec<_>>();
        for id in removed {
            self.stereo_elements[id.index()] = None;
            self.remove_stereo_element_from_groups(id);
        }
        self.stereo_bond_marks.retain(|mark| mark.bond != bond);
        self.invalidate_stereo();
    }

    fn remove_stereo_element_from_groups(&mut self, id: StereoElementId) {
        for slot in &mut self.stereo_groups {
            let Some(group) = slot else {
                continue;
            };
            group.members.retain(|member| *member != id);
            if group.members.is_empty() {
                *slot = None;
            }
        }
    }
}

impl Bond {
    fn connects(&self, a: AtomId, b: AtomId) -> bool {
        (self.a == a && self.b == b) || (self.a == b && self.b == a)
    }

    pub(crate) fn other_atom(&self, atom: AtomId) -> AtomId {
        if self.a == atom {
            self.b
        } else {
            self.a
        }
    }
}

pub type Result<T> = std::result::Result<T, MoleculeError>;

/// Fixed-width identifier spaces owned by [`Molecule`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MoleculeIdKind {
    /// Stable atom slots.
    Atom,
    /// Stable bond slots.
    Bond,
    /// Stable conformer slots.
    Conformer,
    /// Stable stereo-element slots.
    StereoElement,
    /// Stable stereo-group slots.
    StereoGroup,
}

impl fmt::Display for MoleculeIdKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Atom => "atom",
            Self::Bond => "bond",
            Self::Conformer => "conformer",
            Self::StereoElement => "stereo element",
            Self::StereoGroup => "stereo group",
        })
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoleculeError {
    InvalidAtomId(AtomId),
    InvalidBondId(BondId),
    InvalidConformerId(ConformerId),
    InvalidStereoElementId(StereoElementId),
    InvalidStereoGroupId(StereoGroupId),
    InvalidStereoReference(&'static str),
    SelfBond(AtomId),
    DuplicateBond {
        a: AtomId,
        b: AtomId,
    },
    /// A new value cannot be represented by the fixed-width ID for `kind`.
    IdentifierCapacityExceeded(MoleculeIdKind),
    UnsupportedFeature(&'static str),
}

impl fmt::Display for MoleculeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(id) => write!(f, "invalid atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid bond id: {id}"),
            Self::InvalidConformerId(id) => write!(f, "invalid conformer id: {id}"),
            Self::InvalidStereoElementId(id) => write!(f, "invalid stereo element id: {id}"),
            Self::InvalidStereoGroupId(id) => write!(f, "invalid stereo group id: {id}"),
            Self::InvalidStereoReference(message) => {
                write!(f, "invalid stereo reference: {message}")
            }
            Self::SelfBond(id) => write!(f, "cannot create a bond from atom {id} to itself"),
            Self::DuplicateBond { a, b } => write!(f, "duplicate bond between {a} and {b}"),
            Self::IdentifierCapacityExceeded(kind) => {
                write!(f, "{kind} identifier capacity exceeded")
            }
            Self::UnsupportedFeature(name) => write!(f, "unsupported feature: {name}"),
        }
    }
}

impl std::error::Error for MoleculeError {}

fn checked_molecule_id<T>(
    length: usize,
    kind: MoleculeIdKind,
    construct: impl FnOnce(u32) -> T,
) -> Result<T> {
    checked_raw_id(length)
        .map(construct)
        .map_err(|_| MoleculeError::IdentifierCapacityExceeded(kind))
}

#[cfg(all(test, target_pointer_width = "64"))]
mod capacity_tests {
    use super::*;

    fn max_slot() -> usize {
        usize::try_from(u64::from(u32::MAX)).expect("64-bit usize")
    }

    fn first_unsupported_slot() -> usize {
        usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize")
    }

    #[test]
    fn every_molecule_identifier_space_checks_its_boundary() {
        assert_eq!(
            checked_molecule_id(max_slot(), MoleculeIdKind::Atom, AtomId::new),
            Ok(AtomId::new(u32::MAX))
        );
        assert_eq!(
            checked_molecule_id(first_unsupported_slot(), MoleculeIdKind::Atom, AtomId::new,),
            Err(MoleculeError::IdentifierCapacityExceeded(
                MoleculeIdKind::Atom
            ))
        );
        assert_eq!(
            checked_molecule_id(first_unsupported_slot(), MoleculeIdKind::Bond, BondId::new,),
            Err(MoleculeError::IdentifierCapacityExceeded(
                MoleculeIdKind::Bond
            ))
        );
        assert_eq!(
            checked_molecule_id(
                first_unsupported_slot(),
                MoleculeIdKind::Conformer,
                ConformerId::new,
            ),
            Err(MoleculeError::IdentifierCapacityExceeded(
                MoleculeIdKind::Conformer
            ))
        );
        assert_eq!(
            checked_molecule_id(
                first_unsupported_slot(),
                MoleculeIdKind::StereoElement,
                StereoElementId::new,
            ),
            Err(MoleculeError::IdentifierCapacityExceeded(
                MoleculeIdKind::StereoElement
            ))
        );
        assert_eq!(
            checked_molecule_id(
                first_unsupported_slot(),
                MoleculeIdKind::StereoGroup,
                StereoGroupId::new,
            ),
            Err(MoleculeError::IdentifierCapacityExceeded(
                MoleculeIdKind::StereoGroup
            ))
        );
    }

    #[test]
    fn capacity_rejection_does_not_mutate_molecule_state() {
        let mut molecule = Molecule::new();
        let before = molecule.clone();
        assert_eq!(
            molecule.add_atom_at_slot(
                Atom::new(Element::from_atomic_number(6).expect("carbon")),
                first_unsupported_slot(),
            ),
            Err(MoleculeError::IdentifierCapacityExceeded(
                MoleculeIdKind::Atom
            ))
        );
        assert_eq!(molecule, before);
    }
}
