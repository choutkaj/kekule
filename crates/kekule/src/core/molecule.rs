use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::{Deref, DerefMut};

use super::*;

/// One completed connected canonical molecular graph.
///
/// Empty and single-atom values are valid boundary cases. Build nontrivial
/// graphs with [`Molecule::builder`] and change topology transactionally with
/// [`Molecule::edit`], both of which enforce canonical representation and
/// connectedness before publication. Source-format stereo marks are resolved
/// into canonical stereo elements during interpretation and are never stored
/// in this payload.
///
/// ```compile_fail
/// use kekule::core::Molecule;
///
/// let molecule = Molecule::default();
/// let _ = molecule.stereo_bond_marks();
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Molecule {
    pub(crate) atoms: Vec<Option<Atom>>,
    pub(crate) bonds: Vec<Option<Bond>>,
    pub(crate) adjacency: Vec<Vec<BondId>>,
    pub(crate) conformers: Vec<Option<Conformer>>,
    pub(crate) stereo_elements: Vec<Option<StereoElement>>,
    pub(crate) stereo_groups: Vec<Option<StereoGroup>>,
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
    hydrogens: HydrogenDeclaration,
}

impl From<&Atom> for AtomChemistry {
    fn from(atom: &Atom) -> Self {
        Self {
            element: atom.element,
            isotope: atom.isotope,
            formal_charge: atom.formal_charge,
            radical: atom.radical,
            hydrogens: atom.hydrogens,
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
    /// This aggregate does not require perception.
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

    pub(crate) fn atom_mut(&mut self, id: AtomId) -> Result<AtomMut<'_>> {
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

    /// Returns mutable generic annotations without changing represented chemistry.
    pub fn atom_props_mut(&mut self, id: AtomId) -> Result<&mut PropMap> {
        Ok(&mut self
            .atoms
            .get_mut(id.index())
            .and_then(Option::as_mut)
            .ok_or(MoleculeError::InvalidAtomId(id))?
            .props)
    }

    /// Inserts a bond into crate-private construction or editing staging.
    pub(crate) fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
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

    pub(crate) fn bond_mut(&mut self, id: BondId) -> Result<BondMut<'_>> {
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

    /// Returns mutable generic annotations without changing represented chemistry.
    pub fn bond_props_mut(&mut self, id: BondId) -> Result<&mut PropMap> {
        Ok(&mut self
            .bonds
            .get_mut(id.index())
            .and_then(Option::as_mut)
            .ok_or(MoleculeError::InvalidBondId(id))?
            .props)
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
        validate_perception_state(self, &state)?;
        self.perception = state;
        Ok(())
    }

    pub fn implicit_hydrogens(&self, atom: AtomId) -> Result<Option<u8>> {
        self.atom(atom)?;
        Ok(self.perception.implicit_hydrogens(atom))
    }

    pub fn atom_is_aromatic(&self, atom: AtomId) -> Result<Option<bool>> {
        self.atom(atom)?;
        Ok(self.perception.atom_is_aromatic(atom))
    }

    pub fn bond_is_aromatic(&self, bond: BondId) -> Result<Option<bool>> {
        self.bond(bond)?;
        Ok(self.perception.bond_is_aromatic(bond))
    }

    pub fn cip_descriptor(&self, element: StereoElementId) -> Result<Option<StereoDescriptor>> {
        self.stereo_element(element)?;
        Ok(self.perception.cip_descriptor(element))
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
    /// Carrier order, double-bond endpoints, and reference carriers are
    /// canonicalized from represented graph state before storage.
    ///
    /// Group membership must be established separately through
    /// [`Self::add_stereo_group`]. A pre-grouped element is rejected before
    /// slot allocation or perception invalidation.
    pub fn add_stereo_element(&mut self, mut element: StereoElement) -> Result<StereoElementId> {
        if element.group.is_some() {
            return Err(MoleculeError::InvalidStereoReference(
                "stereo element group membership must be established through add_stereo_group",
            ));
        }
        self.canonicalize_stereo_element(&mut element);
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

    /// Replaces an element through the same canonical storage boundary as
    /// [`Self::add_stereo_element`].
    pub fn replace_stereo_element(
        &mut self,
        id: StereoElementId,
        mut replacement: StereoElement,
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
        self.canonicalize_stereo_element(&mut replacement);
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

    pub fn invalidate_topology(&mut self) {
        self.perception = PerceptionState::default();
    }

    fn remove_incident_bond(&mut self, atom: AtomId, bond: BondId) {
        if let Some(incident) = self.adjacency.get_mut(atom.index()) {
            incident.retain(|id| *id != bond);
        }
    }

    pub(crate) fn invalidate_stereo(&mut self) {
        self.perception.stereo = None;
    }

    pub(crate) fn install_valence(
        &mut self,
        model: ValenceModel,
        implicit_hydrogens: BTreeMap<AtomId, u8>,
    ) {
        self.perception.valence = Some(ValencePerceptionState {
            model: Some(model),
            implicit_hydrogens,
        });
        self.perception.aromaticity = None;
        self.perception.stereo = None;
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
    }

    pub(crate) fn install_ring_membership(&mut self, membership: RingMembership) {
        self.perception.rings = Some(RingPerceptionState {
            membership,
            basis: None,
        });
    }

    pub(crate) fn install_ring_basis(
        &mut self,
        membership: RingMembership,
        model: RingBasisModel,
        rings: RingSet,
    ) {
        self.perception.rings = Some(RingPerceptionState {
            membership,
            basis: Some(RingBasisState::new(Some(model), rings)),
        });
    }

    pub(crate) fn begin_aromaticity(&mut self, model: AromaticityModel) {
        self.perception.aromaticity = Some(AromaticityPerceptionState {
            model,
            atoms: BTreeSet::new(),
            bonds: BTreeSet::new(),
        });
        self.perception.stereo = None;
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
    }

    pub(crate) fn install_cip_descriptor(
        &mut self,
        element: StereoElementId,
        descriptor: StereoDescriptor,
    ) {
        self.perception
            .stereo
            .get_or_insert_default()
            .cip_descriptors
            .insert(element, descriptor);
    }

    pub(crate) fn replace_stereo_perception(
        &mut self,
        state: Option<StereoPerceptionState>,
    ) -> Option<StereoPerceptionState> {
        std::mem::replace(&mut self.perception.stereo, state)
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

    fn canonicalize_stereo_element(&self, element: &mut StereoElement) {
        match &mut element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                if sort_stereo_carriers(&mut stereo.carriers) {
                    stereo.orientation = stereo.orientation.map(flip_tetrahedral_orientation);
                }
            }
            StereoElementKind::DoubleBond(stereo) => self.canonicalize_double_bond_stereo(stereo),
            StereoElementKind::Axis(stereo) => self.canonicalize_axis_stereo(stereo),
        }
    }

    pub(crate) fn canonicalize_stored_stereo_elements(&mut self) {
        let stored = self
            .stereo_elements
            .iter()
            .enumerate()
            .filter_map(|(index, element)| element.clone().map(|element| (index, element)))
            .collect::<Vec<_>>();
        let mut replacements = Vec::new();
        for (index, mut element) in stored {
            self.canonicalize_stereo_element(&mut element);
            if self.stereo_elements[index].as_ref() != Some(&element) {
                replacements.push((index, element));
            }
        }
        if replacements.is_empty() {
            return;
        }
        for (index, element) in replacements {
            self.stereo_elements[index] = Some(element);
        }
        self.invalidate_stereo();
    }

    fn canonicalize_double_bond_stereo(&self, stereo: &mut DoubleBondStereo) {
        let Ok(bond) = self.bond(stereo.bond) else {
            return;
        };
        if !bond.connects(stereo.left, stereo.right) {
            return;
        }

        let (left, right) = sorted_atom_pair(stereo.left, stereo.right);
        // Exchanging both endpoint labels and their references preserves the
        // together/opposite relation; only changing one reference flips it.
        let (mut left_carrier, mut right_carrier) = if stereo.left == left {
            (stereo.left_carrier, stereo.right_carrier)
        } else {
            (stereo.right_carrier, stereo.left_carrier)
        };
        let mut reference_changes = 0;
        if let Some(canonical) =
            self.canonical_endpoint_reference(left, right, stereo.bond, left_carrier)
        {
            reference_changes += usize::from(canonical != left_carrier);
            left_carrier = canonical;
        }
        if let Some(canonical) =
            self.canonical_endpoint_reference(right, left, stereo.bond, right_carrier)
        {
            reference_changes += usize::from(canonical != right_carrier);
            right_carrier = canonical;
        }

        stereo.left = left;
        stereo.right = right;
        stereo.left_carrier = left_carrier;
        stereo.right_carrier = right_carrier;
        if reference_changes % 2 == 1 {
            stereo.orientation = stereo.orientation.map(flip_double_bond_orientation);
        }
    }

    fn canonical_endpoint_reference(
        &self,
        endpoint: AtomId,
        other_endpoint: AtomId,
        focus: BondId,
        current: StereoCarrier,
    ) -> Option<StereoCarrier> {
        let atom_references = self.atom_stereo_references(endpoint, other_endpoint, focus);
        match current {
            StereoCarrier::Atom(atom) if atom_references.contains(&StereoCarrier::Atom(atom)) => {
                atom_references.first().copied()
            }
            StereoCarrier::ImplicitHydrogen => atom_references
                .first()
                .copied()
                .or(Some(StereoCarrier::ImplicitHydrogen)),
            StereoCarrier::Atom(_) | StereoCarrier::ImplicitLonePair => None,
        }
    }

    fn canonicalize_axis_stereo(&self, stereo: &mut AxisStereo) {
        if stereo.carriers.len() != 2 {
            return;
        }
        let Ok(axis) = self.bond(stereo.axis) else {
            return;
        };
        // Reversing the axis exchanges both endpoint reference vectors as
        // well, leaving the stored handedness unchanged.
        let (left, right) = sorted_atom_pair(axis.a(), axis.b());
        let Some(mut left_carrier) =
            self.axis_reference_on_endpoint(left, right, stereo.axis, &stereo.carriers)
        else {
            return;
        };
        let Some(mut right_carrier) =
            self.axis_reference_on_endpoint(right, left, stereo.axis, &stereo.carriers)
        else {
            return;
        };

        let mut reference_changes = 0;
        if let Some(canonical) =
            self.canonical_endpoint_reference(left, right, stereo.axis, left_carrier)
        {
            reference_changes += usize::from(canonical != left_carrier);
            left_carrier = canonical;
        }
        if let Some(canonical) =
            self.canonical_endpoint_reference(right, left, stereo.axis, right_carrier)
        {
            reference_changes += usize::from(canonical != right_carrier);
            right_carrier = canonical;
        }

        stereo.carriers = vec![left_carrier, right_carrier];
        if reference_changes % 2 == 1 {
            stereo.orientation = stereo.orientation.map(flip_axis_orientation);
        }
    }

    fn axis_reference_on_endpoint(
        &self,
        endpoint: AtomId,
        other_endpoint: AtomId,
        axis: BondId,
        carriers: &[StereoCarrier],
    ) -> Option<StereoCarrier> {
        let references = self.atom_stereo_references(endpoint, other_endpoint, axis);
        let matches = carriers
            .iter()
            .copied()
            .filter(|carrier| references.contains(carrier))
            .collect::<Vec<_>>();
        if matches.len() == 1 {
            Some(matches[0])
        } else {
            None
        }
    }

    fn atom_stereo_references(
        &self,
        endpoint: AtomId,
        other_endpoint: AtomId,
        focus: BondId,
    ) -> Vec<StereoCarrier> {
        let mut references = self
            .incident_bonds(endpoint)
            .ok()
            .into_iter()
            .flatten()
            .filter_map(|(bond_id, bond)| {
                if bond_id == focus {
                    return None;
                }
                let atom = bond.other_atom(endpoint);
                (atom != other_endpoint).then_some(StereoCarrier::Atom(atom))
            })
            .collect::<Vec<_>>();
        references.sort_by_key(|carrier| match carrier {
            StereoCarrier::Atom(atom) => {
                let is_hydrogen = self
                    .atom(*atom)
                    .is_ok_and(|atom| atom.element.atomic_number() == 1);
                (u8::from(is_hydrogen), atom.raw())
            }
            StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => {
                unreachable!("represented atom references contain only atom carriers")
            }
        });
        references.dedup();
        references
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

fn sort_stereo_carriers(carriers: &mut Vec<StereoCarrier>) -> bool {
    let mut indexed = carriers.iter().copied().enumerate().collect::<Vec<_>>();
    indexed.sort_by_key(|(_, carrier)| stereo_carrier_key(carrier));

    // The orientation changes sign exactly for an odd permutation from the
    // caller's carrier order to the canonical order.
    let mut odd_permutation = false;
    for left in 0..indexed.len() {
        for right in (left + 1)..indexed.len() {
            if indexed[left].0 > indexed[right].0 {
                odd_permutation = !odd_permutation;
            }
        }
    }
    *carriers = indexed.into_iter().map(|(_, carrier)| carrier).collect();
    odd_permutation
}

fn stereo_carrier_key(carrier: &StereoCarrier) -> (u8, u32) {
    match carrier {
        StereoCarrier::Atom(atom) => (0, atom.raw()),
        StereoCarrier::ImplicitHydrogen => (1, 0),
        StereoCarrier::ImplicitLonePair => (2, 0),
    }
}

fn sorted_atom_pair(left: AtomId, right: AtomId) -> (AtomId, AtomId) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}

fn flip_tetrahedral_orientation(orientation: TetrahedralOrientation) -> TetrahedralOrientation {
    match orientation {
        TetrahedralOrientation::Clockwise => TetrahedralOrientation::CounterClockwise,
        TetrahedralOrientation::CounterClockwise => TetrahedralOrientation::Clockwise,
    }
}

fn flip_double_bond_orientation(orientation: DoubleBondOrientation) -> DoubleBondOrientation {
    match orientation {
        DoubleBondOrientation::Together => DoubleBondOrientation::Opposite,
        DoubleBondOrientation::Opposite => DoubleBondOrientation::Together,
    }
}

fn flip_axis_orientation(orientation: AxisOrientation) -> AxisOrientation {
    match orientation {
        AxisOrientation::Clockwise => AxisOrientation::CounterClockwise,
        AxisOrientation::CounterClockwise => AxisOrientation::Clockwise,
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
