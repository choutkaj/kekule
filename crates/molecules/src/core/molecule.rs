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
struct ValencePerception {
    model: Option<ValenceModel>,
    implicit_hydrogens: BTreeMap<AtomId, u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingPerception {
    membership: RingMembership,
    rings: Option<RingSet>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AromaticityPerception {
    provenance: AromaticityProvenance,
    atoms: BTreeSet<AtomId>,
    bonds: BTreeSet<BondId>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PerceptionState {
    valence: Option<ValencePerception>,
    rings: Option<RingPerception>,
    aromaticity: Option<AromaticityPerception>,
    cip_descriptors: BTreeMap<StereoElementId, StereoDescriptor>,
}

impl PerceptionState {
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

    pub fn add_atom(&mut self, atom: Atom) -> AtomId {
        let id = AtomId::new(self.atoms.len() as u32);
        self.atoms.push(Some(atom));
        self.adjacency.push(Vec::new());
        self.invalidate_topology();
        id
    }

    pub fn delete_atom(&mut self, id: AtomId) -> Result<Atom> {
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
        self.atoms
            .iter()
            .enumerate()
            .filter_map(|(index, atom)| atom.as_ref().map(|atom| (AtomId::new(index as u32), atom)))
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atoms().map(|(id, _)| id)
    }

    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<BondId> {
        self.atom(a)?;
        self.atom(b)?;
        if a == b {
            return Err(MoleculeError::SelfBond(a));
        }
        if self.bond_between(a, b)?.is_some() {
            return Err(MoleculeError::DuplicateBond { a, b });
        }
        let id = BondId::new(self.bonds.len() as u32);
        self.bonds.push(Some(Bond::new(a, b, order)));
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        self.invalidate_topology();
        Ok(id)
    }

    pub fn delete_bond(&mut self, id: BondId) -> Result<Bond> {
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
        self.bonds
            .iter()
            .enumerate()
            .filter_map(|(index, bond)| bond.as_ref().map(|bond| (BondId::new(index as u32), bond)))
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

    pub fn add_conformer(&mut self, mut conformer: Conformer) -> Result<ConformerId> {
        for (index, position) in conformer.positions.iter().enumerate() {
            if position.is_some() && self.atoms.get(index).and_then(Option::as_ref).is_none() {
                return Err(MoleculeError::InvalidAtomId(AtomId::new(index as u32)));
            }
        }
        if conformer.positions.len() < self.atoms.len() {
            conformer.positions.resize(self.atoms.len(), None);
        }
        let id = ConformerId::new(self.conformers.len() as u32);
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
        self.conformers
            .iter()
            .enumerate()
            .filter_map(|(index, conformer)| {
                conformer
                    .as_ref()
                    .map(|conformer| (ConformerId::new(index as u32), conformer))
            })
    }

    pub fn first_conformer(&self) -> Option<(ConformerId, &Conformer)> {
        self.conformers().next()
    }

    pub fn add_stereo_element(&mut self, element: StereoElement) -> Result<StereoElementId> {
        self.validate_stereo_element_refs(&element)?;
        let id = StereoElementId::new(self.stereo_elements.len() as u32);
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

    pub fn remove_stereo_element(&mut self, id: StereoElementId) -> Result<StereoElement> {
        let element = self
            .stereo_elements
            .get_mut(id.index())
            .and_then(Option::take)
            .ok_or(MoleculeError::InvalidStereoElementId(id))?;
        self.remove_stereo_element_from_groups(id);
        self.invalidate_stereo();
        Ok(element)
    }

    pub fn stereo_elements(&self) -> impl Iterator<Item = (StereoElementId, &StereoElement)> {
        self.stereo_elements
            .iter()
            .enumerate()
            .filter_map(|(index, element)| {
                element
                    .as_ref()
                    .map(|element| (StereoElementId::new(index as u32), element))
            })
    }

    pub fn stereo_element_ids(&self) -> impl Iterator<Item = StereoElementId> + '_ {
        self.stereo_elements().map(|(id, _)| id)
    }

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
        let id = StereoGroupId::new(self.stereo_groups.len() as u32);
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
        self.stereo_groups
            .iter()
            .enumerate()
            .filter_map(|(index, group)| {
                group
                    .as_ref()
                    .map(|group| (StereoGroupId::new(index as u32), group))
            })
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
        for (index, atom) in self.atoms.iter_mut().enumerate() {
            if let Some(atom) = atom {
                atom.implicit_hydrogens =
                    implicit_hydrogens.get(&AtomId::new(index as u32)).copied();
            }
        }
        self.perception.valence = Some(ValencePerception {
            model: Some(model),
            implicit_hydrogens,
        });
        self.perception.aromaticity = None;
        self.perception.cip_descriptors.clear();
    }

    pub(crate) fn set_implicit_hydrogens(&mut self, atom: AtomId, count: u8) {
        self.perception
            .valence
            .get_or_insert_with(|| ValencePerception {
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
        self.perception.rings = Some(RingPerception {
            membership,
            rings: None,
        });
    }

    pub(crate) fn install_rings(&mut self, membership: RingMembership, rings: RingSet) {
        self.perception.rings = Some(RingPerception {
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
        self.perception.aromaticity = Some(AromaticityPerception {
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

    pub(crate) fn without_conformers(mut self) -> Self {
        self.conformers.clear();
        self
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
        for group in self.stereo_groups.iter_mut().flatten() {
            group.members.retain(|member| *member != id);
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
    DuplicateBond { a: AtomId, b: AtomId },
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
            Self::UnsupportedFeature(name) => write!(f, "unsupported feature: {name}"),
        }
    }
}

impl std::error::Error for MoleculeError {}
