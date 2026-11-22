//! Represented stereo: the one validator that explicit stereo edits,
//! publication, and [`MoleculeEditor::validate`] share, canonical carrier
//! storage, and the pruning that graph edits apply.
use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostics::Countable;
use crate::fixed_u32_id::checked_id;
use crate::stereo::StereoViolation;

use super::super::*;
use super::MoleculeEditor;

impl MoleculeEditor {
    /// Adds a represented stereo element after validating it against the draft
    /// graph, and returns its draft ID. The element is stored in canonical
    /// carrier order, and installed CIP descriptors are cleared.
    ///
    /// # Errors
    ///
    /// - [`MoleculeError::InvalidAtomId`] or [`MoleculeError::InvalidBondId`]
    ///   when the element names a missing focus or carrier.
    /// - [`MoleculeError::InvalidStereoReference`] when the element already
    ///   names a group (use [`Self::add_stereo_group`]), its carriers are
    ///   duplicated, non-adjacent or incomplete, its focus bond has the wrong
    ///   order, or another element already describes the same focus.
    /// - [`MoleculeError::CapacityExceeded`] when no stereo element
    ///   ID is left.
    pub fn add_stereo_element(
        &mut self,
        mut element: StereoElement,
    ) -> Result<StereoElementId, MoleculeError> {
        let pending = self.graph.next_stereo_element_id();
        if element.group.is_some() {
            return Err(MoleculeError::InvalidStereoReference(
                StereoViolation::GroupMembershipOutsideGroupEdit { element: pending },
            ));
        }
        self.graph.validate_stereo_element(pending, &element)?;
        self.graph
            .validate_unique_stereo_focus(pending, &element, None)?;
        self.graph.canonicalize_stereo_element(&mut element);
        checked_id::<StereoElementId>(
            self.graph.stereo_elements.slot_count(),
            Countable::StereoElement,
        )
        .map_err(MoleculeError::CapacityExceeded)?;
        let id = self.graph.push_stereo_element(element);
        self.properties.owner_mut().clear();
        self.invalidate_cip();
        Ok(id)
    }

    /// Replaces a stereo element under the same validation as
    /// [`Self::add_stereo_element`] and returns the previous element. Group
    /// membership cannot change here: the replacement must name the same
    /// group as the current element.
    ///
    /// # Errors
    ///
    /// [`MoleculeError::InvalidStereoElementId`] when `id` names no element;
    /// otherwise as [`Self::add_stereo_element`], including
    /// [`MoleculeError::InvalidStereoReference`] for a changed group.
    pub fn replace_stereo_element(
        &mut self,
        id: StereoElementId,
        mut replacement: StereoElement,
    ) -> Result<StereoElement, MoleculeError> {
        let current = self.graph.require_stereo_element(id)?;
        if replacement.group != current.group {
            return Err(MoleculeError::InvalidStereoReference(
                StereoViolation::GroupMembershipOutsideGroupEdit { element: id },
            ));
        }
        self.graph.validate_stereo_element(id, &replacement)?;
        self.graph
            .validate_unique_stereo_focus(id, &replacement, Some(id))?;
        self.graph.canonicalize_stereo_element(&mut replacement);
        let previous = std::mem::replace(
            self.graph
                .stereo_element_mut(id)
                .expect("validated stereo element should remain live"),
            replacement,
        );
        self.properties.owner_mut().clear();
        self.invalidate_cip();
        Ok(previous)
    }

    /// Stores `element` in place of element `id` without validation or
    /// canonicalization, for hydrogen transforms that rename a carrier while
    /// its graph atom is about to change; publication validates the result.
    pub(crate) fn set_stereo_element_unchecked(
        &mut self,
        id: StereoElementId,
        element: StereoElement,
    ) {
        *self
            .graph
            .stereo_element_mut(id)
            .expect("live stereo element") = element;
        self.invalidate_cip();
    }

    /// Removes a stereo element, drops it from its group, and returns it with
    /// no group. The removal is not listed in
    /// [`Self::pruned_stereo_elements`].
    ///
    /// # Errors
    ///
    /// [`MoleculeError::InvalidStereoElementId`] when `id` names no element.
    pub fn remove_stereo_element(
        &mut self,
        id: StereoElementId,
    ) -> Result<StereoElement, MoleculeError> {
        let mut element = self
            .graph
            .take_stereo_element(id)
            .ok_or(MoleculeError::InvalidStereoElementId(id))?;
        self.properties.owner_mut().clear();
        self.invalidate_cip();
        element.group = None;
        Ok(element)
    }

    /// Adds an enhanced stereo group, assigns its members to it, and returns
    /// its draft ID.
    ///
    /// # Errors
    ///
    /// - [`MoleculeError::InvalidStereoReference`] when the group is empty,
    ///   lists a member twice, or a member already belongs to a group.
    /// - [`MoleculeError::InvalidStereoElementId`] when a member names no
    ///   element.
    /// - [`MoleculeError::CapacityExceeded`] when no stereo group ID
    ///   is left.
    pub fn add_stereo_group(&mut self, group: StereoGroup) -> Result<StereoGroupId, MoleculeError> {
        let pending = self.graph.next_stereo_group_id();
        validate_group_members(pending, &group)?;
        for &member in &group.members {
            if let Some(group) = self.graph.require_stereo_element(member)?.group {
                return Err(MoleculeError::InvalidStereoReference(
                    StereoViolation::StereoElementAlreadyGrouped {
                        element: member,
                        group,
                    },
                ));
            }
        }
        checked_id::<StereoGroupId>(
            self.graph.stereo_groups.slot_count(),
            Countable::StereoGroup,
        )
        .map_err(MoleculeError::CapacityExceeded)?;
        let members = group.members.clone();
        let id = self.graph.push_stereo_group(group);
        for member in members {
            self.graph
                .stereo_element_mut(member)
                .expect("validated stereo group member should remain live")
                .group = Some(id);
        }
        self.properties.owner_mut().clear();
        self.invalidate_cip();
        Ok(id)
    }

    /// Removes an enhanced stereo group and returns it. Its member elements
    /// remain, ungrouped.
    ///
    /// # Errors
    ///
    /// [`MoleculeError::InvalidStereoGroupId`] when `id` names no group.
    pub fn remove_stereo_group(&mut self, id: StereoGroupId) -> Result<StereoGroup, MoleculeError> {
        let group = self
            .graph
            .take_stereo_group(id)
            .ok_or(MoleculeError::InvalidStereoGroupId(id))?;
        for &member in &group.members {
            if let Some(element) = self.graph.stereo_element_mut(member) {
                if element.group == Some(id) {
                    element.group = None;
                }
            }
        }
        self.properties.owner_mut().clear();
        self.invalidate_cip();
        Ok(group)
    }

    /// Replaces a relation group without changing its ID. All membership checks
    /// precede mutation, and members already in another group are rejected.
    pub fn replace_stereo_group(
        &mut self,
        id: StereoGroupId,
        replacement: StereoGroup,
    ) -> Result<StereoGroup, MoleculeError> {
        let previous = self.graph.require_stereo_group(id)?.clone();
        validate_group_members(id, &replacement)?;
        for &member in &replacement.members {
            if let Some(group) = self
                .graph
                .require_stereo_element(member)?
                .group
                .filter(|&group| group != id)
            {
                return Err(MoleculeError::InvalidStereoReference(
                    StereoViolation::StereoElementAlreadyGrouped {
                        element: member,
                        group,
                    },
                ));
            }
        }
        if previous == replacement {
            return Ok(previous);
        }
        for member in &previous.members {
            self.graph
                .stereo_element_mut(*member)
                .expect("validated group member")
                .group = None;
        }
        for member in &replacement.members {
            self.graph
                .stereo_element_mut(*member)
                .expect("validated replacement member")
                .group = Some(id);
        }
        *self
            .graph
            .stereo_group_mut(id)
            .expect("validated stereo group") = replacement;
        self.properties.owner_mut().clear();
        self.invalidate_cip();
        Ok(previous)
    }

    /// The four carriers of a tetrahedral assertion read from a source, or the
    /// error inserting an element with this many carriers reports.
    pub(crate) fn tetrahedral_carriers(
        &self,
        center: AtomId,
        carriers: Vec<StereoCarrier>,
    ) -> Result<[StereoCarrier; 4], MoleculeError> {
        let carrier_count = carriers.len();
        carriers.try_into().map_err(|_| {
            MoleculeError::InvalidStereoReference(StereoViolation::InvalidTetrahedralCarrierCount {
                element: self.graph.next_stereo_element_id(),
                center,
                carrier_count,
            })
        })
    }

    /// Restores canonical carrier order on every stored element, as
    /// publication does after graph edits changed the references, and returns
    /// whether any element changed.
    pub(crate) fn canonicalize_stored_stereo_elements(&mut self) -> bool {
        let replacements = self
            .graph
            .stereo_elements()
            .filter_map(|(id, element)| {
                let mut canonical = element.clone();
                self.graph.canonicalize_stereo_element(&mut canonical);
                (&canonical != element).then_some((id, canonical))
            })
            .collect::<Vec<_>>();
        if replacements.is_empty() {
            return false;
        }
        for (id, element) in replacements {
            *self.graph.stereo_element_mut(id).expect("live element") = element;
        }
        self.invalidate_cip();
        true
    }

    /// Removes the stereo elements that a bond between `atoms` (added, or
    /// deleted as `bond`) leaves invalid, recording them as pruned.
    ///
    /// Validity depends only on the neighbors of focus atoms and on
    /// carrier-focus adjacency, so the edit can invalidate only elements whose
    /// focus atoms, axis endpoints, focus bond, or carriers include its
    /// endpoints. Every other element is untouched, valid or not.
    pub(super) fn prune_invalid_stereo_at(&mut self, atoms: [AtomId; 2], bond: BondId) {
        let graph = &self.graph;
        let involves = |element: &StereoElement| {
            element.references_bond(bond)
                || atoms.iter().any(|&atom| element.references_atom(atom))
                || matches!(&element.kind, StereoElementKind::Axis(axis)
                if graph.bond(axis.axis).is_some_and(|axis| {
                    atoms.iter().any(|&atom| axis.a() == atom || axis.b() == atom)
                }))
        };
        let removed = graph
            .stereo_elements()
            .filter_map(|(id, element)| {
                (involves(element) && graph.validate_stereo_element_refs(id, element).is_err())
                    .then_some(id)
            })
            .collect::<Vec<_>>();
        self.remove_pruned_stereo(removed);
    }

    /// Removes and records the stereo elements that reference `atom`.
    pub(super) fn prune_stereo_for_atom(&mut self, atom: AtomId) {
        let removed = self
            .graph
            .stereo_elements()
            .filter_map(|(id, element)| element.references_atom(atom).then_some(id))
            .collect();
        self.remove_pruned_stereo(removed);
    }

    /// Removes and records the stereo elements focused on `bond`.
    pub(super) fn prune_stereo_for_bond(&mut self, bond: BondId) {
        let removed = self
            .graph
            .stereo_elements()
            .filter_map(|(id, element)| element.references_bond(bond).then_some(id))
            .collect();
        self.remove_pruned_stereo(removed);
    }

    fn remove_pruned_stereo(&mut self, removed: Vec<StereoElementId>) {
        for &id in &removed {
            if let Some(element) = self.graph.take_stereo_element(id) {
                self.pruned_values.push(element);
            }
        }
        self.invalidate_cip();
        self.pruned_stereo.extend(removed);
    }
}

impl<L: SlotLiveness> Graph<L> {
    /// Checks `element`, which has or receives ID `id`, completely: its
    /// graph references and the hydrogen declarations its implicit-hydrogen
    /// carriers rely on. Explicit stereo edits and publication use this;
    /// graph-edit pruning checks references alone, because hydrogen
    /// transforms change a carrier and its declaration in separate steps.
    pub(crate) fn validate_stereo_element(
        &self,
        id: StereoElementId,
        element: &StereoElement,
    ) -> Result<(), MoleculeError> {
        self.validate_stereo_element_refs(id, element)?;
        let invalid = |issue| Err(MoleculeError::InvalidStereoReference(issue));
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                let center = stereo.center;
                let carried = stereo
                    .carriers
                    .iter()
                    .filter(|carrier| **carrier == StereoCarrier::ImplicitHydrogen)
                    .count();
                if let ImplicitHydrogens::Fixed(fixed_hydrogens) =
                    self.require_atom(center)?.hydrogens()
                {
                    if usize::from(fixed_hydrogens) != carried {
                        return invalid(StereoViolation::ImplicitHydrogenCarrierMismatch {
                            element: id,
                            atom: center,
                            fixed_hydrogens,
                        });
                    }
                }
            }
            StereoElementKind::DoubleBond(stereo) => {
                for (endpoint, carrier) in [
                    (stereo.left, stereo.left_carrier),
                    (stereo.right, stereo.right_carrier),
                ] {
                    if carrier == StereoCarrier::ImplicitHydrogen
                        && self.require_atom(endpoint)?.hydrogens() == ImplicitHydrogens::Fixed(0)
                    {
                        return invalid(StereoViolation::ImplicitHydrogenCarrierMismatch {
                            element: id,
                            atom: endpoint,
                            fixed_hydrogens: 0,
                        });
                    }
                }
            }
            StereoElementKind::Axis(_) => {}
        }
        Ok(())
    }

    /// Checks the references of `element`, which has or receives ID `id`,
    /// against the represented graph.
    pub(crate) fn validate_stereo_element_refs(
        &self,
        id: StereoElementId,
        element: &StereoElement,
    ) -> Result<(), MoleculeError> {
        let invalid = |issue| Err(MoleculeError::InvalidStereoReference(issue));
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                let center = stereo.center;
                self.require_atom(center)?;
                for (index, &carrier) in stereo.carriers.iter().enumerate() {
                    if stereo.carriers[..index].contains(&carrier) {
                        return invalid(StereoViolation::DuplicateTetrahedralCarrier {
                            element: id,
                            center,
                            carrier,
                        });
                    }
                }
                self.validate_stereo_carriers(&stereo.carriers)?;
                for &carrier in &stereo.carriers {
                    if !self.carrier_is_adjacent(center, None, carrier)? {
                        return invalid(StereoViolation::TetrahedralCarrierNotAdjacent {
                            element: id,
                            center,
                            carrier,
                        });
                    }
                }
                if let Some(neighbor) = self
                    .require_neighbors(center)?
                    .find(|neighbor| !stereo.carriers.contains(&StereoCarrier::Atom(*neighbor)))
                {
                    return invalid(StereoViolation::UnrepresentedTetrahedralNeighbor {
                        element: id,
                        center,
                        neighbor,
                    });
                }
            }
            StereoElementKind::DoubleBond(stereo) => {
                let bond = self.require_bond(stereo.bond)?;
                if bond.order != BondOrder::Double {
                    return invalid(StereoViolation::InvalidDoubleBondOrder {
                        element: id,
                        bond: stereo.bond,
                        order: bond.order,
                    });
                }
                if !bond.connects(stereo.left, stereo.right) {
                    return invalid(StereoViolation::DoubleBondFocusMismatch {
                        element: id,
                        bond: stereo.bond,
                        left: stereo.left,
                        right: stereo.right,
                    });
                }
                for (endpoint, carrier) in [
                    (stereo.left, stereo.left_carrier),
                    (stereo.right, stereo.right_carrier),
                ] {
                    if carrier == StereoCarrier::ImplicitLonePair {
                        return invalid(StereoViolation::UnsupportedDoubleBondCarrier {
                            element: id,
                            endpoint,
                            carrier,
                        });
                    }
                }
                self.validate_stereo_carriers(&[stereo.left_carrier, stereo.right_carrier])?;
                for (endpoint, other, carrier) in [
                    (stereo.left, stereo.right, stereo.left_carrier),
                    (stereo.right, stereo.left, stereo.right_carrier),
                ] {
                    if carrier == StereoCarrier::Atom(other) {
                        return invalid(StereoViolation::DoubleBondCarrierIsFocusAtom {
                            element: id,
                            endpoint,
                            carrier: other,
                        });
                    }
                    if !self.carrier_is_adjacent(endpoint, Some(other), carrier)? {
                        return invalid(StereoViolation::DoubleBondCarrierNotAdjacent {
                            element: id,
                            endpoint,
                            carrier,
                        });
                    }
                }
                for (endpoint, other) in [(stereo.left, stereo.right), (stereo.right, stereo.left)]
                {
                    let substituent_count = self
                        .require_neighbors(endpoint)?
                        .filter(|neighbor| *neighbor != other)
                        .count();
                    if substituent_count > 2 {
                        return invalid(StereoViolation::DoubleBondEndpointOvercoordinated {
                            element: id,
                            endpoint,
                            substituent_count,
                        });
                    }
                }
            }
            StereoElementKind::Axis(stereo) => {
                let axis = self.require_bond(stereo.axis)?;
                self.validate_stereo_carriers(&stereo.carriers)?;
                let mut endpoints = [false; 2];
                for &carrier in &stereo.carriers {
                    let StereoCarrier::Atom(atom) = carrier else {
                        return invalid(StereoViolation::UnsupportedAxisCarrier {
                            element: id,
                            axis: stereo.axis,
                            carrier,
                        });
                    };
                    let adjacent = [
                        atom != axis.b() && self.require_bond_between(axis.a(), atom)?.is_some(),
                        atom != axis.a() && self.require_bond_between(axis.b(), atom)?.is_some(),
                    ];
                    if adjacent[0] == adjacent[1] {
                        return invalid(StereoViolation::AxisCarrierNotAdjacent {
                            element: id,
                            axis: stereo.axis,
                            carrier,
                        });
                    }
                    let endpoint = usize::from(adjacent[1]);
                    if endpoints[endpoint] {
                        return invalid(StereoViolation::AxisCarrierEndpointMismatch {
                            element: id,
                            axis: stereo.axis,
                        });
                    }
                    endpoints[endpoint] = true;
                }
                for (endpoint, other) in [(axis.a(), axis.b()), (axis.b(), axis.a())] {
                    let substituent_count = self
                        .require_neighbors(endpoint)?
                        .filter(|neighbor| *neighbor != other)
                        .count();
                    if substituent_count > 2 {
                        return invalid(StereoViolation::AxisEndpointOvercoordinated {
                            element: id,
                            endpoint,
                            substituent_count,
                        });
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn validate_unique_stereo_focus(
        &self,
        id: StereoElementId,
        element: &StereoElement,
        replacing: Option<StereoElementId>,
    ) -> Result<(), MoleculeError> {
        if let Some((previous, _)) = self.stereo_elements().find(|(stored_id, stored)| {
            Some(*stored_id) != replacing && stored.kind.focus() == element.kind.focus()
        }) {
            return Err(MoleculeError::InvalidStereoReference(
                StereoViolation::DuplicateStereoFocus {
                    element: id,
                    previous,
                },
            ));
        }
        Ok(())
    }

    pub(super) fn canonicalize_stereo_element(&self, element: &mut StereoElement) {
        match &mut element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                if sort_stereo_carriers(&mut stereo.carriers) {
                    stereo.orientation = stereo.orientation.map(TetrahedralOrientation::inverted);
                }
            }
            StereoElementKind::DoubleBond(stereo) => self.canonicalize_double_bond_stereo(stereo),
            StereoElementKind::Axis(stereo) => self.canonicalize_axis_stereo(stereo),
        }
    }

    fn canonicalize_double_bond_stereo(&self, stereo: &mut DoubleBondStereo) {
        let Some(bond) = self.bond(stereo.bond) else {
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
            stereo.orientation = stereo.orientation.map(DoubleBondOrientation::inverted);
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
        let Some(axis) = self.bond(stereo.axis) else {
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
        if let Some(canonical) = self.canonical_axis_reference(left, right, stereo.axis) {
            reference_changes += usize::from(canonical != left_carrier);
            left_carrier = canonical;
        }
        if let Some(canonical) = self.canonical_axis_reference(right, left, stereo.axis) {
            reference_changes += usize::from(canonical != right_carrier);
            right_carrier = canonical;
        }

        stereo.carriers = [left_carrier, right_carrier];
        if reference_changes % 2 == 1 {
            stereo.orientation = stereo.orientation.map(AxisOrientation::inverted);
        }
    }

    fn canonical_axis_reference(
        &self,
        endpoint: AtomId,
        other_endpoint: AtomId,
        axis: BondId,
    ) -> Option<StereoCarrier> {
        self.atom_stereo_references(endpoint, other_endpoint, axis)
            .into_iter()
            .find(|carrier| {
                matches!(carrier, StereoCarrier::Atom(atom)
                if self.bond_between(other_endpoint, *atom).is_none())
            })
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
                    .is_some_and(|atom| atom.element.atomic_number() == 1);
                (u8::from(is_hydrogen), atom.raw())
            }
            StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => {
                unreachable!("represented atom references contain only atom carriers")
            }
        });
        references.dedup();
        references
    }

    fn validate_stereo_carriers(&self, carriers: &[StereoCarrier]) -> Result<(), MoleculeError> {
        for carrier in carriers {
            if let StereoCarrier::Atom(atom) = carrier {
                self.require_atom(*atom)?;
            }
        }
        Ok(())
    }

    /// Whether an atom carrier, other than `excluded`, is bonded to `center`;
    /// non-atom carriers always are.
    fn carrier_is_adjacent(
        &self,
        center: AtomId,
        excluded: Option<AtomId>,
        carrier: StereoCarrier,
    ) -> Result<bool, MoleculeError> {
        Ok(match carrier {
            StereoCarrier::Atom(atom) => {
                Some(atom) != excluded && self.require_bond_between(center, atom)?.is_some()
            }
            _ => true,
        })
    }
}

/// The one validator of represented stereo: every element's references and
/// hydrogen declarations, unique foci, and consistent groups. Each element or
/// group reports at most its first issue; the result lists every failing one
/// in element order, then group order.
pub(crate) fn stereo_violations<L: SlotLiveness>(molecule: &Graph<L>) -> Vec<StereoViolation> {
    let mut issues = Vec::new();
    let mut foci = BTreeMap::new();
    for (element_id, element) in molecule.stereo_elements() {
        if let Some(previous) = foci.insert(element.kind.focus(), element_id) {
            issues.push(StereoViolation::DuplicateStereoFocus {
                element: element_id,
                previous,
            });
            continue;
        }
        if let Err(error) = molecule.validate_stereo_element(element_id, element) {
            issues.push(match error {
                MoleculeError::InvalidStereoReference(issue) => issue,
                MoleculeError::InvalidBondId(bond) => StereoViolation::MissingStereoBond {
                    element: element_id,
                    bond,
                },
                MoleculeError::InvalidAtomId(atom) => StereoViolation::MissingStereoAtom {
                    element: element_id,
                    atom,
                },
                error => unreachable!("stereo references fail only on atoms and bonds: {error}"),
            });
            continue;
        }
        if let Some(group_id) = element.group {
            match molecule.stereo_group(group_id) {
                None => issues.push(StereoViolation::MissingStereoGroup {
                    element: element_id,
                    group: group_id,
                }),
                Some(group) if !group.members.contains(&element_id) => {
                    issues.push(StereoViolation::InconsistentStereoGroupMembership {
                        group: group_id,
                        element: element_id,
                    });
                }
                Some(_) => {}
            }
        }
    }
    for (group_id, group) in molecule.stereo_groups() {
        if let Err(issue) = validate_stereo_group(molecule, group_id, group) {
            issues.push(issue);
        }
    }
    issues
}

fn validate_stereo_group<L: SlotLiveness>(
    molecule: &Graph<L>,
    group_id: StereoGroupId,
    group: &StereoGroup,
) -> Result<(), StereoViolation> {
    if group.members.is_empty() {
        return Err(StereoViolation::EmptyStereoGroup { group: group_id });
    }
    let mut seen = BTreeSet::new();
    for &element_id in &group.members {
        if !seen.insert(element_id) {
            return Err(StereoViolation::DuplicateStereoGroupMember {
                group: group_id,
                element: element_id,
            });
        }
        let element = molecule.stereo_element(element_id).ok_or({
            StereoViolation::MissingStereoGroupMember {
                group: group_id,
                element: element_id,
            }
        })?;
        if element.group != Some(group_id) {
            return Err(StereoViolation::InconsistentStereoGroupMembership {
                group: group_id,
                element: element_id,
            });
        }
    }
    Ok(())
}

/// [`stereo_violations`] as a result for callers that need valid
/// represented stereo.
pub(crate) fn validate_represented_stereo<L: SlotLiveness>(
    molecule: &Graph<L>,
) -> Result<(), crate::algorithms::StereoValidationError> {
    let issues = stereo_violations(molecule);
    if issues.is_empty() {
        Ok(())
    } else {
        Err(crate::algorithms::StereoValidationError { violations: issues })
    }
}

/// Checks that a group's members are present and unique; `group` names it.
pub(crate) fn validate_group_members(
    group_id: StereoGroupId,
    group: &StereoGroup,
) -> Result<(), MoleculeError> {
    if group.members.is_empty() {
        return Err(MoleculeError::InvalidStereoReference(
            StereoViolation::EmptyStereoGroup { group: group_id },
        ));
    }
    let mut seen = BTreeSet::new();
    for &member in &group.members {
        if !seen.insert(member) {
            return Err(MoleculeError::InvalidStereoReference(
                StereoViolation::DuplicateStereoGroupMember {
                    group: group_id,
                    element: member,
                },
            ));
        }
    }
    Ok(())
}

fn sort_stereo_carriers(carriers: &mut [StereoCarrier]) -> bool {
    let mut indexed = carriers.iter().copied().enumerate().collect::<Vec<_>>();
    indexed.sort_by_key(|(_, carrier)| carrier.canonical_order_key());

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
    for (slot, (_, carrier)) in carriers.iter_mut().zip(indexed) {
        *slot = carrier;
    }
    odd_permutation
}

fn sorted_atom_pair(left: AtomId, right: AtomId) -> (AtomId, AtomId) {
    if left <= right {
        (left, right)
    } else {
        (right, left)
    }
}
