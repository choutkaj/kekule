//! Represented stereo through graph edits: carriers that follow a deleted
//! hydrogen or a bond that takes an implicit hydrogen's place, and the
//! elements edits remove, reported by handle.
use super::{EditAtomId, EditBondId, TopologyEditError, TopologyEditor};
use crate::algorithms::collapse_stereo_hydrogen_carrier;
use crate::molecule::BondOrder;
use crate::stereo::{StereoCarrier, StereoElement, StereoElementKind};
use std::collections::BTreeSet;

/// A represented stereo element an edit removed because it no longer
/// describes valid stereo, named by the handles of what it was about.
/// Handles of atoms and bonds the edit deleted name nothing live.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PrunedStereo {
    /// The stereo of a tetrahedral center.
    Tetrahedral {
        /// The center.
        center: EditAtomId,
    },
    /// The stereo of a double bond.
    DoubleBond {
        /// The bond's endpoints.
        atoms: [EditAtomId; 2],
    },
    /// The stereo of an axis.
    Axis {
        /// The axis bond.
        bond: EditBondId,
    },
}

impl TopologyEditor {
    /// The stereo elements edits removed since the draft started or
    /// [`Self::take_pruned_stereo`] last cleared the record, in removal
    /// order: those that adding, deleting, or rewiring bonds, changing an
    /// order, or deleting atoms left invalid. A carrier that can follow an
    /// edit does, and its element is not listed: deleting an explicit
    /// hydrogen that is its atom's only hydrogen hands its place to an
    /// implicit one, and a single bond that takes an atom's only implicit
    /// hydrogen takes its place too.
    pub fn pruned_stereo(&self) -> &[PrunedStereo] {
        &self.pruned_stereo
    }

    /// Returns [`Self::pruned_stereo`] and clears the record, so a caller can
    /// report one edit at a time.
    pub fn take_pruned_stereo(&mut self) -> Vec<PrunedStereo> {
        std::mem::take(&mut self.pruned_stereo)
    }

    /// Records the elements the draft of `group` pruned, while the group still
    /// maps every local ID the edit touched to its handle.
    pub(super) fn record_pruned(&mut self, group: usize) {
        let super::GroupSlot::Live(stored) = &mut self.groups[group] else {
            return;
        };
        let super::GroupChemistry::Draft(draft) = &mut stored.chemistry else {
            return;
        };
        let pruned = draft.take_pruned_stereo();
        for element in pruned {
            let entry = self.pruned_entry(group, &element);
            self.pruned_stereo.push(entry);
        }
    }

    fn pruned_entry(&self, group: usize, element: &StereoElement) -> PrunedStereo {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => PrunedStereo::Tetrahedral {
                center: self.group_atom(group, stereo.center),
            },
            StereoElementKind::DoubleBond(stereo) => PrunedStereo::DoubleBond {
                atoms: [stereo.left, stereo.right].map(|atom| self.group_atom(group, atom)),
            },
            StereoElementKind::Axis(stereo) => PrunedStereo::Axis {
                bond: self.group_bond(group, stereo.axis),
            },
        }
    }

    /// Before `hydrogen` is deleted: when the atom it is singly bonded to has
    /// no implicit hydrogen, is not in `handed`, and gains one, the stereo
    /// carriers naming it at that atom become the implicit hydrogen the
    /// deletion leaves, as hydrogen collapse makes them, and returns the
    /// atom. An inferred count gains one on perception; a fixed count gains
    /// one when the deletion settles it in `one_hydrogen`.
    pub(super) fn hand_carriers_to_implicit_hydrogen(
        &mut self,
        hydrogen: EditAtomId,
        handed: &BTreeSet<EditAtomId>,
        one_hydrogen: &BTreeSet<EditAtomId>,
    ) -> Result<Option<EditAtomId>, TopologyEditError> {
        // An implicit hydrogen is protium, so only a protium carrier can
        // pass its place to one.
        let atom = self.require_atom(hydrogen)?;
        if atom.element.atomic_number() != 1 || atom.isotope.is_some_and(|mass| mass.get() != 1) {
            return Ok(None);
        }
        let bonds = self.incident_bonds(hydrogen)?.collect::<Vec<_>>();
        let [bond] = bonds.as_slice() else {
            return Ok(None);
        };
        let bond = self.require_bond(*bond)?;
        if bond.order != BondOrder::Single {
            return Ok(None);
        }
        let parent = if bond.a == hydrogen { bond.b } else { bond.a };
        let (implicit, _) = self.materializable_hydrogens(parent)?;
        let fixed = !self.require_atom(parent)?.hydrogens.is_inferred();
        if implicit != 0 || handed.contains(&parent) || (fixed && !one_hydrogen.contains(&parent)) {
            return Ok(None);
        }
        let loc = self.atom_location(hydrogen)?;
        let parent_local = self.atom_location(parent)?.local;
        let elements = self
            .molecule(loc.group)
            .stereo_elements()
            .filter(|(_, element)| match &element.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    stereo.center == parent_local
                        && stereo.carriers.contains(&StereoCarrier::Atom(loc.local))
                }
                StereoElementKind::DoubleBond(stereo) => {
                    (stereo.left == parent_local
                        && stereo.left_carrier == StereoCarrier::Atom(loc.local))
                        || (stereo.right == parent_local
                            && stereo.right_carrier == StereoCarrier::Atom(loc.local))
                }
                StereoElementKind::Axis(_) => false,
            })
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        if elements.is_empty() {
            return Ok(None);
        }
        collapse_stereo_hydrogen_carrier(self.draft_mut(loc.group), loc.local, &elements)
            .expect("the elements are live");
        Ok(Some(parent))
    }
}
