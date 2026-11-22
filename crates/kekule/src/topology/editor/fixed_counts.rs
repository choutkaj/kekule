//! Fixed implicit-hydrogen counts that follow the bond valence an edit takes
//! from or gives to an atom, so a graph edit keeps each atom's chemistry: a
//! lost bond becomes a hydrogen and a new bond takes the place of one, as
//! inferred counts do on perception.
use super::{EditAtomId, TopologyEditError, TopologyEditor};
use crate::molecule::{BondOrder, ImplicitHydrogens};
use std::collections::BTreeMap;

/// The bond valence one edit takes from or gives to atoms: units lost are
/// positive and units gained negative.
#[derive(Debug, Default)]
pub(super) struct ValenceChanges(BTreeMap<EditAtomId, isize>);

impl ValenceChanges {
    /// `atom` loses a bond of `order`.
    pub(super) fn lose(&mut self, atom: EditAtomId, order: BondOrder) {
        *self.0.entry(atom).or_default() += isize::from(order.valence());
    }

    /// `atom` gains a bond of `order`.
    pub(super) fn gain(&mut self, atom: EditAtomId, order: BondOrder) {
        *self.0.entry(atom).or_default() -= isize::from(order.valence());
    }

    /// A bond between `a` and `b` changes from `from` to `to`.
    pub(super) fn reorder(&mut self, [a, b]: [EditAtomId; 2], from: BondOrder, to: BondOrder) {
        for atom in [a, b] {
            self.lose(atom, from);
            self.gain(atom, to);
        }
    }
}

impl TopologyEditor {
    /// The fixed counts `changes` give, checked before the edit is applied:
    /// a lost unit of bond valence adds a hydrogen and a gained one removes
    /// one, down to zero, which leaves an atom over its valence for
    /// perception to report. Inferred counts are left to perception.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::HydrogenCountOverflow`] when a count would pass
    /// 255.
    pub(super) fn fixed_counts_after(
        &self,
        changes: &ValenceChanges,
    ) -> Result<Vec<(EditAtomId, u8)>, TopologyEditError> {
        let mut counts = Vec::new();
        for (&atom, &change) in &changes.0 {
            let ImplicitHydrogens::Fixed(count) = self.require_atom(atom)?.hydrogens else {
                continue;
            };
            let target = isize::from(count) + change;
            let target = u8::try_from(target.max(0)).map_err(|_| {
                TopologyEditError::HydrogenCountOverflow {
                    atom,
                    count: target.unsigned_abs(),
                }
            })?;
            if target != count {
                counts.push((atom, target));
            }
        }
        Ok(counts)
    }

    /// Writes counts from [`Self::fixed_counts_after`] once the edit that
    /// changed the bonds is applied.
    pub(super) fn set_fixed_counts(
        &mut self,
        counts: Vec<(EditAtomId, u8)>,
    ) -> Result<(), TopologyEditError> {
        for (atom, count) in counts {
            let loc = self.atom_location(atom)?;
            self.draft_mut(loc.group)
                .require_atom_mut(loc.local)?
                .hydrogens = ImplicitHydrogens::Fixed(count);
        }
        Ok(())
    }
}
