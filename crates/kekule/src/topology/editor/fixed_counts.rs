//! Fixed implicit-hydrogen counts that follow the bond valence an edit takes
//! from or gives to an atom, so a graph edit keeps each atom's chemistry: a
//! lost bond becomes a hydrogen and a new bond takes the place of one, as
//! inferred counts do on perception.
//!
//! A bond to a metal is coordination, not a place a hydrogen takes, so it
//! changes neither endpoint's count. An atom with declared radical electrons
//! keeps its count too. Valence an atom gains past its last hydrogen fills
//! open valence; past the atom's largest allowed valence it is owed: later
//! losses pay it back first, and publication fails while it is owed.
use super::{EditAtomId, TopologyEditError, TopologyEditor};
use crate::algorithms::{allowed_valences, explicit_valence};
use crate::molecule::{BondOrder, ImplicitHydrogens};
use std::collections::BTreeMap;

/// The bond valence one edit takes from or gives to atoms: units lost are
/// positive and units gained negative.
#[derive(Debug, Default)]
pub(super) struct ValenceChanges(BTreeMap<EditAtomId, isize>);

impl ValenceChanges {
    /// The bond between `a` and `b` changes from `from` to `to`, `None` for
    /// no bond. A bond with a metal endpoint changes nothing.
    pub(super) fn bond(
        &mut self,
        editor: &TopologyEditor,
        [a, b]: [EditAtomId; 2],
        from: Option<BondOrder>,
        to: Option<BondOrder>,
    ) -> Result<(), TopologyEditError> {
        let metal = |atom| {
            editor
                .require_atom(atom)
                .map(|atom| atom.element.is_metal())
        };
        if metal(a)? || metal(b)? {
            return Ok(());
        }
        let units = |order: Option<BondOrder>| isize::from(order.map_or(0, BondOrder::valence));
        let change = units(from) - units(to);
        for atom in [a, b] {
            *self.0.entry(atom).or_default() += change;
        }
        Ok(())
    }

    /// Drops `atom`, which the edit deletes.
    pub(super) fn forget(&mut self, atom: EditAtomId) {
        self.0.remove(&atom);
    }
}

/// An atom's fixed count and the hydrogens it owes after an edit.
#[derive(Debug, Clone, Copy)]
pub(super) struct FixedCount {
    atom: EditAtomId,
    count: u8,
    owed: usize,
}

impl TopologyEditor {
    /// The fixed counts `changes` give, checked before the edit is applied:
    /// a lost unit of bond valence first pays back an owed hydrogen, then
    /// adds one, and a gained unit removes one or, with none left, is owed.
    /// Inferred counts are left to perception.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::HydrogenCountOverflow`] when a count would pass
    /// 255.
    pub(super) fn fixed_counts_after(
        &self,
        changes: &ValenceChanges,
    ) -> Result<Vec<FixedCount>, TopologyEditError> {
        let mut counts = Vec::new();
        for (&atom, &change) in &changes.0 {
            let value = self.require_atom(atom)?;
            let ImplicitHydrogens::Fixed(count) = value.hydrogens else {
                continue;
            };
            if change == 0 || value.radical.is_some() {
                continue;
            }
            let owed = self.hydrogen_deficits.get(&atom).copied().unwrap_or(0);
            let total = isize::from(count) + change - owed as isize;
            let fixed = if total < 0 {
                FixedCount {
                    atom,
                    count: 0,
                    owed: total.unsigned_abs(),
                }
            } else {
                let count =
                    u8::try_from(total).map_err(|_| TopologyEditError::HydrogenCountOverflow {
                        atom,
                        count: total.unsigned_abs(),
                    })?;
                FixedCount {
                    atom,
                    count,
                    owed: 0,
                }
            };
            counts.push(fixed);
        }
        Ok(counts)
    }

    /// Writes counts from [`Self::fixed_counts_after`] once the edit that
    /// changed the bonds is applied. What an atom owes is only what takes it
    /// past its largest allowed valence; the rest filled open valence.
    pub(super) fn set_fixed_counts(
        &mut self,
        counts: Vec<FixedCount>,
    ) -> Result<(), TopologyEditError> {
        for FixedCount { atom, count, owed } in counts {
            let loc = self.atom_location(atom)?;
            let draft = self.draft_mut(loc.group);
            draft.require_atom_mut(loc.local)?.hydrogens = ImplicitHydrogens::Fixed(count);
            let owed = if owed == 0 {
                0
            } else {
                let graph = draft.graph();
                let valence = explicit_valence(graph, loc.local) + usize::from(count);
                graph
                    .atom(loc.local)
                    .and_then(allowed_valences)
                    .and_then(|valences| valences.iter().max())
                    .map_or(0, |&largest| {
                        owed.min(valence.saturating_sub(usize::from(largest)))
                    })
            };
            if owed == 0 {
                self.hydrogen_deficits.remove(&atom);
            } else {
                self.hydrogen_deficits.insert(atom, owed);
            }
        }
        Ok(())
    }

    /// Fixes `atom` at `count` hydrogens, owing none.
    pub(super) fn fixed_count(atom: EditAtomId, count: u8) -> FixedCount {
        FixedCount {
            atom,
            count,
            owed: 0,
        }
    }

    /// Fails for the first atom that still owes hydrogens.
    pub(super) fn check_hydrogen_deficits(&self) -> Result<(), TopologyEditError> {
        match self.hydrogen_deficits.iter().next() {
            Some((&atom, &missing)) => Err(TopologyEditError::HydrogenDeficit { atom, missing }),
            None => Ok(()),
        }
    }
}
