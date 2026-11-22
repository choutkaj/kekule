//! Fixed implicit-hydrogen counts that follow the bond valence an edit takes
//! from or gives to an atom, so a graph edit keeps each atom's chemistry: a
//! lost bond becomes a hydrogen and a new bond takes the place of one, as
//! inferred counts do on perception.
//!
//! Every bond counts with its order's valence, so a single bond to a metal
//! trades a hydrogen as any single bond does, and a dative or zero-order
//! bond, which coordination is, trades none. A metal's own count never
//! changes. A loss fills the atom up to the smallest valence its element and
//! charge allow at or above what remains, counting its hydrogens and radical
//! electrons, so a sulfoxide that loses its oxygen is a sulfide. A gain takes
//! a radical electron, then a hydrogen, then fills open valence; past the
//! atom's largest allowed valence it is owed, later losses pay it back
//! first, and publication fails while it is owed.
use super::{EditAtomId, TopologyEditError, TopologyEditor};
use crate::algorithms::{allowed_valences, is_metal};
use crate::molecule::{Atom, AtomRadical, BondOrder, ImplicitHydrogens};
use std::collections::{BTreeMap, BTreeSet};

/// The valence a bond of `order` gives each endpoint.
pub(crate) fn bond_valence(order: BondOrder) -> usize {
    usize::from(order.valence())
}

/// An atom's fixed count, radical, and owed hydrogens after an edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Settled {
    pub(crate) count: u8,
    pub(crate) radical: Option<AtomRadical>,
    pub(crate) owed: usize,
}

fn radical_electrons(atom: &Atom) -> usize {
    atom.radical
        .map_or(0, |radical| usize::from(radical.electron_count()))
}

fn allowed(atom: &Atom) -> impl Iterator<Item = usize> + '_ {
    allowed_valences(atom)
        .unwrap_or(&[])
        .iter()
        .map(|&value| usize::from(value))
}

/// What `atom` owes with valence `total` from bonds, hydrogens, and radical
/// electrons: what lies past its largest allowed valence, at most `owed`.
fn still_owed(atom: &Atom, owed: usize, total: usize) -> usize {
    allowed(atom)
        .max()
        .map_or(0, |largest| owed.min(total.saturating_sub(largest)))
}

/// The fixed count, radical, and owed hydrogens of `atom`, owing `owed`,
/// after its bond valence goes from `before` to `after`; `None` when nothing
/// changes: the count is inferred, the atom is a metal, or the valence
/// stays.
///
/// # Errors
///
/// The count, when it would pass 255.
pub(crate) fn settle(
    atom: &Atom,
    mut owed: usize,
    before: usize,
    after: usize,
) -> Result<Option<Settled>, usize> {
    let ImplicitHydrogens::Fixed(count) = atom.hydrogens else {
        return Ok(None);
    };
    if is_metal(atom.element) || before == after {
        return Ok(None);
    }
    let electrons = radical_electrons(atom);
    let mut count = usize::from(count);
    let mut remaining = electrons;
    if after < before {
        let filled = after + count + electrons;
        if let Some(target) = allowed(atom).filter(|&value| value >= filled).min() {
            count += target - filled;
        }
    } else {
        let mut gain = after - before;
        let taken = gain.min(remaining);
        remaining -= taken;
        gain -= taken;
        let taken = gain.min(count);
        count -= taken;
        gain -= taken;
        owed += gain;
    }
    finish(atom, count, remaining, owed, after)
}

fn finish(
    atom: &Atom,
    count: usize,
    electrons: usize,
    owed: usize,
    valence: usize,
) -> Result<Option<Settled>, usize> {
    let owed = still_owed(atom, owed, valence + count + electrons);
    let radical = if electrons == radical_electrons(atom) {
        atom.radical
    } else {
        AtomRadical::new(
            u8::try_from(electrons).expect("no more electrons than declared"),
            None,
        )
    };
    Ok(Some(Settled {
        count: u8::try_from(count).map_err(|_| count)?,
        radical,
        owed,
    }))
}

/// What `atom`, owing `owed` with bond valence `valence`, still owes once it
/// replaces the atom it was: only valence past the largest its new state
/// allows.
pub(crate) fn owed_after_replace(atom: &Atom, owed: usize, valence: usize) -> usize {
    let ImplicitHydrogens::Fixed(count) = atom.hydrogens else {
        return 0;
    };
    still_owed(
        atom,
        owed,
        valence + usize::from(count) + radical_electrons(atom),
    )
}

/// The bond valence one edit gives to or takes from atoms.
#[derive(Debug, Default)]
pub(super) struct ValenceChanges(BTreeMap<EditAtomId, isize>);

impl ValenceChanges {
    /// The bond between `a` and `b` changes from `from` to `to`, `None` for
    /// no bond.
    pub(super) fn bond(
        &mut self,
        [a, b]: [EditAtomId; 2],
        from: Option<BondOrder>,
        to: Option<BondOrder>,
    ) {
        self.change(a, from, to);
        self.change(b, from, to);
    }

    /// One bond of `atom` changes from `from` to `to`.
    pub(super) fn change(
        &mut self,
        atom: EditAtomId,
        from: Option<BondOrder>,
        to: Option<BondOrder>,
    ) {
        let units = |order: Option<BondOrder>| isize::from(order.map_or(0, BondOrder::valence));
        *self.0.entry(atom).or_default() += units(to) - units(from);
    }
}

/// An atom's fixed count, radical, and owed hydrogens after an edit.
#[derive(Debug, Clone, Copy)]
pub(super) struct FixedCount {
    atom: EditAtomId,
    settled: Settled,
}

impl FixedCount {
    /// The atoms `counts` leave with exactly one implicit hydrogen.
    pub(super) fn with_one_hydrogen(counts: &[Self]) -> BTreeSet<EditAtomId> {
        counts
            .iter()
            .filter(|count| count.settled.count == 1)
            .map(|count| count.atom)
            .collect()
    }
}

impl TopologyEditor {
    /// The bond valence of `atom` in the draft.
    pub(super) fn bond_valence(&self, atom: EditAtomId) -> Result<usize, TopologyEditError> {
        let mut valence = 0;
        for bond in self.incident_bonds(atom)?.collect::<Vec<_>>() {
            valence += bond_valence(self.require_bond(bond)?.order);
        }
        Ok(valence)
    }

    /// The fixed counts `changes` give, checked before the edit is applied;
    /// see [`settle`]. Inferred counts are left to perception.
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
            let before = self.bond_valence(atom)?;
            let after = before.saturating_add_signed(change);
            let owed = self.hydrogen_deficits.get(&atom).copied().unwrap_or(0);
            let settled = settle(value, owed, before, after)
                .map_err(|count| TopologyEditError::HydrogenCountOverflow { atom, count })?;
            if let Some(settled) = settled {
                counts.push(FixedCount { atom, settled });
            }
        }
        Ok(counts)
    }

    /// Writes counts from [`Self::fixed_counts_after`] once the edit that
    /// changed the bonds is applied.
    pub(super) fn set_fixed_counts(
        &mut self,
        counts: Vec<FixedCount>,
    ) -> Result<(), TopologyEditError> {
        for FixedCount { atom, settled } in counts {
            let loc = self.atom_location(atom)?;
            let mut value = self.draft_mut(loc.group).require_atom_mut(loc.local)?;
            value.hydrogens = ImplicitHydrogens::Fixed(settled.count);
            value.radical = settled.radical;
            if settled.owed == 0 {
                self.hydrogen_deficits.remove(&atom);
            } else {
                self.hydrogen_deficits.insert(atom, settled.owed);
            }
        }
        Ok(())
    }

    /// Keeps of what `atom` owes only what its replaced state still owes;
    /// see [`owed_after_replace`].
    pub(super) fn settle_replaced_deficit(
        &mut self,
        atom: EditAtomId,
    ) -> Result<(), TopologyEditError> {
        let Some(&owed) = self.hydrogen_deficits.get(&atom) else {
            return Ok(());
        };
        let owed = owed_after_replace(self.require_atom(atom)?, owed, self.bond_valence(atom)?);
        if owed == 0 {
            self.hydrogen_deficits.remove(&atom);
        } else {
            self.hydrogen_deficits.insert(atom, owed);
        }
        Ok(())
    }

    /// Fails for the first atom that still owes hydrogens.
    pub(super) fn check_hydrogen_deficits(&self) -> Result<(), TopologyEditError> {
        match self.hydrogen_deficits.iter().next() {
            Some((&atom, &missing)) => Err(TopologyEditError::HydrogenDeficit { atom, missing }),
            None => Ok(()),
        }
    }
}
