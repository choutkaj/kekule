//! Fixed implicit-hydrogen counts that follow the bond valence an edit takes
//! from or gives to an atom, so a graph edit keeps each atom's chemistry: a
//! lost bond becomes a hydrogen and a new bond takes the place of one, as
//! inferred counts do on perception.
//!
//! A metal's count never changes. A non-metal's bonds to metals count as
//! covalent unless it donates to them, as an ammine nitrogen does; donor
//! bonds are coordination and change no count. A lost unit of valence first
//! pays back valence the atom's radical electrons and hydrogens over-fill,
//! then adds a hydrogen. A gained unit takes a radical electron, then a
//! hydrogen, then fills open valence; past the atom's largest allowed
//! valence it is owed, later losses pay it back first, and publication
//! fails while it is owed.
use super::{EditAtomId, TopologyEditError, TopologyEditor};
use crate::algorithms::allowed_valences;
use crate::molecule::{Atom, AtomRadical, BondOrder, ImplicitHydrogens};
use std::collections::BTreeMap;

/// An atom's bond valence, split by whether the neighbor is a metal.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct BondValence {
    /// Valence of bonds to non-metals.
    pub(crate) other: usize,
    /// Valence of bonds to metals.
    pub(crate) metal: usize,
}

impl BondValence {
    /// Adds the valence of a bond of `order`, to a metal or not.
    pub(crate) fn add(&mut self, metal: bool, order: BondOrder) {
        let units = usize::from(order.valence());
        if metal {
            self.metal += units;
        } else {
            self.other += units;
        }
    }
}

/// An atom's fixed count, radical, and owed hydrogens after an edit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Settled {
    pub(crate) count: u8,
    pub(crate) radical: Option<AtomRadical>,
    pub(crate) owed: usize,
}

/// Whether a non-metal donates to metals when its other bonds and
/// hydrogens give it the valence `other` and its bonds to metals add
/// `metal`: `other + metal` is not a valence its element and charge allow,
/// but some allowed valence lies from `other` up to, and not including,
/// `other + metal`.
fn donates_to_metals(atom: &Atom, other: usize, metal: usize) -> bool {
    let Some(allowed) = allowed_valences(atom) else {
        return false;
    };
    let allowed = || allowed.iter().map(|value| usize::from(*value));
    !allowed().any(|value| value == other + metal)
        && allowed().any(|value| (other..other + metal).contains(&value))
}

/// The valence of `atom`'s bonds that hydrogens trade against: its bonds to
/// non-metals, and its bonds to metals unless it donates to them.
fn traded_valence(atom: &Atom, count: usize, electrons: usize, bonds: BondValence) -> usize {
    let donor =
        bonds.metal > 0 && donates_to_metals(atom, bonds.other + count + electrons, bonds.metal);
    bonds.other + if donor { 0 } else { bonds.metal }
}

fn radical_electrons(atom: &Atom) -> usize {
    atom.radical
        .map_or(0, |radical| usize::from(radical.electron_count()))
}

fn largest_valence(atom: &Atom) -> Option<usize> {
    allowed_valences(atom)?
        .iter()
        .max()
        .map(|&value| usize::from(value))
}

/// The fixed count, radical, and owed hydrogens of `atom`, owing `owed`,
/// after its bonds go from `before` to `after`; `None` when nothing
/// changes: the count is inferred, the atom is a metal, or the traded
/// valence stays.
///
/// # Errors
///
/// The count, when it would pass 255.
pub(crate) fn settle(
    atom: &Atom,
    owed: usize,
    before: BondValence,
    after: BondValence,
) -> Result<Option<Settled>, usize> {
    let ImplicitHydrogens::Fixed(count) = atom.hydrogens else {
        return Ok(None);
    };
    if atom.element.is_metal() {
        return Ok(None);
    }
    let count = usize::from(count);
    let electrons = radical_electrons(atom);
    let from = traded_valence(atom, count, electrons, before);
    let to = traded_valence(atom, count, electrons, after);
    let mut new_count = count;
    let mut new_electrons = electrons;
    let new_owed = if from > to {
        let loss = from - to;
        let total = from + count + electrons;
        let target = allowed_valences(atom)
            .unwrap_or(&[])
            .iter()
            .map(|&value| usize::from(value))
            .find(|&value| value >= from + count)
            .or_else(|| largest_valence(atom))
            .unwrap_or(total);
        let excess = total.saturating_sub(target);
        new_count += loss - loss.min(excess);
        owed.saturating_sub(loss)
    } else if to > from {
        let mut gain = to - from;
        let taken = gain.min(new_electrons);
        new_electrons -= taken;
        gain -= taken;
        let taken = gain.min(new_count);
        new_count -= taken;
        gain -= taken;
        let total = to + new_count + new_electrons;
        largest_valence(atom).map_or(0, |largest| {
            (owed + gain).min(total.saturating_sub(largest))
        })
    } else {
        return Ok(None);
    };
    let count = u8::try_from(new_count).map_err(|_| new_count)?;
    let radical = if new_electrons == electrons {
        atom.radical
    } else {
        AtomRadical::new(
            u8::try_from(new_electrons).expect("no more electrons than declared"),
            None,
        )
    };
    Ok(Some(Settled {
        count,
        radical,
        owed: new_owed,
    }))
}

/// What `atom`, owing `owed` with bonds `bonds`, still owes once it replaces
/// the atom it was: only valence past the largest its new state allows.
pub(crate) fn owed_after_replace(atom: &Atom, owed: usize, bonds: BondValence) -> usize {
    let ImplicitHydrogens::Fixed(count) = atom.hydrogens else {
        return 0;
    };
    let count = usize::from(count);
    let electrons = radical_electrons(atom);
    let total = traded_valence(atom, count, electrons, bonds) + count + electrons;
    largest_valence(atom).map_or(0, |largest| owed.min(total.saturating_sub(largest)))
}

/// The bond valence one edit gives to or takes from atoms, as changes to
/// their bonds to non-metals and to metals.
#[derive(Debug, Default)]
pub(super) struct ValenceChanges(BTreeMap<EditAtomId, (isize, isize)>);

impl ValenceChanges {
    /// The bond between `a` and `b` changes from `from` to `to`, `None` for
    /// no bond.
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
        let (a_metal, b_metal) = (metal(a)?, metal(b)?);
        self.change(a, b_metal, from, to);
        self.change(b, a_metal, from, to);
        Ok(())
    }

    /// One bond of `atom`, to a metal or not, changes from `from` to `to`.
    pub(super) fn change(
        &mut self,
        atom: EditAtomId,
        metal: bool,
        from: Option<BondOrder>,
        to: Option<BondOrder>,
    ) {
        let units = |order: Option<BondOrder>| isize::from(order.map_or(0, BondOrder::valence));
        let change = units(to) - units(from);
        let entry = self.0.entry(atom).or_default();
        if metal {
            entry.1 += change;
        } else {
            entry.0 += change;
        }
    }

    /// Drops `atom`, which the edit deletes.
    pub(super) fn forget(&mut self, atom: EditAtomId) {
        self.0.remove(&atom);
    }
}

/// An atom's fixed count, radical, and owed hydrogens after an edit.
#[derive(Debug, Clone, Copy)]
pub(super) struct FixedCount {
    atom: EditAtomId,
    settled: Settled,
}

impl TopologyEditor {
    /// The bond valence of `atom` in the draft.
    pub(super) fn bond_valence(&self, atom: EditAtomId) -> Result<BondValence, TopologyEditError> {
        let mut valence = BondValence::default();
        for bond in self.incident_bonds(atom)?.collect::<Vec<_>>() {
            let bond = self.require_bond(bond)?;
            let other = if bond.a == atom { bond.b } else { bond.a };
            valence.add(self.require_atom(other)?.element.is_metal(), bond.order);
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
        for (&atom, &(other, metal)) in &changes.0 {
            let value = self.require_atom(atom)?;
            let before = self.bond_valence(atom)?;
            let after = BondValence {
                other: before.other.saturating_add_signed(other),
                metal: before.metal.saturating_add_signed(metal),
            };
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

    /// Fixes `atom` at `count` hydrogens, keeping its radical and owing
    /// none.
    pub(super) fn fixed_count(
        &self,
        atom: EditAtomId,
        count: u8,
    ) -> Result<FixedCount, TopologyEditError> {
        Ok(FixedCount {
            atom,
            settled: Settled {
                count,
                radical: self.require_atom(atom)?.radical,
                owed: 0,
            },
        })
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
