//! Ring fusion chemistry: the bond of a ring fragment that fuses onto a
//! target bond, and the bond orders of the fused ring.
//!
//! The double bonds of a fragment's Kekulé structure pair its π atoms. A
//! fused ring keeps every fragment atom outside the fused bond as π or not,
//! and pairs the π atoms again: the fragment atoms, and a target atom whose
//! fragment atom is π when it has no multiple bond besides the target bond.
//! The pairing that includes the most target atoms wins, then the one that
//! keeps the target bond's order, then the one that keeps the most fragment
//! double bonds.
use crate::molecule::{AtomId, BondId, BondOrder, Molecule};
use std::collections::{BTreeMap, BTreeSet};

/// At most this many pairings are compared for one fusion.
const MAX_PAIRINGS: usize = 256;

/// One atom of the target bond, as fusion sees it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FusionEnd {
    /// The hydrogens it carries, explicit and implicit.
    pub(crate) hydrogens: usize,
    /// Whether a bond other than the target bond is a multiple bond.
    pub(crate) multiple_elsewhere: bool,
}

/// The bond orders of a fused ring.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FusionOrders {
    /// The order of each fragment bond the fusion adds: between kept
    /// fragment atoms, or from a fused fragment atom to a kept one.
    pub(crate) bonds: BTreeMap<BondId, BondOrder>,
    /// The target bond's order afterwards.
    pub(crate) shared: BondOrder,
    /// The bond valence each target atom gains, which is the number of
    /// hydrogens it gives up.
    pub(crate) gains: [usize; 2],
    /// How good the pairing is: target atoms paired, target bond order
    /// kept, fragment double bonds kept; larger is better.
    pub(crate) rank: (usize, bool, usize),
}

/// The atoms of `molecule` with exactly one double bond and no triple or
/// higher bond.
fn pi_atoms(molecule: &Molecule) -> BTreeSet<AtomId> {
    molecule
        .atom_ids()
        .filter(|&atom| {
            let orders = molecule
                .incident_bonds(atom)
                .into_iter()
                .flatten()
                .map(|(_, bond)| bond.order)
                .collect::<Vec<_>>();
            orders.iter().filter(|o| **o == BondOrder::Double).count() == 1
                && orders.iter().all(|o| o.valence() <= 2)
        })
        .collect()
}

/// The single and double bonds of `molecule` between two of `pi` atoms.
fn rebondable(molecule: &Molecule, pi: &BTreeSet<AtomId>) -> Vec<(BondId, AtomId, AtomId)> {
    molecule
        .bonds()
        .filter(|(_, bond)| {
            matches!(bond.order, BondOrder::Single | BondOrder::Double)
                && pi.contains(&bond.a())
                && pi.contains(&bond.b())
        })
        .map(|(id, bond)| (id, bond.a(), bond.b()))
        .collect()
}

/// Every perfect matching of `atoms` over `edges`, up to `MAX_PAIRINGS`, in
/// a deterministic order; each is the list of matched edges.
fn perfect_matchings(
    atoms: &BTreeSet<AtomId>,
    edges: &[(BondId, AtomId, AtomId)],
) -> Vec<Vec<BondId>> {
    fn search(
        unmatched: &mut BTreeSet<AtomId>,
        edges: &[(BondId, AtomId, AtomId)],
        chosen: &mut Vec<BondId>,
        found: &mut Vec<Vec<BondId>>,
    ) {
        if found.len() >= MAX_PAIRINGS {
            return;
        }
        let Some(&first) = unmatched.iter().next() else {
            found.push(chosen.clone());
            return;
        };
        unmatched.remove(&first);
        for &(bond, a, b) in edges {
            let other = if a == first {
                b
            } else if b == first {
                a
            } else {
                continue;
            };
            if unmatched.remove(&other) {
                chosen.push(bond);
                search(unmatched, edges, chosen, found);
                chosen.pop();
                unmatched.insert(other);
            }
        }
        unmatched.insert(first);
    }
    let edges = edges
        .iter()
        .copied()
        .filter(|(_, a, b)| atoms.contains(a) && atoms.contains(b))
        .collect::<Vec<_>>();
    let mut found = Vec::new();
    search(&mut atoms.clone(), &edges, &mut Vec::new(), &mut found);
    found
}

/// The hydrogens bonded to `atom` in `molecule`.
fn hydrogens_on(molecule: &Molecule, atom: AtomId) -> Vec<AtomId> {
    molecule
        .neighbors(atom)
        .into_iter()
        .flatten()
        .filter(|&neighbor| is_hydrogen(molecule, neighbor))
        .collect()
}

fn is_hydrogen(molecule: &Molecule, atom: AtomId) -> bool {
    molecule
        .atom(atom)
        .is_some_and(|value| value.element.atomic_number() == 1)
}

/// The one heavy neighbor of `atom` other than `other`, when it has
/// exactly one.
pub(crate) fn single_heavy_neighbor(
    molecule: &Molecule,
    atom: AtomId,
    other: AtomId,
) -> Option<AtomId> {
    let mut heavy = molecule
        .neighbors(atom)
        .into_iter()
        .flatten()
        .filter(|&neighbor| neighbor != other && !is_hydrogen(molecule, neighbor));
    let first = heavy.next()?;
    heavy.next().is_none().then_some(first)
}

/// The atoms a fusion through the fragment bond between `fused` leaves out:
/// those two atoms and their hydrogens.
pub(crate) fn fusion_dropped(molecule: &Molecule, fused: [AtomId; 2]) -> BTreeSet<AtomId> {
    fused
        .iter()
        .flat_map(|&atom| std::iter::once(atom).chain(hydrogens_on(molecule, atom)))
        .collect()
}

/// The default fusion bond of a ring fragment whose ring atoms are `ring`:
/// the first bond, in atom order, between two ring carbons that each carry
/// a hydrogen and one other heavy neighbor, and without which the
/// fragment's other π atoms still pair up, so the fused ring keeps its
/// double bonds whatever the target bond's order.
pub(crate) fn default_fusion_bond(molecule: &Molecule, ring: &BTreeSet<AtomId>) -> Option<BondId> {
    let pi = pi_atoms(molecule);
    let edges = rebondable(molecule, &pi);
    let order = molecule.atom_ids().collect::<Vec<_>>();
    let place = |atom: AtomId| order.iter().position(|id| *id == atom);
    let mut bonds = molecule
        .bonds()
        .map(|(id, bond)| {
            let (a, b) = (place(bond.a()), place(bond.b()));
            ((a.min(b), a.max(b)), id, bond.a(), bond.b())
        })
        .collect::<Vec<_>>();
    bonds.sort_by_key(|(key, ..)| *key);
    bonds.into_iter().find_map(|(_, id, a, b)| {
        let fusable = |atom: AtomId, other: AtomId| {
            ring.contains(&atom)
                && molecule
                    .atom(atom)
                    .is_some_and(|value| value.element.atomic_number() == 6)
                && !hydrogens_on(molecule, atom).is_empty()
                && single_heavy_neighbor(molecule, atom, other).is_some()
        };
        let rest = &pi - &BTreeSet::from([a, b]);
        (fusable(a, b) && fusable(b, a) && !perfect_matchings(&rest, &edges).is_empty())
            .then_some(id)
    })
}

/// The bond orders of the ring the fragment `molecule` fuses through its
/// bond between `fused` onto a target bond of order `shared` (single or
/// double), `fused[i]` merging into the target atom `ends[i]`; `None` when
/// no pairing leaves each target atom enough hydrogens to give up.
pub(crate) fn fusion_orders(
    molecule: &Molecule,
    fused: [AtomId; 2],
    shared: BondOrder,
    ends: [FusionEnd; 2],
) -> Option<FusionOrders> {
    let pi = pi_atoms(molecule);
    let edges = rebondable(molecule, &pi);
    let dropped = fusion_dropped(molecule, fused);
    let fused_bond = molecule.bond_between(fused[0], fused[1])?;
    let rest = &pi - &BTreeSet::from(fused);
    let eligible = |end: usize| pi.contains(&fused[end]) && !ends[end].multiple_elsewhere;
    // A double target bond keeps both its atoms paired or neither.
    let subsets: &[&[usize]] = if shared == BondOrder::Double {
        &[&[0, 1], &[]]
    } else {
        &[&[0, 1], &[0], &[1], &[]]
    };
    let mut best: Option<FusionOrders> = None;
    for subset in subsets {
        if !subset.iter().all(|&end| eligible(end)) {
            continue;
        }
        let mut atoms = rest.clone();
        atoms.extend(subset.iter().map(|&end| fused[end]));
        for matching in perfect_matchings(&atoms, &edges) {
            let paired = |bond: BondId| matching.contains(&bond);
            let mut bonds = BTreeMap::new();
            let mut kept = 0;
            for (id, bond) in molecule.bonds() {
                let lost = |atom: AtomId| dropped.contains(&atom) && !fused.contains(&atom);
                if id == fused_bond || lost(bond.a()) || lost(bond.b()) {
                    continue;
                }
                let order = if edges.iter().any(|(edge, ..)| *edge == id) {
                    let order = if paired(id) {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    };
                    kept += usize::from(order == bond.order);
                    order
                } else {
                    bond.order
                };
                bonds.insert(id, order);
            }
            let new_shared = if subset.len() == 2 {
                if paired(fused_bond) {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                }
            } else {
                shared
            };
            let mut gains = [0usize; 2];
            let mut feasible = true;
            for end in 0..2 {
                let added = molecule
                    .incident_bonds(fused[end])
                    .into_iter()
                    .flatten()
                    .filter_map(|(id, _)| bonds.get(&id))
                    .map(|order| usize::from(order.valence()))
                    .sum::<usize>();
                let gain = (added + usize::from(new_shared.valence()))
                    .checked_sub(usize::from(shared.valence()));
                match gain {
                    Some(gain) if gain <= ends[end].hydrogens => gains[end] = gain,
                    _ => feasible = false,
                }
            }
            if !feasible {
                continue;
            }
            let rank = (subset.len(), new_shared == shared, kept);
            if best.as_ref().is_none_or(|best| rank > best.rank) {
                best = Some(FusionOrders {
                    bonds,
                    shared: new_shared,
                    gains,
                    rank,
                });
            }
        }
    }
    best
}

#[cfg(test)]
mod tests;
