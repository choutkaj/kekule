//! Resonance delocalization: formal charges moved along alternating paths.
//!
//! Only charged atoms start a search, and the search reads only the atoms and
//! bonds its paths reach, so neutral regions of a large molecule cost one scan
//! of the formal charges.
use super::super::rings::compute_ring_membership;
use super::super::valence::{allowed_valences_with_charge, explicit_valence};
use super::ResonancePerceptionError;
use crate::diagnostics::{Cancellation, LimitExceeded, Resource, UnknownHydrogenCount};
use crate::molecule::{
    AromaticityPerception, AtomId, BondId, BondOrder, Molecule, ResonanceGroup, RingMembership,
};
use std::borrow::Cow;

const C: u8 = 6;
const N: u8 = 7;
const O: u8 = 8;
const P: u8 = 15;
const S: u8 = 16;
const AS: u8 = 33;
const SE: u8 = 34;

/// Atoms a path may pass through.
const PATH: [u8; 9] = [5, C, N, O, 14, P, S, AS, SE];
/// Atoms that carry a mobile anion.
const ANION: [u8; 5] = [C, N, O, S, SE];
/// Atoms that carry a mobile onium cation.
const ONIUM: [u8; 5] = [N, O, P, S, SE];

/// The work bound grows with the graph, as other default perception bounds do.
const MIN_WORK: usize = 5_000_000;
const WORK_PER_SLOT: usize = 1_000;

/// Bond-order steps along a path.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Step {
    Raise,
    Lower,
}
impl Step {
    const fn other(self) -> Self {
        match self {
            Self::Raise => Self::Lower,
            Self::Lower => Self::Raise,
        }
    }
}

/// Mobile charges; each moves only to an atom that then carries the same kind.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Charge {
    /// A lone-pair anion, such as carboxylate oxygen.
    Anion,
    /// A cation with a complete octet, such as iminium nitrogen.
    Onium,
    /// A carbocation with an empty orbital.
    Vacant,
}
impl Charge {
    const fn first(self) -> Step {
        match self {
            Self::Anion | Self::Vacant => Step::Raise,
            Self::Onium => Step::Lower,
        }
    }
}

#[derive(Clone, Copy)]
struct AtomState {
    number: u8,
    charge: i8,
    valence: usize,
    radical: bool,
}

/// Bond orders the search can change; others are read as `FIXED`.
const FIXED: u8 = 0;
const UNREAD: u8 = u8::MAX;

struct Search<'a> {
    mol: &'a Molecule,
    aromaticity: &'a AromaticityPerception,
    rings: Cow<'a, RingMembership>,
    atoms: Vec<Option<AtomState>>,
    orders: Vec<u8>,
    left: usize,
    limit: usize,
    cancellation: Option<&'a Cancellation>,
}

/// Groups of resonance-delocalized bonds, as installed by
/// [`perceive_resonance`](super::perceive_resonance).
pub(super) fn delocalized_groups<'a>(
    mol: &'a Molecule,
    cancellation: Option<&'a Cancellation>,
) -> Result<Vec<ResonanceGroup>, ResonancePerceptionError> {
    let aromaticity = mol
        .perception()
        .aromaticity()
        .ok_or(ResonancePerceptionError::MissingAromaticity)?;
    let charged: Vec<_> = mol
        .atoms()
        .filter(|(_, atom)| atom.formal_charge != 0)
        .map(|(id, _)| id.index())
        .collect();
    if charged.is_empty() {
        return Ok(Vec::new());
    }
    let rings = match mol.perception().rings() {
        Some(rings) => Cow::Borrowed(rings.membership()),
        None => Cow::Owned(compute_ring_membership(mol)),
    };
    let slots = mol.graph.atom_slot_count() + mol.graph.bond_slot_count();
    let limit = MIN_WORK.max(slots.saturating_mul(WORK_PER_SLOT));
    let mut search = Search {
        mol,
        aromaticity,
        rings,
        atoms: vec![None; mol.graph.atom_slot_count()],
        orders: vec![UNREAD; mol.graph.bond_slot_count()],
        left: limit,
        limit,
        cancellation,
    };
    for &atom in &charged {
        search.normalize_ylides(atom)?;
    }
    let mut delocalized = vec![false; search.orders.len()];
    let mut visited = vec![false; search.atoms.len()];
    for &carrier in &charged {
        if let Some(charge) = search.mobile_charge(carrier)? {
            search.explore(carrier, charge, &mut visited, &mut delocalized)?;
        }
    }
    Ok(search.groups(&delocalized))
}

impl Search<'_> {
    fn allows(number: u8, charge: i8, valence: usize) -> bool {
        allowed_valences_with_charge(number, charge)
            .is_some_and(|allowed| allowed.iter().any(|&v| usize::from(v) == valence))
    }

    fn atom(&mut self, index: usize) -> Result<AtomState, ResonancePerceptionError> {
        if let Some(state) = self.atoms[index] {
            return Ok(state);
        }
        let id = index_atom(index);
        let atom = self.mol.atom(id).expect("adjacency reaches live atoms");
        let hydrogens = self.mol.implicit_hydrogens(id).ok_or(
            ResonancePerceptionError::UnknownHydrogenCount(UnknownHydrogenCount::new(id)),
        )?;
        let state = AtomState {
            number: atom.element.atomic_number(),
            charge: atom.formal_charge,
            valence: explicit_valence(self.mol, id) + hydrogens,
            radical: atom.radical.is_some(),
        };
        self.atoms[index] = Some(state);
        Ok(state)
    }

    /// The order of a bond the search may change, or [`FIXED`] for aromatic,
    /// zero-order, dative and quadruple bonds.
    fn order(&mut self, index: usize) -> u8 {
        if self.orders[index] == UNREAD {
            let id = index_bond(index);
            let bond = self.mol.bond(id).expect("adjacency reaches live bonds");
            self.orders[index] = match bond.order {
                _ if self.aromaticity.bond_is_aromatic(id) => FIXED,
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Zero | BondOrder::Dative | BondOrder::Quadruple => FIXED,
            };
        }
        self.orders[index]
    }

    /// Reads each ylidic X⁺–Y⁻ pair of hypervalent phosphorus, arsenic, sulfur
    /// or selenium and a terminal chalcogen as X=Y, so both drawings of a
    /// phosphate or sulfonate delocalize alike. Terminal partners are taken in
    /// bond ID order.
    fn normalize_ylides(&mut self, center: usize) -> Result<(), ResonancePerceptionError> {
        let state = self.atom(center)?;
        if !matches!(state.number, P | AS | S | SE) || state.charge <= 0 || state.radical {
            return Ok(());
        }
        let mut bonds: Vec<_> = self.mol.graph.adjacency[center]
            .iter()
            .map(|bond| bond.index())
            .collect();
        bonds.sort_unstable();
        for bond in bonds {
            let terminal = self.other(bond, center);
            let (x, y) = (self.atom(center)?, self.atom(terminal)?);
            let pair = self.order(bond) == 1
                && x.charge > 0
                && matches!(y.number, O | S | SE)
                && y.charge == -1
                && !y.radical
                && y.valence == 1
                && Self::allows(x.number, x.charge - 1, x.valence + 1)
                && Self::allows(y.number, 0, 2);
            if pair {
                self.orders[bond] = 2;
                self.atoms[center] = Some(AtomState {
                    charge: x.charge - 1,
                    valence: x.valence + 1,
                    ..x
                });
                self.atoms[terminal] = Some(AtomState {
                    charge: 0,
                    valence: 2,
                    ..y
                });
            }
        }
        Ok(())
    }

    fn other(&self, bond: usize, atom: usize) -> usize {
        let bond = self
            .mol
            .bond(index_bond(bond))
            .expect("adjacency reaches live bonds");
        let (a, b) = (bond.a().index(), bond.b().index());
        if a == atom {
            b
        } else {
            a
        }
    }

    fn mobile_charge(&mut self, index: usize) -> Result<Option<Charge>, ResonancePerceptionError> {
        let atom = self.atom(index)?;
        if atom.radical {
            return Ok(None);
        }
        Ok(match atom.charge {
            -1 if ANION.contains(&atom.number)
                && Self::allows(atom.number, 0, atom.valence + 1) =>
            {
                Some(Charge::Anion)
            }
            1 if ONIUM.contains(&atom.number)
                && Self::allows(atom.number, 0, atom.valence.wrapping_sub(1)) =>
            {
                Some(Charge::Onium)
            }
            1 if atom.number == C && Self::allows(C, 0, atom.valence + 1) => Some(Charge::Vacant),
            _ => None,
        })
    }

    /// Whether an atom reached by `last` can take over `charge`.
    fn accepts(atom: AtomState, charge: Charge, last: Step) -> bool {
        if atom.charge != 0 || atom.radical {
            return false;
        }
        match (charge, last) {
            (Charge::Anion, Step::Lower) => {
                ANION.contains(&atom.number)
                    && Self::allows(atom.number, -1, atom.valence.wrapping_sub(1))
            }
            (Charge::Onium, Step::Raise) => {
                ONIUM.contains(&atom.number) && Self::allows(atom.number, 1, atom.valence + 1)
            }
            (Charge::Vacant, Step::Lower) => {
                atom.number == C && Self::allows(C, 1, atom.valence.wrapping_sub(1))
            }
            _ => false,
        }
    }

    /// Whether `step` can change a bond: raising a ring bond to triple and
    /// lowering a single bond are not.
    fn steps(&mut self, bond: usize, step: Step) -> bool {
        match (step, self.order(bond)) {
            (Step::Raise, 1) => true,
            (Step::Raise, 2) => !self.rings.bond_in_ring(index_bond(bond)),
            (Step::Lower, 2 | 3) => true,
            _ => false,
        }
    }

    fn charge_work(&mut self) -> Result<(), ResonancePerceptionError> {
        self.left = self.left.checked_sub(1).ok_or_else(|| {
            ResonancePerceptionError::LimitExceeded(LimitExceeded::first_beyond(
                Resource::Work,
                self.limit,
            ))
        })?;
        Cancellation::poll_every(self.cancellation, self.left)
            .map_err(ResonancePerceptionError::Cancelled)
    }

    /// Marks every bond of every simple alternating path that moves the
    /// carrier's charge to an atom able to take it over.
    fn explore(
        &mut self,
        carrier: usize,
        charge: Charge,
        visited: &mut [bool],
        delocalized: &mut [bool],
    ) -> Result<(), ResonancePerceptionError> {
        let mol = self.mol;
        visited[carrier] = true;
        let mut frames = vec![(carrier, charge.first(), 0_usize)];
        let mut path: Vec<usize> = Vec::new();
        // The path's first `marked` bonds are already delocalized, so each
        // push is marked at most once.
        let mut marked = 0;
        while let Some(frame) = frames.last_mut() {
            let (atom, step, cursor) = *frame;
            let Some(&bond) = mol.graph.adjacency[atom].get(cursor) else {
                visited[atom] = false;
                frames.pop();
                path.pop();
                marked = marked.min(path.len());
                continue;
            };
            frame.2 += 1;
            self.charge_work()?;
            let bond = bond.index();
            let next = self.other(bond, atom);
            if visited[next] || !self.steps(bond, step) {
                continue;
            }
            let state = self.atom(next)?;
            if !PATH.contains(&state.number) {
                continue;
            }
            path.push(bond);
            if Self::accepts(state, charge, step) {
                for &member in &path[marked..] {
                    delocalized[member] = true;
                }
                marked = path.len();
            }
            visited[next] = true;
            frames.push((next, step.other(), 0));
        }
        Ok(())
    }

    /// Connected components of the delocalized bonds.
    fn groups(&self, delocalized: &[bool]) -> Vec<ResonanceGroup> {
        let mut seen = vec![false; delocalized.len()];
        let mut reached = vec![false; self.atoms.len()];
        let mut groups = Vec::new();
        for root in 0..delocalized.len() {
            if !delocalized[root] || seen[root] {
                continue;
            }
            seen[root] = true;
            let mut bonds = vec![root];
            let mut atoms = Vec::new();
            let mut cursor = 0;
            while let Some(&bond) = bonds.get(cursor) {
                cursor += 1;
                let ends = self
                    .mol
                    .bond(index_bond(bond))
                    .expect("delocalized bonds are live");
                for atom in [ends.a().index(), ends.b().index()] {
                    if std::mem::replace(&mut reached[atom], true) {
                        continue;
                    }
                    atoms.push(atom);
                    for &next in &self.mol.graph.adjacency[atom] {
                        let next = next.index();
                        if delocalized[next] && !seen[next] {
                            seen[next] = true;
                            bonds.push(next);
                        }
                    }
                }
            }
            atoms.sort_unstable();
            bonds.sort_unstable();
            groups.push(ResonanceGroup {
                atoms: atoms.into_iter().map(index_atom).collect(),
                bonds: bonds.into_iter().map(index_bond).collect(),
            });
        }
        groups
    }
}

fn index_atom(index: usize) -> AtomId {
    AtomId::new(u32::try_from(index).expect("atom slot within ID capacity"))
}

fn index_bond(index: usize) -> BondId {
    BondId::new(u32::try_from(index).expect("bond slot within ID capacity"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn long_paths_mark_each_bond_once_and_poll_cancellation() {
        let text = format!("[CH2-]{}", "C=C".repeat(5_000));
        let mut mol = crate::smiles::to_molecules(&text).unwrap().remove(0);
        mol.perceive().unwrap();
        let groups = delocalized_groups(&mol, None).unwrap();
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].bonds.len(), mol.bond_count());
        let cancelled = Cancellation::new(|| true);
        assert!(matches!(
            delocalized_groups(&mol, Some(&cancelled)),
            Err(ResonancePerceptionError::Cancelled(_))
        ));
    }
}
