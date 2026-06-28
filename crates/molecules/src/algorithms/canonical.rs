use std::collections::BTreeMap;

use crate::{Atom, AtomId, AtomRadical, Bond, BondOrder, Molecule};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalAtomRanking {
    ranks: Vec<(AtomId, u32)>,
}

impl CanonicalAtomRanking {
    pub fn iter(&self) -> impl Iterator<Item = (AtomId, u32)> + '_ {
        self.ranks.iter().copied()
    }

    pub fn rank_of(&self, atom: AtomId) -> Option<u32> {
        self.ranks
            .iter()
            .find_map(|(id, rank)| (*id == atom).then_some(*rank))
    }

    pub fn rank_count(&self) -> usize {
        self.ranks
            .iter()
            .map(|(_, rank)| *rank)
            .collect::<std::collections::BTreeSet<_>>()
            .len()
    }
}

pub fn canonical_atom_ranking(mol: &Molecule) -> CanonicalAtomRanking {
    let atoms = mol.atom_ids().collect::<Vec<_>>();
    let mut ranks = initial_ranks(mol, &atoms);

    for _ in 0..atoms.len() {
        let refined = refine_ranks(mol, &atoms, &ranks);
        if refined == ranks {
            break;
        }
        ranks = refined;
    }

    CanonicalAtomRanking {
        ranks: atoms
            .into_iter()
            .map(|atom| (atom, ranks[atom.index()]))
            .collect(),
    }
}

fn initial_ranks(mol: &Molecule, atoms: &[AtomId]) -> Vec<u32> {
    let signatures = atoms
        .iter()
        .map(|atom| {
            let payload = mol
                .atom(*atom)
                .expect("atom_ids should only yield live atoms");
            (*atom, atom_signature(payload, degree(mol, *atom)))
        })
        .collect::<Vec<_>>();
    compress_signatures(mol, signatures)
}

fn refine_ranks(mol: &Molecule, atoms: &[AtomId], ranks: &[u32]) -> Vec<u32> {
    let signatures = atoms
        .iter()
        .map(|atom| {
            let mut neighbors = mol
                .incident_bonds(*atom)
                .expect("atom_ids should only yield live atoms")
                .map(|(_, bond)| {
                    let neighbor = bond.other_atom(*atom);
                    (bond_code(bond), ranks[neighbor.index()])
                })
                .collect::<Vec<_>>();
            neighbors.sort_unstable();
            (*atom, (ranks[atom.index()], neighbors))
        })
        .collect::<Vec<_>>();
    compress_signatures(mol, signatures)
}

fn compress_signatures<T: Clone + Ord>(mol: &Molecule, signatures: Vec<(AtomId, T)>) -> Vec<u32> {
    let mut ordered = signatures
        .iter()
        .map(|(_, signature)| signature.clone())
        .collect::<Vec<_>>();
    ordered.sort_unstable();
    ordered.dedup();

    let rank_by_signature = ordered
        .into_iter()
        .enumerate()
        .map(|(rank, signature)| (signature, rank as u32))
        .collect::<BTreeMap<_, _>>();

    let mut ranks = vec![u32::MAX; mol.atoms.len()];
    for (atom, signature) in signatures {
        ranks[atom.index()] = rank_by_signature[&signature];
    }
    ranks
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct AtomSignature {
    atomic_number: u8,
    isotope: u16,
    formal_charge: i8,
    radical: u8,
    explicit_hydrogens: u8,
    implicit_hydrogens: u8,
    no_implicit_hydrogens: bool,
    aromatic: bool,
    atom_map: u32,
    degree: u16,
}

fn atom_signature(atom: &Atom, degree: usize) -> AtomSignature {
    AtomSignature {
        atomic_number: atom.element.atomic_number(),
        isotope: atom.isotope.unwrap_or(0),
        formal_charge: atom.formal_charge,
        radical: radical_code(atom.radical),
        explicit_hydrogens: atom.explicit_hydrogens,
        implicit_hydrogens: atom.implicit_hydrogens.unwrap_or(u8::MAX),
        no_implicit_hydrogens: atom.no_implicit_hydrogens,
        aromatic: atom.aromatic,
        atom_map: atom.atom_map.unwrap_or(0),
        degree: degree as u16,
    }
}

fn degree(mol: &Molecule, atom: AtomId) -> usize {
    mol.incident_bonds(atom)
        .expect("atom_ids should only yield live atoms")
        .count()
}

fn radical_code(radical: Option<AtomRadical>) -> u8 {
    match radical {
        None => 0,
        Some(AtomRadical::Singlet) => 1,
        Some(AtomRadical::Doublet) => 2,
        Some(AtomRadical::Triplet) => 3,
    }
}

fn bond_code(bond: &Bond) -> (u8, bool) {
    let order = match bond.order {
        BondOrder::Zero => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::Dative => 6,
    };
    (order, bond.aromatic)
}
