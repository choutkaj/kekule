use std::collections::{BTreeMap, VecDeque};

use super::rings::compute_ring_membership;
use crate::core::{Atom, AtomId, Bond, BondId, BondOrder, Molecule};

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

    // Ordinary color refinement cannot distinguish every inequivalent atom in
    // highly symmetric fused and bridged ring systems. RDKit handles the same
    // limitation with a second, ring-topology-aware partition pass. Use a
    // rooted breadth-first signature of the cyclic subgraph here: level widths
    // capture global ring distances and revisit multiplicities capture where
    // paths reconverge. These are graph invariants, not atom-order tie breakers,
    // so genuinely symmetry-equivalent atoms retain the same rank.
    if should_refine_ring_topology(mol, &atoms, &ranks) {
        let ring_topology = ring_topology_signatures(mol, &atoms);
        for _ in 0..atoms.len() {
            let refined = refine_ranks_with_ring_topology(mol, &atoms, &ranks, &ring_topology);
            if refined == ranks {
                break;
            }
            ranks = refined;
        }
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
            (
                *atom,
                atom_signature(mol, *atom, payload, degree(mol, *atom)),
            )
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
                .map(|(bond_id, bond)| {
                    let neighbor = bond.other_atom(*atom);
                    (bond_code(mol, bond_id, bond), ranks[neighbor.index()])
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
    total_hydrogens: usize,
    atom_map: u32,
    degree: usize,
}

fn atom_signature(mol: &Molecule, atom_id: AtomId, atom: &Atom, degree: usize) -> AtomSignature {
    AtomSignature {
        atomic_number: atom.element.atomic_number(),
        isotope: atom.isotope.unwrap_or(0),
        formal_charge: atom.formal_charge,
        total_hydrogens: usize::from(atom.explicit_hydrogens)
            + usize::from(mol.implicit_hydrogens(atom_id).ok().flatten().unwrap_or(0)),
        atom_map: atom.atom_map.unwrap_or(0),
        degree,
    }
}

fn degree(mol: &Molecule, atom: AtomId) -> usize {
    mol.incident_bonds(atom)
        .expect("atom_ids should only yield live atoms")
        .count()
}

fn bond_code(mol: &Molecule, bond_id: BondId, bond: &Bond) -> (u8, bool) {
    if mol.bond_is_aromatic(bond_id).ok().flatten() == Some(true) {
        return (5, true);
    }
    let order = match bond.order {
        BondOrder::Zero => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Aromatic => 5,
        BondOrder::Dative => 6,
    };
    (order, false)
}

fn refine_ranks_with_ring_topology(
    mol: &Molecule,
    atoms: &[AtomId],
    ranks: &[u32],
    topology: &[RingTopologySignature],
) -> Vec<u32> {
    let signatures = atoms
        .iter()
        .map(|atom| {
            let mut neighbors = mol
                .incident_bonds(*atom)
                .expect("atom_ids should only yield live atoms")
                .map(|(bond_id, bond)| {
                    let neighbor = bond.other_atom(*atom);
                    (bond_code(mol, bond_id, bond), ranks[neighbor.index()])
                })
                .collect::<Vec<_>>();
            neighbors.sort_unstable();
            (
                *atom,
                (
                    ranks[atom.index()],
                    topology[atom.index()].clone(),
                    neighbors,
                ),
            )
        })
        .collect::<Vec<_>>();
    compress_signatures(mol, signatures)
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct RingTopologySignature {
    level_widths: Vec<usize>,
    revisit_multiplicities: Vec<Vec<usize>>,
}

fn should_refine_ring_topology(mol: &Molecule, atoms: &[AtomId], ranks: &[u32]) -> bool {
    let (membership, _) = compute_ring_membership(mol);
    let ring_atoms = atoms
        .iter()
        .copied()
        .filter(|atom| membership.atom_in_ring(*atom))
        .collect::<Vec<_>>();
    if ring_atoms.is_empty() {
        return false;
    }

    let mut class_sizes = BTreeMap::<u32, usize>::new();
    for atom in atoms {
        *class_sizes.entry(ranks[atom.index()]).or_default() += 1;
    }
    let symmetric_ring_atoms = ring_atoms
        .iter()
        .filter(|atom| class_sizes[&ranks[atom.index()]] > 2)
        .count();
    if symmetric_ring_atoms * 2 <= ring_atoms.len() {
        return false;
    }

    let mut computed_ring_set = None;
    let ring_set = if let Some(ring_set) = mol.ring_set() {
        ring_set
    } else {
        let mut staged = mol.clone();
        let Ok(ring_set) = super::rings::perceive_ring_set(&mut staged) else {
            return false;
        };
        computed_ring_set.insert(ring_set)
    };
    let mut ring_counts = vec![0usize; mol.atoms.len()];
    for ring in ring_set.rings() {
        for atom in &ring.atoms {
            ring_counts[atom.index()] += 1;
        }
    }
    ring_atoms
        .iter()
        .any(|atom| ring_counts[atom.index()] > 1 && class_sizes[&ranks[atom.index()]] > 1)
}

fn ring_topology_signatures(mol: &Molecule, atoms: &[AtomId]) -> Vec<RingTopologySignature> {
    let (membership, _) = compute_ring_membership(mol);
    let mut signatures = vec![RingTopologySignature::default(); mol.atoms.len()];
    for &root in atoms {
        if membership.atom_in_ring(root) {
            signatures[root.index()] = ring_topology_signature(mol, root, &membership);
        }
    }
    signatures
}

fn ring_topology_signature(
    mol: &Molecule,
    root: AtomId,
    membership: &super::RingMembership,
) -> RingTopologySignature {
    let mut signature = RingTopologySignature::default();
    let mut visited = vec![false; mol.atoms.len()];
    let mut frontier = VecDeque::from([root]);
    visited[root.index()] = true;

    while !frontier.is_empty() {
        let current_level = frontier.drain(..).collect::<Vec<_>>();
        let mut current_flags = vec![false; mol.atoms.len()];
        for atom in &current_level {
            current_flags[atom.index()] = true;
        }

        let mut next_level = Vec::<AtomId>::new();
        for atom in current_level {
            // Substituents are counted in the preceding level but ring
            // topology traversal stops when it leaves the cyclic subgraph.
            if !membership.atom_in_ring(atom) {
                continue;
            }
            for (_, bond) in mol.incident_bonds(atom).expect("ring atom should be live") {
                let neighbor = bond.other_atom(atom);
                if !visited[neighbor.index()] {
                    visited[neighbor.index()] = true;
                    next_level.push(neighbor);
                }
            }
        }

        let mut next_flags = vec![false; mol.atoms.len()];
        for atom in &next_level {
            next_flags[atom.index()] = true;
        }
        let mut revisits = vec![0usize; mol.atoms.len()];
        for atom in &next_level {
            for (_, bond) in mol
                .incident_bonds(*atom)
                .expect("ring-neighbor atom should be live")
            {
                let neighbor = bond.other_atom(*atom);
                if current_flags[neighbor.index()] || next_flags[neighbor.index()] {
                    revisits[neighbor.index()] += 1;
                }
            }
        }
        let mut revisit_multiplicities = revisits
            .into_iter()
            .filter(|count| *count != 0)
            .collect::<Vec<_>>();
        revisit_multiplicities.sort_unstable();

        signature.level_widths.push(next_level.len());
        signature
            .revisit_multiplicities
            .push(revisit_multiplicities);
        frontier.extend(next_level);
    }

    signature
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::Element;

    #[test]
    fn atom_signature_preserves_large_graph_degree() {
        let mut mol = Molecule::new();
        let atom = mol.add_atom(Atom::new(Element::from_symbol("C").expect("carbon")));
        let degree = usize::from(u16::MAX) + 1;

        let signature = atom_signature(&mol, atom, mol.atom(atom).expect("atom"), degree);

        assert_eq!(signature.degree, degree);
    }

    #[test]
    fn ring_topology_refines_regular_graphs_without_breaking_true_ties() {
        let mut mol = Molecule::new();
        let atoms = (0..8)
            .map(|_| mol.add_atom(Atom::new(Element::from_symbol("C").expect("carbon"))))
            .collect::<Vec<_>>();

        // Two K4 graphs with one edge removed from each, joined across the
        // removed edges. Every vertex has degree three, so ordinary 1-WL
        // refinement leaves one class even though the graph has two distinct
        // vertex orbits.
        for (a, b) in [
            (0, 2),
            (0, 3),
            (1, 2),
            (1, 3),
            (2, 3),
            (4, 6),
            (4, 7),
            (5, 6),
            (5, 7),
            (6, 7),
            (0, 4),
            (1, 5),
        ] {
            mol.add_bond(atoms[a], atoms[b], BondOrder::Single)
                .expect("regular test graph bond");
        }

        let mut ordinary = initial_ranks(&mol, &atoms);
        for _ in 0..atoms.len() {
            let refined = refine_ranks(&mol, &atoms, &ordinary);
            if refined == ordinary {
                break;
            }
            ordinary = refined;
        }
        assert_eq!(
            ordinary
                .iter()
                .copied()
                .collect::<std::collections::BTreeSet<_>>()
                .len(),
            1
        );

        let ranking = canonical_atom_ranking(&mol);
        assert_eq!(ranking.rank_count(), 2);
        let mut class_sizes = BTreeMap::<u32, usize>::new();
        for (_, rank) in ranking.iter() {
            *class_sizes.entry(rank).or_default() += 1;
        }
        assert_eq!(class_sizes.into_values().collect::<Vec<_>>(), vec![4, 4]);
    }
}
