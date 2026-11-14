use super::*;
use crate::molecule::Element;

#[test]
fn atom_signature_preserves_large_graph_degree() {
    let mut mol = crate::molecule::MoleculeEditor::new();
    let atom = mol
        .add_atom(Atom::new(Element::from_symbol("C").expect("carbon")))
        .expect("atom identifier capacity");
    let degree = usize::from(u16::MAX) + 1;

    let mol = mol.into_unchecked();
    let signature = atom_signature(&mol, atom, mol.atom(atom).expect("atom"), degree);

    assert_eq!(signature.degree, degree);
}

#[test]
fn ring_topology_refines_regular_graphs_without_breaking_true_ties() {
    let mut mol = crate::molecule::MoleculeEditor::new();
    let atoms = (0..8)
        .map(|_| {
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("carbon")))
                .expect("atom identifier capacity")
        })
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

    let mol = mol.into_unchecked();
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
