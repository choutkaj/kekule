use super::*;

pub(super) fn carbon() -> Atom {
    Atom::new(Element::from_symbol("C").expect("carbon should be available"))
}

pub(super) fn oxygen() -> Atom {
    Atom::new(Element::from_symbol("O").expect("oxygen should be available"))
}

pub(super) fn ring_molecule(
    symbols: &[&str],
    orders: &[BondOrder],
) -> (Molecule, Vec<AtomId>, Vec<BondId>) {
    assert_eq!(symbols.len(), orders.len());
    let mut mol = Molecule::new();
    let atoms = symbols
        .iter()
        .map(|symbol| {
            mol.add_atom(Atom::new(
                Element::from_symbol(symbol).expect("test element should be available"),
            ))
        })
        .collect::<Vec<_>>();
    let mut bonds = Vec::new();
    for index in 0..atoms.len() {
        let next = (index + 1) % atoms.len();
        bonds.push(
            mol.add_bond(atoms[index], atoms[next], orders[index])
                .expect("ring bond should be valid"),
        );
    }
    (mol, atoms, bonds)
}

pub(super) fn sorted_atom_ids(ids: impl IntoIterator<Item = AtomId>) -> Vec<AtomId> {
    let mut ids = ids.into_iter().collect::<Vec<_>>();
    ids.sort();
    ids
}

pub(super) fn sorted_bond_ids(ids: impl IntoIterator<Item = BondId>) -> Vec<BondId> {
    let mut ids = ids.into_iter().collect::<Vec<_>>();
    ids.sort();
    ids
}

pub(super) fn deterministic_text_mutations(seed: &str) -> Vec<String> {
    let mut mutations = vec![String::new(), seed.to_owned()];
    for index in 0..=seed.len().min(128) {
        for inserted in ["\0", "\n", "%", "[", "]", "é"] {
            let mut value = seed.to_owned();
            value.insert_str(index, inserted);
            mutations.push(value);
        }
        if index < seed.len() {
            let mut removed = seed.to_owned();
            removed.remove(index);
            mutations.push(removed);

            let mut replaced = seed.to_owned();
            replaced.replace_range(index..index + 1, "\u{7f}");
            mutations.push(replaced);
        }
    }
    mutations
}

pub(super) fn mark_all_fresh(mol: &mut Molecule) {
    mol.perception.valence = ComputedState::Fresh;
    mol.perception.rings = ComputedState::Fresh;
    mol.perception.aromaticity = ComputedState::Fresh;
    mol.perception.stereo = ComputedState::Fresh;
}

pub(super) fn assert_all_stale(mol: &Molecule) {
    assert_eq!(mol.perception().valence, ComputedState::Stale);
    assert_eq!(mol.perception().rings, ComputedState::Stale);
    assert_eq!(mol.perception().aromaticity, ComputedState::Stale);
    assert_eq!(mol.perception().stereo, ComputedState::Stale);
}

mod bio;
mod chemistry;
mod core_payload;
mod graph;
mod perception;
mod ring_limits;
mod smiles;
mod v2000;
