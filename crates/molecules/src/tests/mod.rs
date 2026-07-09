use crate::bio::*;
use crate::core::*;
use crate::perception::{
    aromaticity as aromaticity_api, aromaticity::*, rings as rings_api, rings::*,
    stereo as stereo_api, stereo::*, valence as valence_api, valence::*,
};
use crate::sdf::*;
use crate::small::*;
use crate::smiles::*;
use crate::{
    bio as bio_api, canon, molfile, perception as perception_api, sdf, smiles as smiles_api,
};

pub(super) fn carbon() -> Atom {
    Atom::new(Element::from_symbol("C").expect("carbon should be available"))
}

pub(super) fn oxygen() -> Atom {
    Atom::new(Element::from_symbol("O").expect("oxygen should be available"))
}

pub(super) fn element_atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).expect("test element should be available"))
}

pub(super) fn charged_atom(symbol: &str, formal_charge: i8) -> Atom {
    let mut atom = element_atom(symbol);
    atom.formal_charge = formal_charge;
    atom
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

pub(super) fn implicit_h_wedge_geometry_molblock() -> &'static str {
    "\
implicit H geometry wedge
molecules

  4  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000   -1.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  6  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  END
"
}

pub(super) fn rdkit_rp6306_atrop_molblock() -> &'static str {
    include_str!("../../../../validation/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop1.mol")
}

pub(super) fn rdkit_rp6306_atrop3_molblock() -> &'static str {
    include_str!("../../../../validation/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop3.mol")
}

pub(super) fn rdkit_rp6306_atrop4_molblock() -> &'static str {
    include_str!("../../../../validation/corpora/smoke/data/rdkit_atropisomers/RP-6306_atrop4.mol")
}

pub(super) fn rdkit_bms986142_atrop4_molblock() -> &'static str {
    include_str!(
        "../../../../validation/corpora/smoke/data/rdkit_atropisomers/BMS-986142_atrop4.mol"
    )
}

pub(super) fn rdkit_bms986142_atrop5_molblock() -> &'static str {
    include_str!(
        "../../../../validation/corpora/smoke/data/rdkit_atropisomers/BMS-986142_atrop5.mol"
    )
}

pub(super) fn rdkit_jdq443_atrop1_molblock() -> &'static str {
    include_str!("../../../../validation/corpora/smoke/data/rdkit_atropisomers/JDQ443_atrop1.mol")
}

pub(super) fn rdkit_zm374979_atrop1_molblock() -> &'static str {
    include_str!("../../../../validation/corpora/smoke/data/rdkit_atropisomers/ZM374979_atrop1.mol")
}

pub(super) fn rdkit_zm374979_atrop2_molblock() -> &'static str {
    include_str!("../../../../validation/corpora/smoke/data/rdkit_atropisomers/ZM374979_atrop2.mol")
}

mod bio;
mod canonical;
mod chemistry;
mod cip;
mod core_payload;
mod graph;
mod perception;
mod public_api;
mod ring_limits;
mod smiles;
mod v2000;
mod v3000;
