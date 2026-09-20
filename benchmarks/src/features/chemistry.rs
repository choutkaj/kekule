use kekule::{
    core::{Atom, AtomId, AtomRadical, BondOrder, Molecule},
    perception::valence,
};
use serde_json::{json, Value};

pub(crate) fn atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "atom_map": atom.atom_map,
        "spin_multiplicity": atom.radical.and_then(AtomRadical::spin_multiplicity),
        "radical_electrons": atom.radical.map(AtomRadical::electron_count).unwrap_or(0),
        "aromatic": mol.atom_is_aromatic(id).expect("live atom"),
    })
}

pub(crate) fn basic_atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| basic_atom_json(mol, id, atom))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "atom_map": atom.atom_map,
        "aromatic": mol.atom_is_aromatic(id).expect("live atom"),
    })
}

pub(crate) fn valence_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "implicit_hydrogens": mol.implicit_hydrogens(id).expect("live atom"),
        "explicit_valence": valence::represented_valence(mol, id).expect("live atom"),
    })
}

pub(crate) fn bond_order_json(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Zero => "ZERO",
        BondOrder::Single => "SINGLE",
        BondOrder::Double => "DOUBLE",
        BondOrder::Triple => "TRIPLE",
        BondOrder::Quadruple => "QUADRUPLE",
        BondOrder::Dative => "DATIVE",
    }
}
