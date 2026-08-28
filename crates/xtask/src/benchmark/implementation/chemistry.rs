use crate::*;

pub(crate) fn atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| atom_json(mol, id, atom))
        .collect::<Vec<_>>()
}

pub(crate) fn atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "atom_map": atom.atom_map,
        "radical": atom.radical.map(radical_json),
        "unpaired_electrons": atom.radical.map(AtomRadical::unpaired_electron_count).unwrap_or(0),
        "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
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
        "aromatic": mol.atom_is_aromatic(id).ok().flatten().unwrap_or(false),
    })
}

pub(crate) fn valence_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "explicit_hydrogens": atom.hydrogens.explicit_count(),
        "implicit_hydrogens": mol.implicit_hydrogens(id).ok().flatten().unwrap_or(0),
        "explicit_valence": explicit_valence_json(mol, id) + atom.hydrogens.explicit_count(),
    })
}

pub(crate) fn explicit_valence_json(mol: &Molecule, atom: AtomId) -> u8 {
    let atom_record = mol.atom(atom).ok();
    let bonds = mol
        .incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let has_non_aromatic_bond = bonds
        .iter()
        .any(|(id, _)| mol.bond_is_aromatic(*id).ok().flatten() != Some(true));
    let has_non_aromatic_multiple_bond = bonds.iter().any(|(id, bond)| {
        mol.bond_is_aromatic(*id).ok().flatten() != Some(true)
            && matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            )
    });
    let has_marked_aromatic_high_order_bond = bonds.iter().any(|(id, bond)| {
        mol.bond_is_aromatic(*id).ok().flatten() == Some(true)
            && matches!(bond.order, BondOrder::Triple | BondOrder::Quadruple)
    });
    let aromatic_bond_count = bonds
        .iter()
        .filter(|(id, _)| mol.bond_is_aromatic(*id).ok().flatten() == Some(true))
        .count();
    // The RDKit semantic record treats a pyrrolic donor H as explicit after
    // RDKit's prepared state. Kekule keeps an inferred H in Perception, so derive
    // the comparable bond-valence contribution without rewriting the atom.
    let has_aromatic_nitrogen_hydrogen = atom_record.is_some_and(|atom_record| {
        atom_record.element.symbol() == "N"
            && atom_record.formal_charge == 0
            && mol.atom_is_aromatic(atom).ok().flatten() == Some(true)
            && (atom_record.hydrogens.explicit_count() > 0
                || mol.implicit_hydrogens(atom).ok().flatten() == Some(1))
    });
    let doubled: u8 = bonds
        .into_iter()
        .map(|(id, bond)| {
            if mol.bond_is_aromatic(id).ok().flatten() == Some(true) {
                if has_marked_aromatic_high_order_bond {
                    return match bond.order {
                        BondOrder::Triple => 6,
                        BondOrder::Quadruple => 8,
                        _ => 2,
                    };
                }
                return aromatic_bond_valence_twice(
                    atom_record,
                    mol.atom_is_aromatic(atom).ok().flatten() == Some(true),
                    has_non_aromatic_bond,
                    has_non_aromatic_multiple_bond,
                    aromatic_bond_count,
                    has_aromatic_nitrogen_hydrogen,
                );
            }
            match bond.order {
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Single => 2,
                BondOrder::Double => 4,
                BondOrder::Triple => 6,
                BondOrder::Quadruple => 8,
            }
        })
        .sum();
    doubled / 2
}

fn aromatic_bond_valence_twice(
    atom: Option<&Atom>,
    atom_aromatic: bool,
    has_non_aromatic_bond: bool,
    has_non_aromatic_multiple_bond: bool,
    aromatic_bond_count: usize,
    has_aromatic_nitrogen_hydrogen: bool,
) -> u8 {
    let Some(atom) = atom else {
        return 2;
    };
    if atom_aromatic && has_non_aromatic_multiple_bond {
        return 2;
    }
    match atom.element.symbol() {
        "C" if atom.formal_charge < 0
            && (atom.hydrogens.explicit_count() > 0
                || has_non_aromatic_bond
                || aromatic_bond_count >= 3) =>
        {
            2
        }
        "P" | "As" | "Sb"
            if atom.formal_charge == 0
                && atom.hydrogens.explicit_count() == 0
                && (has_non_aromatic_bond || aromatic_bond_count >= 3) =>
        {
            2
        }
        "O" | "S" | "Se" | "Te"
            if atom.formal_charge == 0 && atom.hydrogens.explicit_count() == 0 =>
        {
            2
        }
        "N" if atom.formal_charge < 0 => 2,
        "N" if atom.formal_charge == 0 && has_aromatic_nitrogen_hydrogen => 2,
        "N" if atom.formal_charge == 0 && has_non_aromatic_bond => 2,
        "N" if atom.formal_charge == 0 && aromatic_bond_count >= 3 => 2,
        _ => 3,
    }
}

pub(crate) fn bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| bond_json(mol, id, bond))
        .collect::<Vec<_>>()
}

pub(crate) fn bond_json(mol: &Molecule, id: BondId, bond: &Bond) -> Value {
    json!({
        "index": id.raw(),
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false),
        "stereo": "STEREONONE",
        "bond_direction": "NONE",
    })
}

pub(crate) fn basic_bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| basic_bond_json(mol, id, bond))
        .collect::<Vec<_>>()
}

pub(crate) fn basic_bond_json(mol: &Molecule, id: BondId, bond: &Bond) -> Value {
    json!({
        "index": id.raw(),
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": mol.bond_is_aromatic(id).ok().flatten().unwrap_or(false),
        "stereo": "STEREONONE",
    })
}

pub(crate) fn radical_json(radical: AtomRadical) -> &'static str {
    match radical {
        AtomRadical::Singlet => "SINGLET",
        AtomRadical::Doublet => "DOUBLET",
        AtomRadical::Triplet => "TRIPLET",
        AtomRadical::Quartet => "QUARTET",
        AtomRadical::Quintet => "QUINTET",
    }
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
