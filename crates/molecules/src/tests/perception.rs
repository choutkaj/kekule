use super::*;

#[test]
fn ring_membership_empty_and_linear_molecules_have_no_rings() {
    let mut empty = Molecule::new();
    let empty_membership = perceive_ring_membership(&mut empty);
    assert!(empty_membership.ring_atom_ids().next().is_none());
    assert!(empty_membership.ring_bond_ids().next().is_none());

    let mut chain = Molecule::new();
    let a = chain.add_atom(carbon());
    let b = chain.add_atom(carbon());
    let c = chain.add_atom(carbon());
    let ab = chain
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    let bc = chain
        .add_bond(b, c, BondOrder::Single)
        .expect("bond should be valid");
    let chain_membership = perceive_ring_membership(&mut chain);

    assert!(!chain_membership.atom_in_ring(a));
    assert!(!chain_membership.atom_in_ring(b));
    assert!(!chain_membership.bond_in_ring(ab));
    assert!(!chain_membership.bond_in_ring(bc));
    assert_eq!(chain.perception().rings, ComputedState::Fresh);
}

#[test]
fn ring_membership_marks_triangle_atoms_and_bonds() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(carbon());
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");

    let membership = perceive_ring_membership(&mut mol);

    assert_eq!(sorted_atom_ids(membership.ring_atom_ids()), vec![a, b, c]);
    assert_eq!(
        sorted_bond_ids(membership.ring_bond_ids()),
        vec![ab, bc, ca]
    );
}

#[test]
fn ring_membership_excludes_tail_from_ring() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(carbon());
    let tail = mol.add_atom(oxygen());
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
    let tail_bond = mol.add_bond(c, tail, BondOrder::Single).expect("bond");

    let membership = perceive_ring_membership(&mut mol);

    assert_eq!(sorted_atom_ids(membership.ring_atom_ids()), vec![a, b, c]);
    assert_eq!(
        sorted_bond_ids(membership.ring_bond_ids()),
        vec![ab, bc, ca]
    );
    assert!(!membership.atom_in_ring(tail));
    assert!(!membership.bond_in_ring(tail_bond));
}

#[test]
fn ring_membership_handles_fused_and_disconnected_components() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(carbon());
    let d = mol.add_atom(carbon());
    let isolated_a = mol.add_atom(oxygen());
    let isolated_b = mol.add_atom(oxygen());
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
    let cd = mol.add_bond(c, d, BondOrder::Single).expect("bond");
    let da = mol.add_bond(d, a, BondOrder::Single).expect("bond");
    let bridge = mol
        .add_bond(isolated_a, isolated_b, BondOrder::Single)
        .expect("bond");

    let membership = perceive_ring_membership(&mut mol);

    assert_eq!(
        sorted_atom_ids(membership.ring_atom_ids()),
        vec![a, b, c, d]
    );
    assert_eq!(
        sorted_bond_ids(membership.ring_bond_ids()),
        vec![ab, bc, ca, cd, da]
    );
    assert!(!membership.bond_in_ring(bridge));
}

#[test]
fn ring_membership_ignores_deleted_bonds_and_becomes_stale_after_mutation() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(carbon());
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
    mol.delete_bond(ca).expect("bond should delete");

    let membership = perceive_ring_membership(&mut mol);
    assert!(!membership.bond_in_ring(ab));
    assert!(!membership.bond_in_ring(bc));
    assert!(!membership.bond_in_ring(ca));

    mol.add_bond(c, a, BondOrder::Single).expect("bond");
    assert_eq!(mol.perception().rings, ComputedState::Stale);
    assert!(mol.ring_membership().is_none());
    assert!(mol.ring_set().is_none());
}

#[test]
fn aromaticity_marks_benzene_like_ring() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("benzene should be supported");

    assert_eq!(mol.perception().aromaticity, ComputedState::Fresh);
    assert!(atoms
        .iter()
        .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_evaluates_larger_simple_rings_like_rdkit() {
    let alternating_ten = [
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
    ];
    let (mut ten_member, ten_atoms, ten_bonds) = ring_molecule(&["C"; 10], &alternating_ten);

    perceive_aromaticity(&mut ten_member, AromaticityModel::RdkitLike)
        .expect("10 pi-electron annulene-like ring should be supported");

    assert!(ten_atoms
        .iter()
        .all(|atom| ten_member.atom(*atom).expect("atom exists").aromatic));
    assert!(ten_bonds
        .iter()
        .all(|bond| ten_member.bond(*bond).expect("bond exists").aromatic));

    let alternating_twelve = [
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
    ];
    let (mut twelve_member, twelve_atoms, twelve_bonds) =
        ring_molecule(&["C"; 12], &alternating_twelve);

    perceive_aromaticity(&mut twelve_member, AromaticityModel::RdkitLike)
        .expect("12 pi-electron annulene-like ring should be supported");

    assert!(twelve_atoms
        .iter()
        .all(|atom| !twelve_member.atom(*atom).expect("atom exists").aromatic));
    assert!(twelve_bonds
        .iter()
        .all(|bond| !twelve_member.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_leaves_cyclohexane_and_cyclobutadiene_non_aromatic() {
    let (mut cyclohexane, atoms, bonds) =
        ring_molecule(&["C", "C", "C", "C", "C", "C"], &[BondOrder::Single; 6]);
    perceive_aromaticity(&mut cyclohexane, AromaticityModel::RdkitLike)
        .expect("cyclohexane should be supported");
    assert!(atoms
        .iter()
        .all(|atom| !cyclohexane.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !cyclohexane.bond(*bond).expect("bond exists").aromatic));

    let (mut cyclobutadiene, atoms, bonds) = ring_molecule(
        &["C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    perceive_aromaticity(&mut cyclobutadiene, AromaticityModel::RdkitLike)
        .expect("cyclobutadiene should be supported");
    assert!(atoms
        .iter()
        .all(|atom| !cyclobutadiene.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !cyclobutadiene.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_supports_heteroaromatic_ring() {
    let (mut furan_like, atoms, bonds) = ring_molecule(
        &["O", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );

    perceive_aromaticity(&mut furan_like, AromaticityModel::RdkitLike)
        .expect("furan-like ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| furan_like.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| furan_like.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_supports_explicit_nitrogen_lone_pair_donor_ring() {
    let (mut pyrrole_like, atoms, bonds) = ring_molecule(
        &["N", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    {
        let mut nitrogen = pyrrole_like
            .atom_mut(atoms[0])
            .expect("ring nitrogen should exist");
        nitrogen.explicit_hydrogens = 1;
        nitrogen.implicit_hydrogens = Some(0);
        nitrogen.no_implicit_hydrogens = true;
    }

    perceive_aromaticity(&mut pyrrole_like, AromaticityModel::RdkitLike)
        .expect("pyrrole-like ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| pyrrole_like.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| pyrrole_like.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_supports_phosphorus_lone_pair_donor_ring() {
    let (mut phosphole_like, atoms, bonds) = ring_molecule(
        &["P", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );

    perceive_aromaticity(&mut phosphole_like, AromaticityModel::RdkitLike)
        .expect("phosphole-like ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| phosphole_like.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| phosphole_like.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_rejects_ring_atom_above_rdkit_default_valence() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["P", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    let methyl = mol.add_atom(carbon());
    mol.add_bond(atoms[0], methyl, BondOrder::Single)
        .expect("phosphorus substituent bond");

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("hypervalent phosphorus ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
    assert!(!mol.atom(methyl).expect("substituent exists").aromatic);
}

#[test]
fn aromaticity_applies_rdkit_radical_candidate_rules() {
    let (mut neutral_carbon_radical, atoms, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    neutral_carbon_radical
        .atom_mut(atoms[0])
        .expect("ring atom exists")
        .radical = Some(AtomRadical::Doublet);

    perceive_aromaticity(&mut neutral_carbon_radical, AromaticityModel::RdkitLike)
        .expect("neutral carbon radical ring should be supported");

    assert!(atoms.iter().all(|atom| neutral_carbon_radical
        .atom(*atom)
        .expect("atom exists")
        .aromatic));

    let (mut oxygen_radical, atoms, _) = ring_molecule(
        &["O", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    oxygen_radical
        .atom_mut(atoms[0])
        .expect("ring atom exists")
        .radical = Some(AtomRadical::Doublet);

    perceive_aromaticity(&mut oxygen_radical, AromaticityModel::RdkitLike)
        .expect("heteroatom radical ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !oxygen_radical.atom(*atom).expect("atom exists").aromatic));

    let (mut charged_carbon_radical, atoms, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    {
        let mut atom = charged_carbon_radical
            .atom_mut(atoms[0])
            .expect("ring atom exists");
        atom.formal_charge = 1;
        atom.radical = Some(AtomRadical::Doublet);
    }

    perceive_aromaticity(&mut charged_carbon_radical, AromaticityModel::RdkitLike)
        .expect("charged carbon radical ring should be supported");

    assert!(atoms.iter().all(|atom| !charged_carbon_radical
        .atom(*atom)
        .expect("atom exists")
        .aromatic));
}

#[test]
fn aromaticity_rejects_tetracoordinate_ring_atom_candidate() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["N", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    mol.atom_mut(atoms[0])
        .expect("ring atom exists")
        .formal_charge = 1;
    let methyl_a = mol.add_atom(carbon());
    let methyl_b = mol.add_atom(carbon());
    mol.add_bond(atoms[0], methyl_a, BondOrder::Single)
        .expect("first substituent bond");
    mol.add_bond(atoms[0], methyl_b, BondOrder::Single)
        .expect("second substituent bond");

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("tetracoordinate ring atom should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_rejects_protonated_saturated_ring_nitrogen_donor() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["N", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    {
        let mut nitrogen = mol.atom_mut(atoms[0]).expect("ring atom exists");
        nitrogen.formal_charge = 1;
        nitrogen.explicit_hydrogens = 1;
        nitrogen.implicit_hydrogens = Some(0);
        nitrogen.no_implicit_hydrogens = true;
    }

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("protonated saturated ring nitrogen should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_accepts_cyclopropenyl_cation_two_electron_ring() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["C", "C", "C"],
        &[BondOrder::Single, BondOrder::Double, BondOrder::Single],
    );
    {
        let mut cation = mol.atom_mut(atoms[0]).expect("ring atom exists");
        cation.formal_charge = 1;
        cation.explicit_hydrogens = 1;
        cation.implicit_hydrogens = Some(0);
        cation.no_implicit_hydrogens = true;
    }

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclopropenyl cation should be supported");

    assert!(atoms
        .iter()
        .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_requires_every_atom_to_be_candidate_before_huckel_count() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    {
        let mut saturated = mol.atom_mut(atoms[0]).expect("ring atom exists");
        saturated.explicit_hydrogens = 2;
        saturated.implicit_hydrogens = Some(0);
        saturated.no_implicit_hydrogens = true;
    }

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("over-valent candidate rejection should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_marks_azulene_fused_perimeter_but_not_shared_bond() {
    let mut mol = Molecule::new();
    let atoms = (0..10).map(|_| mol.add_atom(carbon())).collect::<Vec<_>>();
    let orders = [
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Single,
        BondOrder::Double,
    ];
    let mut perimeter_bonds = Vec::new();
    for index in 0..7 {
        perimeter_bonds.push(
            mol.add_bond(atoms[index], atoms[index + 1], orders[index])
                .expect("perimeter bond"),
        );
    }
    let shared = mol
        .add_bond(atoms[7], atoms[3], BondOrder::Single)
        .expect("fused shared bond");
    perimeter_bonds.push(
        mol.add_bond(atoms[7], atoms[8], BondOrder::Single)
            .expect("perimeter bond"),
    );
    perimeter_bonds.push(
        mol.add_bond(atoms[8], atoms[9], BondOrder::Double)
            .expect("perimeter bond"),
    );
    perimeter_bonds.push(
        mol.add_bond(atoms[9], atoms[0], BondOrder::Single)
            .expect("perimeter bond"),
    );

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("azulene-like fused system should be supported");

    assert!(atoms
        .iter()
        .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
    assert!(perimeter_bonds
        .iter()
        .all(|bond| mol.bond(*bond).expect("bond exists").aromatic));
    assert!(!mol.bond(shared).expect("shared bond exists").aromatic);
}

#[test]
fn aromaticity_preserves_anionic_carbon_donor_with_explicit_hydrogen_bond() {
    let (mut mol, atoms, _) = ring_molecule(
        &["C", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    for atom_id in &atoms {
        mol.atom_mut(*atom_id)
            .expect("ring atom exists")
            .formal_charge = -1;
    }
    let hydrogen = mol.add_atom(Atom::new(
        Element::from_symbol("H").expect("hydrogen should be available"),
    ));
    mol.add_bond(atoms[0], hydrogen, BondOrder::Single)
        .expect("explicit hydrogen bond should be valid");

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclopentadienyl anion should be supported");

    assert!(atoms
        .iter()
        .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
    assert!(!mol.atom(hydrogen).expect("hydrogen exists").aromatic);
}

#[test]
fn aromaticity_rejects_neutral_saturated_carbon_in_conjugated_ring() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["C", "C", "C", "C", "C"],
        &[
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclopentadiene should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_uses_ring_membership_not_acyclic_double_bonds() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(carbon());
    mol.add_bond(a, b, BondOrder::Double).expect("bond");
    mol.add_bond(b, c, BondOrder::Single).expect("bond");

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("acyclic molecule should be supported");

    assert!(!mol.atom(a).expect("atom exists").aromatic);
    assert!(!mol.bond(BondId::new(0)).expect("bond exists").aromatic);
}

#[test]
fn aromaticity_clears_existing_flags_before_assignment() {
    let (mut mol, atoms, bonds) =
        ring_molecule(&["C", "C", "C", "C", "C", "C"], &[BondOrder::Single; 6]);
    for atom in &atoms {
        mol.atom_mut(*atom).expect("atom exists").aromatic = true;
    }
    for bond in &bonds {
        mol.bond_mut(*bond).expect("bond exists").aromatic = true;
    }

    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclohexane should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !mol.atom(*atom).expect("atom exists").aromatic));
    assert!(bonds
        .iter()
        .all(|bond| !mol.bond(*bond).expect("bond exists").aromatic));
}

#[test]
fn aromaticity_becomes_stale_after_topology_mutation() {
    let (mut mol, atoms, _) = ring_molecule(
        &["C", "C", "C", "C", "C", "C"],
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
        ],
    );
    perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("benzene should be supported");

    mol.add_atom(oxygen());
    assert_eq!(mol.perception().aromaticity, ComputedState::Stale);
    assert!(atoms
        .iter()
        .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
}
