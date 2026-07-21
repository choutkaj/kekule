use super::*;

#[test]
fn ring_membership_empty_and_linear_molecules_have_no_rings() {
    let mut empty = Molecule::new();
    let empty_membership = rings_api::perceive_ring_membership(&mut empty);
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
    let chain_membership = rings_api::perceive_ring_membership(&mut chain);

    assert!(!chain_membership.atom_in_ring(a));
    assert!(!chain_membership.atom_in_ring(b));
    assert!(!chain_membership.bond_in_ring(ab));
    assert!(!chain_membership.bond_in_ring(bc));
    assert!(chain.perception().has_rings());
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

    let membership = rings_api::perceive_ring_membership(&mut mol);

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

    let membership = rings_api::perceive_ring_membership(&mut mol);

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

    let membership = rings_api::perceive_ring_membership(&mut mol);

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

    let membership = rings_api::perceive_ring_membership(&mut mol);
    assert!(!membership.bond_in_ring(ab));
    assert!(!membership.bond_in_ring(bc));
    assert!(!membership.bond_in_ring(ca));

    mol.add_bond(c, a, BondOrder::Single).expect("bond");
    assert!(!mol.perception().has_rings());
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("benzene should be supported");

    assert!(mol.perception().has_aromaticity());
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

    aromaticity_api::perceive_aromaticity(&mut ten_member, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut twelve_member, AromaticityModel::RdkitLike)
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
    aromaticity_api::perceive_aromaticity(&mut cyclohexane, AromaticityModel::RdkitLike)
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
    aromaticity_api::perceive_aromaticity(&mut cyclobutadiene, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut furan_like, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut pyrrole_like, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut phosphole_like, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut neutral_carbon_radical, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut oxygen_radical, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut charged_carbon_radical, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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
fn aromaticity_keeps_aromatic_heteroring_bond_shared_with_saturated_ring() {
    let mut mol = Molecule::new();
    let c0 = mol.add_atom(carbon());
    let c1 = mol.add_atom(carbon());
    let c2 = mol.add_atom(carbon());
    let n3 = mol.add_atom(Atom::new(
        Element::from_symbol("N").expect("nitrogen should be available"),
    ));
    let n4 = mol.add_atom(Atom::new(
        Element::from_symbol("N").expect("nitrogen should be available"),
    ));
    let saturated_a = mol.add_atom(carbon());
    let saturated_b = mol.add_atom(carbon());
    let saturated_c = mol.add_atom(carbon());

    let aromatic_bonds = [
        mol.add_bond(c0, c1, BondOrder::Double)
            .expect("aromatic ring bond"),
        mol.add_bond(c1, c2, BondOrder::Single)
            .expect("shared fused bond"),
        mol.add_bond(c2, n3, BondOrder::Double)
            .expect("aromatic ring bond"),
        mol.add_bond(n3, n4, BondOrder::Single)
            .expect("aromatic ring bond"),
        mol.add_bond(n4, c0, BondOrder::Single)
            .expect("aromatic ring bond"),
    ];
    let saturated_bonds = [
        mol.add_bond(c1, saturated_a, BondOrder::Single)
            .expect("saturated ring bond"),
        mol.add_bond(saturated_a, saturated_b, BondOrder::Single)
            .expect("saturated ring bond"),
        mol.add_bond(saturated_b, saturated_c, BondOrder::Single)
            .expect("saturated ring bond"),
        mol.add_bond(saturated_c, c2, BondOrder::Single)
            .expect("saturated ring bond"),
    ];

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("fused heteroaromatic ring should be supported");

    for bond_id in aromatic_bonds {
        assert!(
            mol.bond(bond_id).expect("aromatic bond exists").aromatic,
            "aromatic ring bond {bond_id} should be aromatic"
        );
    }
    for bond_id in saturated_bonds {
        assert!(
            !mol.bond(bond_id).expect("saturated bond exists").aromatic,
            "saturated fused-neighbor bond {bond_id} should stay aliphatic"
        );
    }
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
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
    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("benzene should be supported");

    mol.add_atom(oxygen());
    assert!(!mol.perception().has_aromaticity());
    assert!(atoms
        .iter()
        .all(|atom| mol.atom(*atom).expect("atom exists").aromatic));
}

#[test]
fn stereo_validation_reports_invalid_local_elements_without_mutating() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon());
    let a = mol.add_atom(oxygen());
    let b = mol.add_atom(element_atom("N"));
    mol.add_bond(center, a, BondOrder::Single).expect("bond");
    mark_all_fresh(&mut mol);
    let element = mol
        .add_stereo_element(StereoElement {
            kind: StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(a),
                    StereoCarrier::Atom(b),
                ],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            specifiedness: StereoSpecifiedness::Unknown,
            source: StereoSource::User,
            group: None,
            descriptor: None,
        })
        .expect("stereo element");
    mark_all_fresh(&mut mol);

    let report = stereo_api::validate_stereo(&mol);

    assert!(mol.stereo_elements().next().is_some());
    assert!(report
        .issues
        .contains(&StereoPerceptionIssue::InvalidTetrahedralCarrierCount {
            element,
            center,
            carrier_count: 3,
        }));
    assert!(report
        .issues
        .contains(&StereoPerceptionIssue::DuplicateTetrahedralCarrier {
            element,
            center,
            carrier: StereoCarrier::Atom(a),
        }));
    assert!(report
        .issues
        .contains(&StereoPerceptionIssue::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier: StereoCarrier::Atom(b),
        }));
    assert!(
        mol.stereo_element(element).expect("element").specifiedness == StereoSpecifiedness::Unknown
    );
}

#[test]
fn stereo_candidates_use_sanitized_hydrogen_state_without_cip_assignment() {
    let mut molecule = read_smiles("CC(F)(Cl)Br").expect("smiles should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("molecule should sanitize");

    let report = stereo_api::validate_stereo(molecule.graph());

    assert!(report.is_ok());
    assert!(report.candidates.iter().any(|candidate| matches!(
        candidate,
        StereoCandidate::Tetrahedral { center, carriers }
            if *center == AtomId::new(1)
                && carriers.len() == 4
                && !carriers.contains(&StereoCarrier::ImplicitHydrogen)
    )));
    assert!(molecule.graph().stereo_elements().next().is_none());
}

#[test]
fn stereo_perception_assembles_paired_directional_marks_into_double_bond_element() {
    let mut molecule = read_smiles("C/C=C\\F").expect("directional smiles should parse");
    perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )
    .expect("molecule should sanitize");

    let report = stereo_api::perceive_stereo(molecule.graph_mut());

    assert!(report.is_ok());
    assert_eq!(report.created_elements.len(), 1);
    assert!(molecule.graph().stereo_elements().next().is_some());
    let element = molecule
        .graph()
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    match &element.kind {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.bond, BondId::new(1));
            assert_eq!(stereo.left, AtomId::new(1));
            assert_eq!(stereo.right, AtomId::new(2));
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(AtomId::new(0)));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(AtomId::new(3)));
            assert_eq!(stereo.orientation, DoubleBondOrientation::Together);
        }
        other => panic!("expected double-bond stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_skips_small_ring_double_bond_stereo_boundary() {
    let mut cyclohexene = read_smiles(r"C1/C=C\CCC1").expect("marked cyclohexene parses");
    perception_api::sanitize_with_options(
        &mut cyclohexene,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )
    .expect("marked cyclohexene sanitizes without stereo perception");
    let report = stereo_api::perceive_stereo(cyclohexene.graph_mut_raw());

    assert!(report.created_elements.is_empty());
    assert!(cyclohexene.graph().stereo_elements().next().is_none());

    let mut cyclooctene = read_smiles(r"C1/C=C\CCCCC1").expect("marked cyclooctene parses");
    perception_api::sanitize_with_options(
        &mut cyclooctene,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )
    .expect("marked cyclooctene sanitizes without stereo perception");
    let report = stereo_api::perceive_stereo(cyclooctene.graph_mut_raw());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = cyclooctene
        .graph()
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert!(matches!(element.kind, StereoElementKind::DoubleBond(_)));
}

#[test]
fn stereo_perception_assembles_molfile_wedge_into_tetrahedral_element() {
    let input = "\
wedge
molecular

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0
   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 I   0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
";
    let mut molecule = read_molfile(input).expect("wedge molfile should parse");

    let report = stereo_api::perceive_stereo(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element.specifiedness, StereoSpecifiedness::Specified);
    assert_eq!(element.source, StereoSource::MolfileV2000);
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, AtomId::new(0));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(1)),
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::Atom(AtomId::new(4)),
                ]
            );
            assert_eq!(stereo.orientation, TetrahedralOrientation::CounterClockwise);
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_uses_virtual_implicit_h_for_molfile_wedge_geometry() {
    let mut molecule = read_molfile(implicit_h_wedge_geometry_molblock())
        .expect("implicit-H wedge molfile should parse");
    perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )
    .expect("implicit-H wedge molfile should sanitize");

    let report = stereo_api::perceive_stereo(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, AtomId::new(0));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(1)),
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(3)),
                    StereoCarrier::ImplicitHydrogen,
                ]
            );
            assert_eq!(stereo.orientation, TetrahedralOrientation::CounterClockwise);
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_assembles_wedge_either_as_explicit_unknown() {
    let (mut mol, center, carriers, marked_bond) = tetrahedral_marked_graph();
    mol.set_stereo_bond_mark(StereoBondMark {
        bond: marked_bond,
        kind: StereoBondMarkKind::WedgeEither,
        source: StereoSource::MolfileV2000,
    })
    .expect("wedge mark");

    let report = stereo_api::perceive_stereo(&mut mol);

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element.specifiedness, StereoSpecifiedness::Unknown);
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, center);
            assert_eq!(stereo.carriers[0], StereoCarrier::Atom(carriers[0]));
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_reports_ambiguous_tetrahedral_wedge_marks() {
    let (mut mol, center, _carriers, first_bond) = tetrahedral_marked_graph();
    let second_bond = BondId::new(1);
    mol.set_stereo_bond_mark(StereoBondMark {
        bond: first_bond,
        kind: StereoBondMarkKind::WedgeUp,
        source: StereoSource::MolfileV2000,
    })
    .expect("first wedge mark");
    mol.set_stereo_bond_mark(StereoBondMark {
        bond: second_bond,
        kind: StereoBondMarkKind::WedgeDown,
        source: StereoSource::MolfileV2000,
    })
    .expect("second wedge mark");

    let report = stereo_api::perceive_stereo(&mut mol);

    assert!(report
        .issues
        .contains(&StereoPerceptionIssue::AmbiguousTetrahedralWedgeMarks {
            center,
            mark_count: 2,
        }));
    assert!(report.created_elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn stereo_perception_assigns_tetrahedral_from_3d_coordinates() {
    let (mut mol, center, carriers, _) = tetrahedral_marked_graph();
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            center,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            carriers[0],
            crate::units::Quantity::new(Point3::new(1.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            carriers[1],
            crate::units::Quantity::new(Point3::new(0.0, 1.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            carriers[2],
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 1.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            carriers[3],
            crate::units::Quantity::new(Point3::new(0.0, 0.0, -1.0), crate::units::ANGSTROM),
        )
        .unwrap();
    mol.add_conformer(conformer).expect("valid conformer");

    let report = stereo_api::perceive_stereo(&mut mol);

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element.source, StereoSource::Coordinates3D);
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, center);
            assert_eq!(
                stereo.carriers,
                carriers
                    .iter()
                    .copied()
                    .map(StereoCarrier::Atom)
                    .collect::<Vec<_>>()
            );
            assert_eq!(stereo.orientation, TetrahedralOrientation::Clockwise);
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_assigns_double_bond_from_2d_coordinates() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon());
    let right = mol.add_atom(carbon());
    let left_carrier = mol.add_atom(element_atom("F"));
    let right_carrier = mol.add_atom(element_atom("Cl"));
    let double_bond = mol.add_bond(left, right, BondOrder::Double).expect("bond");
    mol.add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    mol.add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            left,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            right,
            crate::units::Quantity::new(Point3::new(1.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            left_carrier,
            crate::units::Quantity::new(Point3::new(0.0, 1.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            right_carrier,
            crate::units::Quantity::new(Point3::new(1.0, -1.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    mol.add_conformer(conformer).expect("valid conformer");

    let report = stereo_api::perceive_stereo(&mut mol);

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element.source, StereoSource::Coordinates2D);
    match &element.kind {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.bond, double_bond);
            assert_eq!(stereo.left, left);
            assert_eq!(stereo.right, right);
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(left_carrier));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(right_carrier));
            assert_eq!(stereo.orientation, DoubleBondOrientation::Opposite);
        }
        other => panic!("expected double-bond stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_assigns_axis_from_3d_coordinates() {
    let (mut mol, axis) = coordinate_axis_graph(true);

    let report = stereo_api::perceive_stereo_with_options(
        &mut mol,
        StereoPerceptionOptions {
            assign_coordinate_axes: true,
            ..StereoPerceptionOptions::default()
        },
    );

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element.source, StereoSource::Coordinates3D);
    match &element.kind {
        StereoElementKind::Axis(stereo) => {
            assert_eq!(stereo.axis, axis);
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(2)),
                    StereoCarrier::Atom(AtomId::new(4)),
                ]
            );
            assert_eq!(stereo.orientation, AxisOrientation::Clockwise);
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_skips_coordinate_axis_without_3d_handedness() {
    let (mut mol, _axis) = coordinate_axis_graph(false);

    let report = stereo_api::perceive_stereo_with_options(
        &mut mol,
        StereoPerceptionOptions {
            assign_coordinate_axes: true,
            ..StereoPerceptionOptions::default()
        },
    );

    assert!(report.is_ok(), "{:?}", report.issues);
    assert!(report.created_elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn stereo_perception_leaves_coordinate_axes_opt_in_by_default() {
    let (mut mol, _axis) = coordinate_axis_graph(true);

    let report = stereo_api::perceive_stereo(&mut mol);

    assert!(report.is_ok(), "{:?}", report.issues);
    assert!(report.created_elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn stereo_perception_reports_unassembled_marks_and_preserves_absence() {
    let mut marked = Molecule::new();
    let a = marked.add_atom(carbon());
    let b = marked.add_atom(carbon());
    let bond = marked.add_bond(a, b, BondOrder::Single).expect("bond");
    marked
        .set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("mark");

    let marked_report = stereo_api::perceive_stereo(&mut marked);
    assert!(marked_report.issues.contains(
        &StereoPerceptionIssue::UnassembledTetrahedralBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeEither,
        }
    ));
    assert!(marked.stereo_elements().next().is_none());

    let mut unsupported = Molecule::new();
    let c = unsupported.add_atom(carbon());
    let d = unsupported.add_atom(carbon());
    let double_bond = unsupported.add_bond(c, d, BondOrder::Double).expect("bond");
    unsupported
        .set_stereo_bond_mark(StereoBondMark {
            bond: double_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("double bond either mark");
    let unsupported_report = stereo_api::perceive_stereo(&mut unsupported);
    assert!(unsupported_report.issues.contains(
        &StereoPerceptionIssue::UnsupportedSourceBondMark {
            bond: double_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
        }
    ));

    let mut unknown = Molecule::new();
    let left = unknown.add_atom(carbon());
    let right = unknown.add_atom(carbon());
    let left_carrier = unknown.add_atom(carbon());
    let right_carrier = unknown.add_atom(carbon());
    let unknown_bond = unknown
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    unknown
        .add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    unknown
        .add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    unknown
        .set_stereo_bond_mark(StereoBondMark {
            bond: unknown_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("double bond either mark");

    let unknown_report = stereo_api::perceive_stereo(&mut unknown);
    assert!(unknown_report.is_ok(), "{:?}", unknown_report.issues);
    assert_eq!(unknown_report.created_elements.len(), 1);
    let (_, element) = unknown.stereo_elements().next().expect("unknown element");
    assert_eq!(element.specifiedness, StereoSpecifiedness::Unknown);
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo) if stereo.bond == unknown_bond
    ));

    let mut absent = Molecule::new();
    let x = absent.add_atom(carbon());
    let y = absent.add_atom(carbon());
    absent.add_bond(x, y, BondOrder::Single).expect("bond");
    let absent_report = stereo_api::perceive_stereo(&mut absent);
    assert!(absent_report.is_ok());
    assert!(absent.stereo_elements().next().is_none());
    assert!(absent.stereo_bond_marks().next().is_none());
}

#[test]
fn stereo_validation_accepts_structural_axis_elements() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon());
    let right = mol.add_atom(carbon());
    let left_carrier = mol.add_atom(element_atom("I"));
    let right_carrier = mol.add_atom(element_atom("Br"));
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    mol.add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    mol.add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let valid_axis = mol
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Axis(AxisStereo {
                axis,
                carriers: vec![
                    StereoCarrier::Atom(left_carrier),
                    StereoCarrier::Atom(right_carrier),
                ],
                orientation: AxisOrientation::CounterClockwise,
            }),
            StereoSource::User,
        ))
        .expect("axis element");

    let report = stereo_api::validate_stereo(&mol);

    assert!(report.is_ok(), "{:?}", report.issues);

    mol.remove_stereo_element(valid_axis)
        .expect("remove valid axis");
    let invalid_axis = mol
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Axis(AxisStereo {
                axis,
                carriers: vec![StereoCarrier::Atom(left_carrier)],
                orientation: AxisOrientation::CounterClockwise,
            }),
            StereoSource::User,
        ))
        .expect("invalid axis element refs are still structurally present");

    let report = stereo_api::validate_stereo(&mol);

    assert_eq!(
        report.issues,
        vec![StereoPerceptionIssue::InvalidAxisCarrierCount {
            element: invalid_axis,
            axis,
            carrier_count: 1,
        }]
    );
}

#[test]
fn stereo_perception_assembles_molfile_atropisomeric_axis() {
    let mut molecule =
        read_molfile(rdkit_rp6306_atrop_molblock()).expect("RDKit atropisomer fixture parses");

    let report = stereo_api::perceive_stereo(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_elements[0])
        .expect("created axis element");
    assert_eq!(element.source, StereoSource::MolfileV2000);
    match &element.kind {
        StereoElementKind::Axis(stereo) => {
            assert_eq!(stereo.axis, BondId::new(3));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(6)),
                    StereoCarrier::Atom(AtomId::new(11)),
                ]
            );
            assert_eq!(stereo.orientation, AxisOrientation::Clockwise);
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_prefers_exocyclic_molfile_atropisomeric_axis() {
    let mut molecule = read_molfile(rdkit_rp6306_atrop3_molblock())
        .expect("RDKit alternate atropisomer fixture parses");

    let report = stereo_api::perceive_stereo(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_elements[0])
        .expect("created axis element");
    match &element.kind {
        StereoElementKind::Axis(stereo) => {
            assert_eq!(stereo.axis, BondId::new(3));
            assert_eq!(
                stereo.carriers,
                vec![
                    StereoCarrier::Atom(AtomId::new(6)),
                    StereoCarrier::Atom(AtomId::new(11)),
                ]
            );
            assert_eq!(stereo.orientation, AxisOrientation::Clockwise);
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn stereo_perception_consumes_redundant_molfile_atrop_wedges_before_tetrahedral_marks() {
    let mut molecule = read_molfile(rdkit_bms986142_atrop5_molblock())
        .expect("RDKit redundant atropisomer wedge fixture parses");
    perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_stereo: false,
            ..SanitizeOptions::default()
        },
    )
    .expect("fixture prepares before stereo perception");

    let report = stereo_api::perceive_stereo(molecule.graph_mut());

    assert!(report.is_ok(), "{:?}", report.issues);
    assert_eq!(report.created_elements.len(), 2);
    assert!(molecule.graph().stereo_elements().any(|(_, element)| {
        matches!(&element.kind, StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(10))
    }));
    assert!(molecule.graph().stereo_elements().any(|(_, element)| {
        matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(8))
    }));
}

#[test]
fn stereo_perception_assembles_molfile_atrop_axis_with_one_exocyclic_sp2_endpoint() {
    for fixture in [
        rdkit_zm374979_atrop1_molblock(),
        rdkit_zm374979_atrop2_molblock(),
    ] {
        let mut molecule =
            read_molfile(fixture).expect("RDKit one-ring-endpoint atropisomer fixture parses");
        perception_api::sanitize_with_options(
            &mut molecule,
            SanitizeOptions {
                perceive_stereo: false,
                ..SanitizeOptions::default()
            },
        )
        .expect("fixture prepares before stereo perception");

        let report = stereo_api::perceive_stereo(molecule.graph_mut());

        assert!(report.is_ok(), "{:?}", report.issues);
        assert_eq!(report.created_elements.len(), 2);
        assert!(molecule.graph().stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(3))
        }));
        assert!(molecule.graph().stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(33))
        }));
    }
}

#[test]
fn stereo_perception_assembles_ring_internal_molfile_atrop_axis() {
    for fixture in [
        rdkit_macrocycle8_ortho_wedge_molblock(),
        rdkit_macrocycle8_ortho_hash_molblock(),
    ] {
        let mut molecule =
            read_molfile(fixture).expect("RDKit macrocyclic atropisomer fixture parses");
        perception_api::sanitize_with_options(
            &mut molecule,
            SanitizeOptions {
                perceive_stereo: false,
                ..SanitizeOptions::default()
            },
        )
        .expect("fixture prepares before stereo perception");

        let report = stereo_api::perceive_stereo(molecule.graph_mut());

        assert!(report.is_ok(), "{:?}", report.issues);
        assert_eq!(report.created_elements.len(), 1);
        assert!(molecule.graph().stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(15))
        }));
    }
}

fn tetrahedral_marked_graph() -> (Molecule, AtomId, Vec<AtomId>, BondId) {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon());
    let carriers = ["F", "Cl", "Br", "I"]
        .into_iter()
        .map(element_atom)
        .map(|atom| mol.add_atom(atom))
        .collect::<Vec<_>>();
    let mut bonds = Vec::new();
    for carrier in &carriers {
        bonds.push(
            mol.add_bond(center, *carrier, BondOrder::Single)
                .expect("tetrahedral carrier bond"),
        );
    }
    (mol, center, carriers, bonds[0])
}
