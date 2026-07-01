use super::*;

#[test]
fn element_from_atomic_number_accepts_periodic_table_bounds() {
    assert_eq!(
        Element::from_atomic_number(1)
            .expect("hydrogen exists")
            .symbol(),
        "H"
    );
    assert_eq!(
        Element::from_atomic_number(118)
            .expect("oganesson exists")
            .symbol(),
        "Og"
    );
}

#[test]
fn element_from_atomic_number_rejects_out_of_range_values() {
    assert_eq!(Element::from_atomic_number(0), None);
    assert_eq!(Element::from_atomic_number(119), None);
}

#[test]
fn element_from_symbol_is_canonical_and_case_sensitive() {
    assert_eq!(
        Element::from_symbol("C")
            .expect("carbon exists")
            .atomic_number(),
        6
    );
    assert_eq!(
        Element::from_symbol("Cl")
            .expect("chlorine exists")
            .atomic_number(),
        17
    );
    assert_eq!(
        Element::from_symbol("Og")
            .expect("oganesson exists")
            .atomic_number(),
        118
    );
    assert_eq!(Element::from_symbol("CL"), None);
    assert_eq!(Element::from_symbol("Xx"), None);
    assert_eq!(Element::from_symbol("?"), None);
}

#[test]
fn element_symbol_and_display_are_canonical() {
    let iron = Element::from_atomic_number(26).expect("iron exists");

    assert_eq!(iron.symbol(), "Fe");
    assert_eq!(iron.to_string(), "Fe");
}

#[test]
fn atom_new_sets_chemically_general_defaults() {
    let atom = carbon();

    assert_eq!(atom.element.symbol(), "C");
    assert_eq!(atom.isotope, None);
    assert_eq!(atom.formal_charge, 0);
    assert_eq!(atom.radical, None);
    assert_eq!(atom.explicit_hydrogens, 0);
    assert_eq!(atom.implicit_hydrogens, None);
    assert!(!atom.no_implicit_hydrogens);
    assert!(!atom.aromatic);
    assert_eq!(atom.chiral, None);
    assert_eq!(atom.atom_map, None);
    assert!(atom.props.is_empty());
}

#[test]
fn atom_payload_fields_can_be_set_and_read() {
    let mut atom = carbon();
    atom.isotope = Some(13);
    atom.formal_charge = -1;
    atom.radical = Some(AtomRadical::Doublet);
    atom.explicit_hydrogens = 3;
    atom.implicit_hydrogens = Some(1);
    atom.no_implicit_hydrogens = true;
    atom.aromatic = true;
    atom.chiral = Some(AtomStereo::TetrahedralClockwise);
    atom.atom_map = Some(7);
    atom.props
        .insert("label".to_owned(), PropValue::String("alpha".to_owned()));

    assert_eq!(atom.isotope, Some(13));
    assert_eq!(atom.formal_charge, -1);
    assert_eq!(atom.radical, Some(AtomRadical::Doublet));
    assert_eq!(atom.explicit_hydrogens, 3);
    assert_eq!(atom.implicit_hydrogens, Some(1));
    assert!(atom.no_implicit_hydrogens);
    assert!(atom.aromatic);
    assert_eq!(atom.chiral, Some(AtomStereo::TetrahedralClockwise));
    assert_eq!(atom.atom_map, Some(7));
    assert_eq!(
        atom.props.get("label"),
        Some(&PropValue::String("alpha".to_owned()))
    );
}

#[test]
fn radical_multiplicity_reports_unpaired_electrons() {
    assert_eq!(AtomRadical::Singlet.unpaired_electron_count(), 0);
    assert_eq!(AtomRadical::Doublet.unpaired_electron_count(), 1);
    assert_eq!(AtomRadical::Triplet.unpaired_electron_count(), 2);
}

#[test]
fn bond_new_sets_endpoints_order_and_aromatic_default() {
    let a = AtomId::new(3);
    let b = AtomId::new(4);
    let single = Bond::new(a, b, BondOrder::Single);
    let aromatic = Bond::new(a, b, BondOrder::Aromatic);

    assert_eq!(single.a(), a);
    assert_eq!(single.b(), b);
    assert_eq!(single.endpoints(), (a, b));
    assert_eq!(single.order, BondOrder::Single);
    assert!(!single.aromatic);
    assert_eq!(single.stereo, None);
    assert!(single.props.is_empty());
    assert!(aromatic.aromatic);
}

#[test]
fn bond_payload_fields_can_be_set_and_read() {
    let mut bond = Bond::new(AtomId::new(1), AtomId::new(2), BondOrder::Dative);
    bond.stereo = Some(BondStereo::Up);
    bond.props
        .insert("score".to_owned(), PropValue::Float(1.25));

    assert_eq!(bond.order, BondOrder::Dative);
    assert_eq!(bond.stereo, Some(BondStereo::Up));
    assert_eq!(bond.props.get("score"), Some(&PropValue::Float(1.25)));
}

#[test]
fn prop_value_equality_covers_all_initial_variants() {
    assert_eq!(
        PropValue::String("value".to_owned()),
        PropValue::String("value".to_owned())
    );
    assert_eq!(PropValue::Int(42), PropValue::Int(42));
    assert_eq!(PropValue::Float(2.5), PropValue::Float(2.5));
    assert_eq!(PropValue::Bool(true), PropValue::Bool(true));
}

#[test]
fn mutable_payload_access_invalidates_fresh_perception() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(oxygen());
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    mark_all_fresh(&mut mol);
    mol.atom_mut(a).expect("atom exists").formal_charge = 1;
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    mol.bond_mut(bond).expect("bond exists").order = BondOrder::Double;
    assert_all_stale(&mol);
}

#[test]
fn perception_owned_chemistry_edits_invalidate_dependent_state() {
    let mut methane = Molecule::new();
    methane.add_atom(carbon());
    mark_all_fresh(&mut methane);

    let report = valence_api::perceive_valence(&mut methane, ValenceModel::RdkitLike);

    assert!(report.is_ok());
    assert_eq!(methane.perception().valence, ComputedState::Fresh);
    assert_eq!(methane.perception().rings, ComputedState::Fresh);
    assert_eq!(methane.perception().aromaticity, ComputedState::Stale);
    assert_eq!(methane.perception().stereo, ComputedState::Stale);

    let (mut benzene, _, _) = ring_molecule(
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
    mark_all_fresh(&mut benzene);

    aromaticity_api::perceive_aromaticity(&mut benzene, AromaticityModel::RdkitLike)
        .expect("benzene should be supported");

    assert_eq!(benzene.perception().valence, ComputedState::Fresh);
    assert_eq!(benzene.perception().rings, ComputedState::Fresh);
    assert_eq!(benzene.perception().aromaticity, ComputedState::Fresh);
    assert_eq!(benzene.perception().stereo, ComputedState::Stale);
}
