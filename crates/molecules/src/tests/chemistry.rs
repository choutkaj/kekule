use super::*;

#[test]
fn valence_and_sanitization_are_explicit() {
    let mut small =
        smiles_api::read_str_with_options("CCO", SmilesParseOptions).expect("smiles should parse");
    assert_eq!(small.graph().perception().valence, ComputedState::Absent);

    let report = perception_api::sanitize_with_options(&mut small, SanitizeOptions::default())
        .expect("ethanol should sanitize");

    assert!(report.valence.expect("valence report").is_ok());
    assert_eq!(small.graph().perception().valence, ComputedState::Fresh);
    assert_eq!(small.graph().perception().rings, ComputedState::Fresh);
    assert_eq!(
        small
            .graph()
            .atom(AtomId::new(2))
            .expect("oxygen")
            .implicit_hydrogens,
        Some(1)
    );
}

#[test]
fn sanitize_options_do_not_leave_skipped_passes_fresh() {
    let mut baseline = smiles_api::read_str_with_options("C1=CC=CC=C1", SmilesParseOptions)
        .expect("benzene should parse");
    perception_api::sanitize_with_options(&mut baseline, SanitizeOptions::default())
        .expect("benzene should sanitize");

    for mask in 0..8 {
        let options = SanitizeOptions {
            perceive_valence: mask & 1 != 0,
            perceive_rings: mask & 2 != 0,
            perceive_aromaticity: mask & 4 != 0,
        };
        let mut molecule = baseline.clone();
        perception_api::sanitize_with_options(&mut molecule, options)
            .unwrap_or_else(|error| panic!("options {mask:03b} should succeed: {error}"));

        assert_eq!(
            molecule.graph().perception().valence,
            if options.perceive_valence {
                ComputedState::Fresh
            } else {
                ComputedState::Stale
            },
            "valence state for options {mask:03b}"
        );
        assert_eq!(
            molecule.graph().perception().rings,
            if options.perceive_rings {
                ComputedState::Fresh
            } else {
                ComputedState::Stale
            },
            "ring state for options {mask:03b}"
        );
        assert_eq!(
            molecule.graph().perception().aromaticity,
            if options.perceive_aromaticity {
                ComputedState::Fresh
            } else {
                ComputedState::Stale
            },
            "aromaticity state for options {mask:03b}"
        );
        assert_eq!(
            molecule.graph().ring_set().is_some(),
            options.perceive_rings,
            "ring cache exposure for options {mask:03b}"
        );
    }
}

#[test]
fn failed_valence_sanitization_is_transactional() {
    let mut mol = Molecule::new();
    let carbon = mol.add_atom(carbon());
    for _ in 0..5 {
        let hydrogen = mol.add_atom(Atom::new(Element::from_symbol("H").expect("hydrogen")));
        mol.add_bond(carbon, hydrogen, BondOrder::Single)
            .expect("bond");
    }
    rings_api::perceive_ring_set(&mut mol).expect("ring perception should succeed");
    let mut molecule = SmallMolecule::from_graph(mol);
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect_err("pentavalent carbon should fail valence");

    assert!(matches!(error, SanitizeError::Valence(_)));
    assert_eq!(molecule, before);
}

#[test]
fn failed_aromaticity_sanitization_is_transactional() {
    let (mol, _, _) = ring_molecule(&["Si", "C", "C", "C", "C", "C"], &[BondOrder::Aromatic; 6]);
    let mut molecule = SmallMolecule::from_graph(mol);
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect_err("unsupported explicitly aromatic silicon should fail");

    assert!(matches!(error, SanitizeError::Aromaticity(_)));
    assert_eq!(molecule, before);
}

#[test]
fn successful_sanitization_is_idempotent() {
    let mut molecule =
        smiles_api::read_str_with_options("CCO", SmilesParseOptions).expect("ethanol should parse");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("first sanitize should succeed");
    let once = molecule.clone();

    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("second sanitize should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn sanitize_cleanup_invalidates_preexisting_perception() {
    let mut mol = Molecule::new();
    let chlorine = mol.add_atom(Atom::new(Element::from_symbol("Cl").expect("chlorine")));
    let oxygen = mol.add_atom(oxygen());
    mol.add_bond(chlorine, oxygen, BondOrder::Double)
        .expect("bond");
    mark_all_fresh(&mut mol);
    let mut molecule = SmallMolecule::from_graph(mol);

    perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_valence: false,
            perceive_rings: false,
            perceive_aromaticity: false,
        },
    )
    .expect("cleanup-only sanitize should succeed");

    assert_all_stale(molecule.graph());
    assert_eq!(
        molecule
            .graph()
            .atom(chlorine)
            .expect("chlorine")
            .formal_charge,
        1
    );
    assert_eq!(
        molecule.graph().atom(oxygen).expect("oxygen").formal_charge,
        -1
    );
    assert_eq!(
        molecule.graph().bond(BondId::new(0)).expect("bond").order,
        BondOrder::Single
    );
}

#[test]
fn valence_reports_excess_common_valence() {
    let mut mol = Molecule::new();
    let c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("C")));
    for _ in 0..5 {
        let h = mol.add_atom(Atom::new(Element::from_symbol("H").expect("H")));
        mol.add_bond(c, h, BondOrder::Single).expect("bond");
    }

    let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

    assert_eq!(report.issues.len(), 1);
    assert!(!report.is_ok());
}

#[test]
fn valence_supports_simple_pubchem_main_group_ions_and_salts() {
    for (symbol, charge, expected_implicit_hydrogens) in [
        ("H", 1, 0),
        ("H", -1, 0),
        ("Rb", 1, 0),
        ("Cs", 1, 0),
        ("Be", 2, 0),
        ("Al", 3, 0),
        ("Ga", 3, 0),
        ("Tl", 1, 0),
        ("U", 2, 0),
        ("Pb", 2, 0),
        ("S", -2, 0),
        ("Se", -2, 0),
    ] {
        let mut mol = Molecule::new();
        let atom_id = mol.add_atom(charged_atom(symbol, charge));

        let report = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);

        assert!(report.is_ok(), "{symbol}{charge:+} should be supported");
        assert_eq!(
            mol.atom(atom_id).expect("atom").implicit_hydrogens,
            Some(expected_implicit_hydrogens),
            "{symbol}{charge:+} implicit hydrogens"
        );
    }

    let mut covalent_aluminum = Molecule::new();
    let aluminum = covalent_aluminum.add_atom(element_atom("Al"));
    for _ in 0..3 {
        let chlorine = covalent_aluminum.add_atom(element_atom("Cl"));
        covalent_aluminum
            .add_bond(aluminum, chlorine, BondOrder::Single)
            .expect("bond");
    }

    let report = valence_api::perceive_valence(&mut covalent_aluminum, ValenceModel::RdkitLike);

    assert!(
        report.is_ok(),
        "neutral trivalent aluminum should be supported"
    );
    assert_eq!(
        covalent_aluminum
            .atom(aluminum)
            .expect("aluminum")
            .implicit_hydrogens,
        Some(0)
    );

    let mut neutral_magnesium = Molecule::new();
    let magnesium = neutral_magnesium.add_atom(element_atom("Mg"));
    for _ in 0..2 {
        let chlorine = neutral_magnesium.add_atom(element_atom("Cl"));
        neutral_magnesium
            .add_bond(magnesium, chlorine, BondOrder::Single)
            .expect("bond");
    }

    let report = valence_api::perceive_valence(&mut neutral_magnesium, ValenceModel::RdkitLike);

    assert!(
        report.is_ok(),
        "neutral divalent magnesium should be supported"
    );
    assert_eq!(
        neutral_magnesium
            .atom(magnesium)
            .expect("magnesium")
            .implicit_hydrogens,
        Some(0)
    );
}
