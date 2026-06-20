use super::*;

#[test]
fn valence_and_sanitization_are_explicit() {
    let mut small = read_smiles_str("CCO", SmilesParseOptions).expect("smiles should parse");
    assert_eq!(small.mol.perception().valence, ComputedState::Absent);

    let report = sanitize_small_molecule(&mut small, SanitizeOptions::default())
        .expect("ethanol should sanitize");

    assert!(report.valence.expect("valence report").is_ok());
    assert_eq!(small.mol.perception().valence, ComputedState::Fresh);
    assert_eq!(small.mol.perception().rings, ComputedState::Fresh);
    assert_eq!(
        small
            .mol
            .atom(AtomId::new(2))
            .expect("oxygen")
            .implicit_hydrogens,
        Some(1)
    );
}

#[test]
fn sanitize_options_do_not_leave_skipped_passes_fresh() {
    let mut baseline =
        read_smiles_str("C1=CC=CC=C1", SmilesParseOptions).expect("benzene should parse");
    sanitize_small_molecule(&mut baseline, SanitizeOptions::default())
        .expect("benzene should sanitize");

    for mask in 0..8 {
        let options = SanitizeOptions {
            perceive_valence: mask & 1 != 0,
            perceive_rings: mask & 2 != 0,
            perceive_aromaticity: mask & 4 != 0,
        };
        let mut molecule = baseline.clone();
        sanitize_small_molecule(&mut molecule, options)
            .unwrap_or_else(|error| panic!("options {mask:03b} should succeed: {error}"));

        assert_eq!(
            molecule.mol.perception().valence,
            if options.perceive_valence {
                ComputedState::Fresh
            } else {
                ComputedState::Stale
            },
            "valence state for options {mask:03b}"
        );
        assert_eq!(
            molecule.mol.perception().rings,
            if options.perceive_rings {
                ComputedState::Fresh
            } else {
                ComputedState::Stale
            },
            "ring state for options {mask:03b}"
        );
        assert_eq!(
            molecule.mol.perception().aromaticity,
            if options.perceive_aromaticity {
                ComputedState::Fresh
            } else {
                ComputedState::Stale
            },
            "aromaticity state for options {mask:03b}"
        );
        assert_eq!(
            molecule.mol.ring_set().is_some(),
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
    perceive_ring_set(&mut mol).expect("ring perception should succeed");
    let mut molecule = SmallMolecule { mol };
    let before = molecule.clone();

    let error = sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect_err("pentavalent carbon should fail valence");

    assert!(matches!(error, SanitizeError::Valence(_)));
    assert_eq!(molecule, before);
}

#[test]
fn failed_aromaticity_sanitization_is_transactional() {
    let (mol, _, _) = ring_molecule(&["Si", "C", "C", "C", "C", "C"], &[BondOrder::Aromatic; 6]);
    let mut molecule = SmallMolecule { mol };
    let before = molecule.clone();

    let error = sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect_err("unsupported explicitly aromatic silicon should fail");

    assert!(matches!(error, SanitizeError::Aromaticity(_)));
    assert_eq!(molecule, before);
}

#[test]
fn successful_sanitization_is_idempotent() {
    let mut molecule = read_smiles_str("CCO", SmilesParseOptions).expect("ethanol should parse");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("first sanitize should succeed");
    let once = molecule.clone();

    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
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
    let mut molecule = SmallMolecule { mol };

    sanitize_small_molecule(
        &mut molecule,
        SanitizeOptions {
            perceive_valence: false,
            perceive_rings: false,
            perceive_aromaticity: false,
        },
    )
    .expect("cleanup-only sanitize should succeed");

    assert_all_stale(&molecule.mol);
    assert_eq!(
        molecule.mol.atom(chlorine).expect("chlorine").formal_charge,
        1
    );
    assert_eq!(molecule.mol.atom(oxygen).expect("oxygen").formal_charge, -1);
    assert_eq!(
        molecule.mol.bond(BondId::new(0)).expect("bond").order,
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

    let report = perceive_valence(&mut mol, ValenceModel::RdkitLike);

    assert_eq!(report.issues.len(), 1);
    assert!(!report.is_ok());
}
