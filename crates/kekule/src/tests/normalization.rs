use super::*;

fn oxo_halide(oxo_count: usize) -> (SmallMolecule, AtomId, Vec<AtomId>, Vec<BondId>) {
    let mut molecule = Molecule::new();
    let chlorine = molecule
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let hydroxyl = molecule
        .add_atom(oxygen())
        .expect("atom identifier capacity");
    molecule
        .add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("hydroxyl bond");

    let mut oxygens = Vec::with_capacity(oxo_count);
    let mut bonds = Vec::with_capacity(oxo_count);
    for _ in 0..oxo_count {
        let oxygen = molecule
            .add_atom(oxygen())
            .expect("atom identifier capacity");
        let bond = molecule
            .add_bond(chlorine, oxygen, BondOrder::Double)
            .expect("oxo bond");
        oxygens.push(oxygen);
        bonds.push(bond);
    }
    (
        SmallMolecule::from_graph(molecule),
        chlorine,
        oxygens,
        bonds,
    )
}

#[test]
fn normalization_rewrites_hypervalent_oxo_halide_representation() {
    let (mut molecule, chlorine, oxygens, bonds) = oxo_halide(2);

    normalization_api::normalize(molecule.graph_mut()).expect("normalization should succeed");

    assert_eq!(molecule.graph().atom(chlorine).unwrap().formal_charge, 2);
    for oxygen in oxygens {
        assert_eq!(molecule.graph().atom(oxygen).unwrap().formal_charge, -1);
    }
    for bond in bonds {
        assert_eq!(
            molecule.graph().bond(bond).unwrap().order,
            BondOrder::Single
        );
    }
    assert_eq!(molecule.graph().formal_charge(), 0);
}

#[test]
fn normalization_is_idempotent() {
    let (mut molecule, ..) = oxo_halide(2);

    molecule
        .normalize()
        .expect("first normalization should succeed");
    let once = molecule.clone();
    molecule
        .normalize()
        .expect("second normalization should succeed");

    assert_eq!(molecule, once);
}

#[test]
fn successful_normalization_clears_installed_perception() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");
    mark_all_fresh(molecule.graph_mut());
    assert_ne!(molecule.graph().perception(), &PerceptionState::default());

    molecule.normalize().expect("normalization should succeed");

    assert_eq!(molecule.graph().perception(), &PerceptionState::default());
    assert_all_stale(molecule.graph());
}

#[test]
fn failed_normalization_is_transactional() {
    let (mut molecule, chlorine, ..) = oxo_halide(128);
    mark_all_fresh(molecule.graph_mut());
    let before = molecule.clone();

    let error = molecule
        .normalize()
        .expect_err("unrepresentable formal charge should fail");

    assert_eq!(
        error,
        crate::normalization::NormalizationError::FormalChargeOutOfRange {
            atom: chlorine,
            charge: 128,
        }
    );
    assert_eq!(molecule, before);
}

#[test]
fn sanitizer_delegates_normalization_failure_transactionally() {
    let (mut molecule, chlorine, ..) = oxo_halide(128);
    let before = molecule.clone();

    let error = perception_api::sanitize_with_options(
        &mut molecule,
        SanitizeOptions {
            perceive_valence: false,
            perceive_rings: false,
            perceive_aromaticity: false,
            perceive_stereo: false,
        },
    )
    .expect_err("normalization failure should fail sanitization");

    assert_eq!(
        error,
        SanitizeError::Normalization(
            crate::normalization::NormalizationError::FormalChargeOutOfRange {
                atom: chlorine,
                charge: 128,
            }
        )
    );
    assert_eq!(molecule, before);
}
