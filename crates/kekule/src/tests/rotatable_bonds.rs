use super::*;

use crate::rotatable_bonds::{self, RotatableBondOptions};

fn detected(smiles: &str) -> Vec<BondId> {
    detected_with_options(smiles, RotatableBondOptions::STRICT)
}

fn detected_with_options(smiles: &str, options: RotatableBondOptions) -> Vec<BondId> {
    let molecule = read_smiles(smiles).expect("rotatable-bond test SMILES should interpret");
    rotatable_bonds::detect(&molecule, options)
        .bond_ids()
        .to_vec()
}

#[test]
fn rdkit_strict_handles_empty_single_and_linear_molecules() {
    let empty = crate::core::MoleculeEditor::new();
    assert!(rotatable_bonds::detect(&empty, RotatableBondOptions::STRICT).is_empty());

    let mut builder = crate::core::MoleculeEditor::new();
    builder.add_atom(carbon()).expect("single carbon");
    let single = builder.finish().expect("single atom is connected");
    assert!(rotatable_bonds::detect(&single, RotatableBondOptions::STRICT).is_empty());

    assert!(detected("CC").is_empty());
    assert!(detected("CCC").is_empty());
    assert_eq!(detected("CCCC"), vec![BondId::new(1)]);
}

#[test]
fn result_is_self_describing_ordered_and_searchable() {
    let molecule = read_smiles("CCCCCC").expect("hexane should interpret");
    let result = rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT);

    assert_eq!(result.options(), RotatableBondOptions::STRICT);
    assert_eq!(
        result.bond_ids(),
        &[BondId::new(1), BondId::new(2), BondId::new(3),]
    );
    assert_eq!(result.len(), 3);
    assert!(!result.is_empty());
    assert!(result.contains(BondId::new(2)));
    assert!(!result.contains(BondId::new(0)));
    assert!(!result.contains(BondId::new(999)));
}

#[test]
fn rings_are_excluded_but_strict_ring_linkages_are_retained() {
    assert!(detected("C1CCCCC1").is_empty());
    assert_eq!(detected("C1CCCCC1-C1CCCCC1").len(), 1);
    assert_eq!(detected("c1ccccc1-c1ccccc1").len(), 1);
}

#[test]
fn strict_resonance_exclusions_keep_only_unrestricted_neighboring_axes() {
    assert!(detected("CC(=O)NC").is_empty());
    assert!(detected("CC(=O)OC").is_empty());
    assert!(detected("CC(=S)SC").is_empty());
    assert_eq!(detected("CC(=O)OCC"), vec![BondId::new(3)]);
    assert_eq!(detected("CC(=O)NCC"), vec![BondId::new(3)]);
    assert!(detected("CC(=[NH2+])NC").is_empty());
}

#[test]
fn localized_aromatic_ring_bonds_do_not_create_false_resonance_exclusions() {
    let molecule = read_smiles("CC1=NC(=NC(=N1)NC(C)C)NCC(C)C").expect("valid aminopyrimidine");
    let graph = molecule;
    let expected_endpoints = [(3, 11), (5, 7), (7, 8), (11, 12), (12, 13)];
    let mut expected = expected_endpoints.map(|(left, right)| {
        graph
            .bond_between(AtomId::new(left), AtomId::new(right))
            .expect("valid atom IDs")
            .expect("expected bond")
    });
    expected.sort_unstable();

    assert_eq!(
        rotatable_bonds::detect(&graph, RotatableBondOptions::STRICT).bond_ids(),
        &expected
    );
}

#[test]
fn localized_five_member_aromatic_rings_do_not_hide_exocyclic_axes() {
    let molecule = read_smiles("CCN1C=CN=C1[N+](=O)[O-]").expect("valid nitroimidazole");
    let graph = molecule;
    let mut expected = [(1, 2), (6, 7)].map(|(left, right)| {
        graph
            .bond_between(AtomId::new(left), AtomId::new(right))
            .expect("valid atom IDs")
            .expect("expected bond")
    });
    expected.sort_unstable();

    assert_eq!(
        rotatable_bonds::detect(&graph, RotatableBondOptions::STRICT).bond_ids(),
        &expected
    );
}

#[test]
fn strict_terminal_triple_and_symmetric_group_exclusions_are_applied() {
    assert!(detected("CCC#N").is_empty());
    assert!(detected("FC(F)(F)CC").is_empty());
    assert!(detected("ClC(Cl)(Cl)CC").is_empty());
    assert!(detected("BrC(Br)(Br)CC").is_empty());
    assert!(detected("CC(C)(C)CC").is_empty());
    assert_eq!(detected("CC(C)CC").len(), 1);
}

#[test]
fn unsupported_focus_bond_orders_are_never_rotatable() {
    for order in [
        BondOrder::Zero,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
        BondOrder::Dative,
    ] {
        let mut builder = crate::core::MoleculeEditor::new();
        let atoms = (0..4)
            .map(|_| builder.add_atom(carbon()).expect("carbon atom"))
            .collect::<Vec<_>>();
        builder
            .add_bond(atoms[0], atoms[1], BondOrder::Single)
            .expect("left terminal bond");
        builder
            .add_bond(atoms[1], atoms[2], order)
            .expect("focus bond");
        builder
            .add_bond(atoms[2], atoms[3], BondOrder::Single)
            .expect("right terminal bond");
        let molecule = builder.finish().expect("chain should be connected");

        assert!(
            rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT).is_empty(),
            "{order:?} focus bond should not be rotatable"
        );
    }
}

#[test]
fn detection_is_hydrogen_invariant_and_survives_tombstoned_slots() {
    let mut molecule = read_smiles("CCOCC").expect("ether should interpret");
    let expected = detected("CCOCC");
    assert_eq!(expected, vec![BondId::new(1), BondId::new(2)]);

    molecule.perceive().expect("ether should perceive");
    molecule
        .add_hydrogens()
        .expect("ether hydrogens should materialize");
    assert_eq!(
        rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT).bond_ids(),
        expected
    );

    molecule
        .perceive()
        .expect("materialized ether should reperceive");
    molecule
        .remove_hydrogens()
        .expect("ether hydrogens should collapse");
    assert_eq!(
        rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT).bond_ids(),
        expected
    );
}

#[test]
fn detection_reuses_or_computes_rings_without_mutating_perception() {
    let mut molecule = read_smiles("c1ccccc1-CCCC").expect("phenylbutane should interpret");
    let before = molecule.perception().clone();
    let detached = rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT);
    assert_eq!(molecule.perception(), &before);
    assert!(!molecule.perception().has_rings());

    crate::perception::rings::perceive_ring_membership(&mut molecule);
    assert!(molecule.perception().has_rings());
    let installed = rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT);
    assert_eq!(installed, detached);
}

#[test]
fn general_enables_every_optional_category() {
    assert_eq!(
        RotatableBondOptions::GENERAL,
        RotatableBondOptions {
            include_terminal_bonds: true,
            include_resonance_restricted_bonds: true,
            include_symmetric_bonds: true,
            include_ring_bonds: true,
        }
    );
    assert_eq!(
        detected_with_options("CC", RotatableBondOptions::GENERAL),
        vec![BondId::new(0)]
    );
    assert_eq!(
        detected_with_options("CC(=O)NC", RotatableBondOptions::GENERAL),
        vec![BondId::new(0), BondId::new(2), BondId::new(3)]
    );
    assert_eq!(
        detected_with_options("C1CCCCC1", RotatableBondOptions::GENERAL),
        (0..6).map(BondId::new).collect::<Vec<_>>()
    );
}

#[test]
fn each_option_independently_relaxes_strict_detection() {
    assert_eq!(
        detected_with_options(
            "CC",
            RotatableBondOptions {
                include_terminal_bonds: true,
                ..RotatableBondOptions::STRICT
            },
        ),
        vec![BondId::new(0)]
    );
    assert_eq!(
        detected_with_options(
            "CC(=O)NC",
            RotatableBondOptions {
                include_resonance_restricted_bonds: true,
                ..RotatableBondOptions::STRICT
            },
        ),
        vec![BondId::new(2)]
    );
    assert_eq!(
        detected_with_options(
            "FC(F)(F)CC",
            RotatableBondOptions {
                include_symmetric_bonds: true,
                ..RotatableBondOptions::STRICT
            },
        ),
        vec![BondId::new(3)]
    );
    assert_eq!(
        detected_with_options(
            "C1CCCCC1",
            RotatableBondOptions {
                include_ring_bonds: true,
                ..RotatableBondOptions::STRICT
            },
        ),
        (0..6).map(BondId::new).collect::<Vec<_>>()
    );
}

#[test]
fn resonance_filter_also_applies_when_ring_bonds_are_enabled() {
    let molecule = read_smiles("O=C1NCCC1").expect("lactam should interpret");
    let amide = molecule
        .bond_between(AtomId::new(1), AtomId::new(2))
        .expect("valid atom IDs")
        .expect("cyclic amide bond");
    let rings_without_restricted = rotatable_bonds::detect(
        &molecule,
        RotatableBondOptions {
            include_ring_bonds: true,
            ..RotatableBondOptions::STRICT
        },
    );

    assert!(!rings_without_restricted.contains(amide));
    assert!(rotatable_bonds::detect(&molecule, RotatableBondOptions::GENERAL).contains(amide));
}

#[test]
fn general_still_excludes_hydrogen_axes_and_unsupported_orders() {
    let molecule = read_smiles("[H]C").expect("methane fragment should interpret");
    assert!(rotatable_bonds::detect(&molecule, RotatableBondOptions::GENERAL).is_empty());

    let molecule = read_smiles("CC#N").expect("acetonitrile should interpret");
    assert!(rotatable_bonds::detect(&molecule, RotatableBondOptions::GENERAL).is_empty());
}
