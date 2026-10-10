use super::*;
use crate::properties::{PropertyKey, PropertyValue};

// Defensive algorithm validation must also handle invalid internal staging;
// public insertion rejects these malformed assertions before allocating a slot.
fn aromatic_atom(molecule: &Molecule, atom: AtomId) -> bool {
    molecule.atom_is_aromatic(atom).expect("atom exists") == Some(true)
}

fn aromatic_bond(molecule: &Molecule, bond: BondId) -> bool {
    molecule.bond_is_aromatic(bond).expect("bond exists") == Some(true)
}

fn fully_perceived_aromatic_stereo_fixture() -> Molecule {
    let mut molecule =
        read_smiles("c1ccccc1[C@H](F)Cl").expect("aromatic stereo fixture should parse");
    perceive(&mut molecule).expect("default perception should succeed");
    let report =
        stereo_api::assign_cip_descriptors(&mut molecule).expect("CIP assignment should succeed");
    assert!(!report.assigned.is_empty());
    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_rings());
    assert!(molecule.perception().has_aromaticity());
    assert!(molecule.perception().has_stereo());
    molecule
}

#[test]
fn ring_membership_empty_and_linear_molecules_have_no_rings() {
    let mut empty = crate::core::MoleculeEditor::new();
    let empty_membership = rings_api::perceive_ring_membership(empty.working_mut());
    assert!(empty_membership.ring_atom_ids().next().is_none());
    assert!(empty_membership.ring_bond_ids().next().is_none());

    let mut chain = crate::core::MoleculeEditor::new();
    let a = chain.add_atom(carbon()).expect("atom identifier capacity");
    let b = chain.add_atom(carbon()).expect("atom identifier capacity");
    let c = chain.add_atom(carbon()).expect("atom identifier capacity");
    let ab = chain
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    let bc = chain
        .add_bond(b, c, BondOrder::Single)
        .expect("bond should be valid");
    let chain_membership = rings_api::perceive_ring_membership(chain.working_mut());

    assert!(!chain_membership.atom_in_ring(a));
    assert!(!chain_membership.atom_in_ring(b));
    assert!(!chain_membership.bond_in_ring(ab));
    assert!(!chain_membership.bond_in_ring(bc));
    assert!(chain.perception().has_rings());
}

#[test]
fn ring_membership_reperception_preserves_valence_and_clears_downstream_sections() {
    let mut molecule = fully_perceived_aromatic_stereo_fixture();
    let valence = molecule
        .perception()
        .valence_state()
        .expect("installed valence")
        .clone();

    let membership = rings_api::perceive_ring_membership(&mut molecule);

    let perception = molecule.perception();
    assert_eq!(perception.valence_state(), Some(&valence));
    let rings = perception.ring_state().expect("membership installed");
    assert_eq!(rings.membership(), &membership);
    assert!(rings.basis().is_none());
    assert!(!perception.has_aromaticity());
    assert!(!perception.has_stereo());
}

#[test]
fn ring_basis_reperception_preserves_valence_and_clears_downstream_sections() {
    let mut molecule = fully_perceived_aromatic_stereo_fixture();
    let valence = molecule
        .perception()
        .valence_state()
        .expect("installed valence")
        .clone();

    let ring_set =
        rings_api::perceive_ring_set(&mut molecule).expect("ring basis perception should succeed");

    let perception = molecule.perception();
    assert_eq!(perception.valence_state(), Some(&valence));
    assert_eq!(perception.ring_set(), Some(&ring_set));
    assert_eq!(
        perception.ring_basis_model(),
        Some(RingBasisModel::FiguerasSssrLike)
    );
    assert!(!perception.has_aromaticity());
    assert!(!perception.has_stereo());
}

#[test]
fn implicit_hydrogen_update_preserves_rings_and_rebuilds_downstream_sections() {
    let mut molecule = fully_perceived_aromatic_stereo_fixture();
    let rings = molecule
        .perception()
        .ring_state()
        .expect("installed rings")
        .clone();
    let atom = AtomId::new(0);
    assert_eq!(molecule.perception().inferred_hydrogens(atom), Some(1));

    molecule.set_inferred_hydrogens(atom, 0);

    let perception = molecule.perception();
    assert_eq!(perception.inferred_hydrogens(atom), Some(0));
    assert_eq!(perception.ring_state(), Some(&rings));
    assert!(!perception.has_aromaticity());
    assert!(!perception.has_stereo());

    aromaticity_api::perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
        .expect("aromaticity should rebuild from retained valence and rings");
    assert!(molecule.perception().has_aromaticity());
    assert!(!molecule.perception().has_stereo());
    stereo_api::assign_cip_descriptors(&mut molecule)
        .expect("CIP should rebuild after aromaticity");
    assert!(molecule.perception().has_stereo());
    assert!(molecule.perception().has_cip_descriptors());
    assert_eq!(molecule.perception().ring_state(), Some(&rings));
}

#[test]
fn ring_membership_marks_triangle_atoms_and_bonds() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");

    let membership = rings_api::perceive_ring_membership(mol.working_mut());

    assert_eq!(sorted_atom_ids(membership.ring_atom_ids()), vec![a, b, c]);
    assert_eq!(
        sorted_bond_ids(membership.ring_bond_ids()),
        vec![ab, bc, ca]
    );
}

#[test]
fn ring_membership_excludes_tail_from_ring() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let tail = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
    let tail_bond = mol.add_bond(c, tail, BondOrder::Single).expect("bond");

    let membership = rings_api::perceive_ring_membership(mol.working_mut());

    assert_eq!(sorted_atom_ids(membership.ring_atom_ids()), vec![a, b, c]);
    assert_eq!(
        sorted_bond_ids(membership.ring_bond_ids()),
        vec![ab, bc, ca]
    );
    assert!(!membership.atom_in_ring(tail));
    assert!(!membership.bond_in_ring(tail_bond));
}

#[test]
fn ring_membership_handles_fused_rings_with_an_acyclic_tail() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let d = mol.add_atom(carbon()).expect("atom identifier capacity");
    let tail_a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let tail_b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
    let cd = mol.add_bond(c, d, BondOrder::Single).expect("bond");
    let da = mol.add_bond(d, a, BondOrder::Single).expect("bond");
    let linker = mol
        .add_bond(d, tail_a, BondOrder::Single)
        .expect("tail linker");
    let bridge = mol
        .add_bond(tail_a, tail_b, BondOrder::Single)
        .expect("bond");

    let membership = rings_api::perceive_ring_membership(mol.working_mut());

    assert_eq!(
        sorted_atom_ids(membership.ring_atom_ids()),
        vec![a, b, c, d]
    );
    assert_eq!(
        sorted_bond_ids(membership.ring_bond_ids()),
        vec![ab, bc, ca, cd, da]
    );
    assert!(!membership.bond_in_ring(linker));
    assert!(!membership.bond_in_ring(bridge));
}

#[test]
fn ring_membership_ignores_deleted_bonds_and_becomes_stale_after_mutation() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let ab = mol.add_bond(a, b, BondOrder::Single).expect("bond");
    let bc = mol.add_bond(b, c, BondOrder::Single).expect("bond");
    let ca = mol.add_bond(c, a, BondOrder::Single).expect("bond");
    mol.delete_bond(ca).expect("bond should delete");

    let membership = rings_api::perceive_ring_membership(mol.working_mut());
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
    assert!(atoms.iter().all(|atom| aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| aromatic_bond(&mol, *bond)));
}

#[test]
fn discrete_chemical_perception_changes_only_perception_state() {
    let mut molecule =
        read_smiles("F[C@](Cl)(Br)c1cc[nH]c1").expect("heteroaromatic stereo fixture should parse");

    let atom_ids = molecule.atom_ids().collect::<Vec<_>>();
    let annotated_atom = atom_ids[0];
    let annotated_bond = molecule.bond_ids().next().expect("fixture bond");
    molecule
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("perception_purity_fixture").unwrap(),
            PropertyValue::String("molecule property".to_owned()),
        )
        .unwrap();
    molecule
        .properties_mut()
        .atoms_mut()
        .set_value(
            PropertyKey::new("atom_note").unwrap(),
            annotated_atom,
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    molecule
        .properties_mut()
        .bonds_mut()
        .set_value(
            PropertyKey::new("bond_note").unwrap(),
            annotated_bond,
            Some(PropertyValue::Int(7)),
        )
        .unwrap();

    let stereo_element = molecule
        .stereo_element_ids()
        .next()
        .expect("direct SMILES stereo element");
    molecule
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![stereo_element],
        })
        .expect("valid absolute stereo group");

    assert_eq!(molecule.perception(), &Perception::default());
    let represented_before = represented_molecule_snapshot(&molecule);

    molecule.perceive().expect("default perception");

    assert_eq!(represented_molecule_snapshot(&molecule), represented_before);
    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_rings());
    assert!(molecule.perception().has_aromaticity());
    assert_eq!(molecule.stereo_elements().count(), 1);
    assert_eq!(molecule.stereo_groups().count(), 1);
}

#[test]
fn molecule_perception_queries_read_the_installed_state_directly() {
    let mut molecule =
        read_smiles("F[C@](Cl)(Br)c1cc[nH]c1").expect("stereo aromatic fixture should parse");
    perceive(&mut molecule).expect("fixture should perceive");
    stereo_api::assign_cip_descriptors(&mut molecule).expect("CIP assignment");

    let graph = molecule;
    for atom in graph.atom_ids() {
        assert_eq!(
            graph.perception().inferred_hydrogens(atom),
            graph.perception().inferred_hydrogens(atom)
        );
        assert_eq!(
            graph.atom_is_aromatic(atom).expect("live atom"),
            graph.perception().atom_is_aromatic(atom)
        );
    }
    for bond in graph.bond_ids() {
        assert_eq!(
            graph.bond_is_aromatic(bond).expect("live bond"),
            graph.perception().bond_is_aromatic(bond)
        );
    }
    for element in graph.stereo_element_ids() {
        assert_eq!(
            graph.cip_descriptor(element).expect("live stereo element"),
            graph.perception().cip_descriptor(element)
        );
    }
}

#[test]
fn default_perception_accepts_aromatic_source_localized_by_interpretation() {
    let mut molecule = read_smiles("c1ccccc1").expect("benzene should parse");

    molecule
        .perceive()
        .expect("localized aromatic source should perceive directly");

    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_rings());
    assert!(molecule.perception().has_aromaticity());
}

#[test]
fn default_perception_rolls_back_when_ring_perception_fails_after_valence() {
    const ATOM_COUNT: usize = 4_097;

    let mut molecule = crate::core::MoleculeEditor::new();
    let atoms = (0..ATOM_COUNT)
        .map(|_| {
            molecule
                .add_atom(carbon())
                .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for index in 0..ATOM_COUNT {
        molecule
            .add_bond(
                atoms[index],
                atoms[(index + 1) % ATOM_COUNT],
                BondOrder::Single,
            )
            .expect("large ring bond");
    }
    rings_api::perceive_ring_membership(molecule.working_mut());
    let original = molecule.clone();

    let error = molecule
        .working_mut()
        .perceive()
        .expect_err("default ring cycle-size limit must fail");

    assert!(matches!(
        error,
        perception_api::PerceptionError::Rings(RingPerceptionError::ResourceLimit {
            resource: "cycle size",
            observed: ATOM_COUNT,
            limit: 4_096,
        })
    ));
    assert_eq!(molecule, original);
}

#[test]
fn default_perceive_installs_only_valence_rings_and_aromaticity() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");

    molecule.perceive().expect("ethanol should perceive");

    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_rings());
    assert!(molecule.perception().has_aromaticity());
    assert_eq!(
        molecule.perception().ring_basis_model(),
        Some(RingBasisModel::FiguerasSssrLike)
    );
    assert!(!molecule.perception().has_stereo());
}

#[test]
fn interpretation_canonicalizes_source_stereo_before_perception() {
    let document = smiles_api::parse_str("F/C=C/c1ccccc1").expect("SMILES parses");
    let interpretation = smiles_api::interpret(&document).expect("SMILES interprets");
    let (mut molecule, report) = interpretation.into_parts().expect("one component");

    assert_eq!(report.created_stereo_elements().len(), 1);
    assert!(molecule
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert_eq!(molecule.perception(), &Perception::default());

    let represented_before = molecule.clone();
    molecule.perceive().expect("default perception succeeds");
    assert_eq!(
        molecule.atoms().collect::<Vec<_>>(),
        represented_before.atoms().collect::<Vec<_>>()
    );
    assert_eq!(
        molecule.bonds().collect::<Vec<_>>(),
        represented_before.bonds().collect::<Vec<_>>()
    );
    assert_eq!(
        molecule.stereo_elements().collect::<Vec<_>>(),
        represented_before.stereo_elements().collect::<Vec<_>>()
    );
    assert!(molecule.perception().has_valence());
    assert!(molecule.perception().has_rings());
    assert!(molecule.perception().has_aromaticity());
    assert!(!molecule.perception().has_stereo());
}

#[test]
fn perceive_does_not_infer_or_materialize_coordinate_stereo() {
    let (graph, _center, _carriers, _) = tetrahedral_marked_graph();
    let positions = test_positions(vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Point3::new(0.0, 0.0, -1.0),
    ]);
    let mut molecule = graph.finish().expect("coordinate fixture publishes");

    molecule
        .perceive()
        .expect("ordinary workflow should succeed");

    assert!(molecule.stereo_elements().next().is_none());
    assert_eq!(
        stereo_api::infer_coordinate_stereo(&molecule, &positions)
            .expect("separate coordinate inference")
            .elements
            .len(),
        1
    );
}

#[test]
fn perceive_rolls_back_failure_without_rewriting_canonical_representation() {
    let mut graph = crate::core::MoleculeEditor::new();
    let chlorine = graph
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let oxo = graph.add_atom(oxygen()).expect("atom identifier capacity");
    let hydroxyl = graph.add_atom(oxygen()).expect("atom identifier capacity");
    let carbon = graph.add_atom(carbon()).expect("atom identifier capacity");
    graph
        .add_bond(chlorine, oxo, BondOrder::Double)
        .expect("oxo bond");
    graph
        .add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("hydroxyl bond");
    graph
        .add_bond(chlorine, carbon, BondOrder::Single)
        .expect("connecting bond");
    for symbol in ["F", "F", "F", "F"] {
        let fluorine = graph
            .add_atom(element_atom(symbol))
            .expect("atom identifier capacity");
        graph
            .add_bond(carbon, fluorine, BondOrder::Single)
            .expect("pentavalent carbon bond");
    }
    let mut molecule = graph;
    molecule
        .working_mut()
        .canonicalize_fixture()
        .expect("fixture canonicalization should succeed");
    let before = molecule.clone();

    let error = molecule
        .working_mut()
        .perceive()
        .expect_err("default perception should reject pentavalent carbon");

    assert!(matches!(error, perception_api::PerceptionError::Valence(_)));
    assert_eq!(molecule, before);
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
        .all(|atom| aromatic_atom(&ten_member, *atom)));
    assert!(ten_bonds
        .iter()
        .all(|bond| aromatic_bond(&ten_member, *bond)));

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
        .all(|atom| !aromatic_atom(&twelve_member, *atom)));
    assert!(twelve_bonds
        .iter()
        .all(|bond| !aromatic_bond(&twelve_member, *bond)));
}

#[test]
fn aromaticity_leaves_cyclohexane_and_cyclobutadiene_non_aromatic() {
    let (mut cyclohexane, atoms, bonds) =
        ring_molecule(&["C", "C", "C", "C", "C", "C"], &[BondOrder::Single; 6]);
    aromaticity_api::perceive_aromaticity(&mut cyclohexane, AromaticityModel::RdkitLike)
        .expect("cyclohexane should be supported");
    assert!(atoms.iter().all(|atom| !aromatic_atom(&cyclohexane, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&cyclohexane, *bond)));

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
        .all(|atom| !aromatic_atom(&cyclobutadiene, *atom)));
    assert!(bonds
        .iter()
        .all(|bond| !aromatic_bond(&cyclobutadiene, *bond)));
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

    assert!(atoms.iter().all(|atom| aromatic_atom(&furan_like, *atom)));
    assert!(bonds.iter().all(|bond| aromatic_bond(&furan_like, *bond)));
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
        nitrogen.hydrogens = ImplicitHydrogens::Fixed(1);
    }
    pyrrole_like.set_inferred_hydrogens(atoms[0], 0);

    aromaticity_api::perceive_aromaticity(&mut pyrrole_like, AromaticityModel::RdkitLike)
        .expect("pyrrole-like ring should be supported");

    assert!(atoms.iter().all(|atom| aromatic_atom(&pyrrole_like, *atom)));
    assert!(bonds.iter().all(|bond| aromatic_bond(&pyrrole_like, *bond)));
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
        .all(|atom| aromatic_atom(&phosphole_like, *atom)));
    assert!(bonds
        .iter()
        .all(|bond| aromatic_bond(&phosphole_like, *bond)));
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
    let methyl = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(atoms[0], methyl, BondOrder::Single)
        .expect("phosphorus substituent bond");

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("hypervalent phosphorus ring should be supported");

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
    assert!(!aromatic_atom(&mol, methyl));
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
        .radical = AtomRadical::new(1, Some(2));

    aromaticity_api::perceive_aromaticity(&mut neutral_carbon_radical, AromaticityModel::RdkitLike)
        .expect("neutral carbon radical ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| aromatic_atom(&neutral_carbon_radical, *atom)));

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
        .radical = AtomRadical::new(1, Some(2));

    aromaticity_api::perceive_aromaticity(&mut oxygen_radical, AromaticityModel::RdkitLike)
        .expect("heteroatom radical ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !aromatic_atom(&oxygen_radical, *atom)));

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
        atom.radical = AtomRadical::new(1, Some(2));
    }

    aromaticity_api::perceive_aromaticity(&mut charged_carbon_radical, AromaticityModel::RdkitLike)
        .expect("charged carbon radical ring should be supported");

    assert!(atoms
        .iter()
        .all(|atom| !aromatic_atom(&charged_carbon_radical, *atom)));
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
    let methyl_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let methyl_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(atoms[0], methyl_a, BondOrder::Single)
        .expect("first substituent bond");
    mol.add_bond(atoms[0], methyl_b, BondOrder::Single)
        .expect("second substituent bond");

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("tetracoordinate ring atom should be supported");

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
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
        nitrogen.hydrogens = ImplicitHydrogens::Fixed(1);
    }
    mol.set_inferred_hydrogens(atoms[0], 0);

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("protonated saturated ring nitrogen should be supported");

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
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
        cation.hydrogens = ImplicitHydrogens::Fixed(1);
    }
    mol.set_inferred_hydrogens(atoms[0], 0);

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclopropenyl cation should be supported");

    assert!(atoms.iter().all(|atom| aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| aromatic_bond(&mol, *bond)));
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
        saturated.hydrogens = ImplicitHydrogens::Fixed(2);
    }
    mol.set_inferred_hydrogens(atoms[0], 0);

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("over-valent candidate rejection should be supported");

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
}

#[test]
fn aromaticity_marks_azulene_fused_perimeter_but_not_shared_bond() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atoms = (0..10)
        .map(|_| mol.add_atom(carbon()).expect("atom identifier capacity"))
        .collect::<Vec<_>>();
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

    aromaticity_api::perceive_aromaticity(mol.working_mut(), AromaticityModel::RdkitLike)
        .expect("azulene-like fused system should be supported");

    assert!(atoms.iter().all(|atom| aromatic_atom(mol.working(), *atom)));
    assert!(perimeter_bonds
        .iter()
        .all(|bond| aromatic_bond(mol.working(), *bond)));
    assert!(!aromatic_bond(mol.working(), shared));
}

#[test]
fn aromaticity_keeps_aromatic_heteroring_bond_shared_with_saturated_ring() {
    let mut mol = crate::core::MoleculeEditor::new();
    let c0 = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c1 = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c2 = mol.add_atom(carbon()).expect("atom identifier capacity");
    let n3 = mol
        .add_atom(Atom::new(
            Element::from_symbol("N").expect("nitrogen should be available"),
        ))
        .expect("atom identifier capacity");
    let n4 = mol
        .add_atom(Atom::new(
            Element::from_symbol("N").expect("nitrogen should be available"),
        ))
        .expect("atom identifier capacity");
    let saturated_a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let saturated_b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let saturated_c = mol.add_atom(carbon()).expect("atom identifier capacity");

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

    aromaticity_api::perceive_aromaticity(mol.working_mut(), AromaticityModel::RdkitLike)
        .expect("fused heteroaromatic ring should be supported");

    for bond_id in aromatic_bonds {
        assert!(
            aromatic_bond(mol.working(), bond_id),
            "aromatic ring bond {bond_id} should be aromatic"
        );
    }
    for bond_id in saturated_bonds {
        assert!(
            !aromatic_bond(mol.working(), bond_id),
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
    let hydrogen = mol
        .add_atom(Atom::new(
            Element::from_symbol("H").expect("hydrogen should be available"),
        ))
        .expect("atom identifier capacity");
    mol.add_bond(atoms[0], hydrogen, BondOrder::Single)
        .expect("explicit hydrogen bond should be valid");

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclopentadienyl anion should be supported");

    assert!(atoms.iter().all(|atom| aromatic_atom(&mol, *atom)));
    assert!(!aromatic_atom(&mol, hydrogen));
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

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
}

#[test]
fn aromaticity_uses_ring_membership_not_acyclic_double_bonds() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(a, b, BondOrder::Double).expect("bond");
    mol.add_bond(b, c, BondOrder::Single).expect("bond");

    aromaticity_api::perceive_aromaticity(mol.working_mut(), AromaticityModel::RdkitLike)
        .expect("acyclic molecule should be supported");

    assert!(!aromatic_atom(mol.working(), a));
    assert!(!aromatic_bond(mol.working(), BondId::new(0)));
}

#[test]
fn aromaticity_clears_existing_flags_before_assignment() {
    let (mut mol, atoms, bonds) =
        ring_molecule(&["C", "C", "C", "C", "C", "C"], &[BondOrder::Single; 6]);
    mol.begin_aromaticity(AromaticityModel::RdkitLike);
    for atom in &atoms {
        mol.set_atom_aromatic(*atom, true);
    }
    for bond in &bonds {
        mol.set_bond_aromatic(*bond, true);
    }

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("cyclohexane should be supported");

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
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

    mol.add_atom(oxygen()).expect("atom identifier capacity");
    assert!(!mol.perception().has_aromaticity());
    assert!(atoms
        .iter()
        .all(|atom| mol.atom_is_aromatic(*atom).expect("atom exists").is_none()));
}
