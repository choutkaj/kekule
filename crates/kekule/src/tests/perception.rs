use super::*;

fn aromatic_atom(molecule: &Molecule, atom: AtomId) -> bool {
    molecule.atom_is_aromatic(atom).expect("atom exists") == Some(true)
}

fn aromatic_bond(molecule: &Molecule, bond: BondId) -> bool {
    molecule.bond_is_aromatic(bond).expect("bond exists") == Some(true)
}

fn fully_perceived_aromatic_stereo_fixture() -> SmallMolecule {
    let mut molecule =
        read_smiles("c1ccccc1[C@H](F)Cl").expect("aromatic stereo fixture should parse");
    perceive(&mut molecule).expect("default perception should succeed");
    let report = stereo_api::assign_cip_descriptors(molecule.graph_mut())
        .expect("CIP assignment should succeed");
    assert!(!report.assigned.is_empty());
    assert!(molecule.graph().perception().has_valence());
    assert!(molecule.graph().perception().has_rings());
    assert!(molecule.graph().perception().has_aromaticity());
    assert!(molecule.graph().perception().has_stereo());
    molecule
}

#[test]
fn ring_membership_empty_and_linear_molecules_have_no_rings() {
    let mut empty = Molecule::new();
    let empty_membership = rings_api::perceive_ring_membership(&mut empty);
    assert!(empty_membership.ring_atom_ids().next().is_none());
    assert!(empty_membership.ring_bond_ids().next().is_none());

    let mut chain = Molecule::new();
    let a = chain.add_atom(carbon()).expect("atom identifier capacity");
    let b = chain.add_atom(carbon()).expect("atom identifier capacity");
    let c = chain.add_atom(carbon()).expect("atom identifier capacity");
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
fn ring_membership_reperception_preserves_valence_and_clears_downstream_sections() {
    let mut molecule = fully_perceived_aromatic_stereo_fixture();
    let valence = molecule
        .graph()
        .perception()
        .valence_state()
        .expect("installed valence")
        .clone();

    let membership = rings_api::perceive_ring_membership(molecule.graph_mut());

    let perception = molecule.graph().perception();
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
        .graph()
        .perception()
        .valence_state()
        .expect("installed valence")
        .clone();

    let ring_set = rings_api::perceive_ring_set(molecule.graph_mut())
        .expect("ring basis perception should succeed");

    let perception = molecule.graph().perception();
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
        .graph()
        .perception()
        .ring_state()
        .expect("installed rings")
        .clone();
    let atom = AtomId::new(0);
    assert_eq!(molecule.graph().implicit_hydrogens(atom), Ok(Some(1)));

    molecule.graph_mut().set_implicit_hydrogens(atom, 0);

    let perception = molecule.graph().perception();
    assert_eq!(perception.implicit_hydrogens(atom), Some(0));
    assert_eq!(perception.ring_state(), Some(&rings));
    assert!(!perception.has_aromaticity());
    assert!(!perception.has_stereo());

    aromaticity_api::perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
        .expect("aromaticity should rebuild from retained valence and rings");
    assert!(molecule.graph().perception().has_aromaticity());
    assert!(!molecule.graph().perception().has_stereo());
    stereo_api::assign_cip_descriptors(molecule.graph_mut())
        .expect("CIP should rebuild after aromaticity");
    assert!(molecule.graph().perception().has_stereo());
    assert!(molecule.graph().perception().has_cip_descriptors());
    assert_eq!(molecule.graph().perception().ring_state(), Some(&rings));
}

#[test]
fn ring_membership_marks_triangle_atoms_and_bonds() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
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
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let tail = mol.add_atom(oxygen()).expect("atom identifier capacity");
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
fn ring_membership_handles_fused_rings_with_an_acyclic_tail() {
    let mut mol = Molecule::new();
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

    let membership = rings_api::perceive_ring_membership(&mut mol);

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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
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
    assert!(atoms.iter().all(|atom| aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| aromatic_bond(&mol, *bond)));
}

#[test]
fn discrete_chemical_perception_changes_only_perception_state() {
    let mut molecule =
        read_smiles("F[C@](Cl)(Br)c1cc[nH]c1").expect("heteroaromatic stereo fixture should parse");

    let atom_ids = molecule.graph().atom_ids().collect::<Vec<_>>();
    let annotated_atom = atom_ids[0];
    let annotated_bond = molecule.graph().bond_ids().next().expect("fixture bond");
    molecule.graph_mut().props_mut().insert(
        "perception_purity_fixture".to_owned(),
        PropValue::String("molecule property".to_owned()),
    );
    molecule
        .graph_mut()
        .atom_props_mut(annotated_atom)
        .expect("fixture atom")
        .insert("atom_note".to_owned(), PropValue::Bool(true));
    molecule
        .graph_mut()
        .bond_props_mut(annotated_bond)
        .expect("fixture bond")
        .insert("bond_note".to_owned(), PropValue::Int(7));

    let stereo_element = molecule
        .graph()
        .stereo_element_ids()
        .next()
        .expect("direct SMILES stereo element");
    molecule
        .graph_mut()
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![stereo_element],
        })
        .expect("valid absolute stereo group");

    let mut conformer = Conformer::new(crate::units::ANGSTROM).expect("angstrom conformer");
    conformer.props_mut().insert(
        "conformer_note".to_owned(),
        PropValue::String("source coordinates".to_owned()),
    );
    for (index, atom) in atom_ids.iter().copied().enumerate() {
        conformer
            .set_position(
                atom,
                crate::units::Quantity::new(
                    Point3::new(index as f64, (index % 3) as f64, 0.25),
                    crate::units::ANGSTROM,
                ),
            )
            .expect("finite fixture coordinate");
    }
    molecule
        .graph_mut()
        .add_conformer(conformer)
        .expect("valid complete conformer");

    assert_eq!(molecule.graph().perception(), &PerceptionState::default());
    let represented_before = represented_molecule_snapshot(molecule.graph());

    perception_api::perceive(molecule.graph_mut()).expect("default perception");

    assert_eq!(
        represented_molecule_snapshot(molecule.graph()),
        represented_before
    );
    assert!(molecule.graph().perception().has_valence());
    assert!(molecule.graph().perception().has_rings());
    assert!(molecule.graph().perception().has_aromaticity());
    assert_eq!(molecule.graph().stereo_elements().count(), 1);
    assert_eq!(molecule.graph().stereo_groups().count(), 1);
}

#[test]
fn molecule_perception_queries_read_the_installed_state_directly() {
    let mut molecule =
        read_smiles("F[C@](Cl)(Br)c1cc[nH]c1").expect("stereo aromatic fixture should parse");
    perceive(&mut molecule).expect("fixture should perceive");
    stereo_api::assign_cip_descriptors(molecule.graph_mut()).expect("CIP assignment");

    let graph = molecule.graph();
    for atom in graph.atom_ids() {
        assert_eq!(
            graph.implicit_hydrogens(atom).expect("live atom"),
            graph.perception().implicit_hydrogens(atom)
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

    perception_api::perceive(molecule.graph_mut())
        .expect("localized aromatic source should perceive directly");

    assert!(molecule.graph().perception().has_valence());
    assert!(molecule.graph().perception().has_rings());
    assert!(molecule.graph().perception().has_aromaticity());
}

#[test]
fn default_perception_rolls_back_when_ring_perception_fails_after_valence() {
    const ATOM_COUNT: usize = 4_097;

    let mut molecule = SmallMolecule::default();
    let atoms = (0..ATOM_COUNT)
        .map(|_| {
            molecule
                .graph_mut()
                .add_atom(carbon())
                .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for index in 0..ATOM_COUNT {
        molecule
            .graph_mut()
            .add_bond(
                atoms[index],
                atoms[(index + 1) % ATOM_COUNT],
                BondOrder::Single,
            )
            .expect("large ring bond");
    }
    rings_api::perceive_ring_membership(molecule.graph_mut());
    let original = molecule.clone();

    let error = perception_api::perceive(molecule.graph_mut())
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

    assert!(molecule.graph().perception().has_valence());
    assert!(molecule.graph().perception().has_rings());
    assert!(molecule.graph().perception().has_aromaticity());
    assert_eq!(
        molecule.graph().perception().ring_basis_model(),
        Some(RingBasisModel::FiguerasSssrLike)
    );
    assert!(!molecule.graph().perception().has_stereo());
}

#[test]
fn interpretation_canonicalizes_source_stereo_before_perception() {
    let document = smiles_api::parse_str("F/C=C/c1ccccc1").expect("SMILES parses");
    let interpretation = smiles_api::interpret(&document).expect("SMILES interprets");
    let (mut molecule, report) = interpretation.into_parts().expect("one component");

    assert_eq!(report.created_stereo_elements().len(), 1);
    assert!(molecule
        .graph()
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert_eq!(molecule.graph().perception(), &PerceptionState::default());

    let represented_before = molecule.graph().clone();
    molecule.perceive().expect("default perception succeeds");
    assert_eq!(
        molecule.graph().atoms().collect::<Vec<_>>(),
        represented_before.atoms().collect::<Vec<_>>()
    );
    assert_eq!(
        molecule.graph().bonds().collect::<Vec<_>>(),
        represented_before.bonds().collect::<Vec<_>>()
    );
    assert_eq!(
        molecule.graph().stereo_elements().collect::<Vec<_>>(),
        represented_before.stereo_elements().collect::<Vec<_>>()
    );
    assert!(molecule.graph().perception().has_valence());
    assert!(molecule.graph().perception().has_rings());
    assert!(molecule.graph().perception().has_aromaticity());
    assert!(!molecule.graph().perception().has_stereo());
}

#[test]
fn perceive_does_not_infer_or_materialize_coordinate_stereo() {
    let (mut graph, center, carriers, _) = tetrahedral_marked_graph();
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    for (atom, point) in
        std::iter::once((center, Point3::new(0.0, 0.0, 0.0))).chain(carriers.iter().copied().zip([
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, -1.0),
        ]))
    {
        conformer
            .set_position(
                atom,
                crate::units::Quantity::new(point, crate::units::ANGSTROM),
            )
            .unwrap();
    }
    graph.add_conformer(conformer).expect("complete conformer");
    let mut molecule = SmallMolecule::from_graph(graph);

    molecule
        .perceive()
        .expect("ordinary workflow should succeed");

    assert!(molecule.graph().stereo_elements().next().is_none());
    assert_eq!(
        stereo_api::infer_coordinate_stereo(molecule.graph())
            .expect("separate coordinate inference")
            .elements
            .len(),
        1
    );
}

#[test]
fn perceive_rolls_back_failure_without_rewriting_canonical_representation() {
    let mut graph = Molecule::new();
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
    let mut molecule = SmallMolecule::from_graph(graph);
    molecule
        .canonicalize_fixture()
        .expect("fixture canonicalization should succeed");
    let before = molecule.clone();

    let error = molecule
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
        nitrogen.hydrogens = HydrogenDeclaration::Fixed(1);
    }
    pyrrole_like.set_implicit_hydrogens(atoms[0], 0);

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
        .radical = Some(AtomRadical::Doublet);

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
        .radical = Some(AtomRadical::Doublet);

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
        atom.radical = Some(AtomRadical::Doublet);
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
        nitrogen.hydrogens = HydrogenDeclaration::Fixed(1);
    }
    mol.set_implicit_hydrogens(atoms[0], 0);

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
        cation.hydrogens = HydrogenDeclaration::Fixed(1);
    }
    mol.set_implicit_hydrogens(atoms[0], 0);

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
        saturated.hydrogens = HydrogenDeclaration::Fixed(2);
    }
    mol.set_implicit_hydrogens(atoms[0], 0);

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("over-valent candidate rejection should be supported");

    assert!(atoms.iter().all(|atom| !aromatic_atom(&mol, *atom)));
    assert!(bonds.iter().all(|bond| !aromatic_bond(&mol, *bond)));
}

#[test]
fn aromaticity_marks_azulene_fused_perimeter_but_not_shared_bond() {
    let mut mol = Molecule::new();
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("azulene-like fused system should be supported");

    assert!(atoms.iter().all(|atom| aromatic_atom(&mol, *atom)));
    assert!(perimeter_bonds
        .iter()
        .all(|bond| aromatic_bond(&mol, *bond)));
    assert!(!aromatic_bond(&mol, shared));
}

#[test]
fn aromaticity_keeps_aromatic_heteroring_bond_shared_with_saturated_ring() {
    let mut mol = Molecule::new();
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

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("fused heteroaromatic ring should be supported");

    for bond_id in aromatic_bonds {
        assert!(
            aromatic_bond(&mol, bond_id),
            "aromatic ring bond {bond_id} should be aromatic"
        );
    }
    for bond_id in saturated_bonds {
        assert!(
            !aromatic_bond(&mol, bond_id),
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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(a, b, BondOrder::Double).expect("bond");
    mol.add_bond(b, c, BondOrder::Single).expect("bond");

    aromaticity_api::perceive_aromaticity(&mut mol, AromaticityModel::RdkitLike)
        .expect("acyclic molecule should be supported");

    assert!(!aromatic_atom(&mol, a));
    assert!(!aromatic_bond(&mol, BondId::new(0)));
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

#[test]
fn stereo_validation_reports_invalid_local_elements_without_mutating() {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let a = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let b = mol
        .add_atom(element_atom("N"))
        .expect("atom identifier capacity");
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
                orientation: None,
            }),
            group: None,
        })
        .expect("stereo element");
    mark_all_fresh(&mut mol);

    let error = stereo_api::validate_stereo(&mol).expect_err("invalid stored stereo");

    assert!(mol.stereo_elements().next().is_some());
    assert!(error
        .issues
        .contains(&StereoValidationIssue::InvalidTetrahedralCarrierCount {
            element,
            center,
            carrier_count: 3,
        }));
    assert!(error
        .issues
        .contains(&StereoValidationIssue::DuplicateTetrahedralCarrier {
            element,
            center,
            carrier: StereoCarrier::Atom(a),
        }));
    assert!(error
        .issues
        .contains(&StereoValidationIssue::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier: StereoCarrier::Atom(b),
        }));
    assert!(mol
        .stereo_element(element)
        .expect("element")
        .is_explicitly_unknown());
}

#[test]
fn stereo_validation_checks_implicit_carrier_form_without_perception_state() {
    let mut tetrahedral = Molecule::new();
    let center = tetrahedral
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let mut atom_carriers = Vec::new();
    for symbol in ["F", "Cl", "Br"] {
        let carrier = tetrahedral
            .add_atom(element_atom(symbol))
            .expect("atom identifier capacity");
        tetrahedral
            .add_bond(center, carrier, BondOrder::Single)
            .expect("carrier bond");
        atom_carriers.push(StereoCarrier::Atom(carrier));
    }
    let mut carriers = atom_carriers.clone();
    carriers.push(StereoCarrier::ImplicitHydrogen);
    let hydrogen_element = tetrahedral
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers,
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("tetrahedral stereo element");

    stereo_api::validate_stereo(&tetrahedral)
        .expect("implicit hydrogen availability is chemically interpretive");

    tetrahedral
        .remove_stereo_element(hydrogen_element)
        .expect("remove hydrogen-carrier element");
    atom_carriers.push(StereoCarrier::ImplicitLonePair);
    tetrahedral
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: atom_carriers,
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("tetrahedral stereo element");
    stereo_api::validate_stereo(&tetrahedral)
        .expect("implicit lone-pair availability is chemically interpretive");

    let mut double_bond = Molecule::new();
    let left = double_bond
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right = double_bond
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let bond = double_bond
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let double_element = double_bond
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond,
                left,
                right,
                left_carrier: StereoCarrier::ImplicitHydrogen,
                right_carrier: StereoCarrier::ImplicitLonePair,
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .expect("double-bond stereo element");

    let error = stereo_api::validate_stereo(&double_bond)
        .expect_err("unavailable double-bond carriers should be reported");
    assert_eq!(
        error.issues,
        vec![StereoValidationIssue::UnsupportedDoubleBondCarrier {
            element: double_element,
            endpoint: right,
            carrier: StereoCarrier::ImplicitLonePair,
        }]
    );

    let mut axis = Molecule::new();
    let axis_left = axis.add_atom(carbon()).expect("atom identifier capacity");
    let axis_right = axis.add_atom(carbon()).expect("atom identifier capacity");
    let axis_bond = axis
        .add_bond(axis_left, axis_right, BondOrder::Single)
        .expect("axis bond");
    let axis_element = axis
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis: axis_bond,
            carriers: vec![
                StereoCarrier::ImplicitHydrogen,
                StereoCarrier::ImplicitLonePair,
            ],
            orientation: Some(AxisOrientation::Clockwise),
        })))
        .expect("axis stereo element");

    let error = stereo_api::validate_stereo(&axis)
        .expect_err("implicit axis carriers should be unsupported");
    assert_eq!(
        error.issues,
        vec![
            StereoValidationIssue::UnsupportedAxisCarrier {
                element: axis_element,
                axis: axis_bond,
                carrier: StereoCarrier::ImplicitHydrogen,
            },
            StereoValidationIssue::UnsupportedAxisCarrier {
                element: axis_element,
                axis: axis_bond,
                carrier: StereoCarrier::ImplicitLonePair,
            },
        ]
    );
}

#[test]
fn stereo_candidates_use_normalized_and_perceived_hydrogen_state_without_cip_assignment() {
    let mut molecule = read_smiles("CC(F)(Cl)Br").expect("smiles should parse");
    perceive(&mut molecule).expect("molecule should perceive");

    stereo_api::validate_stereo(molecule.graph()).expect("stored stereo should be valid");
    let candidates = stereo_api::detect_stereo_candidates(molecule.graph());

    assert!(candidates.iter().any(|candidate| matches!(
        candidate,
        StereoCandidate::Tetrahedral { center, carriers }
            if *center == AtomId::new(1)
                && carriers.len() == 4
                && !carriers.contains(&StereoCarrier::ImplicitHydrogen)
    )));
    assert!(molecule.graph().stereo_elements().next().is_none());
}

#[test]
fn interpretation_assembles_paired_directional_marks_into_double_bond_element() {
    let (molecule, report) =
        read_smiles_with_report("C/C=C\\F").expect("directional smiles should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    assert!(molecule.graph().stereo_elements().next().is_some());
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    match &element.kind {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.bond, BondId::new(1));
            assert_eq!(stereo.left, AtomId::new(1));
            assert_eq!(stereo.right, AtomId::new(2));
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(AtomId::new(0)));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(AtomId::new(3)));
            assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Together));
        }
        other => panic!("expected double-bond stereo, found {other:?}"),
    }
}

#[test]
fn equivalent_smiles_direction_tokens_publish_equivalent_canonical_stereo() {
    for (first, second) in [("C/C=C/C", r"C\C=C\C"), (r"C/C=C\C", r"C\C=C/C")] {
        let first = read_smiles(first).expect("first directional spelling should interpret");
        let second = read_smiles(second).expect("second directional spelling should interpret");
        let first = first
            .graph()
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();
        let second = second
            .graph()
            .stereo_elements()
            .map(|(_, element)| element.clone())
            .collect::<Vec<_>>();

        assert_eq!(first, second);
        assert_eq!(first.len(), 1);
    }
}

#[test]
fn alternate_directional_source_carriers_publish_identical_double_bond_stereo() {
    let source_graph = || {
        let mut molecule = Molecule::new();
        let left = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        let right = molecule
            .add_atom(carbon())
            .expect("atom identifier capacity");
        let left_reference = molecule
            .add_atom(element_atom("F"))
            .expect("atom identifier capacity");
        let left_alternative = molecule
            .add_atom(element_atom("Cl"))
            .expect("atom identifier capacity");
        let right_reference = molecule
            .add_atom(element_atom("Br"))
            .expect("atom identifier capacity");
        let right_alternative = molecule
            .add_atom(element_atom("I"))
            .expect("atom identifier capacity");
        molecule
            .add_bond(left, right, BondOrder::Double)
            .expect("double bond");
        let left_reference_bond = molecule
            .add_bond(left, left_reference, BondOrder::Single)
            .expect("left reference bond");
        let left_alternative_bond = molecule
            .add_bond(left, left_alternative, BondOrder::Single)
            .expect("left alternative bond");
        let right_reference_bond = molecule
            .add_bond(right, right_reference, BondOrder::Single)
            .expect("right reference bond");
        let right_alternative_bond = molecule
            .add_bond(right, right_alternative, BondOrder::Single)
            .expect("right alternative bond");
        (
            molecule,
            left,
            right,
            left_reference_bond,
            left_alternative_bond,
            right_reference_bond,
            right_alternative_bond,
        )
    };

    let canonicalize = |mut molecule: Molecule, marks: &[SourceStereoBondMark]| {
        let report = canonicalize_molecule_for_publication(&mut molecule, marks)
            .expect("paired directional source marks should canonicalize");
        assert_eq!(report.created_stereo_elements.len(), 1);
        molecule
            .stereo_element(report.created_stereo_elements[0])
            .expect("created double-bond element")
            .clone()
    };

    let (molecule, left, right, left_reference, _, right_reference, _) = source_graph();
    let expected = canonicalize(
        molecule,
        &[
            SourceStereoBondMark {
                bond: left_reference,
                from: left,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
            SourceStereoBondMark {
                bond: right_reference,
                from: right,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
        ],
    );

    let (molecule, left, right, _, left_alternative, right_reference, _) = source_graph();
    let alternate_left = canonicalize(
        molecule,
        &[
            SourceStereoBondMark {
                bond: left_alternative,
                from: left,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
            SourceStereoBondMark {
                bond: right_reference,
                from: right,
                kind: SourceStereoBondMarkKind::DirectionalDown,
            },
        ],
    );

    let (molecule, left, right, _, left_alternative, _, right_alternative) = source_graph();
    let both_alternatives = canonicalize(
        molecule,
        &[
            SourceStereoBondMark {
                bond: left_alternative,
                from: left,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
            SourceStereoBondMark {
                bond: right_alternative,
                from: right,
                kind: SourceStereoBondMarkKind::DirectionalUp,
            },
        ],
    );

    assert_eq!(alternate_left, expected);
    assert_eq!(both_alternatives, expected);
}

#[test]
fn smiles_ring_direction_preserves_the_textual_origin_endpoint() {
    let marked_when_opened =
        read_smiles(r"F/C=C/1CCCCC1").expect("opening ring direction should interpret");
    let marked_when_closed =
        read_smiles(r"F/C=C1CCCCC\1").expect("closing ring direction should interpret");
    let opened_stereo = marked_when_opened
        .graph()
        .stereo_elements()
        .map(|(_, element)| element.clone())
        .collect::<Vec<_>>();
    let closed_stereo = marked_when_closed
        .graph()
        .stereo_elements()
        .map(|(_, element)| element.clone())
        .collect::<Vec<_>>();

    assert_eq!(opened_stereo.len(), 1);
    assert_eq!(closed_stereo, opened_stereo);
}

#[test]
fn interpretation_enforces_small_ring_double_bond_boundary() {
    let document = smiles_api::parse_str(r"C1/C=C\CCC1").expect("source syntax should parse");
    let error = smiles_api::interpret(&document)
        .expect_err("excluded small-ring directional marks must reject interpretation");
    assert_eq!(error.offset(), 2);
    assert!(error.message().contains("UnpairedDirectionalBondMark"));

    let (cyclooctene, report) =
        read_smiles_with_report(r"C1/C=C\CCCCC1").expect("marked cyclooctene interprets");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = cyclooctene
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    assert!(matches!(element.kind, StereoElementKind::DoubleBond(_)));
}

#[test]
fn interpretation_assembles_molfile_wedge_into_tetrahedral_element() {
    let input = "\
wedge
kekule

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
    let (molecule, report) =
        read_molfile_with_report(input).expect("wedge molfile should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
        .expect("created stereo element");
    assert!(element.is_specified());
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
            assert_eq!(
                stereo.orientation,
                Some(TetrahedralOrientation::CounterClockwise)
            );
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn canonical_tetrahedral_stereo_is_identical_across_smiles_molfile_and_manual_sources() {
    let smiles = read_smiles("F[C@](Cl)(Br)I").expect("tetrahedral SMILES should interpret");
    let expected = smiles
        .graph()
        .stereo_elements()
        .next()
        .expect("SMILES should create canonical tetrahedral stereo")
        .1
        .clone();

    for written in [
        molfile::write_v2000(&smiles).expect("canonical stereo should project to V2000"),
        molfile::write_v3000(&smiles).expect("canonical stereo should project to V3000"),
    ] {
        let interpreted = read_molfile(&written).expect("projected Molfile should interpret");
        let actual = interpreted
            .graph()
            .stereo_elements()
            .next()
            .expect("Molfile should recreate canonical tetrahedral stereo")
            .1;
        assert_eq!(actual, &expected);
    }

    let mut manual = Molecule::new();
    let fluorine = manual
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let center = manual.add_atom(carbon()).expect("atom identifier capacity");
    let chlorine = manual
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    let bromine = manual
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let iodine = manual
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    for carrier in [fluorine, chlorine, bromine, iodine] {
        manual
            .add_bond(center, carrier, BondOrder::Single)
            .expect("tetrahedral carrier bond");
    }
    let manual_id = manual
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(bromine),
                    StereoCarrier::Atom(iodine),
                ],
                orientation: Some(TetrahedralOrientation::Clockwise),
            },
        )))
        .expect("manual canonical stereo element");
    assert_eq!(manual.stereo_element(manual_id).unwrap(), &expected);
}

#[test]
fn interpretation_uses_source_declared_h_for_molfile_wedge_geometry() {
    let (molecule, report) = read_molfile_with_report(implicit_h_wedge_geometry_molblock())
        .expect("implicit-H wedge molfile should interpret");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
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
            assert_eq!(
                stereo.orientation,
                Some(TetrahedralOrientation::CounterClockwise)
            );
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn normalization_assembles_wedge_either_as_explicit_unknown() {
    let (mut mol, center, carriers, marked_bond) = tetrahedral_marked_graph();
    let source_stereo = [SourceStereoBondMark {
        bond: marked_bond,
        from: center,
        kind: SourceStereoBondMarkKind::WedgeEither,
    }];

    let report = canonicalize_molecule_for_publication(&mut mol, &source_stereo)
        .expect("wedge/either should assemble as unknown stereo");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_stereo_elements[0])
        .expect("created stereo element");
    assert!(element.is_explicitly_unknown());
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assert_eq!(stereo.center, center);
            assert_eq!(stereo.carriers[0], StereoCarrier::Atom(carriers[0]));
            assert_eq!(stereo.orientation, None);
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }
}

#[test]
fn alternate_tetrahedral_wedge_carriers_publish_identical_canonical_stereo() {
    let canonicalize = |kind, bond| {
        let (mut molecule, center, _, _) = tetrahedral_marked_graph();
        let report = canonicalize_molecule_for_publication(
            &mut molecule,
            &[SourceStereoBondMark {
                bond,
                from: center,
                kind,
            }],
        )
        .expect("tetrahedral source wedge should canonicalize");
        assert_eq!(report.created_stereo_elements.len(), 1);
        molecule
            .stereo_element(report.created_stereo_elements[0])
            .expect("created tetrahedral element")
            .clone()
    };

    let wedge_on_first = canonicalize(SourceStereoBondMarkKind::WedgeUp, BondId::new(0));
    let wedge_on_second = canonicalize(SourceStereoBondMarkKind::WedgeDown, BondId::new(1));
    assert_eq!(wedge_on_second, wedge_on_first);

    let unknown_on_first = canonicalize(SourceStereoBondMarkKind::WedgeEither, BondId::new(0));
    let unknown_on_second = canonicalize(SourceStereoBondMarkKind::WedgeEither, BondId::new(1));
    assert_eq!(unknown_on_second, unknown_on_first);
    assert!(unknown_on_first.is_explicitly_unknown());
}

#[test]
fn source_stereo_rejects_an_origin_outside_the_marked_bond() {
    let mut molecule = Molecule::new();
    let a = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let b = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let outside = molecule
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let bond = molecule
        .add_bond(a, b, BondOrder::Single)
        .expect("marked bond");
    let source_stereo = [SourceStereoBondMark {
        bond,
        from: outside,
        kind: SourceStereoBondMarkKind::WedgeUp,
    }];

    let error = canonicalize_molecule_for_publication(&mut molecule, &source_stereo)
        .expect_err("the marked origin must be an endpoint of its bond");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues == vec![SourceStereoNormalizationIssue::InvalidSourceBondMarkEndpoint {
                bond,
                from: outside,
            }]
    ));
}

#[test]
fn normalization_reports_ambiguous_tetrahedral_wedge_marks() {
    let (mut mol, center, _carriers, first_bond) = tetrahedral_marked_graph();
    let second_bond = BondId::new(1);
    let source_stereo = [
        SourceStereoBondMark {
            bond: first_bond,
            from: center,
            kind: SourceStereoBondMarkKind::WedgeUp,
        },
        SourceStereoBondMark {
            bond: second_bond,
            from: center,
            kind: SourceStereoBondMarkKind::WedgeDown,
        },
    ];

    let report = canonicalize_molecule_for_publication(&mut mol, &source_stereo)
        .expect("ambiguous wedges should warn without failing");

    assert!(report
        .warnings
        .contains(&NormalizationWarning::AmbiguousTetrahedralWedgeMarks {
            center,
            mark_count: 2,
        }));
    assert!(report.created_stereo_elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn coordinate_stereo_inference_is_read_only_and_materializes_tetrahedral_stereo() {
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
    mark_all_fresh(&mut mol);
    let before = mol.clone();

    let inferred = stereo_api::infer_coordinate_stereo(&mol)
        .expect("3D tetrahedral stereo should be inferred");
    assert_eq!(mol, before);
    assert_eq!(inferred.elements.len(), 1);
    let proposed = &inferred.elements[0];
    match &proposed.kind {
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
            assert_eq!(stereo.orientation, Some(TetrahedralOrientation::Clockwise));
        }
        other => panic!("expected tetrahedral stereo, found {other:?}"),
    }

    let report = stereo_api::materialize_coordinate_stereo(&mut mol)
        .expect("3D tetrahedral stereo should materialize");
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element, proposed);
}

#[test]
fn coordinate_stereo_inference_is_read_only_and_materializes_double_bond_stereo() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_carrier = mol
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    let right_carrier = mol
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
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
    let before = mol.clone();

    let inferred = stereo_api::infer_coordinate_stereo(&mol)
        .expect("2D double-bond stereo should be inferred");
    assert_eq!(mol, before);
    assert_eq!(inferred.elements.len(), 1);
    let proposed = &inferred.elements[0];
    match &proposed.kind {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.bond, double_bond);
            assert_eq!(stereo.left, left);
            assert_eq!(stereo.right, right);
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(left_carrier));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(right_carrier));
            assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
        }
        other => panic!("expected double-bond stereo, found {other:?}"),
    }

    let report = stereo_api::materialize_coordinate_stereo(&mut mol)
        .expect("2D double-bond stereo should materialize");
    assert_eq!(report.created_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_elements[0])
        .expect("created stereo element");
    assert_eq!(element, proposed);
}

#[test]
fn coordinate_stereo_inference_assigns_axis_only_when_requested() {
    let (mol, axis) = coordinate_axis_graph(true);
    let before = mol.clone();

    let default = stereo_api::infer_coordinate_stereo(&mol)
        .expect("default coordinate-stereo inference should succeed");
    assert!(default.elements.is_empty());
    let inferred = stereo_api::infer_coordinate_stereo_with_options(
        &mol,
        CoordinateStereoOptions { infer_axes: true },
    )
    .expect("3D axis stereo should be inferred");
    assert_eq!(mol, before);
    assert_eq!(inferred.elements.len(), 1);
    let element = &inferred.elements[0];
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
            assert_eq!(stereo.orientation, Some(AxisOrientation::Clockwise));
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn coordinate_stereo_inference_skips_axis_without_3d_handedness() {
    let (mol, _axis) = coordinate_axis_graph(false);

    let result = stereo_api::infer_coordinate_stereo_with_options(
        &mol,
        CoordinateStereoOptions { infer_axes: true },
    )
    .expect("flat coordinates should be a successful non-assignment");
    assert!(result.elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
}

#[test]
fn coordinate_stereo_does_not_duplicate_existing_represented_stereo() {
    let (mut mol, center, carriers, _) = tetrahedral_marked_graph();
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    for (atom, point) in
        std::iter::once((center, Point3::new(0.0, 0.0, 0.0))).chain(carriers.iter().copied().zip([
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, -1.0),
        ]))
    {
        conformer
            .set_position(
                atom,
                crate::units::Quantity::new(point, crate::units::ANGSTROM),
            )
            .unwrap();
    }
    mol.add_conformer(conformer).expect("valid conformer");
    mol.add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
        TetrahedralStereo {
            center,
            carriers: carriers.iter().copied().map(StereoCarrier::Atom).collect(),
            orientation: Some(TetrahedralOrientation::CounterClockwise),
        },
    )))
    .expect("represented source stereo");

    let inferred = stereo_api::infer_coordinate_stereo(&mol)
        .expect("represented stereo should make coordinate inference a no-op");
    assert!(inferred.elements.is_empty());
    let report = stereo_api::materialize_coordinate_stereo(&mut mol)
        .expect("represented stereo should make materialization a no-op");
    assert!(report.created_elements.is_empty());
    assert_eq!(mol.stereo_elements().count(), 1);
}

#[test]
fn coordinate_stereo_materialization_is_transactional_on_invalid_representation() {
    let (mut mol, center, carriers, _) = tetrahedral_marked_graph();
    mol.add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
        TetrahedralStereo {
            center,
            carriers: vec![
                StereoCarrier::Atom(carriers[0]),
                StereoCarrier::Atom(carriers[0]),
                StereoCarrier::Atom(carriers[1]),
            ],
            orientation: Some(TetrahedralOrientation::Clockwise),
        },
    )))
    .expect("reference-valid but structurally invalid stereo");
    let before = mol.clone();

    let error = stereo_api::materialize_coordinate_stereo(&mut mol)
        .expect_err("invalid represented stereo must reject materialization");

    assert!(matches!(error, CoordinateStereoError::InvalidStereo(_)));
    assert_eq!(mol, before);
}

#[test]
fn invalid_source_stereo_reports_an_issue_without_publishing_a_placeholder_element() {
    let mut marked = Molecule::new();
    let a = marked.add_atom(carbon()).expect("atom identifier capacity");
    let b = marked.add_atom(carbon()).expect("atom identifier capacity");
    let bond = marked.add_bond(a, b, BondOrder::Single).expect("bond");
    let marked_source = [SourceStereoBondMark {
        bond,
        from: a,
        kind: SourceStereoBondMarkKind::WedgeEither,
    }];

    let coordinate_result = stereo_api::infer_coordinate_stereo(&marked)
        .expect("coordinate inference is independent of detached source marks");
    assert!(coordinate_result.elements.is_empty());
    let marked_before = marked.clone();

    let marked_error = canonicalize_molecule_for_publication(&mut marked, &marked_source)
        .expect_err("unassembled tetrahedral mark should fail");
    assert!(matches!(
        marked_error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnassembledTetrahedralBondMark {
            bond,
            kind: SourceStereoBondMarkKind::WedgeEither,
        })
    ));
    assert!(marked.stereo_elements().next().is_none());
    assert_eq!(marked, marked_before);

    let mut unsupported = Molecule::new();
    let c = unsupported
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let d = unsupported
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let double_bond = unsupported.add_bond(c, d, BondOrder::Double).expect("bond");
    let unsupported_source = [SourceStereoBondMark {
        bond: double_bond,
        from: c,
        kind: SourceStereoBondMarkKind::DoubleBondEither,
    }];
    let unsupported_error =
        canonicalize_molecule_for_publication(&mut unsupported, &unsupported_source)
            .expect_err("unsupported double-bond mark should fail");
    assert!(matches!(
        unsupported_error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnsupportedSourceBondMark {
            bond: double_bond,
            kind: SourceStereoBondMarkKind::DoubleBondEither,
        })
    ));

    let mut unknown = Molecule::new();
    let left = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let left_carrier = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let right_carrier = unknown
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let unknown_bond = unknown
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    unknown
        .add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    unknown
        .add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let unknown_source = [SourceStereoBondMark {
        bond: unknown_bond,
        from: left,
        kind: SourceStereoBondMarkKind::DoubleBondEither,
    }];

    let unknown_report = canonicalize_molecule_for_publication(&mut unknown, &unknown_source)
        .expect("double-bond either should assemble as unknown stereo");
    assert_eq!(unknown_report.created_stereo_elements.len(), 1);
    let (_, element) = unknown.stereo_elements().next().expect("unknown element");
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo)
            if stereo.bond == unknown_bond && stereo.orientation.is_none()
    ));

    let mut absent = Molecule::new();
    let x = absent.add_atom(carbon()).expect("atom identifier capacity");
    let y = absent.add_atom(carbon()).expect("atom identifier capacity");
    absent.add_bond(x, y, BondOrder::Single).expect("bond");
    let absent_report = canonicalize_molecule_for_publication(&mut absent, &[])
        .expect("unmarked molecule should normalize");
    assert!(absent_report.created_stereo_elements.is_empty());
    assert!(absent.stereo_elements().next().is_none());
}

#[test]
fn failed_source_stereo_canonicalization_reports_the_unpaired_mark() {
    let mut molecule = read_smiles("F[C@](Cl)(Br)I").expect("stereo SMILES should parse");
    perceive(&mut molecule).expect("stored stereo should prepare");
    let marked_bond = molecule.graph().bond_ids().next().expect("single bond");
    let source_stereo = [SourceStereoBondMark {
        bond: marked_bond,
        from: molecule.graph().bond(marked_bond).expect("marked bond").a(),
        kind: SourceStereoBondMarkKind::DirectionalUp,
    }];
    let cip = stereo_api::assign_cip_descriptors(molecule.graph_mut())
        .expect("CIP assignment should succeed");
    assert_eq!(cip.assigned.len(), 1);
    stereo_api::validate_stereo(molecule.graph()).expect("stored-state validation should succeed");
    let error = molecule
        .canonicalize_fixture_with_source_stereo(&source_stereo)
        .expect_err("unpaired directional mark should fail canonicalization");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnpairedDirectionalBondMark {
                bond: marked_bond,
            })
    ));
}

#[test]
fn stereo_validation_accepts_structural_axis_elements() {
    let mut mol = Molecule::new();
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left_carrier = mol
        .add_atom(element_atom("I"))
        .expect("atom identifier capacity");
    let right_carrier = mol
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    let axis = mol.add_bond(left, right, BondOrder::Single).expect("axis");
    mol.add_bond(left, left_carrier, BondOrder::Single)
        .expect("left carrier");
    mol.add_bond(right, right_carrier, BondOrder::Single)
        .expect("right carrier");
    let valid_axis = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![
                StereoCarrier::Atom(left_carrier),
                StereoCarrier::Atom(right_carrier),
            ],
            orientation: Some(AxisOrientation::CounterClockwise),
        })))
        .expect("axis element");

    stereo_api::validate_stereo(&mol).expect("axis should be structurally valid");

    mol.remove_stereo_element(valid_axis)
        .expect("remove valid axis");
    let invalid_axis = mol
        .add_stereo_element(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![StereoCarrier::Atom(left_carrier)],
            orientation: Some(AxisOrientation::CounterClockwise),
        })))
        .expect("invalid axis element refs are still structurally present");

    let error = stereo_api::validate_stereo(&mol).expect_err("axis should be invalid");

    assert_eq!(
        error.issues,
        vec![StereoValidationIssue::InvalidAxisCarrierCount {
            element: invalid_axis,
            axis,
            carrier_count: 1,
        }]
    );
}

#[test]
fn interpretation_assembles_molfile_atropisomeric_axis() {
    let (molecule, report) = read_molfile_with_report(rdkit_rp6306_atrop_molblock())
        .expect("RDKit atropisomer fixture interprets");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
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
            assert_eq!(stereo.orientation, Some(AxisOrientation::Clockwise));
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn molfile_writers_project_tetrahedral_stereo_independent_of_bond_endpoint_storage() {
    let molecule = canonical_tetrahedral_molecule();
    let reversed = reverse_bond_endpoint_storage(&molecule);
    let expected = molecule
        .graph()
        .stereo_elements()
        .next()
        .expect("canonical tetrahedral element")
        .1
        .clone();

    for molecule in [&molecule, &reversed] {
        let before = molecule.clone();
        for written in [
            molfile::write_v2000(molecule).expect("V2000 projects tetrahedral stereo"),
            molfile::write_v3000(molecule).expect("V3000 projects tetrahedral stereo"),
        ] {
            let (reparsed, report) =
                read_molfile_with_report(&written).expect("projected tetrahedral stereo reparses");
            assert_eq!(report.created_stereo_elements().len(), 1);
            let actual = reparsed
                .graph()
                .stereo_element(report.created_stereo_elements()[0])
                .expect("reparsed tetrahedral element");
            assert_eq!(actual.kind, expected.kind);
        }
        assert_eq!(*molecule, before);
    }
}

#[test]
fn molfile_writers_project_canonical_axis_stereo_without_mutating_the_molecule() {
    let (molecule, report) = read_molfile_with_report(rdkit_rp6306_atrop_molblock())
        .expect("RDKit atropisomer fixture interprets");
    let expected = molecule
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
        .expect("canonical axis element")
        .clone();
    let axis = match &expected.kind {
        StereoElementKind::Axis(stereo) => stereo.axis,
        other => panic!("expected canonical axis stereo, found {other:?}"),
    };
    let reversed = reverse_bond_endpoint_storage_except(&molecule, &[axis]);

    for molecule in [&molecule, &reversed] {
        let before = molecule.clone();
        for written in [
            molfile::write_v2000(molecule).expect("V2000 projects axis stereo"),
            molfile::write_v3000(molecule).expect("V3000 projects axis stereo"),
        ] {
            let (reparsed, report) =
                read_molfile_with_report(&written).expect("projected axis stereo reparses");
            assert_eq!(report.created_stereo_elements().len(), 1);
            let actual = reparsed
                .graph()
                .stereo_element(report.created_stereo_elements()[0])
                .expect("reparsed axis element");
            assert_eq!(actual.kind, expected.kind);
        }
        assert_eq!(*molecule, before);
    }
}

#[test]
fn interpretation_prefers_exocyclic_molfile_atropisomeric_axis() {
    let (molecule, report) = read_molfile_with_report(rdkit_rp6306_atrop3_molblock())
        .expect("RDKit alternate atropisomer fixture interprets");
    assert_eq!(report.created_stereo_elements().len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements()[0])
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
            assert_eq!(stereo.orientation, Some(AxisOrientation::Clockwise));
        }
        other => panic!("expected axis stereo, found {other:?}"),
    }
}

#[test]
fn interpretation_rejects_atrop_fixture_with_an_omitted_tetrahedral_carrier() {
    let error = read_molfile(rdkit_bms986142_atrop5_molblock())
        .expect_err("source stereo without its carrier declaration must not publish");
    assert!(error.to_string().contains("UnassembledTetrahedralBondMark"));
}

#[test]
fn interpretation_rejects_one_ring_endpoint_atrop_fixtures_with_omitted_carriers() {
    for fixture in [
        rdkit_zm374979_atrop1_molblock(),
        rdkit_zm374979_atrop2_molblock(),
    ] {
        let error = read_molfile(fixture)
            .expect_err("omitted tetrahedral carrier must reject interpretation");
        assert!(error.to_string().contains("UnassembledTetrahedralBondMark"));
    }
}

#[test]
fn interpretation_assembles_ring_internal_molfile_atrop_axis() {
    for fixture in [
        rdkit_macrocycle8_ortho_wedge_molblock(),
        rdkit_macrocycle8_ortho_hash_molblock(),
    ] {
        let (molecule, report) = read_molfile_with_report(fixture)
            .expect("RDKit macrocyclic atropisomer fixture interprets");
        assert_eq!(report.created_stereo_elements().len(), 1);
        assert!(molecule.graph().stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(15))
        }));
    }
}

fn tetrahedral_marked_graph() -> (Molecule, AtomId, Vec<AtomId>, BondId) {
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let carriers = ["F", "Cl", "Br", "I"]
        .into_iter()
        .map(element_atom)
        .map(|atom| mol.add_atom(atom).expect("atom identifier capacity"))
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

fn canonical_tetrahedral_molecule() -> SmallMolecule {
    let (mut molecule, center, carriers, _) = tetrahedral_marked_graph();
    molecule
        .add_stereo_element(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers: carriers.into_iter().map(StereoCarrier::Atom).collect(),
                orientation: Some(TetrahedralOrientation::CounterClockwise),
            },
        )))
        .expect("canonical tetrahedral element");
    SmallMolecule::from_graph(molecule)
}

fn reverse_bond_endpoint_storage(molecule: &SmallMolecule) -> SmallMolecule {
    reverse_bond_endpoint_storage_except(molecule, &[])
}

fn reverse_bond_endpoint_storage_except(
    molecule: &SmallMolecule,
    excluded: &[BondId],
) -> SmallMolecule {
    let mut reversed = molecule.clone();
    for (index, bond) in reversed.graph.bonds.iter_mut().enumerate() {
        if excluded.contains(&BondId::new(index as u32)) {
            continue;
        }
        let Some(bond) = bond else {
            continue;
        };
        std::mem::swap(&mut bond.a, &mut bond.b);
    }
    reversed
}
