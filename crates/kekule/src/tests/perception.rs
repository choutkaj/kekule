use super::*;

fn aromatic_atom(molecule: &Molecule, atom: AtomId) -> bool {
    molecule.atom_is_aromatic(atom).expect("atom exists") == Some(true)
}

fn aromatic_bond(molecule: &Molecule, bond: BondId) -> bool {
    molecule.bond_is_aromatic(bond).expect("bond exists") == Some(true)
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
    molecule
        .normalize()
        .expect("heteroaromatic stereo fixture should normalize");

    let atom_ids = molecule.graph().atom_ids().collect::<Vec<_>>();
    let annotated_atom = atom_ids[0];
    let annotated_bond = molecule.graph().bond_ids().next().expect("fixture bond");
    molecule.graph_mut().props_mut().insert(
        "perception_purity_fixture".to_owned(),
        PropValue::String("molecule property".to_owned()),
    );
    molecule
        .graph_mut()
        .atom_mut(annotated_atom)
        .expect("fixture atom")
        .props
        .insert("atom_note".to_owned(), PropValue::Bool(true));
    molecule
        .graph_mut()
        .bond_mut(annotated_bond)
        .expect("fixture bond")
        .props
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
    assert!(molecule.graph().stereo_bond_marks().next().is_none());
}

#[test]
fn molecule_perception_queries_read_the_installed_state_directly() {
    let mut molecule =
        read_smiles("F[C@](Cl)(Br)c1cc[nH]c1").expect("stereo aromatic fixture should parse");
    normalize_and_perceive(&mut molecule).expect("fixture should normalize_and_perceive");
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
fn normalize_and_perceive_succeeds_for_simple_molecule() {
    let mut molecule = read_smiles("CCO").expect("ethanol should parse");

    let report = molecule
        .normalize_and_perceive()
        .expect("ethanol should normalize and perceive");

    assert!(report.warnings.is_empty());
    assert!(report.created_stereo_elements.is_empty());
    assert!(molecule.graph().perception().has_valence());
    assert!(molecule.graph().perception().has_rings());
    assert!(molecule.graph().perception().has_aromaticity());
}

#[test]
fn normalize_and_perceive_matches_explicit_aromatic_source_stereo_workflow() {
    let mut combined =
        read_smiles("F/C=C/c1ccccc1").expect("aromatic directional SMILES should parse");
    let mut explicit = combined.clone();

    let combined_report = combined
        .normalize_and_perceive()
        .expect("combined workflow should succeed");
    let explicit_report = explicit.normalize().expect("explicit normalization");
    explicit.perceive().expect("explicit default perception");

    assert_eq!(combined_report, explicit_report);
    assert_eq!(combined, explicit);
    assert_eq!(combined_report.created_stereo_elements.len(), 1);
    assert!(combined.graph().stereo_bond_marks().next().is_none());
    assert!(combined
        .graph()
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert!(combined.graph().perception().has_valence());
    assert!(combined.graph().perception().has_rings());
    assert!(combined.graph().perception().has_aromaticity());
    assert!(!combined.graph().perception().has_cip_descriptors());
}

#[test]
fn normalize_and_perceive_does_not_infer_or_materialize_coordinate_stereo() {
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
        .normalize_and_perceive()
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
fn normalize_and_perceive_rolls_back_normalization_failure() {
    let mut graph = Molecule::new();
    let a = graph.add_atom(carbon()).expect("atom identifier capacity");
    let b = graph.add_atom(carbon()).expect("atom identifier capacity");
    let bond = graph.add_bond(a, b, BondOrder::Single).expect("bond");
    graph
        .set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("source mark");
    let mut molecule = SmallMolecule::from_graph(graph);
    let before = molecule.clone();

    let error = molecule
        .normalize_and_perceive()
        .expect_err("invalid source stereo must fail normalization");

    assert!(matches!(
        error,
        NormalizeAndPerceiveError::Normalization(NormalizationError::SourceStereo(_))
    ));
    assert_eq!(molecule, before);
}

#[test]
fn normalize_and_perceive_rolls_back_perception_after_effective_normalization() {
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
    let before = molecule.clone();
    let mut normalized = molecule.clone();
    normalized
        .normalize()
        .expect("representation normalization should succeed");
    assert_ne!(normalized, before);
    assert!(normalized.perceive().is_err());

    let error = molecule
        .normalize_and_perceive()
        .expect_err("default perception should reject pentavalent carbon");

    assert!(matches!(error, NormalizeAndPerceiveError::Perception(_)));
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
        nitrogen.explicit_hydrogens = 1;
        nitrogen.no_implicit_hydrogens = true;
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
        nitrogen.explicit_hydrogens = 1;
        nitrogen.no_implicit_hydrogens = true;
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
        cation.explicit_hydrogens = 1;
        cation.no_implicit_hydrogens = true;
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
        saturated.explicit_hydrogens = 2;
        saturated.no_implicit_hydrogens = true;
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
                orientation: TetrahedralOrientation::Clockwise,
            }),
            specifiedness: StereoSpecifiedness::Unknown,
            source: StereoSource::User,
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
    assert!(
        mol.stereo_element(element).expect("element").specifiedness == StereoSpecifiedness::Unknown
    );
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
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers,
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("tetrahedral stereo element");

    stereo_api::validate_stereo(&tetrahedral)
        .expect("implicit hydrogen availability is chemically interpretive");

    tetrahedral
        .remove_stereo_element(hydrogen_element)
        .expect("remove hydrogen-carrier element");
    atom_carriers.push(StereoCarrier::ImplicitLonePair);
    tetrahedral
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: atom_carriers,
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
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
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond,
                left,
                right,
                left_carrier: StereoCarrier::ImplicitHydrogen,
                right_carrier: StereoCarrier::ImplicitLonePair,
                orientation: DoubleBondOrientation::Together,
            }),
            StereoSource::User,
        ))
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
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Axis(AxisStereo {
                axis: axis_bond,
                carriers: vec![
                    StereoCarrier::ImplicitHydrogen,
                    StereoCarrier::ImplicitLonePair,
                ],
                orientation: AxisOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
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
    normalize_and_perceive(&mut molecule).expect("molecule should normalize_and_perceive");

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
fn normalization_assembles_paired_directional_marks_into_double_bond_element() {
    let mut molecule = read_smiles("C/C=C\\F").expect("directional smiles should parse");
    let report = molecule
        .normalize()
        .expect("directional marks should assemble");
    assert_eq!(report.created_stereo_elements.len(), 1);
    assert!(molecule.graph().stereo_elements().next().is_some());
    assert!(molecule.graph().stereo_bond_marks().next().is_none());
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements[0])
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
fn source_stereo_normalization_keeps_small_ring_double_bond_boundary() {
    let mut cyclohexene = read_smiles(r"C1/C=C\CCC1").expect("marked cyclohexene parses");
    let error = cyclohexene
        .normalize()
        .expect_err("excluded small-ring directional marks remain unpaired");
    let NormalizationError::SourceStereo(SourceStereoNormalizationError { issues }) = error else {
        panic!("expected source-stereo normalization error");
    };
    assert_eq!(issues.len(), 2);
    assert!(issues.iter().all(|issue| matches!(
        issue,
        SourceStereoNormalizationIssue::UnpairedDirectionalBondMark { .. }
    )));
    assert!(cyclohexene.graph().stereo_elements().next().is_none());

    let mut cyclooctene = read_smiles(r"C1/C=C\CCCCC1").expect("marked cyclooctene parses");
    let report = cyclooctene
        .normalize()
        .expect("cyclooctene stereo should assemble");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = cyclooctene
        .graph()
        .stereo_element(report.created_stereo_elements[0])
        .expect("created stereo element");
    assert!(matches!(element.kind, StereoElementKind::DoubleBond(_)));
}

#[test]
fn normalization_assembles_molfile_wedge_into_tetrahedral_element() {
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
    let mut molecule = read_molfile(input).expect("wedge molfile should parse");

    let report = molecule.normalize().expect("Molfile wedge should assemble");
    assert_eq!(report.created_stereo_elements.len(), 1);
    assert!(molecule.graph().stereo_bond_marks().next().is_none());
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements[0])
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
fn normalization_uses_source_declared_h_for_molfile_wedge_geometry() {
    let mut molecule = read_molfile(implicit_h_wedge_geometry_molblock())
        .expect("implicit-H wedge molfile should parse");
    let report = molecule
        .normalize()
        .expect("implicit-H wedge should assemble");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements[0])
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
fn normalization_assembles_wedge_either_as_explicit_unknown() {
    let (mut mol, center, carriers, marked_bond) = tetrahedral_marked_graph();
    mol.set_stereo_bond_mark(StereoBondMark {
        bond: marked_bond,
        kind: StereoBondMarkKind::WedgeEither,
        source: StereoSource::MolfileV2000,
    })
    .expect("wedge mark");

    let report = normalization_api::normalize(&mut mol)
        .expect("wedge/either should assemble as unknown stereo");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = mol
        .stereo_element(report.created_stereo_elements[0])
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
fn normalization_reports_ambiguous_tetrahedral_wedge_marks() {
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

    let report = normalization_api::normalize(&mut mol)
        .expect("ambiguous wedges should warn without failing");

    assert!(report
        .warnings
        .contains(&NormalizationWarning::AmbiguousTetrahedralWedgeMarks {
            center,
            mark_count: 2,
        }));
    assert!(report.created_stereo_elements.is_empty());
    assert!(mol.stereo_elements().next().is_none());
    assert!(mol.stereo_bond_marks().next().is_none());
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
    assert_eq!(proposed.source, StereoSource::Coordinates3D);
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
            assert_eq!(stereo.orientation, TetrahedralOrientation::Clockwise);
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
    assert_eq!(proposed.source, StereoSource::Coordinates2D);
    match &proposed.kind {
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
    mol.add_stereo_element(StereoElement::specified(
        StereoElementKind::Tetrahedral(TetrahedralStereo {
            center,
            carriers: carriers.iter().copied().map(StereoCarrier::Atom).collect(),
            orientation: TetrahedralOrientation::CounterClockwise,
        }),
        StereoSource::MolfileV2000,
    ))
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
    mol.add_stereo_element(StereoElement::specified(
        StereoElementKind::Tetrahedral(TetrahedralStereo {
            center,
            carriers: vec![
                StereoCarrier::Atom(carriers[0]),
                StereoCarrier::Atom(carriers[0]),
                StereoCarrier::Atom(carriers[1]),
            ],
            orientation: TetrahedralOrientation::Clockwise,
        }),
        StereoSource::User,
    ))
    .expect("reference-valid but structurally invalid stereo");
    let before = mol.clone();

    let error = stereo_api::materialize_coordinate_stereo(&mut mol)
        .expect_err("invalid represented stereo must reject materialization");

    assert!(matches!(error, CoordinateStereoError::InvalidStereo(_)));
    assert_eq!(mol, before);
}

#[test]
fn normalization_reports_unassembled_marks_and_preserves_absence() {
    let mut marked = Molecule::new();
    let a = marked.add_atom(carbon()).expect("atom identifier capacity");
    let b = marked.add_atom(carbon()).expect("atom identifier capacity");
    let bond = marked.add_bond(a, b, BondOrder::Single).expect("bond");
    marked
        .set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("mark");

    let coordinate_result = stereo_api::infer_coordinate_stereo(&marked)
        .expect("coordinate inference must ignore source marks");
    assert!(coordinate_result.elements.is_empty());
    assert!(marked.stereo_bond_mark(bond).is_some());

    let marked_error = normalization_api::normalize(&mut marked)
        .expect_err("unassembled tetrahedral mark should fail");
    assert!(matches!(
        marked_error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnassembledTetrahedralBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeEither,
        })
    ));
    assert!(marked.stereo_elements().next().is_none());
    assert!(marked.stereo_bond_mark(bond).is_some());

    let mut unsupported = Molecule::new();
    let c = unsupported
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let d = unsupported
        .add_atom(carbon())
        .expect("atom identifier capacity");
    let double_bond = unsupported.add_bond(c, d, BondOrder::Double).expect("bond");
    unsupported
        .set_stereo_bond_mark(StereoBondMark {
            bond: double_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("double bond either mark");
    let unsupported_error = normalization_api::normalize(&mut unsupported)
        .expect_err("unsupported double-bond mark should fail");
    assert!(matches!(
        unsupported_error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnsupportedSourceBondMark {
            bond: double_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
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
    unknown
        .set_stereo_bond_mark(StereoBondMark {
            bond: unknown_bond,
            kind: StereoBondMarkKind::DoubleBondEither,
            source: StereoSource::MolfileV2000,
        })
        .expect("double bond either mark");

    let unknown_report = normalization_api::normalize(&mut unknown)
        .expect("double-bond either should assemble as unknown stereo");
    assert_eq!(unknown_report.created_stereo_elements.len(), 1);
    assert!(unknown.stereo_bond_mark(unknown_bond).is_none());
    let (_, element) = unknown.stereo_elements().next().expect("unknown element");
    assert_eq!(element.specifiedness, StereoSpecifiedness::Unknown);
    assert!(matches!(
        &element.kind,
        StereoElementKind::DoubleBond(stereo) if stereo.bond == unknown_bond
    ));

    let mut absent = Molecule::new();
    let x = absent.add_atom(carbon()).expect("atom identifier capacity");
    let y = absent.add_atom(carbon()).expect("atom identifier capacity");
    absent.add_bond(x, y, BondOrder::Single).expect("bond");
    let absent_report =
        normalization_api::normalize(&mut absent).expect("unmarked molecule should normalize");
    assert!(absent_report.created_stereo_elements.is_empty());
    assert!(absent.stereo_elements().next().is_none());
    assert!(absent.stereo_bond_marks().next().is_none());
}

#[test]
fn failed_source_stereo_normalization_preserves_complete_original_state() {
    let mut molecule = read_smiles("F[C@](Cl)(Br)I").expect("stereo SMILES should parse");
    normalize_and_perceive(&mut molecule).expect("stored stereo should prepare");
    let marked_bond = molecule.graph().bond_ids().next().expect("single bond");
    molecule
        .graph_mut()
        .set_stereo_bond_mark(StereoBondMark {
            bond: marked_bond,
            kind: StereoBondMarkKind::DirectionalUp,
            source: StereoSource::Smiles,
        })
        .expect("directional source mark");
    let cip = stereo_api::assign_cip_descriptors(molecule.graph_mut())
        .expect("CIP assignment should succeed");
    assert_eq!(cip.assigned.len(), 1);
    stereo_api::validate_stereo(molecule.graph())
        .expect("stored-state validation must ignore source marks");
    let before = molecule.clone();

    let error = molecule
        .normalize()
        .expect_err("unpaired directional mark should fail normalization");

    assert!(matches!(
        error,
        NormalizationError::SourceStereo(SourceStereoNormalizationError { issues })
            if issues.contains(&SourceStereoNormalizationIssue::UnpairedDirectionalBondMark {
                bond: marked_bond,
            })
    ));
    assert_eq!(molecule, before);
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

    stereo_api::validate_stereo(&mol).expect("axis should be structurally valid");

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
fn normalization_assembles_molfile_atropisomeric_axis() {
    let mut molecule =
        read_molfile(rdkit_rp6306_atrop_molblock()).expect("RDKit atropisomer fixture parses");

    let report = molecule
        .normalize()
        .expect("Molfile atrop axis should assemble");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements[0])
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
fn normalization_prefers_exocyclic_molfile_atropisomeric_axis() {
    let mut molecule = read_molfile(rdkit_rp6306_atrop3_molblock())
        .expect("RDKit alternate atropisomer fixture parses");

    let report = molecule
        .normalize()
        .expect("preferred exocyclic Molfile atrop axis should assemble");
    assert_eq!(report.created_stereo_elements.len(), 1);
    let element = molecule
        .graph()
        .stereo_element(report.created_stereo_elements[0])
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
fn normalization_consumes_redundant_molfile_atrop_wedges_before_tetrahedral_marks() {
    let mut molecule = read_molfile(rdkit_bms986142_atrop5_molblock())
        .expect("RDKit redundant atropisomer wedge fixture parses");
    // The external fixture omits this tetrahedral carrier declaration. Keep
    // normalization model-free by supplying the represented H explicitly.
    declare_explicit_fixture_hydrogen(&mut molecule, AtomId::new(10));
    let report = molecule
        .normalize()
        .expect("redundant Molfile atrop wedges should assemble");
    assert_eq!(report.created_stereo_elements.len(), 2);
    assert!(molecule.graph().stereo_bond_marks().next().is_none());
    assert!(molecule.graph().stereo_elements().any(|(_, element)| {
        matches!(&element.kind, StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(10))
    }));
    assert!(molecule.graph().stereo_elements().any(|(_, element)| {
        matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(8))
    }));
}

#[test]
fn normalization_assembles_molfile_atrop_axis_with_one_exocyclic_sp2_endpoint() {
    for fixture in [
        rdkit_zm374979_atrop1_molblock(),
        rdkit_zm374979_atrop2_molblock(),
    ] {
        let mut molecule =
            read_molfile(fixture).expect("RDKit one-ring-endpoint atropisomer fixture parses");
        declare_explicit_fixture_hydrogen(&mut molecule, AtomId::new(3));
        let report = molecule
            .normalize()
            .expect("one-ring-endpoint Molfile atrop axis should assemble");
        assert_eq!(report.created_stereo_elements.len(), 2);
        assert!(molecule.graph().stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Tetrahedral(stereo) if stereo.center == AtomId::new(3))
        }));
        assert!(molecule.graph().stereo_elements().any(|(_, element)| {
            matches!(&element.kind, StereoElementKind::Axis(stereo) if stereo.axis == BondId::new(33))
        }));
    }
}

#[test]
fn normalization_assembles_ring_internal_molfile_atrop_axis() {
    for fixture in [
        rdkit_macrocycle8_ortho_wedge_molblock(),
        rdkit_macrocycle8_ortho_hash_molblock(),
    ] {
        let mut molecule =
            read_molfile(fixture).expect("RDKit macrocyclic atropisomer fixture parses");
        let report = molecule
            .normalize()
            .expect("ring-internal Molfile atrop axis should assemble");
        assert_eq!(report.created_stereo_elements.len(), 1);
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
