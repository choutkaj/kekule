use super::*;
use crate::hydrogens::{
    AddHydrogensOptions, AddedHydrogenOrigin, HydrogenNormalizationError, RetainedHydrogen,
    RetainedHydrogenReason,
};

#[test]
fn add_hydrogens_materializes_perceived_counts_and_invalidates_perception() {
    let mut molecule = SmallMolecule::from_smiles_sanitized("C").expect("sanitized methane");
    let carbon = molecule.graph().atom_ids().next().expect("carbon");
    assert_eq!(molecule.graph().implicit_hydrogens(carbon), Ok(Some(4)));

    let report = molecule.add_hydrogens().expect("materialize hydrogens");

    assert_eq!(report.added.len(), 4);
    assert!(report
        .added
        .iter()
        .all(|entry| entry.parent == carbon && entry.origin == AddedHydrogenOrigin::Implicit));
    assert_eq!(molecule.atom_count(), 5);
    assert_eq!(molecule.bond_count(), 4);
    assert!(!molecule.graph().perception().has_valence());
    for entry in &report.added {
        assert_eq!(
            molecule
                .graph()
                .atom(entry.hydrogen)
                .expect("added hydrogen")
                .element
                .symbol(),
            "H"
        );
        assert_eq!(
            molecule
                .graph()
                .neighbors(entry.hydrogen)
                .expect("hydrogen neighbor")
                .collect::<Vec<_>>(),
            vec![carbon]
        );
    }
}

#[test]
fn add_hydrogens_is_transactional_for_missing_perception_and_resource_limits() {
    let mut unsanitized = SmallMolecule::from_smiles("C").expect("methane");
    let original = unsanitized.clone();
    assert_eq!(
        unsanitized.add_hydrogens(),
        Err(HydrogenNormalizationError::MissingValencePerception)
    );
    assert_eq!(unsanitized, original);

    let mut sanitized = SmallMolecule::from_smiles_sanitized("C").expect("methane");
    let original = sanitized.clone();
    let options = AddHydrogensOptions {
        max_added_hydrogens: 3,
        ..AddHydrogensOptions::default()
    };
    assert_eq!(
        sanitized.add_hydrogens_with_options(options),
        Err(HydrogenNormalizationError::ResourceLimit {
            requested_hydrogens: 4,
            limit: 3,
        })
    );
    assert_eq!(sanitized, original);
}

#[test]
fn explicit_only_materializes_bracket_counts_without_implicit_hydrogens() {
    let mut molecule = SmallMolecule::from_smiles_sanitized("[CH3]").expect("methyl radical");
    let carbon = molecule.graph().atom_ids().next().expect("carbon");
    let report = molecule
        .add_hydrogens_with_options(AddHydrogensOptions {
            explicit_only: true,
            ..AddHydrogensOptions::default()
        })
        .expect("materialize explicit count");

    assert_eq!(report.added.len(), 3);
    assert!(report
        .added
        .iter()
        .all(|entry| entry.origin == AddedHydrogenOrigin::ExplicitCount));
    assert_eq!(
        molecule
            .graph()
            .atom(carbon)
            .expect("carbon")
            .explicit_hydrogens,
        0
    );
}

#[test]
fn add_and_remove_hydrogens_round_trip_methane_semantics() {
    let mut molecule = SmallMolecule::from_smiles_sanitized("C").expect("methane");
    let carbon = molecule.graph().atom_ids().next().expect("carbon");
    let added = molecule.add_hydrogens().expect("add hydrogens");
    molecule.sanitize().expect("resanitize explicit methane");

    let removed = molecule.remove_hydrogens().expect("remove hydrogens");

    assert_eq!(removed.removed.len(), 4);
    assert!(removed.retained.is_empty());
    assert_eq!(molecule.atom_count(), 1);
    assert_eq!(molecule.bond_count(), 0);
    assert_eq!(removed.adjustments.len(), 1);
    assert_eq!(removed.adjustments[0].parent, carbon);
    assert_eq!(removed.adjustments[0].explicit_hydrogens, 0);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 4);
    assert!(!molecule.graph().perception().has_valence());
    assert!(added
        .added
        .iter()
        .all(|entry| molecule.graph().atom(entry.hydrogen).is_err()));
    molecule.sanitize().expect("resanitize collapsed methane");
    assert_eq!(molecule.to_canonical_smiles().expect("canonical"), "C");
}

#[test]
fn remove_hydrogens_preserves_aromatic_bracket_hydrogen_counts() {
    let mut molecule = SmallMolecule::from_smiles_sanitized("c1cc[nH]c1").expect("pyrrole");
    let nitrogen = molecule
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "N").then_some(id))
        .expect("nitrogen");
    let added = molecule
        .add_hydrogens_with_options(AddHydrogensOptions {
            explicit_only: true,
            ..AddHydrogensOptions::default()
        })
        .expect("materialize bracket hydrogen");
    assert_eq!(added.added.len(), 1);
    assert_eq!(added.added[0].parent, nitrogen);
    molecule.sanitize().expect("resanitize explicit pyrrole");

    let removed = molecule.remove_hydrogens().expect("collapse hydrogen");

    assert_eq!(removed.removed.len(), 1);
    assert_eq!(removed.adjustments[0].parent, nitrogen);
    assert_eq!(removed.adjustments[0].explicit_hydrogens, 1);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 0);
    assert_eq!(
        molecule
            .graph()
            .atom(nitrogen)
            .expect("nitrogen")
            .explicit_hydrogens,
        1
    );
}

#[test]
fn hydrogen_materialization_and_collapse_preserve_tetrahedral_stereo_carriers() {
    let mut molecule =
        SmallMolecule::from_smiles_sanitized("F[C@H](Cl)Br").expect("chiral molecule");
    let (element_id, before) = molecule
        .graph()
        .stereo_elements()
        .next()
        .map(|(id, element)| (id, element.clone()))
        .expect("tetrahedral stereo");
    let center = match &before.kind {
        StereoElementKind::Tetrahedral(stereo) => stereo.center,
        _ => panic!("expected tetrahedral stereo"),
    };

    let added = molecule.add_hydrogens().expect("materialize hydrogen");
    let hydrogen = added
        .added
        .iter()
        .find(|entry| entry.parent == center)
        .expect("center hydrogen")
        .hydrogen;
    match &molecule
        .graph()
        .stereo_element(element_id)
        .expect("stereo after addition")
        .kind
    {
        StereoElementKind::Tetrahedral(stereo) => {
            assert!(stereo.carriers.contains(&StereoCarrier::Atom(hydrogen)));
        }
        _ => panic!("expected tetrahedral stereo"),
    }
    molecule.sanitize().expect("resanitize explicit hydrogen");

    let removed = molecule.remove_hydrogens().expect("collapse hydrogen");
    assert_eq!(removed.adjustments[0].explicit_hydrogens, 1);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 0);
    match &molecule
        .graph()
        .stereo_element(element_id)
        .expect("stereo after removal")
        .kind
    {
        StereoElementKind::Tetrahedral(stereo) => {
            assert!(stereo.carriers.contains(&StereoCarrier::ImplicitHydrogen));
        }
        _ => panic!("expected tetrahedral stereo"),
    }
}

#[test]
fn added_hydrogens_have_explicitly_missing_conformer_positions() {
    let mut molecule = SmallMolecule::from_smiles_sanitized("C").expect("methane");
    let carbon = molecule.graph().atom_ids().next().expect("carbon");
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            carbon,
            crate::units::Quantity::new(Point3::new(1.0, 2.0, 3.0), crate::units::ANGSTROM),
        )
        .unwrap();
    let conformer_id = molecule
        .graph_mut()
        .add_conformer(conformer)
        .expect("conformer");

    let report = molecule.add_hydrogens().expect("materialize hydrogens");
    let conformer = molecule
        .graph()
        .conformer(conformer_id)
        .expect("conformer remains");
    assert_eq!(
        conformer.position(carbon),
        Some(crate::units::Quantity::new(
            Point3::new(1.0, 2.0, 3.0),
            crate::units::ANGSTROM,
        ))
    );
    assert!(report
        .added
        .iter()
        .all(|entry| conformer.position(entry.hydrogen).is_none()));
}

#[test]
fn remove_hydrogens_reports_lossy_hydrogens_as_retained() {
    let mut graph = Molecule::new();
    let first_carbon = graph.add_atom(carbon());

    let mut isotope = element_atom("H");
    isotope.isotope = Some(2);
    let isotope = graph.add_atom(isotope);
    graph
        .add_bond(first_carbon, isotope, BondOrder::Single)
        .expect("isotope bond");

    let second_carbon = graph.add_atom(carbon());
    let mut mapped = element_atom("H");
    mapped.atom_map = Some(7);
    let mapped = graph.add_atom(mapped);
    graph
        .add_bond(second_carbon, mapped, BondOrder::Single)
        .expect("mapped bond");

    let third_carbon = graph.add_atom(carbon());
    let property_hydrogen = graph.add_atom(element_atom("H"));
    graph
        .atom_mut(property_hydrogen)
        .expect("property hydrogen")
        .props
        .insert("source".into(), PropValue::String("kept".into()));
    graph
        .add_bond(third_carbon, property_hydrogen, BondOrder::Single)
        .expect("property bond");

    let _ = valence_api::perceive_valence(&mut graph, ValenceModel::RdkitLike);
    let mut molecule = SmallMolecule::from_graph(graph);
    let report = molecule.remove_hydrogens().expect("conservative removal");

    assert!(report.removed.is_empty());
    assert_eq!(
        report
            .retained
            .iter()
            .map(|entry| (entry.hydrogen, entry.reason))
            .collect::<Vec<_>>(),
        vec![
            (isotope, RetainedHydrogenReason::Isotopic),
            (mapped, RetainedHydrogenReason::Mapped),
            (property_hydrogen, RetainedHydrogenReason::AtomProperties),
        ]
    );
}

#[test]
fn remove_hydrogens_is_transactional_when_encoded_count_overflows() {
    let mut graph = Molecule::new();
    let mut parent = carbon();
    parent.explicit_hydrogens = u8::MAX;
    parent.no_implicit_hydrogens = true;
    let parent = graph.add_atom(parent);
    let hydrogen = graph.add_atom(element_atom("H"));
    graph
        .add_bond(parent, hydrogen, BondOrder::Single)
        .expect("hydrogen bond");
    let _ = valence_api::perceive_valence(&mut graph, ValenceModel::RdkitLike);
    let mut molecule = SmallMolecule::from_graph(graph);
    let original = molecule.clone();

    assert_eq!(
        molecule.remove_hydrogens(),
        Err(HydrogenNormalizationError::HydrogenCountOverflow {
            atom: parent,
            count: 256,
        })
    );
    assert_eq!(molecule, original);
}

#[test]
fn remove_hydrogens_preserves_double_bond_stereo_carriers() {
    let mut graph = Molecule::new();
    let left = graph.add_atom(carbon());
    let right = graph.add_atom(carbon());
    let double_bond = graph
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let hydrogen = graph.add_atom(element_atom("H"));
    graph
        .add_bond(left, hydrogen, BondOrder::Single)
        .expect("hydrogen bond");
    let fluorine = graph.add_atom(element_atom("F"));
    graph
        .add_bond(left, fluorine, BondOrder::Single)
        .expect("fluorine bond");
    let chlorine = graph.add_atom(element_atom("Cl"));
    graph
        .add_bond(right, chlorine, BondOrder::Single)
        .expect("chlorine bond");
    let bromine = graph.add_atom(element_atom("Br"));
    graph
        .add_bond(right, bromine, BondOrder::Single)
        .expect("bromine bond");
    let _ = valence_api::perceive_valence(&mut graph, ValenceModel::RdkitLike);
    let stereo = graph
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond: double_bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(hydrogen),
                right_carrier: StereoCarrier::Atom(chlorine),
                orientation: DoubleBondOrientation::Opposite,
            }),
            StereoSource::User,
        ))
        .expect("double-bond stereo");
    let mut molecule = SmallMolecule::from_graph(graph);

    let report = molecule.remove_hydrogens().expect("collapse hydrogen");

    assert_eq!(report.removed[0].hydrogen, hydrogen);
    assert_eq!(report.adjustments[0].explicit_hydrogens, 1);
    assert_eq!(report.adjustments[0].implicit_hydrogens, 0);
    match &molecule
        .graph()
        .stereo_element(stereo)
        .expect("stereo survives")
        .kind
    {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.left_carrier, StereoCarrier::ImplicitHydrogen);
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(chlorine));
            assert_eq!(stereo.orientation, DoubleBondOrientation::Opposite);
        }
        _ => panic!("expected double-bond stereo"),
    }
}

#[test]
fn remove_hydrogens_retains_source_marked_bonds() {
    let mut graph = Molecule::new();
    let parent = graph.add_atom(carbon());
    let hydrogen = graph.add_atom(element_atom("H"));
    let bond = graph
        .add_bond(parent, hydrogen, BondOrder::Single)
        .expect("hydrogen bond");
    graph
        .set_stereo_bond_mark(StereoBondMark {
            bond,
            kind: StereoBondMarkKind::WedgeUp,
            source: StereoSource::MolfileV2000,
        })
        .expect("source mark");
    let _ = valence_api::perceive_valence(&mut graph, ValenceModel::RdkitLike);
    let mut molecule = SmallMolecule::from_graph(graph);

    let report = molecule.remove_hydrogens().expect("conservative removal");

    assert!(report.removed.is_empty());
    assert_eq!(
        report.retained,
        vec![RetainedHydrogen {
            hydrogen,
            reason: RetainedHydrogenReason::StereoBondMark,
        }]
    );
    assert!(molecule.graph().atom(hydrogen).is_ok());
    assert!(molecule.graph().stereo_bond_mark(bond).is_some());
}
