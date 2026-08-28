use super::*;
use crate::hydrogens::{
    AddHydrogensOptions, AddedHydrogenOrigin, HydrogenTransformError, RetainedHydrogenReason,
};
use crate::properties::{PropertyKey, PropertyValue};

fn perceived_smiles(input: &str) -> Molecule {
    let mut molecule = read_smiles(input).expect("SMILES should parse");
    perceive(&mut molecule).expect("molecule should perceive");
    molecule
}

#[test]
fn add_hydrogens_materializes_perceived_counts_and_invalidates_perception() {
    let mut molecule = perceived_smiles("C");
    let carbon = molecule.atom_ids().next().expect("carbon");
    assert_eq!(molecule.implicit_hydrogens(carbon), Ok(Some(4)));

    let report = molecule.add_hydrogens().expect("materialize hydrogens");

    assert_eq!(report.added.len(), 4);
    assert!(report
        .added
        .iter()
        .all(|entry| entry.parent == carbon && entry.origin == AddedHydrogenOrigin::Implicit));
    assert_eq!(molecule.atom_count(), 5);
    assert_eq!(molecule.bond_count(), 4);
    assert!(!molecule.perception().has_valence());
    for entry in &report.added {
        assert_eq!(
            molecule
                .atom(entry.hydrogen)
                .expect("added hydrogen")
                .element
                .symbol(),
            "H"
        );
        assert_eq!(
            molecule
                .neighbors(entry.hydrogen)
                .expect("hydrogen neighbor")
                .collect::<Vec<_>>(),
            vec![carbon]
        );
    }
}

#[test]
fn add_hydrogens_is_transactional_for_missing_perception_and_resource_limits() {
    let mut unperceived = read_smiles("C").expect("methane");
    let original = unperceived.clone();
    assert_eq!(
        unperceived.add_hydrogens(),
        Err(HydrogenTransformError::MissingValencePerception)
    );
    assert_eq!(unperceived, original);

    let mut perceived = perceived_smiles("C");
    let original = perceived.clone();
    let options = AddHydrogensOptions {
        max_added_hydrogens: 3,
        ..AddHydrogensOptions::default()
    };
    assert_eq!(
        perceived.add_hydrogens_with_options(options),
        Err(HydrogenTransformError::ResourceLimit {
            requested_hydrogens: 4,
            limit: 3,
        })
    );
    assert_eq!(perceived, original);
}

#[test]
fn explicit_only_materializes_bracket_counts_without_implicit_hydrogens() {
    let mut molecule = perceived_smiles("[CH3]");
    let carbon = molecule.atom_ids().next().expect("carbon");
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
        molecule.atom(carbon).expect("carbon").hydrogens,
        HydrogenDeclaration::Fixed(0)
    );

    perceive(&mut molecule).expect("materialized fixed hydrogens perceive");
    molecule
        .remove_hydrogens()
        .expect("fixed graph hydrogens collapse");
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        HydrogenDeclaration::Fixed(3)
    );
}

#[test]
fn materializing_inferred_declaration_preserves_inference_policy() {
    let mut graph = crate::core::MoleculeEditor::new();
    let mut carbon_atom = carbon();
    carbon_atom.hydrogens = HydrogenDeclaration::Infer { explicit: 1 };
    let carbon = graph.add_atom(carbon_atom).expect("carbon");
    let graph = graph.finish().expect("single atom graph");
    let mut molecule = graph;
    perceive(&mut molecule).expect("represented-plus-inferred carbon perceives");
    assert_eq!(molecule.implicit_hydrogens(carbon), Ok(Some(3)));

    let report = molecule
        .add_hydrogens_with_options(AddHydrogensOptions {
            explicit_only: true,
            ..AddHydrogensOptions::default()
        })
        .expect("represented hydrogen materializes");
    assert_eq!(report.added.len(), 1);
    assert_eq!(report.added[0].origin, AddedHydrogenOrigin::ExplicitCount);
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        HydrogenDeclaration::Infer { explicit: 0 }
    );

    perceive(&mut molecule).expect("materialized inference policy perceives");
    assert_eq!(molecule.implicit_hydrogens(carbon), Ok(Some(3)));
    molecule
        .remove_hydrogens()
        .expect("materialized hydrogen collapses");
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        HydrogenDeclaration::Infer { explicit: 0 }
    );
    perceive(&mut molecule).expect("collapsed inference policy perceives");
    assert_eq!(molecule.implicit_hydrogens(carbon), Ok(Some(4)));
}

#[test]
fn add_and_remove_hydrogens_round_trip_methane_semantics() {
    let mut molecule = perceived_smiles("C");
    let carbon = molecule.atom_ids().next().expect("carbon");
    let added = molecule.add_hydrogens().expect("add hydrogens");
    perceive(&mut molecule).expect("re-perceive explicit methane");

    let removed = molecule.remove_hydrogens().expect("remove hydrogens");

    assert_eq!(removed.removed.len(), 4);
    assert!(removed.retained.is_empty());
    assert_eq!(molecule.atom_count(), 1);
    assert_eq!(molecule.bond_count(), 0);
    assert_eq!(removed.adjustments.len(), 1);
    assert_eq!(removed.adjustments[0].parent, carbon);
    assert_eq!(removed.adjustments[0].explicit_hydrogens, 0);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 4);
    assert!(!molecule.perception().has_valence());
    assert!(added
        .added
        .iter()
        .all(|entry| molecule.atom(entry.hydrogen).is_err()));
    perceive(&mut molecule).expect("re-perceive collapsed methane");
    assert_eq!(
        smiles_api::write_canonical(&molecule).expect("canonical"),
        "C"
    );
}

#[test]
fn remove_and_add_hydrogens_round_trip_graph_methane() {
    let mut molecule = perceived_smiles("[H]C([H])([H])[H]");
    let carbon = molecule
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "C").then_some(id))
        .expect("carbon");

    let removed = molecule.remove_hydrogens().expect("collapse graph methane");
    assert_eq!(removed.removed.len(), 4);
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        HydrogenDeclaration::Infer { explicit: 0 }
    );

    perceive(&mut molecule).expect("collapsed methane perceives");
    let added = molecule.add_hydrogens().expect("methane materializes");
    assert_eq!(added.added.len(), 4);
    assert_eq!(molecule.atom_count(), 5);
    assert_eq!(molecule.bond_count(), 4);
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        HydrogenDeclaration::Infer { explicit: 0 }
    );
}

#[test]
fn remove_hydrogens_preserves_aromatic_bracket_hydrogen_counts() {
    let mut molecule = perceived_smiles("c1cc[nH]c1");
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
    perceive(&mut molecule).expect("re-perceive explicit pyrrole");

    let removed = molecule.remove_hydrogens().expect("collapse hydrogen");

    assert_eq!(removed.removed.len(), 1);
    assert_eq!(removed.adjustments[0].parent, nitrogen);
    assert_eq!(removed.adjustments[0].explicit_hydrogens, 1);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 0);
    assert_eq!(
        molecule
            .atom(nitrogen)
            .expect("nitrogen")
            .hydrogens
            .explicit_count(),
        1
    );
}

#[test]
fn hydrogen_materialization_and_collapse_preserve_tetrahedral_stereo_carriers() {
    let mut molecule = perceived_smiles("F[C@H](Cl)Br");
    let (element_id, before) = molecule
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
        .stereo_element(element_id)
        .expect("stereo after addition")
        .kind
    {
        StereoElementKind::Tetrahedral(stereo) => {
            assert!(stereo.carriers.contains(&StereoCarrier::Atom(hydrogen)));
        }
        _ => panic!("expected tetrahedral stereo"),
    }
    perceive(&mut molecule).expect("re-perceive explicit hydrogen");

    let removed = molecule.remove_hydrogens().expect("collapse hydrogen");
    assert_eq!(removed.adjustments[0].explicit_hydrogens, 1);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 0);
    match &molecule
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
fn remove_hydrogens_reports_lossy_hydrogens_as_retained() {
    let mut graph = crate::core::MoleculeEditor::new();
    let first_carbon = graph.add_atom(carbon()).expect("atom identifier capacity");

    let mut isotope = element_atom("H");
    isotope.isotope = Some(2);
    let isotope = graph.add_atom(isotope).expect("atom identifier capacity");
    graph
        .add_bond(first_carbon, isotope, BondOrder::Single)
        .expect("isotope bond");

    let second_carbon = graph.add_atom(carbon()).expect("atom identifier capacity");
    let mut mapped = element_atom("H");
    mapped.atom_map = Some(7);
    let mapped = graph.add_atom(mapped).expect("atom identifier capacity");
    graph
        .add_bond(second_carbon, mapped, BondOrder::Single)
        .expect("mapped bond");

    let third_carbon = graph.add_atom(carbon()).expect("atom identifier capacity");
    let property_hydrogen = graph
        .add_atom(element_atom("H"))
        .expect("atom identifier capacity");
    graph
        .set_atom_property(
            property_hydrogen,
            PropertyKey::new("source").unwrap(),
            Some(PropertyValue::String("kept".into())),
        )
        .unwrap();
    graph
        .add_bond(third_carbon, property_hydrogen, BondOrder::Single)
        .expect("property bond");
    graph
        .add_bond(first_carbon, second_carbon, BondOrder::Single)
        .expect("first carbon link");
    graph
        .add_bond(second_carbon, third_carbon, BondOrder::Single)
        .expect("second carbon link");

    let _ = valence_api::perceive_valence(&mut graph, ValenceModel::RdkitLike);
    let mut molecule = graph;
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
    let mut graph = crate::core::MoleculeEditor::new();
    let mut parent = carbon();
    parent.hydrogens = HydrogenDeclaration::Fixed(u8::MAX);
    let parent = graph.add_atom(parent).expect("atom identifier capacity");
    let hydrogen = graph
        .add_atom(element_atom("H"))
        .expect("atom identifier capacity");
    graph
        .add_bond(parent, hydrogen, BondOrder::Single)
        .expect("hydrogen bond");
    valence_api::perceive_valence_with_options(
        &mut graph,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect("permissive valence perception");
    let mut molecule = graph;
    let original = molecule.clone();

    assert_eq!(
        molecule.remove_hydrogens(),
        Err(HydrogenTransformError::HydrogenCountOverflow {
            atom: parent,
            count: 256,
        })
    );
    assert_eq!(molecule, original);
}

#[test]
fn remove_hydrogens_preserves_double_bond_stereo_carriers() {
    let mut graph = crate::core::MoleculeEditor::new();
    let left = graph.add_atom(carbon()).expect("atom identifier capacity");
    let right = graph.add_atom(carbon()).expect("atom identifier capacity");
    let double_bond = graph
        .add_bond(left, right, BondOrder::Double)
        .expect("double bond");
    let hydrogen = graph
        .add_atom(element_atom("H"))
        .expect("atom identifier capacity");
    graph
        .add_bond(left, hydrogen, BondOrder::Single)
        .expect("hydrogen bond");
    let fluorine = graph
        .add_atom(element_atom("F"))
        .expect("atom identifier capacity");
    graph
        .add_bond(left, fluorine, BondOrder::Single)
        .expect("fluorine bond");
    let chlorine = graph
        .add_atom(element_atom("Cl"))
        .expect("atom identifier capacity");
    graph
        .add_bond(right, chlorine, BondOrder::Single)
        .expect("chlorine bond");
    let bromine = graph
        .add_atom(element_atom("Br"))
        .expect("atom identifier capacity");
    graph
        .add_bond(right, bromine, BondOrder::Single)
        .expect("bromine bond");
    let _ = valence_api::perceive_valence(&mut graph, ValenceModel::RdkitLike);
    let stereo = graph
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(hydrogen),
                right_carrier: StereoCarrier::Atom(chlorine),
                orientation: Some(DoubleBondOrientation::Opposite),
            },
        )))
        .expect("double-bond stereo");
    let mut molecule = graph;

    let report = molecule.remove_hydrogens().expect("collapse hydrogen");

    assert_eq!(report.removed[0].hydrogen, hydrogen);
    assert_eq!(report.adjustments[0].explicit_hydrogens, 0);
    assert_eq!(report.adjustments[0].implicit_hydrogens, 1);
    match &molecule
        .stereo_element(stereo)
        .expect("stereo survives")
        .kind
    {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(stereo.left_carrier, StereoCarrier::Atom(fluorine));
            assert_eq!(stereo.right_carrier, StereoCarrier::Atom(chlorine));
            assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Together));
        }
        _ => panic!("expected double-bond stereo"),
    }
}
