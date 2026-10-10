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
fn hydrogen_collapse_preserves_a_double_bond_with_only_hydrogen_references() {
    let mut molecule = perceived_smiles("[H]/N=N/[H]");
    let id = molecule.stereo_element_ids().next().unwrap();
    let mut editor = molecule.into_editor();
    let group = editor
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Relative,
            members: vec![id],
        })
        .unwrap();
    molecule = editor.finish().unwrap();
    perceive(&mut molecule).unwrap();
    let before = molecule.stereo_elements().next().unwrap().1.clone();
    assert_eq!(molecule.remove_hydrogens().unwrap().removed.len(), 2);
    let element = molecule
        .stereo_elements()
        .next()
        .expect("diazene stereo survives hydrogen collapse")
        .1;
    let StereoElementKind::DoubleBond(stereo) = &element.kind else {
        unreachable!()
    };
    assert_eq!(stereo.left_carrier, StereoCarrier::ImplicitHydrogen);
    assert_eq!(stereo.right_carrier, StereoCarrier::ImplicitHydrogen);
    let StereoElementKind::DoubleBond(original) = &before.kind else {
        unreachable!()
    };
    assert_eq!(stereo.orientation, original.orientation);
    assert_eq!(element.group, Some(group));
    assert_eq!(molecule.stereo_group(group).unwrap().members, vec![id]);

    perceive(&mut molecule).unwrap();
    molecule.add_hydrogens().unwrap();
    let restored = molecule.stereo_elements().next().unwrap().1;
    let StereoElementKind::DoubleBond(stereo) = &restored.kind else {
        unreachable!()
    };
    assert!(matches!(stereo.left_carrier, StereoCarrier::Atom(_)));
    assert!(matches!(stereo.right_carrier, StereoCarrier::Atom(_)));
    assert_eq!(stereo.orientation, original.orientation);
}

#[test]
fn hydrogen_collapse_keeps_all_tetrahedral_centers_complete_during_remapping() {
    let mut molecule = perceived_smiles("F[C@H](Cl)[C@H](F)Cl");
    let expected = stereo_api::assign_cip_descriptors(&mut molecule)
        .unwrap()
        .assigned;
    molecule.add_hydrogens().unwrap();
    let ids = molecule.stereo_element_ids().collect::<Vec<_>>();
    assert_eq!(ids.len(), 2);
    let mut editor = molecule.into_editor();
    let group = editor
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Relative,
            members: ids.clone(),
        })
        .unwrap();
    molecule = editor.finish().unwrap();
    perceive(&mut molecule).unwrap();
    assert_eq!(molecule.remove_hydrogens().unwrap().removed.len(), 2);
    assert_eq!(molecule.stereo_element_ids().collect::<Vec<_>>(), ids);
    assert_eq!(molecule.stereo_group(group).unwrap().members, ids);
    for (_, element) in molecule.stereo_elements() {
        let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
            unreachable!()
        };
        assert!(stereo.carriers.contains(&StereoCarrier::ImplicitHydrogen));
    }
    perceive(&mut molecule).unwrap();
    assert_eq!(
        stereo_api::assign_cip_descriptors(&mut molecule)
            .unwrap()
            .assigned,
        expected
    );
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
        .all(|entry| entry.parent == carbon && entry.origin == AddedHydrogenOrigin::Inferred));
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
fn added_hydrogens_use_ordinary_atom_defaults_and_do_not_expand_again() {
    for source in ["C", "N", "O", "[NH4+]", "c1cc[nH]c1", "[CH3]"] {
        for fixed_only in [false, true] {
            let mut molecule = perceived_smiles(source);
            let declarations = molecule
                .atoms()
                .map(|(id, atom)| (id, atom.hydrogens))
                .collect::<Vec<_>>();
            let options = AddHydrogensOptions {
                fixed_only,
                ..Default::default()
            };
            let added = molecule.add_hydrogens_with_options(options).unwrap();
            for entry in &added.added {
                let atom = molecule.atom(entry.hydrogen).unwrap();
                assert_eq!(atom.hydrogens, ImplicitHydrogens::default(), "{source}");
                assert_eq!(molecule.neighbors(entry.hydrogen).unwrap().count(), 1);
            }
            for (id, declaration) in declarations {
                let expected = match declaration {
                    ImplicitHydrogens::Fixed(_) => ImplicitHydrogens::Fixed(0),
                    ImplicitHydrogens::Inferred => ImplicitHydrogens::Inferred,
                };
                assert_eq!(molecule.atom(id).unwrap().hydrogens, expected, "{source}");
            }
            perceive(&mut molecule).unwrap();
            for entry in &added.added {
                assert_eq!(molecule.implicit_hydrogens(entry.hydrogen), Ok(Some(0)));
            }
            let atom_count = molecule.atom_count();
            let bond_count = molecule.bond_count();
            assert!(molecule
                .add_hydrogens_with_options(options)
                .unwrap()
                .added
                .is_empty());
            assert_eq!(molecule.atom_count(), atom_count);
            assert_eq!(molecule.bond_count(), bond_count);
        }
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
fn fixed_only_materializes_bracket_counts_without_inferred_hydrogens() {
    let mut molecule = perceived_smiles("[CH3]");
    let carbon = molecule.atom_ids().next().expect("carbon");
    let report = molecule
        .add_hydrogens_with_options(AddHydrogensOptions {
            fixed_only: true,
            ..AddHydrogensOptions::default()
        })
        .expect("materialize fixed count");

    assert_eq!(report.added.len(), 3);
    assert!(report
        .added
        .iter()
        .all(|entry| entry.origin == AddedHydrogenOrigin::Fixed));
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        ImplicitHydrogens::Fixed(0)
    );

    perceive(&mut molecule).expect("materialized fixed hydrogens perceive");
    molecule
        .remove_hydrogens()
        .expect("fixed graph hydrogens collapse");
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        ImplicitHydrogens::Fixed(3)
    );
}

#[test]
fn fixed_only_leaves_inferred_counts_implicit_without_perception() {
    let mut molecule = read_smiles("[CH3]C").expect("ethane parses");
    let (fixed, inferred) = (AtomId::new(0), AtomId::new(1));
    assert!(!molecule.perception().has_valence());

    let report = molecule
        .add_hydrogens_with_options(AddHydrogensOptions {
            fixed_only: true,
            ..AddHydrogensOptions::default()
        })
        .expect("fixed count materializes without perception");

    assert_eq!(
        report
            .added
            .iter()
            .map(|entry| (entry.parent, entry.origin))
            .collect::<Vec<_>>(),
        vec![(fixed, AddedHydrogenOrigin::Fixed); 3]
    );
    assert_eq!(
        molecule.atom(fixed).unwrap().hydrogens,
        ImplicitHydrogens::Fixed(0)
    );
    assert_eq!(
        molecule.atom(inferred).unwrap().hydrogens,
        ImplicitHydrogens::Inferred
    );
    perceive(&mut molecule).expect("partially materialized ethane perceives");
    assert_eq!(molecule.explicit_hydrogens(fixed), Ok(3));
    assert_eq!(molecule.implicit_hydrogens(fixed), Ok(Some(0)));
    assert_eq!(molecule.explicit_hydrogens(inferred), Ok(0));
    assert_eq!(molecule.implicit_hydrogens(inferred), Ok(Some(3)));
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
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 4);
    assert_eq!(
        removed.adjustments[0].hydrogens,
        ImplicitHydrogens::Inferred
    );
    assert!(!molecule.perception().has_valence());
    assert!(added
        .added
        .iter()
        .all(|entry| molecule.atom(entry.hydrogen).is_err()));
    perceive(&mut molecule).expect("re-perceive collapsed methane");
    assert_eq!(
        smiles_api::write(&molecule, smiles_api::SmilesWriteOptions::canonical())
            .expect("canonical"),
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
    let carbon = removed
        .correspondence
        .atom(carbon)
        .expect("carbon survives");
    assert_eq!(carbon, AtomId::new(0));
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        ImplicitHydrogens::Inferred
    );

    perceive(&mut molecule).expect("collapsed methane perceives");
    let added = molecule.add_hydrogens().expect("methane materializes");
    assert_eq!(added.added.len(), 4);
    assert_eq!(molecule.atom_count(), 5);
    assert_eq!(molecule.bond_count(), 4);
    assert_eq!(
        molecule.atom(carbon).expect("carbon").hydrogens,
        ImplicitHydrogens::Inferred
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
            fixed_only: true,
            ..AddHydrogensOptions::default()
        })
        .expect("materialize bracket hydrogen");
    assert_eq!(added.added.len(), 1);
    assert_eq!(added.added[0].parent, nitrogen);
    perceive(&mut molecule).expect("re-perceive explicit pyrrole");

    let removed = molecule.remove_hydrogens().expect("collapse hydrogen");

    assert_eq!(removed.removed.len(), 1);
    assert_eq!(removed.adjustments[0].parent, nitrogen);
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 1);
    assert_eq!(
        removed.adjustments[0].hydrogens,
        ImplicitHydrogens::Fixed(1)
    );
    assert_eq!(
        molecule.atom(nitrogen).expect("nitrogen").hydrogens,
        ImplicitHydrogens::Fixed(1)
    );
}

#[test]
fn hydrogen_collapse_fixes_counts_that_inference_cannot_reproduce() {
    // Inference assigns no hydrogens to these metals and only the lowest
    // allowed valence to S and P (SH2, PH3), so keeping inference would drop
    // hydrogens. Each expected count is the input's own graph-hydrogen count.
    for (symbol, count) in [("Ir", 1), ("Mo", 2), ("S", 4), ("P", 5)] {
        let mut editor = MoleculeEditor::new();
        let parent = editor.add_atom(element_atom(symbol)).unwrap();
        for _ in 0..count {
            let hydrogen = editor.add_atom(element_atom("H")).unwrap();
            editor
                .add_bond(parent, hydrogen, BondOrder::Single)
                .unwrap();
        }
        let mut molecule = editor.finish().unwrap();
        perceive(&mut molecule).unwrap();
        assert_eq!(molecule.implicit_hydrogens(parent), Ok(Some(0)), "{symbol}");
        let report = molecule.remove_hydrogens().unwrap();
        assert_eq!(report.removed.len(), count, "{symbol}");
        assert_eq!(
            report.adjustments[0].hydrogens,
            ImplicitHydrogens::Fixed(count as u8),
            "{symbol}"
        );
        assert_eq!(molecule.atom_count(), 1);
        assert_eq!(
            molecule.atom(parent).unwrap().hydrogens,
            ImplicitHydrogens::Fixed(count as u8)
        );
        perceive(&mut molecule).unwrap();
        assert_eq!(molecule.implicit_hydrogens(parent), Ok(Some(count)));
        assert_eq!(molecule.add_hydrogens().unwrap().added.len(), count);
        assert_eq!(molecule.atom_count(), count + 1);
        assert_eq!(molecule.bond_count(), count);
    }
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
    assert_eq!(removed.adjustments[0].implicit_hydrogens, 1);
    assert_eq!(
        removed.adjustments[0].hydrogens,
        ImplicitHydrogens::Fixed(1)
    );
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
        .properties_mut()
        .atoms_mut()
        .set_value(
            PropertyKey::new("source").unwrap(),
            property_hydrogen,
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

    let _ = valence_api::perceive_valence(graph.working_mut(), ValenceModel::RdkitLike);
    let mut molecule = graph;
    let report = molecule
        .working_mut()
        .remove_hydrogens()
        .expect("conservative removal");

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
    parent.hydrogens = ImplicitHydrogens::Fixed(u8::MAX);
    let parent = graph.add_atom(parent).expect("atom identifier capacity");
    let hydrogen = graph
        .add_atom(element_atom("H"))
        .expect("atom identifier capacity");
    graph
        .add_bond(parent, hydrogen, BondOrder::Single)
        .expect("hydrogen bond");
    valence_api::perceive_valence_with_options(
        graph.working_mut(),
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    )
    .expect("permissive valence perception");
    let mut molecule = graph;
    let original = molecule.clone();

    assert_eq!(
        molecule.working_mut().remove_hydrogens(),
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
    let _ = valence_api::perceive_valence(graph.working_mut(), ValenceModel::RdkitLike);
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

    let report = molecule
        .working_mut()
        .remove_hydrogens()
        .expect("collapse hydrogen");

    assert_eq!(report.removed[0].hydrogen, hydrogen);
    assert_eq!(report.adjustments[0].implicit_hydrogens, 1);
    assert_eq!(report.adjustments[0].hydrogens, ImplicitHydrogens::Inferred);
    // Atoms after the removed hydrogen are renumbered densely.
    let ids = &report.correspondence;
    assert_eq!(ids.atom(hydrogen), None);
    assert_eq!(ids.atom(fluorine), Some(hydrogen));
    match &molecule
        .stereo_element(ids.stereo_element(stereo).unwrap())
        .expect("stereo survives")
        .kind
    {
        StereoElementKind::DoubleBond(stereo) => {
            assert_eq!(
                stereo.left_carrier,
                StereoCarrier::Atom(ids.atom(fluorine).unwrap())
            );
            assert_eq!(
                stereo.right_carrier,
                StereoCarrier::Atom(ids.atom(chlorine).unwrap())
            );
            assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Together));
        }
        _ => panic!("expected double-bond stereo"),
    }
}

#[test]
fn hydrogen_collapse_retains_a_double_bond_reference_beside_another_hydrogen() {
    // On a terminal =CH2 the reference is one of two equivalent hydrogens.
    // An implicit reference there names neither: CIP could not rank it and
    // materialization rejected it, so the reference stays a graph atom.
    // A toy molecule belongs here, not in a fixture; the external JDQ443 and
    // Sotorasib 3D records cover this through the hydrogen round-trip
    // invariant. Reference: RDKit 2026.03.3 FindPotentialStereo and
    // rdCIPLabeler find no stereo on this C=C before or after collapse
    // (F/C([H])=C(/[H])[H], F/C([H])=C/[H] and F/C=C/[H]).
    for sibling_is_explicit in [true, false] {
        let mut editor = MoleculeEditor::new();
        let fluorine = editor.add_atom(element_atom("F")).unwrap();
        let left = editor.add_atom(carbon()).unwrap();
        let right = editor.add_atom(carbon()).unwrap();
        let left_hydrogen = editor.add_atom(element_atom("H")).unwrap();
        let reference = editor.add_atom(element_atom("H")).unwrap();
        editor.add_bond(fluorine, left, BondOrder::Single).unwrap();
        let bond = editor.add_bond(left, right, BondOrder::Double).unwrap();
        editor
            .add_bond(left, left_hydrogen, BondOrder::Single)
            .unwrap();
        editor
            .add_bond(right, reference, BondOrder::Single)
            .unwrap();
        if sibling_is_explicit {
            let sibling = editor.add_atom(element_atom("H")).unwrap();
            editor.add_bond(right, sibling, BondOrder::Single).unwrap();
        }
        editor
            .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
                DoubleBondStereo {
                    bond,
                    left,
                    right,
                    left_carrier: StereoCarrier::Atom(fluorine),
                    right_carrier: StereoCarrier::Atom(reference),
                    orientation: Some(DoubleBondOrientation::Opposite),
                },
            )))
            .unwrap();
        let mut molecule = editor.finish().unwrap();
        perceive(&mut molecule).unwrap();
        let id = molecule.stereo_element_ids().next().unwrap();
        let nonstereogenic = Ok(CipAssignmentReport {
            assigned: Vec::new(),
            skipped: vec![CipSkipped {
                element: id,
                reason: CipSkippedReason::NotStereogenic,
            }],
        });
        assert_eq!(
            stereo_api::assign_cip_descriptors(&mut molecule),
            nonstereogenic
        );

        let report = molecule.remove_hydrogens().unwrap();

        assert_eq!(
            report
                .retained
                .iter()
                .map(|entry| (entry.hydrogen, entry.reason))
                .collect::<Vec<_>>(),
            vec![(reference, RetainedHydrogenReason::UnsupportedStereoRole)],
            "sibling_is_explicit: {sibling_is_explicit}"
        );
        let ids = &report.correspondence;
        let StereoElementKind::DoubleBond(stereo) = &molecule
            .stereo_element(ids.stereo_element(id).unwrap())
            .unwrap()
            .kind
        else {
            unreachable!()
        };
        assert_eq!(
            stereo.right_carrier,
            StereoCarrier::Atom(ids.atom(reference).unwrap())
        );
        assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
        perceive(&mut molecule).unwrap();
        assert_eq!(
            molecule.implicit_hydrogens(ids.atom(right).unwrap()),
            Ok(Some(1))
        );
        assert_eq!(
            stereo_api::assign_cip_descriptors(&mut molecule),
            nonstereogenic
        );
        assert_eq!(molecule.add_hydrogens().unwrap().added.len(), 2);
    }
}

#[test]
fn hydrogen_collapse_requires_only_the_affected_parent_counts() {
    let mut fixed = crate::tests::read_smiles("[H][CH2]C").unwrap();
    assert!(!fixed.perception().has_valence());
    let report = fixed.remove_hydrogens().unwrap();
    assert_eq!(report.removed.len(), 1);
    let parent = report
        .correspondence
        .atom(report.removed[0].parent)
        .unwrap();
    assert_eq!(
        fixed.atom(parent).unwrap().hydrogens,
        ImplicitHydrogens::Fixed(3)
    );
    let mut unknown = crate::tests::read_smiles("[H]C").unwrap();
    let before = unknown.clone();
    assert_eq!(
        unknown.remove_hydrogens(),
        Err(HydrogenTransformError::MissingValencePerception)
    );
    assert_eq!(unknown, before);
    let mut unchanged = crate::tests::read_smiles("[H][H]").unwrap();
    let before = unchanged.clone();
    assert!(unchanged.remove_hydrogens().unwrap().removed.is_empty());
    assert_eq!(unchanged, before);
}
