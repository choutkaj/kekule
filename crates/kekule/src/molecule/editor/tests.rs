use super::publication::{validate_graph, GraphValidationError};
use super::*;
use crate::stereo::StereoViolation;
use crate::stereo::{
    DoubleBondOrientation, DoubleBondStereo, StereoCarrier, StereoElement, StereoElementKind,
};

fn carbon() -> Atom {
    Atom::new(Element::from_atomic_number(6).expect("carbon"))
}

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).expect("fixture element"))
}

#[test]
fn graph_validation_checks_property_domains() {
    for bonds in [false, true] {
        let mut editor = MoleculeEditor::new();
        editor.add_atom(carbon()).unwrap();
        if bonds {
            editor.properties.bonds_mut().resize_missing(1);
        } else {
            editor.properties.atoms_mut().resize_missing(0);
        }
        let expected = if bonds {
            GraphValidationError::BondPropertySlotCount {
                expected: 0,
                actual: 1,
            }
        } else {
            GraphValidationError::AtomPropertySlotCount {
                expected: 1,
                actual: 0,
            }
        };
        assert_eq!(validate_graph(&editor), Err(expected));
    }
}

#[test]
fn stereo_publication_issues_have_diagnostic_display_messages() {
    let error =
        MoleculePublicationError::InvalidStereo(StereoViolation::DuplicateStereoGroupMember {
            group: StereoGroupId::new(2),
            element: StereoElementId::new(7),
        });
    assert_eq!(
        error.to_string(),
        "Stereo group 2 lists stereo element 7 twice."
    );
}

#[test]
fn builder_rejects_disconnected_final_graph() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    builder.add_atom(carbon()).expect("first atom");
    builder.add_atom(carbon()).expect("second atom");
    assert_eq!(
        builder.finish(),
        Err(MoleculePublicationError::Disconnected { components: 2 })
    );
}

#[test]
fn builder_accepts_connected_graph() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let a = builder.add_atom(carbon()).expect("first atom");
    let b = builder.add_atom(carbon()).expect("second atom");
    builder
        .add_bond(a, b, BondOrder::Single)
        .expect("connecting bond");
    assert_eq!(
        builder
            .finish()
            .expect("connected molecule")
            .connected_components()
            .len(),
        1
    );
}

#[test]
fn publication_rejects_a_stereo_focus_changed_by_internal_graph_code() {
    let molecule = crate::smiles::to_molecules("F/C=C/F")
        .unwrap()
        .pop()
        .unwrap();
    let (element, stereo) = molecule.stereo_elements().next().unwrap();
    let StereoElementKind::DoubleBond(stereo) = &stereo.kind else {
        panic!("expected double-bond stereo");
    };
    let bond = stereo.bond;
    let mut editor = molecule.into_editor();
    // Internal algorithms can touch graph storage directly; publication
    // must still reject an assertion invalidated outside the public setter.
    editor.graph.bond_mut(bond).unwrap().order = BondOrder::Single;
    let error = editor.try_finish().unwrap_err();
    assert_eq!(
        error.error(),
        &MoleculePublicationError::InvalidStereo(StereoViolation::InvalidDoubleBondOrder {
            element,
            bond,
            order: BondOrder::Single,
        })
    );
    let mut repaired = error.into_draft();
    repaired.remove_stereo_element(element).unwrap();
    assert_eq!(repaired.finish().unwrap().stereo_elements().count(), 0);
}

#[test]
fn publication_rejects_nonadjacent_stereo_carriers_from_internal_graph_code() {
    let molecule = crate::smiles::to_molecules("F[C@]1(Cl)CCNC1")
        .unwrap()
        .pop()
        .unwrap();
    let (element, original) = molecule.stereo_elements().next().unwrap();
    let mut invalid = original.clone();
    let StereoElementKind::Tetrahedral(stereo) = &mut invalid.kind else {
        panic!("expected tetrahedral stereo");
    };
    let center = stereo.center;
    stereo.carriers[0] = StereoCarrier::Atom(center);
    let mut editor = molecule.into_editor();
    *editor.graph.stereo_element_mut(element).unwrap() = invalid;
    let failure = editor.try_finish().unwrap_err();
    assert_eq!(
        failure.error(),
        &MoleculePublicationError::InvalidStereo(StereoViolation::TetrahedralCarrierNotAdjacent {
            element,
            center,
            carrier: StereoCarrier::Atom(center),
        })
    );
    let mut repaired = failure.into_draft();
    repaired.remove_stereo_element(element).unwrap();
    repaired.finish().unwrap();
}

#[test]
fn canonical_bond_order_rewriting_prunes_obsolete_stereo_focus() {
    let mut editor = MoleculeEditor::new();
    let chlorine = editor.add_atom(atom("Cl")).unwrap();
    let oxo = editor.add_atom(atom("O")).unwrap();
    let hydroxyl = editor.add_atom(atom("O")).unwrap();
    let bond = editor.add_bond(chlorine, oxo, BondOrder::Double).unwrap();
    editor
        .add_bond(chlorine, hydroxyl, BondOrder::Single)
        .unwrap();
    editor
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond,
                left: chlorine,
                right: oxo,
                left_carrier: StereoCarrier::Atom(hydroxyl),
                right_carrier: StereoCarrier::ImplicitHydrogen,
                orientation: None,
            },
        )))
        .unwrap();
    assert_eq!(editor.normalize_oxohalogens(), [chlorine]);
    let molecule = editor.finish().unwrap();
    assert_eq!(molecule.bond(bond).unwrap().order, BondOrder::Single);
    assert_eq!(molecule.stereo_elements().count(), 0);
}

#[test]
fn explicit_normalization_publishes_canonical_hypervalent_representation() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
    let oxo = builder.add_atom(atom("O")).expect("oxo oxygen");
    let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
    builder
        .add_bond(chlorine, oxo, BondOrder::Double)
        .expect("source-convention oxo bond");
    builder
        .add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("hydroxyl bond");

    let before = builder.clone();
    assert_eq!(
        builder.validate(),
        Err(MoleculePublicationError::UnnormalizedOxohalogen { atom: chlorine })
    );
    assert_eq!(builder, before, "validation never rewrites chemistry");
    let error = builder.try_finish().unwrap_err();
    assert_eq!(
        error.error(),
        &MoleculePublicationError::UnnormalizedOxohalogen { atom: chlorine }
    );
    let mut builder = error.into_draft();
    assert_eq!(builder, before);
    assert_eq!(builder.normalize_oxohalogens(), [chlorine]);
    assert!(
        builder.normalize_oxohalogens().is_empty(),
        "normalization is idempotent"
    );
    let molecule = builder.finish().expect("canonical publication");

    assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
    assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
    let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap();
    assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
    assert_eq!(molecule.perception(), &Perception::default());
}

#[test]
fn failed_editor_commit_leaves_original_unchanged() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let a = builder.add_atom(carbon()).expect("first atom");
    let b = builder.add_atom(carbon()).expect("second atom");
    let bond = builder
        .add_bond(a, b, BondOrder::Single)
        .expect("connecting bond");
    let molecule = builder.finish().expect("connected molecule");
    let before = molecule.clone();

    let mut editor = molecule.to_editor();
    editor.delete_bond(bond).expect("delete in working copy");
    assert_eq!(
        editor.finish(),
        Err(MoleculePublicationError::Disconnected { components: 2 })
    );
    assert_eq!(molecule, before);
}

#[test]
fn editor_can_pass_through_temporarily_disconnected_state() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let a = builder.add_atom(carbon()).expect("first atom");
    let b = builder.add_atom(carbon()).expect("second atom");
    let c = builder.add_atom(carbon()).expect("third atom");
    let ab = builder
        .add_bond(a, b, BondOrder::Single)
        .expect("first bond");
    builder
        .add_bond(b, c, BondOrder::Single)
        .expect("second bond");
    let molecule = builder.finish().expect("connected molecule");

    let mut editor = molecule.to_editor();
    editor.delete_bond(ab).expect("temporary disconnect");
    editor
        .add_bond(a, c, BondOrder::Single)
        .expect("reconnect through another edge");
    let molecule = editor.finish().expect("connected final graph");
    assert_eq!(molecule.connected_components().len(), 1);
}

#[test]
fn editor_rejects_unnormalized_oxohalogen_edits_until_normalized() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let chlorine = builder.add_atom(atom("Cl")).expect("chlorine");
    let hydroxyl = builder.add_atom(atom("O")).expect("hydroxyl oxygen");
    builder
        .add_bond(chlorine, hydroxyl, BondOrder::Single)
        .expect("initial bond");
    let molecule = builder.finish().expect("initial molecule");

    let mut editor = molecule.to_editor();
    let oxo = editor.add_atom(atom("O")).expect("oxo oxygen");
    editor
        .add_bond(chlorine, oxo, BondOrder::Double)
        .expect("source-convention oxo bond");
    assert!(matches!(
        editor.validate(),
        Err(MoleculePublicationError::UnnormalizedOxohalogen { .. })
    ));
    editor.normalize_oxohalogens();
    let molecule = editor.finish().expect("canonical edit publication");

    assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
    assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
    let oxo_bond = molecule.bond_between(chlorine, oxo).unwrap();
    assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
}

#[test]
fn explicit_normalization_covers_direct_bond_order_changes() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let chlorine = builder.add_atom(atom("Cl")).unwrap();
    let anchor = builder.add_atom(atom("O")).unwrap();
    let oxo = builder.add_atom(atom("O")).unwrap();
    builder
        .add_bond(chlorine, anchor, BondOrder::Single)
        .unwrap();
    let oxo_bond = builder.add_bond(chlorine, oxo, BondOrder::Single).unwrap();
    let molecule = builder.finish().unwrap();

    let mut editor = molecule.to_editor();
    editor
        .bond_mut(oxo_bond)
        .unwrap()
        .set_order(BondOrder::Double);
    assert_eq!(editor.normalize_oxohalogens(), [chlorine]);
    let molecule = editor.finish().expect("canonical edit publication");

    assert_eq!(molecule.atom(chlorine).unwrap().formal_charge, 1);
    assert_eq!(molecule.atom(oxo).unwrap().formal_charge, -1);
    assert_eq!(molecule.bond(oxo_bond).unwrap().order, BondOrder::Single);
}

#[test]
fn editor_recanonicalizes_stereo_references_after_topology_changes() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let left = builder.add_atom(carbon()).unwrap();
    let right = builder.add_atom(carbon()).unwrap();
    let hydrogen = builder.add_atom(atom("H")).unwrap();
    let chlorine = builder.add_atom(atom("Cl")).unwrap();
    let bromine = builder.add_atom(atom("Br")).unwrap();
    let double_bond = builder.add_bond(left, right, BondOrder::Double).unwrap();
    builder.add_bond(left, hydrogen, BondOrder::Single).unwrap();
    builder
        .add_bond(right, chlorine, BondOrder::Single)
        .unwrap();
    builder.add_bond(right, bromine, BondOrder::Single).unwrap();
    let molecule = builder.finish().unwrap();
    let mut editor = molecule.to_editor();
    let element = editor
        .add_stereo_element(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: double_bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(hydrogen),
                right_carrier: StereoCarrier::Atom(chlorine),
                orientation: Some(DoubleBondOrientation::Together),
            },
        )))
        .unwrap();
    let molecule = editor.finish().unwrap();

    let mut editor = molecule.to_editor();
    let fluorine = editor.add_atom(atom("F")).unwrap();
    editor.add_bond(left, fluorine, BondOrder::Single).unwrap();
    let molecule = editor.finish().expect("canonical edit publication");

    let StereoElementKind::DoubleBond(stereo) = &molecule.stereo_element(element).unwrap().kind
    else {
        panic!("expected double-bond stereo");
    };
    assert_eq!(stereo.left_carrier, StereoCarrier::Atom(fluorine));
    assert_eq!(stereo.right_carrier, StereoCarrier::Atom(chlorine));
    assert_eq!(stereo.orientation, Some(DoubleBondOrientation::Opposite));
}

#[test]
fn editor_publication_leaves_out_of_model_oxohalogens_for_valence_validation() {
    let mut builder = crate::molecule::MoleculeEditor::new();
    let chlorine = builder.add_atom(atom("Cl")).unwrap();
    let anchor = builder.add_atom(atom("O")).unwrap();
    builder
        .add_bond(chlorine, anchor, BondOrder::Single)
        .unwrap();
    let mut oxo_bonds = Vec::new();
    for _ in 0..128 {
        let oxygen = builder.add_atom(atom("O")).unwrap();
        oxo_bonds.push(
            builder
                .add_bond(chlorine, oxygen, BondOrder::Single)
                .unwrap(),
        );
    }
    let molecule = builder.finish().unwrap();
    let before = molecule.clone();

    let mut editor = molecule.to_editor();
    for &bond in &oxo_bonds {
        editor.bond_mut(bond).unwrap().set_order(BondOrder::Double);
    }
    let mut edited = editor
        .finish()
        .expect("publication does not validate valence");
    assert_eq!(edited.atom(chlorine).unwrap().formal_charge, 0);
    assert!(oxo_bonds
        .iter()
        .all(|&bond| edited.bond(bond).unwrap().order == BondOrder::Double));
    assert!(edited.perceive().is_err(), "valence must reject this graph");
    assert_eq!(molecule, before);
}

#[cfg(target_pointer_width = "64")]
mod capacity;
