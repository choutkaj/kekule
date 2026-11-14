use super::*;
use crate::properties::{PropertyKey, PropertyValue};

#[test]
fn generic_properties_neither_invalidate_perception_nor_affect_equality() {
    let mut editor = MoleculeEditor::new();
    let atom = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let mut molecule = editor.finish().unwrap();
    molecule.perceive().unwrap();
    let perceived = molecule.perception().clone();
    let represented = molecule.clone();

    let mut properties = molecule.properties_mut();
    properties
        .owner_mut()
        .insert(
            PropertyKey::new("source").unwrap(),
            PropertyValue::String("generated".into()),
        )
        .unwrap();
    properties
        .atoms_mut()
        .set_value(
            PropertyKey::new("selected").unwrap(),
            atom,
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();

    assert_eq!(molecule.perception(), &perceived);
    assert_eq!(molecule, represented);
    molecule.clear_perception();
    assert_eq!(molecule, represented);
}
