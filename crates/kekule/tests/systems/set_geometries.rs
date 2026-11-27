//! Geometries of a structure set: bond rows bound on entry, the editor a set
//! lends, and the set's own properties beside each geometry's tables and row.

use std::sync::Arc;

use kekule::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use kekule::math::Point3;
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{
    PropertyColumn, PropertyError, PropertyKey, PropertyValue, PropertyValueRef,
};
use kekule::structure::{
    DenseArrayError, Geometry, GeometryError, GeometryIndex, Positions, Structure, StructureSet,
};
use kekule::topology::{
    AtomSelection, InstanceAtomId, TopologyAtomIndex, TopologyBondIndex, TopologyBuilder,
};
use kekule::units::{Quantity, NANOMETER};

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);
const BOND0: TopologyBondIndex = TopologyBondIndex::new(0);

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

/// Carbon monoxide at the origin and 1 nm along x, and its carbon.
fn carbon_monoxide() -> (Structure, InstanceAtomId) {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(editor.finish().unwrap()).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
        NANOMETER,
    ))
    .unwrap();
    (
        Structure::new(topology, positions).unwrap(),
        InstanceAtomId::new(instance, carbon),
    )
}

#[test]
fn detached_geometries_get_bond_rows_only_when_bound() {
    let (structure, _) = carbon_monoxide();
    let label = key("label");
    let mut member = Geometry::new(structure.positions().clone());
    assert_eq!(member.properties().atoms().len(), 2);
    assert!(member.properties().bonds().is_empty());
    // Detached payloads cannot carry bond annotations yet.
    assert!(matches!(
        member.properties_mut().bonds_mut().set_value(
            label.clone(),
            BOND0,
            Some(PropertyValue::Int(1))
        ),
        Err(PropertyError::InvalidIndex(payload)) if payload == IndexOutOfRange::new(Countable::Row, 0, 0)
    ));
    member
        .properties_mut()
        .owner_mut()
        .insert(label.clone(), PropertyValue::Int(7))
        .unwrap();
    let mut set = StructureSet::from_geometries(structure.shared_topology(), [member]).unwrap();
    let stored = set.get(0).unwrap();
    assert_eq!(stored.properties().atoms().len(), 2);
    assert_eq!(stored.properties().bonds().len(), 1);
    assert!(stored.properties().owner().is_empty());
    assert_eq!(
        stored.row().unwrap().get(&label),
        Some(PropertyValueRef::Int(7))
    );

    let mut editor = set.get_mut(0).unwrap();
    let before = editor.properties().clone();
    assert!(editor
        .bond_properties_mut()
        .insert(label.clone(), PropertyColumn::Int(vec![Some(8); 2]))
        .is_err());
    assert_eq!(editor.properties(), &before);
    editor
        .bond_properties_mut()
        .insert(label.clone(), PropertyColumn::Int(vec![Some(8)]))
        .unwrap();
    assert_eq!(
        editor.properties().bonds().value(&label, BOND0).unwrap(),
        Some(PropertyValue::Int(8))
    );
}

#[test]
fn geometry_editor_preserves_dimensions_and_validated_values() {
    let (structure, _) = carbon_monoxide();
    let mut set = StructureSet::from_structures([structure]).unwrap();
    {
        let mut member = set.get_mut(0).unwrap();
        let properties = member.properties().clone();
        assert!(member
            .atom_properties_mut()
            .insert(key("short"), PropertyColumn::Int(vec![Some(1)]))
            .is_err());
        assert!(member
            .bond_properties_mut()
            .insert(key("long"), PropertyColumn::Int(vec![Some(1); 2]))
            .is_err());
        assert_eq!(member.properties(), &properties);
        assert_eq!(
            member.set_positions(Quantity::new(vec![Point3::origin(); 3], NANOMETER)),
            Err(GeometryError::Array(Box::new(
                DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 3))
            )))
        );
        assert_eq!(member.positions().len(), 2);
    }
    let member = set.get(0).unwrap();
    assert_eq!(member.as_structure_view().positions().len(), 2);
    assert_eq!(member.to_structure().properties().bonds().len(), 1);
}

#[test]
fn set_and_geometry_properties_are_separate() {
    let (mut structure, atom) = carbon_monoxide();
    let mut geometry = structure.geometry_mut();
    let mut properties = geometry.properties_mut();
    properties
        .atoms_mut()
        .set_value(key("member_atom"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    properties
        .bonds_mut()
        .set_value(key("member_bond"), BOND0, Some(PropertyValue::Int(2)))
        .unwrap();
    let mut set = StructureSet::from_structures([structure]).unwrap();
    set.properties_mut()
        .insert(key("collection"), PropertyValue::Bool(true))
        .unwrap();
    set.columns_mut()
        .set_value(
            key("member"),
            GeometryIndex::new(0),
            Some(PropertyValue::Int(1)),
        )
        .unwrap();
    assert!(!set.properties().is_empty());
    let member = set.get(0).unwrap();
    assert!(member.properties().owner().is_empty());
    assert_eq!(
        member.row().unwrap().get(&key("member")),
        Some(PropertyValueRef::Int(1))
    );
    assert_eq!(
        member
            .properties()
            .atoms()
            .value(&key("member_atom"), ATOM0)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        member
            .properties()
            .bonds()
            .value(&key("member_bond"), BOND0)
            .unwrap(),
        Some(PropertyValue::Int(2))
    );

    let cloned = set.clone();
    assert_eq!(cloned.properties(), set.properties());
    let selection = AtomSelection::from_atoms(set.shared_topology(), [atom]).unwrap();
    let subset = set.subset(&selection).unwrap();
    assert!(subset.properties().is_empty());
    assert_eq!(subset.columns(), set.columns());
}
