//! Atom views of a structure: the topology's identity beside the geometry's
//! state.

use std::sync::Arc;

use kekule::math::Point3;
use kekule::molecule::{Atom, AtomId, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Positions, Structure};
use kekule::topology::{InstanceAtomId, InstanceBondId, TopologyAtomIndex, TopologyBuilder};
use kekule::units::{Quantity, NANOMETER, SQUARE_ANGSTROM};

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

/// Carbon monoxide at the origin and 1 nm along x, its carbon, and its bond.
fn carbon_monoxide() -> (Structure, InstanceAtomId, InstanceBondId) {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let bond = editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
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
        InstanceBondId::new(instance, bond),
    )
}

#[test]
fn structure_atom_views_join_topology_identity_with_geometry_state() {
    let (mut structure, carbon, _) = carbon_monoxide();
    structure
        .geometry_mut()
        .set_occupancy(ATOM0, Some(0.5))
        .unwrap();
    structure
        .geometry_mut()
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    structure
        .geometry_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(key("charge"), ATOM0, Some(PropertyValue::Int(-1)))
        .unwrap();

    let atom = structure.atom(carbon).unwrap();
    assert_eq!(atom.id(), carbon);
    assert_eq!(atom.index(), ATOM0);
    assert_eq!(atom.element().symbol(), "C");
    assert_eq!(atom.position().into_value(), Point3::origin());
    assert_eq!(atom.occupancy(), Some(0.5));
    let b_factor = atom.b_factor().unwrap().value_in(SQUARE_ANGSTROM).unwrap();
    assert!((b_factor - 20.0).abs() < 1.0e-12);
    assert_eq!(
        atom.geometry_property(&key("charge")),
        Some(PropertyValue::Int(-1))
    );
    // Static annotations stay on the topology view.
    assert_eq!(atom.property(&key("charge")), None);

    let oxygen = structure.atom_at(TopologyAtomIndex::new(1)).unwrap();
    assert_eq!(oxygen.element().symbol(), "O");
    assert_eq!(oxygen.occupancy(), None);
    assert_eq!(oxygen.b_factor(), None);
    assert_eq!(
        oxygen.neighbors().map(|atom| atom.id()).collect::<Vec<_>>(),
        [carbon]
    );
    assert_eq!(
        structure
            .atoms()
            .map(|atom| atom.index())
            .collect::<Vec<_>>(),
        [ATOM0, TopologyAtomIndex::new(1)]
    );
    assert!(structure.atom_at(TopologyAtomIndex::new(2)).is_none());
    let missing = InstanceAtomId::new(carbon.instance(), AtomId::new(9));
    assert!(structure.atom(missing).is_none());
    let view = structure.as_structure_view();
    assert_eq!(view.atom(carbon).unwrap().position(), atom.position());
}
