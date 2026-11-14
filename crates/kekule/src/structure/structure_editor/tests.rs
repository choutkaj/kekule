use super::*;
use crate::properties::{PropertyColumn, PropertyKey, PropertyValue, RawPropertyTable};

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn int_column(table: &RawPropertyTable, key: &PropertyKey) -> *const Option<i64> {
    match table.get(key).unwrap() {
        PropertyColumn::Int(values) => values.as_ptr(),
        _ => panic!("expected an integer column"),
    }
}

// Structure batches stage only their own column, as topology batches do.
#[test]
fn structure_batches_keep_unrelated_column_allocations() {
    let topology = crate::smiles::to_topology("CCC").unwrap();
    let atoms = topology.atom_count();
    let mut editor = Structure::new(topology, Positions::zeros(atoms))
        .unwrap()
        .into_editor();
    let atoms = editor.topology().atom_ids().collect::<Vec<_>>();
    let bonds = editor.topology().bond_ids().collect::<Vec<_>>();
    let untouched = key("untouched");
    let mut properties = editor.properties_mut();
    properties
        .atoms_mut()
        .set_value(untouched.clone(), atoms[0], Some(PropertyValue::Int(1)))
        .unwrap();
    properties
        .bonds_mut()
        .set_value(untouched.clone(), bonds[0], Some(PropertyValue::Int(2)))
        .unwrap();
    let atom_column = int_column(editor.slots.properties.atoms().raw(), &untouched);
    let bond_column = int_column(editor.slots.properties.bonds().raw(), &untouched);
    let mut properties = editor.properties_mut();
    properties
        .atoms_mut()
        .set_values(
            key("edited"),
            [
                (atoms[0], Some(PropertyValue::Int(3))),
                (atoms[1], Some(PropertyValue::Int(4))),
            ],
        )
        .unwrap();
    assert!(properties
        .atoms_mut()
        .set_values(
            key("edited"),
            [(atoms[1], Some(PropertyValue::String("bad type".into())))],
        )
        .is_err());
    properties
        .bonds_mut()
        .set_values(key("edited"), [(bonds[1], Some(PropertyValue::Int(5)))])
        .unwrap();
    properties
        .bonds_mut()
        .set_values(key("edited"), [(bonds[1], None)])
        .unwrap();
    assert_eq!(
        int_column(editor.slots.properties.atoms().raw(), &untouched),
        atom_column
    );
    assert_eq!(
        int_column(editor.slots.properties.bonds().raw(), &untouched),
        bond_column
    );
}
