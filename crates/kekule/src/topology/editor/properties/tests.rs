use super::*;
use crate::properties::{PropertyColumn, RawPropertyTable};

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn int_column(table: &RawPropertyTable, key: &PropertyKey) -> *const Option<i64> {
    match table.get(key).unwrap() {
        PropertyColumn::Int(values) => values.as_ptr(),
        _ => panic!("expected an integer column"),
    }
}

// A batch stages and commits only its own column; unrelated columns keep
// their allocations through successful and failed batches.
#[test]
fn property_batches_keep_unrelated_column_allocations() {
    let mut editor = TopologyEditor::from_topology(crate::smiles::to_topology("CCC").unwrap());
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let bonds = editor.bond_ids().collect::<Vec<_>>();
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
    let atom_column = int_column(editor.properties.atoms().raw(), &untouched);
    let bond_column = int_column(editor.properties.bonds().raw(), &untouched);
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
        int_column(editor.properties.atoms().raw(), &untouched),
        atom_column
    );
    assert_eq!(
        int_column(editor.properties.bonds().raw(), &untouched),
        bond_column
    );
}

// The shared table guard keeps handle addressing: live-handle column
// order, the handle's own error, and no-op writes that keep the source.
#[test]
fn draft_table_guards_address_live_handles_and_keep_unchanged_sources() {
    let tag = key("tag");
    let mut seeded = TopologyEditor::from_topology(crate::smiles::to_topology("CCO").unwrap());
    let atoms = seeded.atom_ids().collect::<Vec<_>>();
    seeded
        .properties_mut()
        .atoms_mut()
        .set_values(
            tag.clone(),
            [
                (atoms[0], Some(PropertyValue::Int(1))),
                (atoms[2], Some(PropertyValue::Int(3))),
            ],
        )
        .unwrap();
    let source = Arc::new(seeded.finish().unwrap());

    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let column = editor.properties().atoms().get(&tag).unwrap();
    let mut properties = editor.properties_mut();
    let mut table = properties.atoms_mut();
    assert_eq!(
        table.insert(tag.clone(), column.clone()).unwrap(),
        Some(column)
    );
    table.set_values(tag.clone(), [(atoms[1], None)]).unwrap();
    table.clear_value(tag.clone(), atoms[1]).unwrap();
    assert_eq!(table.remove(&key("absent")), None);
    properties.chains_mut().clear();
    assert!(editor.finish().unwrap().shares_layout(&source));

    let mut editor = TopologyEditor::from_topology(source);
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    editor.delete_atom(atoms[0]).unwrap();
    assert_eq!(
        editor.properties_mut().atoms_mut().set_value(
            tag.clone(),
            atoms[0],
            Some(PropertyValue::Int(0))
        ),
        Err(TopologyEditError::InvalidAtom(atoms[0]))
    );
    assert_eq!(
        editor.properties().atoms().value(&tag, atoms[0]),
        Err(TopologyEditError::InvalidAtom(atoms[0]))
    );
    assert_eq!(
        editor.properties_mut().atoms_mut().get(&tag),
        Some(PropertyColumn::Int(vec![None, Some(3)]))
    );
}
