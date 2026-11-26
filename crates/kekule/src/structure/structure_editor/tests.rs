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

fn kras() -> Structure {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/mmcif/9UOH.cif");
    crate::mmcif::parse_str(&std::fs::read_to_string(path).unwrap(), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_structure()
}

/// Each site's metadata, residue labels, and chain identifiers, owned.
fn site_snapshot(
    structure: &Structure,
) -> Vec<(
    crate::topology::AtomSiteMetadata,
    String,
    String,
    Option<String>,
)> {
    structure
        .topology()
        .atom_sites()
        .map(|site| {
            let (residue, chain) = (site.residue(), site.residue().chain());
            let labels = format!(
                "{:?}",
                (
                    residue.name(),
                    residue.label_comp_id(),
                    residue.author_comp_id(),
                    residue.label_seq_id(),
                    residue.author_seq_id(),
                    residue.insertion_code(),
                )
            );
            (
                site.metadata().clone(),
                labels,
                chain.label_id().to_owned(),
                chain.author_id().map(str::to_owned),
            )
        })
        .collect()
}

// A moved copy whose chains take new identifiers keeps its source's labels
// and metadata, except that its sites take the new chain identifiers where
// they named their chain's. The source keeps its own.
#[test]
fn moved_copies_rename_the_chain_ids_of_their_sites() {
    let source = kras();
    let before = site_snapshot(&source);
    let mut editor = StructureEditor::new();
    let renamed: &crate::topology::ChainIds<'_> = &|label, author| {
        (
            format!("{label}x"),
            author.map(|author| format!("{author}x")),
        )
    };
    editor
        .append_moved_copy(
            &source,
            crate::math::Matrix3::identity(),
            crate::math::Vector3::new(0.0, 0.0, 0.0),
            renamed,
        )
        .unwrap();
    let copy = editor.finish().unwrap();
    assert_eq!(site_snapshot(&source), before);
    let after = site_snapshot(&copy);
    assert_eq!(after.len(), before.len());
    let mut changed = 0;
    for ((metadata, labels, label, author), (copied, copied_labels, copied_label, copied_author)) in
        before.iter().zip(&after)
    {
        let mut expected = metadata.clone();
        if expected.label_asym_id() == Some(label) {
            expected = expected.with_label_asym_id(Some(format!("{label}x")));
        }
        if expected.auth_asym_id().is_some() && expected.auth_asym_id() == author.as_deref() {
            expected =
                expected.with_auth_asym_id(author.as_ref().map(|author| format!("{author}x")));
        }
        changed += usize::from(expected != *metadata);
        assert_eq!(copied, &expected);
        assert_eq!(copied_labels, labels);
        assert_eq!(copied_label, &format!("{label}x"));
        assert_eq!(
            copied_author,
            &author.as_ref().map(|author| format!("{author}x"))
        );
    }
    assert!(changed > 0);
}

// A structure cut from a much larger one keeps keys into the larger one's
// label and metadata tables. Appending it copies the rows its residues and
// sites name, so the draft keeps none of the others; appending a whole
// structure shares its rows.
#[test]
fn appending_a_cut_structure_keeps_only_its_own_rows() {
    let whole = kras();
    let residue = whole.topology().residues().next().unwrap();
    let mut editor = whole.to_editor();
    let others = whole
        .topology()
        .atoms()
        .filter(|atom| atom.residue().is_none_or(|kept| kept.id() != residue.id()))
        .map(|atom| editor.topology().atom_handle(atom.id()).unwrap())
        .collect::<Vec<_>>();
    editor.delete_atoms(others).unwrap();
    let cut = editor.finish().unwrap();
    let (sites, residues) = (
        cut.topology().atom_sites().len(),
        cut.topology().residues().len(),
    );
    let rows = |structure: &Structure| {
        let hierarchy = structure.topology().hierarchy();
        (
            hierarchy.label_rows().len(),
            hierarchy.metadata_rows().len(),
        )
    };
    assert_eq!(residues, 1);
    assert!(
        rows(&cut).1 > 10 * sites,
        "{:?} rows for {sites} sites",
        rows(&cut)
    );

    let mut editor = StructureEditor::new();
    editor.append_structure(&cut).unwrap();
    let appended = editor.finish().unwrap();
    assert_eq!(rows(&appended), (residues, sites));
    assert_eq!(site_snapshot(&appended), site_snapshot(&cut));

    let mut editor = StructureEditor::new();
    editor.append_structure(&whole).unwrap();
    let appended = editor.finish().unwrap();
    assert_eq!(rows(&appended), rows(&whole));
    assert_eq!(site_snapshot(&appended), site_snapshot(&whole));
}
