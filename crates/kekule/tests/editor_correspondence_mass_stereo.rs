use kekule::{
    core::{Atom, BondOrder, Element, HydrogenDeclaration},
    geometry::Point3,
    structure::{Model, ModelEditor, Positions},
    topology::TopologyEditor,
    units::{Quantity, BOHR, DALTON},
};
use std::sync::Arc;
fn atom(z: u8) -> Atom {
    let mut a = Atom::new(Element::from_atomic_number(z).unwrap());
    a.hydrogens = HydrogenDeclaration::Fixed(0);
    a
}
#[test]
fn publication_correspondence_tracks_reordering_deletion_and_bonds() {
    let mut editor = ModelEditor::new();
    let a = editor
        .add_atom(atom(6), Quantity::new(Point3::new(1., 0., 0.), BOHR))
        .unwrap();
    let b = editor
        .add_atom(atom(8), Quantity::new(Point3::new(2., 0., 0.), BOHR))
        .unwrap();
    let c = editor
        .add_atom(atom(7), Quantity::new(Point3::new(3., 0., 0.), BOHR))
        .unwrap();
    let removed = editor
        .add_atom(atom(1), Quantity::new(Point3::new(4., 0., 0.), BOHR))
        .unwrap();
    editor.delete_atom(removed).unwrap();
    let bond = editor.add_bond(a, c, BondOrder::Single).unwrap();
    let (model, map) = editor.finish_with_correspondence().unwrap();
    assert!(map.atom(removed).is_none());
    for (handle, x) in [(a, 1.), (b, 2.), (c, 3.)] {
        let (id, index) = map.atom(handle).unwrap();
        assert_eq!(model.atom_ids()[index.index()], id);
        assert!((model.position(id).unwrap().value_in(BOHR).unwrap().x - x).abs() < 1e-12);
    }
    let (id, index) = map.bond(bond).unwrap();
    assert_eq!(model.bond_ids()[index.index()], id);
    let mut next = model.edit();
    let handle = next.bond_handle(id).unwrap();
    next.delete_bond(handle).unwrap();
    let (split, second) = next.finish_with_correspondence().unwrap();
    assert!(second.bond(handle).is_none());
    assert_eq!(split.topology().instance_count(), 3);
}
#[test]
fn no_op_publication_retains_snapshot_and_foreign_handles_do_not_map() {
    let mut e = ModelEditor::new();
    e.add_atom(atom(2), Quantity::new(Point3::origin(), BOHR))
        .unwrap();
    let model = e.finish().unwrap();
    let source = model.shared_topology();
    let edit = TopologyEditor::from_topology(source.clone());
    let handle = edit.atom_handle(model.atom_ids()[0]).unwrap();
    let (topology, map) = edit.finish_with_correspondence().unwrap();
    assert!(Arc::ptr_eq(&source, &topology));
    assert_eq!(map.atom(handle).unwrap().0, model.atom_ids()[0]);
    let mut foreign = TopologyEditor::new();
    let other = foreign.add_atom(atom(2)).unwrap();
    assert!(map.atom(other).is_none());
}
#[test]
fn atomic_masses_are_neutral_and_do_not_require_a_molecule() {
    let carbon = Element::from_atomic_number(6).unwrap();
    assert_eq!(
        kekule::descriptors::atomic_mass(carbon, Some(12))
            .unwrap()
            .value_in(DALTON)
            .unwrap(),
        12.
    );
    assert!(
        (kekule::descriptors::atomic_mass(carbon, Some(13))
            .unwrap()
            .value_in(DALTON)
            .unwrap()
            - 13.003354835)
            .abs()
            < 1e-7
    );
    assert!(kekule::descriptors::atomic_mass(carbon, Some(999)).is_none());
}
#[test]
fn stereo_replacement_is_occurrence_local_and_rejects_changed_source() {
    let molecule = kekule::smiles::to_molecules("F[C@](Cl)(Br)I")
        .unwrap()
        .pop()
        .unwrap();
    assert_eq!(molecule.stereo_elements().count(), 1);
    let mut builder = Model::builder();
    builder
        .add_molecule(&molecule, &Positions::zeros(5))
        .unwrap();
    builder
        .add_molecule(&molecule, &Positions::zeros(5))
        .unwrap();
    let original = builder.build().unwrap();
    let instances = original
        .topology()
        .molecules()
        .map(|m| m.id())
        .collect::<Vec<_>>();
    let mut editor = original.edit();
    editor
        .replace_source_instance_stereo(instances[0], &[])
        .unwrap();
    let result = editor.finish().unwrap();
    let counts = result
        .topology()
        .molecules()
        .map(|m| m.molecule().stereo_elements().count())
        .collect::<Vec<_>>();
    assert_eq!(counts, vec![0, 1]);
    assert!(original
        .topology()
        .molecules()
        .all(|m| m.molecule().stereo_elements().count() == 1));
    let mut changed = original.edit();
    let id = changed.atom_handle(original.atom_ids()[0]).unwrap();
    let mut a = changed.atom(id).unwrap().clone();
    a.isotope = Some(19);
    changed.replace_atom(id, a).unwrap();
    assert!(changed
        .replace_source_instance_stereo(instances[0], &[])
        .is_err());
}
