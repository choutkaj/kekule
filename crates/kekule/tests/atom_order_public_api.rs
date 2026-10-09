//! Explicit dense atom order: construction, interpretation, editing, subsets.

use std::sync::Arc;

use kekule::core::{Atom, AtomId, BondOrder, Element};
use kekule::geometry::Point3;
use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::smiles;
use kekule::structure::{Model, ModelBuilder, Positions};
use kekule::topology::{
    transform, AtomSelection, InstanceAtomId, MoleculeInstanceId, Topology, TopologyAtomIndex,
    TopologyBuildError, TopologyBuilder,
};
use kekule::units::{Quantity, ANGSTROM};

fn molecule(source: &str) -> kekule::core::Molecule {
    smiles::to_molecules(source).unwrap().remove(0)
}

fn positions(xs: &[f64]) -> Positions {
    Positions::new(Quantity::new(
        xs.iter()
            .map(|&x| Point3::new(x, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap()
}

fn xs(model: &Model) -> Vec<f64> {
    model
        .positions()
        .values()
        .value()
        .iter()
        .map(|point| (point.x * 10.0).round())
        .collect()
}

fn atom(instance: MoleculeInstanceId, raw: u32) -> InstanceAtomId {
    InstanceAtomId::new(instance, AtomId::new(raw))
}

/// Methanol C and O with a water oxygen between them, at x = 0, 1, 2 Å.
fn interleaved() -> (Model, [InstanceAtomId; 3]) {
    let mut builder = ModelBuilder::new();
    let methanol = builder
        .add_molecule(&molecule("CO"), &positions(&[0.0, 2.0]))
        .unwrap();
    let water = builder
        .add_molecule(&molecule("O"), &positions(&[1.0]))
        .unwrap();
    let order = [atom(methanol, 0), atom(water, 0), atom(methanol, 1)];
    builder.set_atom_order(order).unwrap();
    (builder.build().unwrap(), order)
}

#[test]
fn builder_order_is_a_validated_permutation_carrying_atom_rows() {
    let mut builder = TopologyBuilder::new();
    let methanol = builder.add_molecule(&molecule("CO")).unwrap();
    let water = builder.add_molecule(&molecule("O")).unwrap();
    let [carbon, oxygen, water_oxygen] = [atom(methanol, 0), atom(methanol, 1), atom(water, 0)];
    assert_eq!(builder.atom_ids(), &[carbon, oxygen, water_oxygen]);
    let tag = PropertyKey::new("row").unwrap();
    builder
        .atom_properties_mut()
        .insert(
            tag.clone(),
            PropertyColumn::Int(vec![Some(0), Some(1), Some(2)]),
        )
        .unwrap();

    let before = builder.clone();
    for (order, error) in [
        (
            vec![carbon, water_oxygen],
            TopologyBuildError::AtomOrderLengthMismatch {
                expected: 3,
                actual: 2,
            },
        ),
        (
            vec![carbon, water_oxygen, atom(water, 1)],
            TopologyBuildError::InvalidAtomOrderEntry(atom(water, 1)),
        ),
        (
            vec![carbon, water_oxygen, carbon],
            TopologyBuildError::DuplicateAtomOrderEntry(carbon),
        ),
        (
            vec![carbon, water_oxygen, atom(MoleculeInstanceId::new(9), 0)],
            TopologyBuildError::InvalidAtomOrderEntry(atom(MoleculeInstanceId::new(9), 0)),
        ),
    ] {
        assert_eq!(builder.set_atom_order(order), Err(error));
        assert_eq!(builder, before);
    }

    builder
        .set_atom_order([carbon, water_oxygen, oxygen])
        .unwrap();
    // A later instance appends after the explicit order.
    let ion = builder.add_molecule(&molecule("[Na+]")).unwrap();
    let topology = builder.build().unwrap();
    assert_eq!(
        topology.atom_ids(),
        &[carbon, water_oxygen, oxygen, atom(ion, 0)]
    );
    for (dense, &id) in topology.atom_ids().iter().enumerate() {
        assert_eq!(
            topology.atom_index(id),
            Some(TopologyAtomIndex::new(dense as u32))
        );
    }
    assert_eq!(
        topology.atom_properties().get(&tag),
        Some(&PropertyColumn::Int(vec![Some(0), Some(2), Some(1), None]))
    );
    // Bonds stay in instance order regardless of atom order.
    assert_eq!(topology.bond_ids().len(), 1);
    assert_eq!(topology.bond_ids()[0].molecule(), methanol);

    // Append-oriented rebuilding keeps the explicit order.
    let default_order =
        Topology::from_molecules(&[molecule("CO"), molecule("O"), molecule("[Na+]")]).unwrap();
    assert!(!topology.same_layout(&default_order));
    let ammonia = molecule("N");
    let mut extended = topology.into_builder();
    let appended = extended.add_molecule(&ammonia).unwrap();
    let extended = extended.build().unwrap();
    assert_eq!(
        extended.atom_ids(),
        &[
            carbon,
            water_oxygen,
            oxygen,
            atom(ion, 0),
            atom(appended, 0)
        ]
    );
}

#[test]
fn model_builder_order_moves_positions_with_atoms() {
    let (model, order) = interleaved();
    assert_eq!(model.atom_ids(), order);
    assert_eq!(xs(&model), [0.0, 1.0, 2.0]);
    for (atom, x) in order.into_iter().zip([0.0, 1.0, 2.0]) {
        assert_eq!((model.position(atom).unwrap().x * 10.0).round(), x);
    }
}

#[test]
fn smiles_topology_keeps_source_order_across_dot_separated_components() {
    // The ring closure joins the first and last atoms across the dots.
    let topology = smiles::to_topology("C1.O.C1").unwrap();
    assert_eq!(topology.instance_count(), 2);
    let [first, water] = [0, 1].map(MoleculeInstanceId::new);
    assert_eq!(
        topology.atom_ids(),
        &[atom(first, 0), atom(water, 0), atom(first, 1)]
    );
}

#[test]
fn editors_keep_surviving_atoms_in_source_order_and_append_new_atoms() {
    let (model, [carbon, water_oxygen, oxygen]) = interleaved();
    // Editing handles follow dense order, so editor and model positions agree.
    let editor = model.edit();
    assert_eq!(editor.positions(), *model.positions());
    let unchanged = editor.finish().unwrap();
    assert!(Arc::ptr_eq(
        &unchanged.shared_topology(),
        &model.shared_topology()
    ));

    let mut editor = model.edit();
    let [carbon, water_oxygen, oxygen] =
        [carbon, water_oxygen, oxygen].map(|atom| editor.atom_handle(atom).unwrap());
    // Merging the water into methanol keeps every atom in place.
    editor
        .add_bond(oxygen, water_oxygen, BondOrder::Single)
        .unwrap();
    let added = editor
        .add_atom(
            Atom::new(Element::from_symbol("N").unwrap()),
            Quantity::new(Point3::new(3.0, 0.0, 0.0), ANGSTROM),
        )
        .unwrap();
    editor.add_bond(carbon, added, BondOrder::Single).unwrap();
    let (merged, ids) = editor.finish_with_correspondence().unwrap();
    assert_eq!(merged.topology().instance_count(), 1);
    assert_eq!(xs(&merged), [0.0, 1.0, 2.0, 3.0]);
    for (dense, handle) in [carbon, water_oxygen, oxygen, added]
        .into_iter()
        .enumerate()
    {
        assert_eq!(ids.atom(handle).unwrap().1.index(), dense);
    }

    let mut editor = model.edit();
    let water_oxygen = editor.atom_handle(model.atom_ids()[1]).unwrap();
    editor.delete_atom(water_oxygen).unwrap();
    let removed = editor.finish().unwrap();
    assert_eq!(xs(&removed), [0.0, 2.0]);
    assert_eq!(
        removed
            .topology()
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "O"]
    );
}

#[test]
fn appending_a_model_keeps_its_dense_order_and_positions() {
    let (source, _) = interleaved();
    let mut editor = kekule::structure::ModelEditor::new();
    let appended = editor.append_model(&source).unwrap();
    for (dense, &atom) in source.atom_ids().iter().enumerate() {
        assert_eq!(
            editor.position(appended.atom(atom).unwrap()).unwrap(),
            source.position(atom).unwrap(),
            "dense atom {dense}"
        );
    }
    let model = editor.finish().unwrap();
    assert_eq!(xs(&model), [0.0, 1.0, 2.0]);
    assert_eq!(
        model
            .topology()
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "O", "O"]
    );
}

#[test]
fn subsets_and_instance_filters_keep_relative_source_order() {
    let (model, [carbon, water_oxygen, oxygen]) = interleaved();
    let topology = model.shared_topology();

    let selection = AtomSelection::from_atoms(&topology, [oxygen, water_oxygen]).unwrap();
    let subset = topology.subset(&selection).unwrap();
    assert_eq!(
        subset.correspondence().source_atom_indices(),
        &[TopologyAtomIndex::new(1), TopologyAtomIndex::new(2)]
    );
    let sliced = model.slice(&selection).unwrap();
    assert_eq!(xs(&sliced), [1.0, 2.0]);

    let methanol = transform::retain_instances(&topology, [carbon.molecule()]).unwrap();
    assert_eq!(
        methanol
            .atoms()
            .map(|(_, atom)| atom.element.symbol())
            .collect::<Vec<_>>(),
        ["C", "O"]
    );
    let water = transform::remove_instances(&topology, [carbon.molecule()]).unwrap();
    assert_eq!(water.atom_count(), 1);
}
