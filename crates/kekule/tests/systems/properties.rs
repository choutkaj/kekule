use kekule::diagnostics::Countable;
use kekule::diagnostics::LengthMismatch;
use std::sync::Arc;

use kekule::molecule::Molecule;
use kekule::properties::{PropertyColumn, PropertyError, PropertyKey, PropertyValue};
use kekule::structure::{Geometry, GeometryError, Positions, Structure, StructureSet};
use kekule::topology::{
    AtomSelection, AtomSiteMetadata, HierarchyBuilder, InstanceAtomId, Topology,
    TopologyBuildError, TopologyBuilder,
};

fn key() -> PropertyKey {
    PropertyKey::new("tag").unwrap()
}

fn molecule(smiles: &str) -> Molecule {
    kekule::smiles::to_molecules(smiles).unwrap().pop().unwrap()
}

fn builder() -> TopologyBuilder {
    let molecule = molecule("C");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule.clone()).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "UNL", None, None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, molecule.atom_ids().next().unwrap()),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    builder
}

#[test]
fn owner_sized_tables_reject_wrong_columns_and_extend_with_instances() {
    let mut builder = builder();
    builder
        .properties_mut()
        .atoms_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(11)]))
        .unwrap();
    assert!(matches!(
        builder
            .properties_mut()
            .atoms_mut()
            .insert(key(), PropertyColumn::Int(vec![Some(22), Some(33)])),
        Err(PropertyError::LengthMismatch(payload)) if payload == LengthMismatch::new(Countable::Value, 1, 2)
    ));
    let molecule = kekule::smiles::to_molecules("C").unwrap().pop().unwrap();
    builder.add_molecule(molecule.clone()).unwrap();
    let topology = builder.build().unwrap();
    assert_eq!(
        topology.properties().atoms().get(&key()),
        Some(&PropertyColumn::Int(vec![Some(11), None]))
    );
}

#[test]
fn replacing_hierarchy_cannot_silently_truncate_populated_columns() {
    let mut builder = builder();
    builder
        .properties_mut()
        .chains_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(11)]))
        .unwrap();
    *builder.hierarchy_mut() = HierarchyBuilder::new();
    assert_eq!(
        builder.properties_mut().chains_mut().get(&key()),
        Some(&PropertyColumn::Int(vec![Some(11)]))
    );
    assert!(matches!(
        builder.build(),
        Err(TopologyBuildError::Property(error)) if matches!(*error, PropertyError::LengthMismatch(payload) if payload == LengthMismatch::new(Countable::Value, 0, 1))
    ));
}

#[test]
fn topology_only_domains_are_typed_and_stay_out_of_structures() {
    let mut builder = builder();
    builder
        .properties_mut()
        .molecule_instances_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(1)]))
        .unwrap();
    builder
        .properties_mut()
        .chains_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(2)]))
        .unwrap();
    builder
        .properties_mut()
        .residues_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(3)]))
        .unwrap();
    builder
        .properties_mut()
        .atom_sites_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(4)]))
        .unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let properties = topology.properties();
    let instance = topology.instances().next().unwrap().id();
    let chain = topology.chains().next().unwrap().id();
    let residue = topology.residues().next().unwrap().id();
    let site = topology.atom_sites().next().unwrap().id();
    assert_eq!(
        properties
            .molecule_instances()
            .value(&key(), instance)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        properties.chains().value(&key(), chain).unwrap(),
        Some(PropertyValue::Int(2))
    );
    assert_eq!(
        properties.residues().value(&key(), residue).unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert_eq!(
        properties.atom_sites().value(&key(), site).unwrap(),
        Some(PropertyValue::Int(4))
    );
    // Structure annotations have no hierarchy domains at all.
    let structure = Structure::new(Arc::clone(&topology), Positions::zeros(1)).unwrap();
    assert!(structure.properties().is_empty());
}

#[test]
fn structures_reject_mismatched_property_rows_transactionally() {
    let methane = Arc::new(builder().build().unwrap());
    let ethane = Arc::new(Topology::from_molecule(molecule("C=C")).unwrap());
    let mut source = Structure::new(Arc::clone(&ethane), Positions::zeros(2)).unwrap();
    source
        .geometry_mut()
        .properties_mut()
        .atoms_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    source
        .geometry_mut()
        .properties_mut()
        .bonds_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(3)]))
        .unwrap();
    let bound = source.properties().clone();

    let mut structure = Structure::new(Arc::clone(&methane), Positions::zeros(1)).unwrap();
    let original = structure.properties().clone();
    assert_eq!(
        structure.geometry_mut().set_properties(bound.clone()),
        Err(GeometryError::Array(
            kekule::structure::DenseArrayError::CountMismatch(LengthMismatch::new(
                Countable::Row,
                1,
                2
            ))
        ))
    );
    assert_eq!(structure.properties(), &original);

    // A detached geometry has no bond rows yet.
    let mut member = Geometry::new(Positions::zeros(2));
    let detached = member.properties().clone();
    assert_eq!(
        member.set_properties(bound.clone()),
        Err(GeometryError::BondCountMismatch(LengthMismatch::new(
            Countable::Row,
            0,
            1
        )))
    );
    assert_eq!(member.properties(), &detached);

    let mut set = StructureSet::from_geometries(Arc::clone(&ethane), [member]).unwrap();
    assert_eq!(set.get(0).unwrap().properties().bonds().len(), 1);
    set.get_mut(0)
        .unwrap()
        .set_properties(bound.clone())
        .unwrap();
    assert_eq!(set.get(0).unwrap().properties(), &bound);
    assert!(set.get_mut(0).unwrap().set_properties(original).is_err());
    assert_eq!(set.get(0).unwrap().properties(), &bound);
}

#[test]
fn owner_topology_property_edits_keep_layout_without_touching_other_holders() {
    let topology = kekule::smiles::to_topology("CCO").unwrap();
    let selection = AtomSelection::all(Arc::clone(&topology));

    let mut structure = Structure::new(Arc::clone(&topology), Positions::zeros(3)).unwrap();
    structure
        .topology_properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(1))
        .unwrap();
    // A shared topology is replaced by an edited snapshot of the same layout.
    assert!(!Arc::ptr_eq(&topology, &structure.shared_topology()));
    assert!(structure.topology().shares_layout(&topology));
    assert_eq!(
        structure.topology().properties().owner().get(&key()),
        Some(&PropertyValue::Int(1))
    );
    assert_eq!(topology.properties().owner().get(&key()), None);
    selection
        .ensure_compatible(&structure.shared_topology())
        .unwrap();
    // A uniquely owned snapshot is edited in place.
    let unique = Arc::as_ptr(&structure.shared_topology());
    structure
        .topology_properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(2))
        .unwrap();
    assert_eq!(Arc::as_ptr(&structure.shared_topology()), unique);
    assert_eq!(
        structure.topology().properties().owner().get(&key()),
        Some(&PropertyValue::Int(2))
    );

    let set = || {
        StructureSet::from_geometries(Arc::clone(&topology), [Geometry::new(Positions::zeros(3))])
            .unwrap()
    };
    let (mut first, mut second) = (set(), set());
    first
        .topology_properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(3))
        .unwrap();
    second
        .topology_properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(4))
        .unwrap();

    assert!(first.topology().shares_layout(&topology));
    assert!(second.topology().shares_layout(&topology));
    assert_eq!(
        first.topology().properties().owner().get(&key()),
        Some(&PropertyValue::Int(3))
    );
    assert_eq!(
        second.topology().properties().owner().get(&key()),
        Some(&PropertyValue::Int(4))
    );
    assert_eq!(topology.properties().owner().get(&key()), None);
}
