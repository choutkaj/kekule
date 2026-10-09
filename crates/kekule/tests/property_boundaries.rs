use std::sync::Arc;

use kekule::core::Molecule;
use kekule::properties::{PropertyColumn, PropertyError, PropertyKey, PropertyValue};
use kekule::structure::{ConformationError, Ensemble, EnsembleMember, Model, ModelView, Positions};
use kekule::topology::{
    AtomSiteMetadata, Hierarchy, InstanceAtomId, Topology, TopologyBuildError, TopologyBuilder,
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
        Err(PropertyError::LengthMismatch {
            expected: 1,
            actual: 2
        })
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
    *builder.hierarchy_mut() = Hierarchy::new();
    assert_eq!(
        builder.properties_mut().chains_mut().get(&key()),
        Some(&PropertyColumn::Int(vec![Some(11)]))
    );
    assert!(matches!(
        builder.build(),
        Err(TopologyBuildError::Property(error)) if matches!(*error, PropertyError::LengthMismatch { expected: 0, actual: 1 })
    ));
}

#[test]
fn topology_only_domains_are_typed_and_stay_out_of_realizations() {
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
    let instance = topology.molecules().next().unwrap().id();
    let chain = topology.hierarchy().chains().next().unwrap().0;
    let residue = topology.hierarchy().residues().next().unwrap().0;
    let site = topology.hierarchy().atom_sites().next().unwrap().0;
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
    // Realization annotations have no hierarchy domains at all.
    let model = Model::new(Arc::clone(&topology), Positions::zeros(1)).unwrap();
    assert!(model.properties().is_empty());
}

#[test]
fn realizations_reject_mismatched_property_rows_transactionally() {
    let methane = Arc::new(builder().build().unwrap());
    let ethane = Arc::new(Topology::from_molecule(molecule("C=C")).unwrap());
    let mut source = Model::new(Arc::clone(&ethane), Positions::zeros(2)).unwrap();
    source
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    source
        .conformation_mut()
        .properties_mut()
        .bonds_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(3)]))
        .unwrap();
    let bound = source.properties().clone();

    let mut model = Model::new(Arc::clone(&methane), Positions::zeros(1)).unwrap();
    let original = model.properties().clone();
    assert_eq!(
        model.conformation_mut().set_properties(bound.clone()),
        Err(ConformationError::AtomCountMismatch {
            expected: 1,
            actual: 2
        })
    );
    assert_eq!(model.properties(), &original);

    // A detached member has no bond rows yet.
    let mut member = EnsembleMember::new(Positions::zeros(2));
    let detached = member.properties().clone();
    assert_eq!(
        member.conformation_mut().set_properties(bound.clone()),
        Err(ConformationError::BondCountMismatch {
            expected: 0,
            actual: 1
        })
    );
    assert_eq!(member.properties(), &detached);
    assert!(ModelView::new(&ethane, member.conformation()).is_err());

    let mut ensemble = Ensemble::from_items(Arc::clone(&ethane), [member]).unwrap();
    assert_eq!(ensemble.get(0).unwrap().properties().bonds().len(), 1);
    ensemble
        .get_mut(0)
        .unwrap()
        .conformation_mut()
        .set_properties(bound.clone())
        .unwrap();
    assert_eq!(ensemble.get(0).unwrap().properties(), &bound);
    assert!(ensemble
        .get_mut(0)
        .unwrap()
        .conformation_mut()
        .set_properties(original)
        .is_err());
    assert_eq!(ensemble.get(0).unwrap().properties(), &bound);
    assert!(ModelView::new(&ethane, ensemble.get(0).unwrap().conformation()).is_ok());
}
