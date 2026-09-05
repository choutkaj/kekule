use std::sync::Arc;

use kekule::properties::{PropertyColumn, PropertyError, PropertyKey};
use kekule::structure::{Ensemble, EnsembleMember, Model, ModelView, Positions};
use kekule::topology::{
    AtomSiteMetadata, Hierarchy, InstanceAtomId, TopologyBuildError, TopologyBuilder,
};

fn key() -> PropertyKey {
    PropertyKey::new("tag").unwrap()
}

fn builder() -> TopologyBuilder {
    let molecule = kekule::smiles::to_molecules("C").unwrap().pop().unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
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
        .atom_properties_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(11)]))
        .unwrap();
    assert!(matches!(
        builder
            .atom_properties_mut()
            .insert(key(), PropertyColumn::Int(vec![Some(22), Some(33)])),
        Err(PropertyError::LengthMismatch {
            expected: 1,
            actual: 2
        })
    ));
    let molecule = kekule::smiles::to_molecules("C").unwrap().pop().unwrap();
    builder.add_molecule(&molecule).unwrap();
    let topology = builder.build().unwrap();
    assert_eq!(
        topology.atom_properties().get(&key()),
        Some(&PropertyColumn::Int(vec![Some(11), None]))
    );
}

#[test]
fn replacing_hierarchy_cannot_silently_truncate_populated_columns() {
    let mut builder = builder();
    builder
        .chain_properties_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(11)]))
        .unwrap();
    *builder.hierarchy_mut() = Hierarchy::new();
    assert_eq!(
        builder.chain_properties_mut().get(&key()),
        Some(&PropertyColumn::Int(vec![Some(11)]))
    );
    assert!(matches!(
        builder.build(),
        Err(TopologyBuildError::Property(error)) if matches!(*error, PropertyError::LengthMismatch { expected: 0, actual: 1 })
    ));
}

#[test]
fn realizations_reject_each_populated_topology_only_domain_transactionally() {
    for domain in 0..4 {
        let mut builder = builder();
        let mut table = match domain {
            0 => builder.molecule_instance_properties_mut(),
            1 => builder.chain_properties_mut(),
            2 => builder.residue_properties_mut(),
            _ => builder.atom_site_properties_mut(),
        };
        table
            .insert(key(), PropertyColumn::Int(vec![Some(11)]))
            .unwrap();
        let topology = Arc::new(builder.build().unwrap());
        let foreign = topology.properties().clone();
        assert!(matches!(
            foreign.validate_realization_properties(),
            Err(PropertyError::InvalidRealizationDomain(_))
        ));
        let mut model = Model::new(Arc::clone(&topology), Positions::zeros(1)).unwrap();
        let original = model.properties().clone();
        assert!(model.set_properties(foreign.clone()).is_err());
        assert_eq!(model.properties(), &original);
        assert!(ModelView::new(&topology, model.positions(), None, &foreign).is_err());
        let mut member = EnsembleMember::new(Positions::zeros(1));
        assert!(member.set_properties(foreign.clone()).is_err());
        assert_eq!(member.properties(), &original);
        let mut ensemble = Ensemble::from_members(topology, [member]).unwrap();
        assert!(ensemble
            .member_mut(0)
            .unwrap()
            .set_properties(foreign)
            .is_err());
        assert_eq!(ensemble.member(0).unwrap().properties(), &original);
    }
}
