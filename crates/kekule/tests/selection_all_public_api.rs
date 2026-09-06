use std::sync::Arc;

use kekule::{
    smiles,
    topology::{AtomSelection, TopologyBuilder},
};

#[test]
fn all_selects_dense_atoms_across_reused_definitions_and_binds_the_exact_topology() {
    let molecule = smiles::to_molecules("CO").unwrap().pop().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let all = AtomSelection::all(&topology);
    assert_eq!(
        all.indices()
            .iter()
            .map(|id| id.index())
            .collect::<Vec<_>>(),
        [0, 1, 2, 3]
    );
    assert_eq!(
        all,
        AtomSelection::from_atoms(&topology, topology.atom_ids().iter().copied()).unwrap()
    );
    assert!(std::ptr::eq(all.topology(), topology.as_ref()));
    assert!(all.ensure_compatible(&topology).is_ok());
    let mut independent = TopologyBuilder::new();
    let definition = independent.add_molecule_definition(&molecule).unwrap();
    independent.add_instance(definition).unwrap();
    independent.add_instance(definition).unwrap();
    let independent = Arc::new(independent.build().unwrap());
    assert!(topology.same_layout(&independent));
    assert!(all.ensure_compatible(&independent).is_err());
}
