use std::sync::Arc;

use kekule::{
    core::BondOrder,
    query::parse_smarts,
    smiles,
    substructure::{find_topology_substructure_matches_complete, SubstructureMatchOptions},
    topology::{
        AtomSelection, BondSelection, BondSelectionMode, InstanceAtomId, InstanceBondId,
        MoleculeDefinitionId, MoleculeInstanceId, SelectionError, Topology, TopologyAtomIndex,
        TopologyBondIndex, TopologyBuilder,
    },
};

fn topology() -> Arc<Topology> {
    let molecule = smiles::to_molecules("CC=O").unwrap().pop().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    builder
        .add_molecule(&smiles::to_molecules("[Na+]").unwrap()[0])
        .unwrap();
    Arc::new(builder.build().unwrap())
}

fn atoms(top: &Arc<Topology>, mask: u32) -> AtomSelection {
    AtomSelection::from_atoms(
        top,
        top.atom_ids()
            .iter()
            .enumerate()
            .filter(|(index, _)| mask & (1 << index) != 0)
            .map(|(_, id)| *id),
    )
    .unwrap()
}

fn bonds(top: &Arc<Topology>, mask: u32) -> BondSelection {
    BondSelection::from_bonds(
        top,
        top.bond_ids()
            .iter()
            .enumerate()
            .filter(|(index, _)| mask & (1 << index) != 0)
            .map(|(_, id)| *id),
    )
    .unwrap()
}

#[test]
fn set_algebra_matches_independent_bitset_oracle() {
    let top = topology();
    // Exhaust all subsets of four atoms/bonds, including overlapping and
    // disjoint sets, empty operands, and interleaved dense order.
    for left in 0..16 {
        let a = atoms(&top, left);
        let ba = bonds(&top, left);
        assert_eq!(a.complement(), atoms(&top, 0x7f ^ left));
        assert_eq!(ba.complement(), bonds(&top, 0xf ^ left));
        assert_eq!(a.len(), left.count_ones() as usize);
        assert_eq!(ba.len(), left.count_ones() as usize);
        for right in 0..16 {
            let b = atoms(&top, right);
            let bb = bonds(&top, right);
            assert_eq!(a.union(&b).unwrap(), atoms(&top, left | right));
            assert_eq!(a.intersection(&b).unwrap(), atoms(&top, left & right));
            assert_eq!(a.difference(&b).unwrap(), atoms(&top, left & !right));
            assert_eq!(
                a.symmetric_difference(&b).unwrap(),
                atoms(&top, left ^ right)
            );
            assert_eq!(a.is_subset(&b).unwrap(), left & !right == 0);
            assert_eq!(a.is_disjoint(&b).unwrap(), left & right == 0);
            assert_eq!(ba.union(&bb).unwrap(), bonds(&top, left | right));
            assert_eq!(ba.intersection(&bb).unwrap(), bonds(&top, left & right));
            assert_eq!(ba.difference(&bb).unwrap(), bonds(&top, left & !right));
            assert_eq!(
                ba.symmetric_difference(&bb).unwrap(),
                bonds(&top, left ^ right)
            );
            assert_eq!(ba.is_subset(&bb).unwrap(), left & !right == 0);
            assert_eq!(ba.is_disjoint(&bb).unwrap(), left & right == 0);
        }
    }
}

#[test]
fn construction_and_click_mutations_validate_sort_and_deduplicate() {
    let top = topology();
    let a = top.atom_ids();
    let b = top.bond_ids();
    let mut selected = AtomSelection::empty(&top);
    let mut selected_bonds = BondSelection::empty(&top);
    assert!(selected.is_empty());
    assert!(selected_bonds.is_empty());
    for index in [2, 0, 1] {
        assert!(selected.insert(a[index]).unwrap());
        assert!(!selected.insert(a[index]).unwrap());
        assert!(selected_bonds.insert(b[index]).unwrap());
        assert!(!selected_bonds.insert(b[index]).unwrap());
    }
    assert_eq!(selected.atom_ids().collect::<Vec<_>>(), a[..3]);
    assert_eq!(selected_bonds.bond_ids().collect::<Vec<_>>(), b[..3]);
    assert_eq!(selected_bonds.bond_ids().len(), 3);
    assert_eq!(
        BondSelection::from_bonds(&top, [b[2], b[0], b[1], b[2]]).unwrap(),
        selected_bonds
    );
    assert_eq!(
        BondSelection::from_indices(&top, [2, 0, 1, 2].map(TopologyBondIndex::new)).unwrap(),
        selected_bonds
    );
    assert_eq!(
        AtomSelection::from_indices(&top, [2, 0, 1, 2].map(TopologyAtomIndex::new)).unwrap(),
        selected
    );
    assert!(selected.remove(a[1]).unwrap());
    assert!(!selected.remove(a[1]).unwrap());
    assert!(selected.toggle(a[1]).unwrap());
    assert!(!selected.toggle(a[1]).unwrap());
    assert!(selected_bonds.remove(b[1]).unwrap());
    assert!(!selected_bonds.remove(b[1]).unwrap());
    assert!(selected_bonds.toggle(b[1]).unwrap());
    assert!(!selected_bonds.toggle(b[1]).unwrap());
    for (index, id) in a.iter().enumerate() {
        assert_eq!(selected.contains(*id), index == 0 || index == 2);
        assert_eq!(
            selected.contains_index(TopologyAtomIndex::new(index as u32)),
            index == 0 || index == 2
        );
    }
    for (index, id) in b.iter().enumerate() {
        assert_eq!(selected_bonds.contains(*id), index == 0 || index == 2);
        assert_eq!(
            selected_bonds.contains_index(TopologyBondIndex::new(index as u32)),
            index == 0 || index == 2
        );
    }
    let invalid_atom = InstanceAtomId::new(MoleculeInstanceId::new(u32::MAX), a[0].atom());
    let invalid_bond = InstanceBondId::new(MoleculeInstanceId::new(u32::MAX), b[0].bond());
    let before = selected.clone();
    let bonds_before = selected_bonds.clone();
    assert_eq!(
        selected.insert(invalid_atom),
        Err(SelectionError::InvalidAtomId(invalid_atom))
    );
    assert_eq!(
        selected.remove(invalid_atom),
        Err(SelectionError::InvalidAtomId(invalid_atom))
    );
    assert_eq!(
        selected.toggle(invalid_atom),
        Err(SelectionError::InvalidAtomId(invalid_atom))
    );
    assert_eq!(
        selected_bonds.insert(invalid_bond),
        Err(SelectionError::InvalidBondId(invalid_bond))
    );
    assert_eq!(
        selected_bonds.remove(invalid_bond),
        Err(SelectionError::InvalidBondId(invalid_bond))
    );
    assert_eq!(
        selected_bonds.toggle(invalid_bond),
        Err(SelectionError::InvalidBondId(invalid_bond))
    );
    assert_eq!(selected, before);
    assert_eq!(selected_bonds, bonds_before);
    assert!(!selected.contains(invalid_atom));
    assert!(!selected_bonds.contains(invalid_bond));
    assert!(!selected.contains_index(TopologyAtomIndex::new(u32::MAX)));
    assert!(!selected_bonds.contains_index(TopologyBondIndex::new(u32::MAX)));
    assert_eq!(
        BondSelection::from_bonds(&top, [b[0], invalid_bond]),
        Err(SelectionError::InvalidBondId(invalid_bond))
    );
    assert_eq!(
        BondSelection::from_indices(&top, [TopologyBondIndex::new(4)]),
        Err(SelectionError::InvalidBondIndex(TopologyBondIndex::new(4)))
    );
    selected.clear();
    selected_bonds.clear();
    assert_eq!(selected, AtomSelection::empty(&top));
    assert_eq!(selected_bonds, BondSelection::empty(&top));
    assert!(Arc::ptr_eq(&selected.shared_topology(), &top));
    assert!(Arc::ptr_eq(&selected_bonds.shared_topology(), &top));
    assert!(std::ptr::eq(selected_bonds.topology(), top.as_ref()));
}

#[test]
fn distinct_snapshots_are_rejected_even_when_empty_and_layout_equal() {
    let top = topology();
    let other = topology();
    for selected in [AtomSelection::empty(&top), AtomSelection::all(&top)] {
        let foreign = AtomSelection::empty(&other);
        assert_ne!(selected, foreign);
        for result in [
            selected.union(&foreign),
            selected.intersection(&foreign),
            selected.difference(&foreign),
            selected.symmetric_difference(&foreign),
        ] {
            assert_eq!(result, Err(SelectionError::TopologyMismatch));
        }
        assert_eq!(
            selected.is_subset(&foreign),
            Err(SelectionError::TopologyMismatch)
        );
        assert_eq!(
            selected.is_disjoint(&foreign),
            Err(SelectionError::TopologyMismatch)
        );
    }
    for selected in [BondSelection::empty(&top), BondSelection::all(&top)] {
        let foreign = BondSelection::empty(&other);
        assert_ne!(selected, foreign);
        for result in [
            selected.union(&foreign),
            selected.intersection(&foreign),
            selected.difference(&foreign),
            selected.symmetric_difference(&foreign),
        ] {
            assert_eq!(result, Err(SelectionError::TopologyMismatch));
        }
        assert_eq!(
            selected.is_subset(&foreign),
            Err(SelectionError::TopologyMismatch)
        );
        assert_eq!(
            selected.is_disjoint(&foreign),
            Err(SelectionError::TopologyMismatch)
        );
        assert_eq!(
            selected.ensure_compatible(&other),
            Err(SelectionError::TopologyMismatch)
        );
    }
}

#[test]
fn bond_endpoint_rules_keep_occurrences_distinct_and_cover_boundary_edges() {
    let top = topology();
    // Only the first occurrence: C-C selected, O unselected.
    let selected = atoms(&top, 0b11);
    assert_eq!(
        selected.to_bonds(BondSelectionMode::Internal),
        bonds(&top, 0b1)
    );
    assert_eq!(
        selected.to_bonds(BondSelectionMode::Incident),
        bonds(&top, 0b11)
    );
    assert_eq!(
        selected.to_bonds(BondSelectionMode::Boundary),
        bonds(&top, 0b10)
    );
    assert_eq!(bonds(&top, 0b1010).to_atoms(), atoms(&top, 0b110110));
    assert_eq!(BondSelection::all(&top).to_atoms(), atoms(&top, 0b111111));
    assert_eq!(
        BondSelection::empty(&top).to_atoms(),
        AtomSelection::empty(&top)
    );
    for mode in [
        BondSelectionMode::Internal,
        BondSelectionMode::Incident,
        BondSelectionMode::Boundary,
    ] {
        assert_eq!(
            AtomSelection::empty(&top).to_bonds(mode),
            BondSelection::empty(&top)
        );
    }
    assert_eq!(
        AtomSelection::all(&top).to_bonds(BondSelectionMode::Internal),
        BondSelection::all(&top)
    );
    assert!(AtomSelection::all(&top)
        .to_bonds(BondSelectionMode::Boundary)
        .is_empty());
    let instance = top.atom_ids()[3].molecule();
    assert_eq!(
        BondSelection::for_instances(&top, [instance, instance]).unwrap(),
        bonds(&top, 0b1100)
    );
    let definition = top.instance(instance).unwrap().definition();
    assert_eq!(
        BondSelection::for_definitions(&top, [definition]).unwrap(),
        BondSelection::all(&top)
    );
    assert_eq!(
        BondSelection::for_instances(&top, [MoleculeInstanceId::new(u32::MAX)]),
        Err(SelectionError::InvalidMoleculeInstanceId(
            MoleculeInstanceId::new(u32::MAX)
        ))
    );
    assert_eq!(
        BondSelection::for_definitions(&top, [MoleculeDefinitionId::new(u32::MAX)]),
        Err(SelectionError::InvalidMoleculeDefinitionId(
            MoleculeDefinitionId::new(u32::MAX)
        ))
    );
}

#[test]
fn bond_to_atom_round_trip_may_add_unselected_ring_bonds() {
    let top = Arc::new(smiles::to_topology("C1CC1").unwrap());
    let selected = BondSelection::from_bonds(&top, top.bond_ids()[..2].iter().copied()).unwrap();
    assert_eq!(selected.len(), 2);
    assert_eq!(selected.to_atoms().len(), 3);
    assert_eq!(
        selected
            .to_atoms()
            .to_bonds(BondSelectionMode::Internal)
            .len(),
        3
    );
    let isolated = Arc::new(smiles::to_topology("[Na+].[Cl-]").unwrap());
    assert!(BondSelection::all(&isolated).is_empty());
    assert!(BondSelection::all(&isolated).complement().is_empty());
    assert!(AtomSelection::all(&isolated)
        .to_bonds(BondSelectionMode::Incident)
        .is_empty());
}

#[test]
fn predicates_and_graph_expansion_respect_membership_and_instance_boundaries() {
    let top = topology();
    let oxygen = AtomSelection::from_predicate(&top, |_, atom| atom.element.symbol() == "O");
    assert_eq!(oxygen, atoms(&top, 0b100100));
    assert_eq!(
        oxygen.filter(|id, _| id.molecule() == top.atom_ids()[0].molecule()),
        atoms(&top, 0b100)
    );
    let doubles = BondSelection::from_predicate(&top, |_, bond| bond.order == BondOrder::Double);
    assert_eq!(doubles, bonds(&top, 0b1010));
    assert_eq!(
        doubles.filter(|id, _| id.molecule() == top.atom_ids()[0].molecule()),
        bonds(&top, 0b10)
    );
    assert_eq!(atoms(&top, 1).expand_bonded(0), atoms(&top, 1));
    assert_eq!(atoms(&top, 1).expand_bonded(1), atoms(&top, 0b11));
    assert_eq!(atoms(&top, 1).expand_bonded(2), atoms(&top, 0b111));
    assert_eq!(atoms(&top, 1).expand_bonded(usize::MAX), atoms(&top, 0b111));
    assert_eq!(
        atoms(&top, 0b1000000).expand_bonded(usize::MAX),
        atoms(&top, 0b1000000)
    );
    assert_eq!(
        atoms(&top, 0b1000001).expand_to_instances(),
        atoms(&top, 0b1000111)
    );
    let empty = AtomSelection::empty(&top);
    assert_eq!(empty.expand_bonded(usize::MAX), empty);
    assert_eq!(empty.expand_to_instances(), empty);
    assert_eq!(
        empty.filter(|_, _| panic!("empty predicate must not run")),
        empty
    );
    assert_eq!(
        BondSelection::empty(&top).filter(|_, _| panic!("empty predicate must not run")),
        BondSelection::empty(&top)
    );
}

#[test]
fn topology_query_match_conversion_checks_provenance_before_selecting() {
    let top = Arc::new(topology().perceived().unwrap());
    let query = parse_smarts("[#6]=[#8]").unwrap();
    let matches = find_topology_substructure_matches_complete(
        &top,
        &query,
        SubstructureMatchOptions::default(),
    )
    .unwrap();
    assert_eq!(
        AtomSelection::from_topology_query_matches(&top, &matches).unwrap(),
        atoms(&top, 0b110110)
    );
    assert_eq!(
        AtomSelection::from_topology_query_matches(&top, &[]).unwrap(),
        AtomSelection::empty(&top)
    );
    let foreign = topology();
    assert_eq!(
        AtomSelection::from_topology_query_matches(&foreign, &matches),
        Err(SelectionError::TopologyMismatch)
    );
}
