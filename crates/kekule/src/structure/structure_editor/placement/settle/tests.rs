use crate::structure::{Positions, Structure};

#[test]
fn a_ring_side_that_met_itself_is_weighed_whole() {
    // Butylcyclohexane without hydrogens: at the ring atom that carries the
    // chain, the ring's five other atoms outweigh the chain's four, though
    // the ring's two branches finish first.
    let topology = crate::smiles::to_topology("C1CCCCC1CCCC").unwrap();
    let atoms = topology.atom_count();
    let editor = Structure::new(topology, Positions::zeros(atoms))
        .unwrap()
        .into_editor();
    let center = editor.topology().atom_ids().nth(5).unwrap();
    let others = editor
        .topology()
        .neighbors(center)
        .unwrap()
        .collect::<Vec<_>>();
    let (sides, largest) = editor.neighbor_sides(center, &others).unwrap();
    assert_eq!(sides.len(), 2);
    assert_eq!(sides[largest].neighbors.len(), 2);
    assert_eq!(sides[largest].atoms.len(), 5);
}
