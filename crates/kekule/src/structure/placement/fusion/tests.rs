use super::*;
use crate::structure::placement::fragment;

/// The sorted locants of the fusion bond's atoms of the fragment `name`.
fn fusion_locants(name: &str) -> Option<[String; 2]> {
    let fragment = fragment(name).unwrap();
    let bond = fragment.molecule().bond(fragment.fusion_bond()?).unwrap();
    let mut locants = [bond.a(), bond.b()].map(|atom| fragment.locant(atom).unwrap().to_owned());
    locants.sort();
    Some(locants)
}

fn end(hydrogens: usize) -> FusionEnd {
    FusionEnd {
        hydrogens,
        multiple_elsewhere: false,
    }
}

/// The fragment's molecule and its fusion bond's atoms, first atom first.
fn fusion_of(name: &str) -> (&'static Molecule, [AtomId; 2]) {
    let fragment = fragment(name).unwrap();
    let molecule = fragment.molecule();
    let bond = molecule.bond(fragment.fusion_bond().unwrap()).unwrap();
    (molecule, [bond.a(), bond.b()])
}

fn doubles(orders: &FusionOrders) -> usize {
    orders
        .bonds
        .values()
        .filter(|order| **order == BondOrder::Double)
        .count()
}

#[test]
fn the_fusion_bond_joins_two_ring_carbons_away_from_heteroatoms() {
    for (name, locants) in [
        ("benzene", ["1", "2"]),
        ("cyclohexane", ["1", "2"]),
        ("pyridine", ["2", "3"]),
        ("pyrimidine", ["4", "5"]),
        ("pyrrole", ["2", "3"]),
        ("thiophene", ["2", "3"]),
        ("indole", ["2", "3"]),
    ] {
        assert_eq!(
            fusion_locants(name),
            Some(locants.map(str::to_owned)),
            "{name}"
        );
    }
    assert_eq!(fusion_locants("methyl"), None);
    assert_eq!(fusion_locants("phenyl"), None);
}

#[test]
fn benzene_on_either_kekule_bond_of_benzene_keeps_the_bond_order() {
    let (molecule, fused) = fusion_of("benzene");
    // On a double bond, the bond pairs the fused atoms.
    let orders = fusion_orders(molecule, fused, BondOrder::Double, [end(1), end(1)]).unwrap();
    assert_eq!(orders.shared, BondOrder::Double);
    assert_eq!(orders.gains, [1, 1]);
    assert_eq!(orders.rank.0, 2);
    assert_eq!(doubles(&orders), 2);
    // On a single bond, each bond atom has its double bond elsewhere, and
    // the four new atoms pair among themselves.
    let paired = FusionEnd {
        hydrogens: 1,
        multiple_elsewhere: true,
    };
    let orders = fusion_orders(molecule, fused, BondOrder::Single, [paired, paired]).unwrap();
    assert_eq!(orders.shared, BondOrder::Single);
    assert_eq!(orders.gains, [1, 1]);
    assert_eq!(orders.rank.0, 0);
    assert_eq!(doubles(&orders), 2);
}

#[test]
fn a_saturated_bond_turns_aromatic_when_its_atoms_can_give_two_hydrogens() {
    let (molecule, fused) = fusion_of("benzene");
    let orders = fusion_orders(molecule, fused, BondOrder::Single, [end(2), end(2)]).unwrap();
    // Tetralin: both bond atoms join the aromatic ring.
    assert_eq!(orders.rank.0, 2);
    assert_eq!(orders.gains, [2, 2]);
    assert_eq!(
        doubles(&orders) + usize::from(orders.shared == BondOrder::Double),
        3
    );
}

#[test]
fn one_hydrogen_each_leaves_a_diene_on_a_saturated_bond() {
    let (molecule, fused) = fusion_of("benzene");
    let orders = fusion_orders(molecule, fused, BondOrder::Single, [end(1), end(1)]).unwrap();
    assert_eq!(orders.rank.0, 0);
    assert_eq!(orders.shared, BondOrder::Single);
    assert_eq!(orders.gains, [1, 1]);
    assert_eq!(doubles(&orders), 2);
    // A multiple bond elsewhere at a bond atom keeps it out of the pairing.
    let busy = FusionEnd {
        hydrogens: 2,
        multiple_elsewhere: true,
    };
    let orders = fusion_orders(molecule, fused, BondOrder::Single, [busy, end(2)]).unwrap();
    assert_eq!(orders.rank.0, 0);
    assert_eq!(orders.gains, [1, 1]);
}

#[test]
fn no_hydrogen_to_give_up_means_no_orders() {
    let (molecule, fused) = fusion_of("benzene");
    assert!(fusion_orders(molecule, fused, BondOrder::Double, [end(1), end(0)]).is_none());
    let (molecule, fused) = fusion_of("cyclohexane");
    assert!(fusion_orders(molecule, fused, BondOrder::Single, [end(0), end(2)]).is_none());
}

#[test]
fn a_saturated_ring_keeps_single_bonds() {
    let (molecule, fused) = fusion_of("cyclohexane");
    let orders = fusion_orders(molecule, fused, BondOrder::Single, [end(2), end(2)]).unwrap();
    assert_eq!(orders.shared, BondOrder::Single);
    assert_eq!(orders.gains, [1, 1]);
    assert_eq!(doubles(&orders), 0);
    // Five carbon bonds besides the fusion bond and the eight hydrogens of
    // the four new carbons; the fused atoms' hydrogens are dropped.
    assert_eq!(orders.bonds.len(), 5 + 8);
    assert_eq!(fusion_dropped(molecule, fused).len(), 6);
}
