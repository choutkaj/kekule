use super::*;

fn domains(symbol: &str, charge: i8, substituents: usize, valence: usize) -> (usize, usize) {
    perceive_domains(DomainInput {
        element: Element::from_symbol(symbol).unwrap(),
        formal_charge: charge,
        radical_electrons: 0,
        substituents,
        bond_valence: valence,
        coordinate_bonds: 0,
        conjugated_nitrogen: false,
    })
}

#[test]
fn domain_counts_follow_valence_shell_electron_pairs() {
    assert_eq!(domains("C", 0, 4, 4), (4, 0)); // methane
    assert_eq!(domains("C", 0, 3, 4), (3, 0)); // alkene carbon
    assert_eq!(domains("C", 0, 2, 4), (2, 0)); // alkyne carbon
    assert_eq!(domains("N", 0, 1, 3), (2, 1)); // nitrile nitrogen
    assert_eq!(domains("O", 0, 2, 2), (4, 2)); // water
    assert_eq!(domains("O", 0, 1, 2), (3, 2)); // carbonyl oxygen
    assert_eq!(domains("N", 1, 4, 4), (4, 0)); // ammonium
    assert_eq!(domains("N", 1, 3, 4), (3, 0)); // nitro nitrogen
    assert_eq!(domains("C", 1, 3, 3), (3, 0)); // carbocation
    assert_eq!(domains("C", -1, 3, 3), (4, 1)); // carbanion
    assert_eq!(domains("S", 0, 4, 6), (4, 0)); // sulfonyl
    assert_eq!(domains("S", 0, 3, 4), (4, 1)); // sulfoxide
    assert_eq!(domains("P", 0, 5, 5), (5, 0)); // phosphorus pentachloride
    assert_eq!(domains("S", 0, 6, 6), (6, 0)); // sulfur hexafluoride
    assert_eq!(domains("S", 0, 4, 4), (5, 1)); // sulfur tetrafluoride
    assert_eq!(domains("Xe", 0, 4, 4), (6, 2)); // xenon tetrafluoride
    assert_eq!(domains("Xe", 0, 2, 2), (5, 3)); // xenon difluoride
    assert_eq!(domains("Fe", 0, 6, 6), (6, 0)); // no lone pairs off the main group
}

#[test]
fn radical_pairs_take_a_domain_and_single_radicals_none() {
    let radical = |symbol: &str, electrons: usize, substituents: usize, valence: usize| {
        perceive_domains(DomainInput {
            element: Element::from_symbol(symbol).unwrap(),
            formal_charge: 0,
            radical_electrons: electrons,
            substituents,
            bond_valence: valence,
            coordinate_bonds: 0,
            conjugated_nitrogen: false,
        })
    };
    assert_eq!(radical("C", 2, 2, 2), (3, 1)); // methylene, bent
    assert_eq!(radical("N", 2, 1, 1), (3, 2)); // imidogen
    assert_eq!(radical("C", 1, 3, 3), (3, 0)); // methyl radical, planar
    assert_eq!(radical("N", 1, 2, 2), (3, 1)); // amidogen, bent
    assert_eq!(radical("O", 1, 1, 1), (3, 2)); // hydroxyl
}
