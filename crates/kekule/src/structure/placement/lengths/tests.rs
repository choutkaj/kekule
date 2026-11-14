use super::*;

fn element(symbol: &str) -> Element {
    Element::from_symbol(symbol).unwrap()
}

#[test]
fn table_is_ordered_and_symmetric_lookup_works() {
    for (low, high, _, length) in TABLE {
        assert!(low <= high && *length > 0.5 && *length < 2.5);
    }
    assert_eq!(
        ideal_length_angstrom(element("O"), element("C"), BondOrder::Double),
        Some(1.22)
    );
    assert_eq!(
        ideal_length_angstrom(element("H"), element("C"), BondOrder::Dative),
        Some(1.09)
    );
}

#[test]
fn fallback_shortens_multiple_bonds() {
    let single = ideal_length_angstrom(element("Si"), element("Si"), BondOrder::Single).unwrap();
    let double = ideal_length_angstrom(element("Si"), element("Si"), BondOrder::Double).unwrap();
    let triple = ideal_length_angstrom(element("Si"), element("Si"), BondOrder::Triple).unwrap();
    assert!((single - 2.22).abs() < 1e-9);
    assert!(single > double && double > triple);
}
