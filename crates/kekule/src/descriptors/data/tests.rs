use super::*;

#[test]
fn generated_tables_cover_pinned_reference_rows() {
    assert_eq!(ISOTOPE_MASSES.len(), 3557);
    assert_eq!(
        exact_isotope_mass(Element::from_symbol("C").unwrap(), 13),
        Some(13.003_354_835_34)
    );
    assert_eq!(
        most_abundant_isotope(Element::from_symbol("U").unwrap()),
        Some(238)
    );
    assert_eq!(
        standard_atomic_weight(Element::from_symbol("Tc").unwrap()),
        None
    );
}
