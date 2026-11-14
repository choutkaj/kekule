use super::*;
use crate::units::ANGSTROM;

#[test]
fn carbon_reference_values_match_archived_sources() {
    let carbon = Element::from_atomic_number(6).unwrap();

    assert_eq!(
        carbon.covalent_radius(),
        Some(Quantity::new(0.76, ANGSTROM))
    );
    assert_eq!(
        carbon.dipole_polarizability(),
        Some(Quantity::new(11.3, CUBIC_BOHR))
    );
    assert!(
        (carbon
            .van_der_waals_radius(VanDerWaalsRadiusSource::FreeAtom)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            - 1.91)
            .abs()
            < 0.02
    );
}

#[test]
fn requested_reference_radius_falls_back_to_free_atom_value() {
    let neodymium = Element::from_atomic_number(60).unwrap();
    assert!(REFERENCE_VDW_RADIUS_BOHR[60].is_none());
    assert!(
        (neodymium
            .van_der_waals_radius(VanDerWaalsRadiusSource::Reference)
            .unwrap()
            .value_in(ANGSTROM)
            .unwrap()
            - 2.90)
            .abs()
            < 0.03
    );
}

fn angstrom(element: u8, source: VanDerWaalsRadiusSource) -> Option<f64> {
    Element::from_atomic_number(element)
        .unwrap()
        .van_der_waals_radius(source)
        .map(|radius| radius.value_in(ANGSTROM).unwrap())
}

#[test]
fn bondi_column_follows_bondi_1964_with_mantina_completions() {
    use VanDerWaalsRadiusSource::Bondi;
    // Table I hydrogen, carbon, and bromine; Table XIV lithium, mercury,
    // cadmium, and platinum's X-ray value; a Mantina completion.
    for (element, radius) in [
        (1, 1.20),
        (6, 1.70),
        (35, 1.85),
        (3, 1.82),
        (80, 1.55),
        (48, 1.58),
        (78, 1.72),
        (38, 2.49),
    ] {
        assert_eq!(angstrom(element, Bondi), Some(radius), "element {element}");
    }
    let covered = (1..=118).filter(|&z| angstrom(z, Bondi).is_some()).count();
    assert_eq!(covered, 54);
    // No fallback to another column: iron and the lanthanides have no
    // Bondi radius.
    assert_eq!(angstrom(26, Bondi), None);
    assert_eq!(angstrom(64, Bondi), None);
}

#[test]
fn bondi_mantina_column_is_table_twelve() {
    use VanDerWaalsRadiusSource::{Bondi, BondiMantina};
    let covered = (1..=118)
        .filter(|&z| angstrom(z, BondiMantina).is_some())
        .collect::<Vec<_>>();
    assert_eq!(covered.len(), 44);
    let differing = covered
        .iter()
        .copied()
        .filter(|&z| angstrom(z, BondiMantina) != angstrom(z, Bondi))
        .collect::<Vec<_>>();
    // Hydrogen, lithium, and bromine.
    assert_eq!(differing, [1, 3, 35]);
    assert_eq!(angstrom(1, BondiMantina), Some(1.10));
    assert_eq!(angstrom(29, BondiMantina), None);
}

#[test]
fn alvarez_column_follows_table_one() {
    use VanDerWaalsRadiusSource::Alvarez;
    for (element, radius) in [(1, 1.20), (6, 1.77), (26, 2.44), (54, 2.06), (99, 2.7)] {
        assert_eq!(
            angstrom(element, Alvarez),
            Some(radius),
            "element {element}"
        );
    }
    let missing = (1..=118)
        .filter(|&z| angstrom(z, Alvarez).is_none())
        .collect::<Vec<_>>();
    let expected = [61, 84, 85, 86, 87, 88]
        .into_iter()
        .chain(100..=118)
        .collect::<Vec<_>>();
    assert_eq!(missing, expected);
}

#[test]
fn charry_tkatchenko_conversion_is_unchanged_for_solvation() {
    let oxygen = Element::from_atomic_number(8).unwrap();
    assert_eq!(
        oxygen.van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Reference),
        Some(3.2314 * BOHR_TO_ANGSTROM)
    );
    assert_eq!(
        oxygen.van_der_waals_radius_angstrom(VanDerWaalsRadiusSource::Bondi),
        Some(1.52)
    );
    assert_eq!(
        VanDerWaalsRadiusSource::default(),
        VanDerWaalsRadiusSource::FreeAtom
    );
}

#[test]
fn missing_both_vdw_sources_is_explicit() {
    let livermorium = Element::from_atomic_number(116).unwrap();
    assert_eq!(
        livermorium.van_der_waals_radius(VanDerWaalsRadiusSource::FreeAtom),
        None
    );
}
