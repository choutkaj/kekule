use super::*;

mod reference;

#[test]
fn unit_conversion_and_dimensional_rejection() {
    let bond = parse_quantity(
        "2 * kilocalorie_per_mole * angstrom ** -2",
        CANONICAL_FORCE_CONSTANT_UNIT,
    )
    .unwrap();
    assert!((*bond.value() - 836.8).abs() < 1e-10);
    assert!(
        (*parse_quantity("180 * degree", RADIAN).unwrap().value() - std::f64::consts::PI).abs()
            < 1e-14
    );
    assert!(parse_quantity("1 * angstrom", KILOJOULE_PER_MOLE).is_err());
    assert!(parse_quantity("NaN * angstrom", NANOMETER).is_err());
}
