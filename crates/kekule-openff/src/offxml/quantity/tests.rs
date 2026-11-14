use super::*;

#[test]
fn equivalent_unit_spellings_and_parentheses() {
    for text in [
        "2*kilocalorie_per_mole*angstrom**-2",
        "2 * kilocalorie / (mole * angstrom**2)",
        "2e0*kcal/mol/angstroms**2",
        "2.*(kilocalories/moles)*angstrom**(-2)",
        "+2 * kilocalorie_per_mole / angstrom ** 2",
    ] {
        let value = parse_quantity(text, CANONICAL_FORCE_CONSTANT_UNIT).unwrap();
        assert!((*value.value() - 836.8).abs() < 1e-10, "{text}");
    }
    assert_eq!(
        *parse_quantity("-2**2 * elementary_charge", ELEMENTARY_CHARGE)
            .unwrap()
            .value(),
        -4.0
    );
    assert_eq!(
        *parse_quantity("(-2)**2 * elementary_charge", ELEMENTARY_CHARGE)
            .unwrap()
            .value(),
        4.0
    );
}

#[test]
fn invalid_or_unbounded_expressions_are_rejected() {
    for text in [
        "",
        "1 nm",
        "1*nm+2*nm",
        "1*nm/0",
        "NaN*nm",
        "1e999*nm",
        "1*nm**.5",
        "1*nm**999999999999",
        "1*(nm",
        "1*nm)",
        "1*nm**2**3",
        "1*kelvin",
        "1*__import__('os')",
        "1*nm # comment",
        "1*nm***2",
        "1*nm**33",
        "1*nm/1e-999",
        "1*nm*",
        "1*Å",
    ] {
        assert!(parse_quantity(text, NANOMETER).is_err(), "{text}");
    }
    assert!(parse_quantity(
        &format!("{}nm{}", "(".repeat(34), ")".repeat(34)),
        NANOMETER
    )
    .is_err());
    assert!(parse_quantity(&" ".repeat(4097), NANOMETER).is_err());
    assert!(parse_quantity("1*kJ/mol", NANOMETER).is_err());
}
