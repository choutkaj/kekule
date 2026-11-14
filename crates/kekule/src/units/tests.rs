use super::*;

#[test]
fn converts_scalar_and_collection_quantities() {
    let noncanonical_length = Quantity::new(10.0, ANGSTROM);
    assert!((noncanonical_length.value_in(NANOMETER).unwrap() - 1.0).abs() < 1.0e-12);

    let energy = Quantity::new(1.0, KILOCALORIE_PER_MOLE);
    assert!((energy.value_in(KILOJOULE_PER_MOLE).unwrap() - 4.184).abs() < 1.0e-12);

    let angle = Quantity::new(180.0, DEGREE);
    assert!((angle.value_in(RADIAN).unwrap() - std::f64::consts::PI).abs() < 1.0e-12);

    let length = Quantity::new(1.0, NANOMETER);
    assert_eq!(length.value_in(ANGSTROM).unwrap(), 10.0);

    let positions = Quantity::new(vec![[0.1, 0.2, 0.3]], NANOMETER);
    assert_eq!(positions.value_in(ANGSTROM).unwrap(), vec![[1.0, 2.0, 3.0]]);
}

#[test]
fn composes_and_checks_dimensions() {
    let gradient = KILOJOULE_PER_MOLE.try_div(ANGSTROM).unwrap();
    assert!(gradient.is_compatible(CANONICAL_GRADIENT_UNIT));
    assert_eq!(
        gradient
            .conversion_factor_to(CANONICAL_GRADIENT_UNIT)
            .unwrap(),
        10.0
    );
    assert!(ANGSTROM.conversion_factor_to(KELVIN).is_err());
}

#[test]
fn scalar_arithmetic_preserves_or_composes_units() {
    let sum = Quantity::new(1.0, NANOMETER)
        .try_add(Quantity::new(5.0, ANGSTROM))
        .unwrap();
    assert_eq!(sum, Quantity::new(1.5, NANOMETER));

    let energy = Quantity::new(2.0, CANONICAL_FORCE_CONSTANT_UNIT)
        .try_mul(Quantity::new(0.3, NANOMETER))
        .and_then(|partial| partial.try_mul(Quantity::new(0.3, NANOMETER)))
        .unwrap();
    assert!(energy.unit().is_compatible(CANONICAL_ENERGY_UNIT));
    assert!((energy.value_in(CANONICAL_ENERGY_UNIT).unwrap() - 0.18).abs() < 1.0e-15);
}

#[test]
fn supports_checked_custom_units_and_explicit_equivalence() {
    let picometer = Unit::new(Dimension::LENGTH, 1.0e-12, Some("pm")).unwrap();
    assert_eq!(Quantity::new(100.0, picometer).value_in(ANGSTROM), Ok(1.0));
    assert!(Unit::new(Dimension::LENGTH, 0.0, None).is_err());
    assert!(Quantity::new(1.0, NANOMETER)
        .is_close(
            &Quantity::new(10.0, ANGSTROM),
            1.0e-12,
            Quantity::new(0.0, NANOMETER)
        )
        .unwrap());
}

#[test]
fn ordering_and_equality_compare_the_same_representation() {
    let angstroms = Quantity::new(10.0, ANGSTROM);
    let nanometers = Quantity::new(1.0, NANOMETER);
    assert_ne!(angstroms, nanometers);
    assert_eq!(angstroms.partial_cmp(&nanometers), None);
    assert_eq!(nanometers.partial_cmp(&angstroms), None);
    assert!(angstroms.equivalent_to(&nanometers).unwrap());
    assert_eq!(angstroms.partial_cmp(&angstroms), Some(Ordering::Equal));
    assert!(angstroms < Quantity::new(11.0, ANGSTROM));
    assert!(angstroms > Quantity::new(9.0, ANGSTROM));
    assert_eq!(angstroms.partial_cmp(&Quantity::new(10.0, SECOND)), None);
    assert_eq!(
        angstroms.partial_cmp(&Quantity::new(f64::NAN, ANGSTROM)),
        None
    );
}

#[test]
fn canonical_units_match_the_library_wide_basis() {
    assert_eq!(CANONICAL_LENGTH_UNIT, NANOMETER);
    assert_eq!(CANONICAL_MOLAR_MASS_UNIT, DALTON);
    assert_eq!(CANONICAL_AREA_UNIT, SQUARE_NANOMETER);
    assert_eq!(CANONICAL_TIME_UNIT, PICOSECOND);
    assert_eq!(CANONICAL_ENERGY_UNIT, KILOJOULE_PER_MOLE);
    assert_eq!(CANONICAL_CHARGE_UNIT, ELEMENTARY_CHARGE);
    assert_eq!(CANONICAL_TEMPERATURE_UNIT, KELVIN);
    assert_eq!(CANONICAL_ANGLE_UNIT, RADIAN);
}

#[test]
fn dalton_is_molar_mass_and_composes_coherently() {
    assert_eq!(DALTON.dimension(), Dimension::MOLAR_MASS);
    assert_ne!(DALTON.dimension(), KILOGRAM.dimension());

    let per_square_picosecond = PICOSECOND.try_powi(-2).unwrap();
    let energy = DALTON
        .try_mul(NANOMETER.try_powi(2).unwrap())
        .and_then(|unit| unit.try_mul(per_square_picosecond))
        .unwrap();
    let force = DALTON
        .try_mul(NANOMETER)
        .and_then(|unit| unit.try_mul(per_square_picosecond))
        .unwrap();
    let force_constant = DALTON.try_mul(per_square_picosecond).unwrap();

    assert!(energy.is_compatible(KILOJOULE_PER_MOLE));
    assert!(force.is_compatible(CANONICAL_FORCE_UNIT));
    assert!(force_constant.is_compatible(CANONICAL_FORCE_CONSTANT_UNIT));
    for factor in [
        energy.conversion_factor_to(KILOJOULE_PER_MOLE).unwrap(),
        force.conversion_factor_to(CANONICAL_FORCE_UNIT).unwrap(),
        force_constant
            .conversion_factor_to(CANONICAL_FORCE_CONSTANT_UNIT)
            .unwrap(),
    ] {
        assert!((factor - 1.0).abs() <= 4.0 * f64::EPSILON);
    }
}
