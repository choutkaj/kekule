use super::MOLAR_GAS_CONSTANT;
use crate::units::constants;
use crate::units::{Unit, KILOJOULE};

#[test]
fn the_gas_constant_is_the_si_constant_in_kilojoules() {
    let kilojoule_per_mole_per_kelvin = Unit::new(
        constants::MOLAR_GAS_CONSTANT.unit().dimension(),
        KILOJOULE.scale(),
        None,
    )
    .unwrap();
    let si = constants::MOLAR_GAS_CONSTANT
        .value_in(kilojoule_per_mole_per_kelvin)
        .unwrap();
    assert!(
        (si - MOLAR_GAS_CONSTANT).abs() <= f64::EPSILON * MOLAR_GAS_CONSTANT,
        "{si} != {MOLAR_GAS_CONSTANT}"
    );
}
