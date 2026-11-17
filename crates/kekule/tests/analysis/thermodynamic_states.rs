//! Thermodynamic states: validated temperatures and pressures, in canonical
//! units, and their thermal energy.

use kekule::structure::{StructureSetError, ThermodynamicState};
use kekule::units::{Quantity, Unit, BAR, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

fn kelvin(value: f64) -> Quantity<f64> {
    Quantity::new(value, KELVIN)
}

fn bar(value: f64) -> Quantity<f64> {
    Quantity::new(value, BAR)
}

#[test]
fn a_state_has_a_pressure_exactly_when_it_is_isothermal_isobaric() {
    let canonical = ThermodynamicState::canonical(kelvin(300.0)).unwrap();
    assert_eq!(canonical.temperature(), kelvin(300.0));
    assert_eq!(canonical.pressure(), None);

    let pascal = Unit::new(BAR.dimension(), 1.0, None).unwrap();
    let isobaric =
        ThermodynamicState::isothermal_isobaric(kelvin(310.0), Quantity::new(101_325.0, pascal))
            .unwrap();
    assert_eq!(isobaric.temperature(), kelvin(310.0));
    let pressure = isobaric.pressure().unwrap();
    assert_eq!(pressure.unit(), BAR);
    assert!((*pressure.value() - 1.013_25).abs() <= 2.0 * f64::EPSILON);

    let warmer = isobaric.with_temperature(kelvin(350.0)).unwrap();
    assert_eq!(warmer.temperature(), kelvin(350.0));
    assert_eq!(warmer.pressure(), isobaric.pressure());
    assert_eq!(
        canonical.with_temperature(kelvin(350.0)).unwrap(),
        ThermodynamicState::canonical(kelvin(350.0)).unwrap()
    );
}

#[test]
fn the_thermal_energy_is_rt_in_kilojoules_per_mole() {
    let state = ThermodynamicState::canonical(kelvin(298.15)).unwrap();
    let thermal_energy = state.thermal_energy();
    assert_eq!(thermal_energy.unit(), KILOJOULE_PER_MOLE);
    let expected = 8.314_462_618_153_24e-3 * 298.15;
    assert!((*thermal_energy.value() - expected).abs() <= f64::EPSILON * expected);
}

#[test]
fn temperatures_and_pressures_are_finite_and_positive() {
    // The thermal energy of the smallest positive temperature underflows.
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert_eq!(
            ThermodynamicState::canonical(kelvin(temperature)),
            Err(StructureSetError::InvalidTemperature),
            "{temperature}"
        );
        assert_eq!(
            ThermodynamicState::isothermal_isobaric(kelvin(temperature), bar(1.0)),
            Err(StructureSetError::InvalidTemperature),
            "{temperature}"
        );
    }
    for pressure in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            ThermodynamicState::isothermal_isobaric(kelvin(300.0), bar(pressure)),
            Err(StructureSetError::InvalidPressure),
            "{pressure}"
        );
    }
    let state = ThermodynamicState::canonical(kelvin(300.0)).unwrap();
    assert_eq!(
        state.with_temperature(kelvin(-1.0)),
        Err(StructureSetError::InvalidTemperature)
    );
    let length = Quantity::new(1.0, NANOMETER);
    assert!(matches!(
        ThermodynamicState::canonical(length),
        Err(StructureSetError::Unit(_))
    ));
    assert!(matches!(
        ThermodynamicState::isothermal_isobaric(kelvin(300.0), length),
        Err(StructureSetError::Unit(_))
    ));
}
