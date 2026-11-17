//! Thermodynamic states: the temperature that a weight refers to and, at
//! constant pressure, the pressure.

use crate::structure::StructureSetError;
use crate::units::{
    Quantity, CANONICAL_ENERGY_UNIT, CANONICAL_PRESSURE_UNIT, CANONICAL_TEMPERATURE_UNIT,
};

/// A temperature and, in the isothermal-isobaric ensemble, a pressure.
///
/// A state is what a weight means something against: the state a run
/// sampled, or the state an analysis reweights to. A state with a pressure
/// is isothermal-isobaric (NpT) and one without is canonical (NVT). The
/// pressure is recorded, not reweighted.
///
/// The temperature is finite and positive, and so is its thermal energy `RT`;
/// the pressure is finite and positive. They are stored in
/// [`CANONICAL_TEMPERATURE_UNIT`] and [`CANONICAL_PRESSURE_UNIT`].
///
/// ```
/// use kekule::structure::ThermodynamicState;
/// use kekule::units::{Quantity, BAR, KELVIN, KILOJOULE_PER_MOLE};
///
/// let state = ThermodynamicState::isothermal_isobaric(
///     Quantity::new(300.0, KELVIN),
///     Quantity::new(1.0, BAR),
/// )?;
/// assert!(state.pressure().is_some());
/// let thermal_energy = state.thermal_energy().value_in(KILOJOULE_PER_MOLE)?;
/// assert!((thermal_energy - 2.494_338_785).abs() < 1e-9);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermodynamicState {
    temperature: f64,
    pressure: Option<f64>,
}

impl ThermodynamicState {
    /// The canonical state at `temperature`, without a pressure.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::Unit`] for a temperature in a unit that is not a
    /// temperature, and [`StructureSetError::InvalidTemperature`] for one that
    /// is not finite and positive or whose thermal energy is not.
    pub fn canonical(temperature: Quantity<f64>) -> Result<Self, StructureSetError> {
        Ok(Self {
            temperature: state_temperature(temperature)?,
            pressure: None,
        })
    }

    /// The isothermal-isobaric state at `temperature` and `pressure`.
    ///
    /// # Errors
    ///
    /// The errors of [`ThermodynamicState::canonical`],
    /// [`StructureSetError::Unit`] for a pressure in a unit that is not a
    /// pressure, and [`StructureSetError::InvalidPressure`] for one that is not
    /// finite and positive.
    pub fn isothermal_isobaric(
        temperature: Quantity<f64>,
        pressure: Quantity<f64>,
    ) -> Result<Self, StructureSetError> {
        let temperature = state_temperature(temperature)?;
        let pressure = pressure.value_in(CANONICAL_PRESSURE_UNIT)?;
        if !pressure.is_finite() || pressure <= 0.0 {
            return Err(StructureSetError::InvalidPressure);
        }
        Ok(Self {
            temperature,
            pressure: Some(pressure),
        })
    }

    /// This state at another temperature, with the same pressure or none.
    ///
    /// # Errors
    ///
    /// The errors of [`ThermodynamicState::canonical`].
    pub fn with_temperature(self, temperature: Quantity<f64>) -> Result<Self, StructureSetError> {
        Ok(Self {
            temperature: state_temperature(temperature)?,
            ..self
        })
    }

    /// The temperature, in [`CANONICAL_TEMPERATURE_UNIT`].
    pub fn temperature(&self) -> Quantity<f64> {
        Quantity::new(self.temperature, CANONICAL_TEMPERATURE_UNIT)
    }

    /// The pressure, in [`CANONICAL_PRESSURE_UNIT`], of an isothermal-isobaric
    /// state; `None` for a canonical one.
    pub fn pressure(&self) -> Option<Quantity<f64>> {
        self.pressure
            .map(|pressure| Quantity::new(pressure, CANONICAL_PRESSURE_UNIT))
    }

    /// The thermal energy `RT`, in [`CANONICAL_ENERGY_UNIT`].
    pub fn thermal_energy(&self) -> Quantity<f64> {
        let thermal_energy = thermal_energy(self.temperature)
            .expect("a thermodynamic state's thermal energy is finite and positive");
        Quantity::new(thermal_energy, CANONICAL_ENERGY_UNIT)
    }
}

/// A temperature in [`CANONICAL_TEMPERATURE_UNIT`] that is finite and
/// positive, with a finite and positive thermal energy.
fn state_temperature(temperature: Quantity<f64>) -> Result<f64, StructureSetError> {
    let kelvin = checked_temperature(temperature)?;
    thermal_energy(kelvin)
        .map(|_| kelvin)
        .ok_or(StructureSetError::InvalidTemperature)
}

/// A finite, positive temperature in [`CANONICAL_TEMPERATURE_UNIT`].
fn checked_temperature(temperature: Quantity<f64>) -> Result<f64, StructureSetError> {
    let temperature = temperature.value_in(CANONICAL_TEMPERATURE_UNIT)?;
    if temperature.is_finite() && temperature > 0.0 {
        Ok(temperature)
    } else {
        Err(StructureSetError::InvalidTemperature)
    }
}

/// The molar gas constant `N_A k_B` in kJ/(mol K), exact since the 2019 SI:
/// [`constants::MOLAR_GAS_CONSTANT`](crate::units::constants::MOLAR_GAS_CONSTANT)
/// in kilojoules.
const MOLAR_GAS_CONSTANT: f64 = 8.314_462_618_153_24e-3;

/// `R T` in [`CANONICAL_ENERGY_UNIT`] for a temperature in kelvin, when it is
/// finite and positive.
fn thermal_energy(kelvin: f64) -> Option<f64> {
    let thermal_energy = MOLAR_GAS_CONSTANT * kelvin;
    (thermal_energy.is_finite() && thermal_energy > 0.0).then_some(thermal_energy)
}

#[cfg(test)]
mod tests;
