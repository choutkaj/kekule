//! Runtime physical quantities and composable units.
//!
//! A [`Quantity`] keeps its value in the unit chosen by the caller. Conversion
//! is explicit and applies one scale factor to the complete scalar or
//! collection value. Units use linear scale conversions and integer powers of
//! the base dimensions supported by molecular structure and modelling.

use std::cmp::Ordering;
use std::fmt;
use std::ops::{Deref, DerefMut, Div, Mul, Neg};

const DIMENSION_COUNT: usize = 7;

/// Independent physical dimensions supported by [`Unit`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(usize)]
#[non_exhaustive]
pub enum BaseDimension {
    Length = 0,
    Mass = 1,
    Time = 2,
    Temperature = 3,
    Amount = 4,
    Charge = 5,
    Angle = 6,
}

/// Integer powers of the independent physical dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Dimension {
    exponents: [i32; DIMENSION_COUNT],
}

impl Dimension {
    pub const DIMENSIONLESS: Self = Self::new([0; DIMENSION_COUNT]);
    pub const LENGTH: Self = Self::base(BaseDimension::Length);
    pub const MASS: Self = Self::base(BaseDimension::Mass);
    pub const TIME: Self = Self::base(BaseDimension::Time);
    pub const TEMPERATURE: Self = Self::base(BaseDimension::Temperature);
    pub const AMOUNT: Self = Self::base(BaseDimension::Amount);
    pub const CHARGE: Self = Self::base(BaseDimension::Charge);
    pub const ANGLE: Self = Self::base(BaseDimension::Angle);

    pub const fn new(exponents: [i32; DIMENSION_COUNT]) -> Self {
        Self { exponents }
    }

    pub const fn base(dimension: BaseDimension) -> Self {
        let mut exponents = [0; DIMENSION_COUNT];
        exponents[dimension as usize] = 1;
        Self { exponents }
    }

    pub const fn exponent(self, dimension: BaseDimension) -> i32 {
        self.exponents[dimension as usize]
    }

    pub const fn is_dimensionless(self) -> bool {
        let mut index = 0;
        while index < DIMENSION_COUNT {
            if self.exponents[index] != 0 {
                return false;
            }
            index += 1;
        }
        true
    }

    const fn multiply(self, other: Self) -> Self {
        let mut result = [0; DIMENSION_COUNT];
        let mut index = 0;
        while index < DIMENSION_COUNT {
            result[index] = self.exponents[index] + other.exponents[index];
            index += 1;
        }
        Self::new(result)
    }

    const fn divide(self, other: Self) -> Self {
        let mut result = [0; DIMENSION_COUNT];
        let mut index = 0;
        while index < DIMENSION_COUNT {
            result[index] = self.exponents[index] - other.exponents[index];
            index += 1;
        }
        Self::new(result)
    }

    fn powi(self, power: i32) -> Self {
        Self::new(self.exponents.map(|exponent| exponent * power))
    }
}

impl fmt::Display for Dimension {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_dimensionless() {
            return formatter.write_str("dimensionless");
        }
        let labels = ["L", "M", "T", "Theta", "N", "Q", "A"];
        let mut first = true;
        for (label, exponent) in labels.into_iter().zip(self.exponents) {
            if exponent == 0 {
                continue;
            }
            if !first {
                formatter.write_str(" ")?;
            }
            first = false;
            if exponent == 1 {
                formatter.write_str(label)?;
            } else {
                write!(formatter, "{label}^{exponent}")?;
            }
        }
        Ok(())
    }
}

/// A linear physical unit represented by dimensions and a reference scale.
///
/// `scale` converts a magnitude in this unit to the reference SI-scale unit for
/// its dimensions. Quantity values are not normalized eagerly; the scale is
/// used only for explicit conversions.
#[derive(Debug, Clone, Copy)]
pub struct Unit {
    dimension: Dimension,
    scale: f64,
    symbol: Option<&'static str>,
}

impl Unit {
    /// Creates a custom linear unit.
    ///
    /// The scale converts a value in this unit to the SI-scale reference for
    /// `dimension`. Symbols are static so [`Unit`] remains a small `Copy` value.
    pub fn new(
        dimension: Dimension,
        scale: f64,
        symbol: Option<&'static str>,
    ) -> Result<Self, UnitError> {
        if !scale.is_finite() || scale <= 0.0 {
            return Err(UnitError::InvalidScale(scale));
        }
        Ok(Self {
            dimension,
            scale,
            symbol,
        })
    }

    const fn named(dimension: Dimension, scale: f64, symbol: &'static str) -> Self {
        Self {
            dimension,
            scale,
            symbol: Some(symbol),
        }
    }

    const fn derived(dimension: Dimension, scale: f64) -> Self {
        Self {
            dimension,
            scale,
            symbol: None,
        }
    }

    pub const fn dimension(self) -> Dimension {
        self.dimension
    }

    /// Returns the multiplier from this unit to its SI-scale reference.
    pub const fn scale(self) -> f64 {
        self.scale
    }

    pub const fn is_dimensionless(self) -> bool {
        self.dimension.is_dimensionless()
    }

    pub fn is_compatible(self, other: Self) -> bool {
        self.dimension == other.dimension
    }

    pub const fn symbol(self) -> Option<&'static str> {
        self.symbol
    }

    pub fn conversion_factor_to(self, other: Self) -> Result<f64, UnitError> {
        if !self.is_compatible(other) {
            return Err(UnitError::IncompatibleUnits {
                from: self,
                to: other,
            });
        }
        Ok(self.scale / other.scale)
    }

    pub fn powi(self, power: i32) -> Self {
        Self::derived(self.dimension.powi(power), self.scale.powi(power))
    }
}

impl PartialEq for Unit {
    fn eq(&self, other: &Self) -> bool {
        self.dimension == other.dimension && self.scale == other.scale
    }
}

impl fmt::Display for Unit {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(symbol) = self.symbol {
            formatter.write_str(symbol)
        } else if self.is_dimensionless() && self.scale == 1.0 {
            formatter.write_str("1")
        } else {
            write!(formatter, "{} [{}]", self.scale, self.dimension)
        }
    }
}

impl Mul for Unit {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Self::derived(
            self.dimension.multiply(other.dimension),
            self.scale * other.scale,
        )
    }
}

impl Div for Unit {
    type Output = Self;

    fn div(self, other: Self) -> Self::Output {
        Self::derived(
            self.dimension.divide(other.dimension),
            self.scale / other.scale,
        )
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitError {
    InvalidScale(f64),
    IncompatibleUnits { from: Unit, to: Unit },
}

impl fmt::Display for UnitError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidScale(scale) => {
                write!(
                    formatter,
                    "unit scale must be finite and positive, got {scale}"
                )
            }
            Self::IncompatibleUnits { from, to } => write!(
                formatter,
                "unit {from} with dimension {} is incompatible with unit {to} with dimension {}",
                from.dimension, to.dimension
            ),
        }
    }
}

impl std::error::Error for UnitError {}

/// A value or collection paired with one physical unit.
#[derive(Debug, Clone, Copy)]
pub struct Quantity<T> {
    value: T,
    unit: Unit,
}

impl<T> Quantity<T> {
    pub const fn new(value: T, unit: Unit) -> Self {
        Self { value, unit }
    }

    pub const fn value(&self) -> &T {
        &self.value
    }

    pub fn value_mut(&mut self) -> &mut T {
        &mut self.value
    }

    pub fn into_value(self) -> T {
        self.value
    }

    pub const fn unit(&self) -> Unit {
        self.unit
    }

    pub fn as_ref(&self) -> Quantity<&T> {
        Quantity::new(&self.value, self.unit)
    }

    pub fn map<U>(self, map: impl FnOnce(T) -> U) -> Quantity<U> {
        Quantity::new(map(self.value), self.unit)
    }
}

impl<T> Deref for Quantity<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.value
    }
}

impl<T> DerefMut for Quantity<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.value
    }
}

impl<T> Quantity<T>
where
    T: Clone + ScaleValue,
{
    pub fn to(&self, unit: Unit) -> Result<Self, UnitError> {
        let factor = self.unit.conversion_factor_to(unit)?;
        Ok(Self::new(self.value.clone().scaled(factor), unit))
    }

    pub fn value_in(&self, unit: Unit) -> Result<T, UnitError> {
        self.to(unit).map(Self::into_value)
    }
}

impl<T> Quantity<T>
where
    T: ScaleValue,
{
    pub fn into_unit(self, unit: Unit) -> Result<Self, UnitError> {
        let factor = self.unit.conversion_factor_to(unit)?;
        Ok(Self::new(self.value.scaled(factor), unit))
    }
}

impl<T: PartialEq> PartialEq for Quantity<T> {
    fn eq(&self, other: &Self) -> bool {
        self.unit == other.unit && self.value == other.value
    }
}

impl<T> Quantity<T>
where
    T: Clone + PartialEq + ScaleValue,
{
    /// Compares values after converting `other` to this quantity's unit.
    pub fn equivalent_to(&self, other: &Self) -> Result<bool, UnitError> {
        Ok(self.value == other.value_in(self.unit)?)
    }
}

impl Quantity<f64> {
    pub fn try_add(self, other: Self) -> Result<Self, UnitError> {
        let other_value = other.value_in(self.unit)?;
        Ok(Self::new(self.value + other_value, self.unit))
    }

    pub fn try_sub(self, other: Self) -> Result<Self, UnitError> {
        let other_value = other.value_in(self.unit)?;
        Ok(Self::new(self.value - other_value, self.unit))
    }

    /// Compares compatible scalar quantities with explicit numeric tolerances.
    ///
    /// `absolute_tolerance` is expressed in this quantity's unit.
    pub fn is_close(
        &self,
        other: &Self,
        relative_tolerance: f64,
        absolute_tolerance: f64,
    ) -> Result<bool, UnitError> {
        if !relative_tolerance.is_finite()
            || relative_tolerance < 0.0
            || !absolute_tolerance.is_finite()
            || absolute_tolerance < 0.0
        {
            return Ok(false);
        }
        let other_value = other.value_in(self.unit)?;
        let difference = (self.value - other_value).abs();
        Ok(difference
            <= absolute_tolerance + relative_tolerance * self.value.abs().max(other_value.abs()))
    }
}

impl PartialOrd for Quantity<f64> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if !self.unit.is_compatible(other.unit) {
            return None;
        }
        other
            .value_in(self.unit)
            .ok()
            .and_then(|other_value| self.value.partial_cmp(&other_value))
    }
}

impl Mul<f64> for Quantity<f64> {
    type Output = Self;

    fn mul(self, factor: f64) -> Self::Output {
        Self::new(self.value * factor, self.unit)
    }
}

impl Div<f64> for Quantity<f64> {
    type Output = Self;

    fn div(self, divisor: f64) -> Self::Output {
        Self::new(self.value / divisor, self.unit)
    }
}

impl Mul for Quantity<f64> {
    type Output = Self;

    fn mul(self, other: Self) -> Self::Output {
        Self::new(self.value * other.value, self.unit * other.unit)
    }
}

impl Div for Quantity<f64> {
    type Output = Self;

    fn div(self, other: Self) -> Self::Output {
        Self::new(self.value / other.value, self.unit / other.unit)
    }
}

impl Neg for Quantity<f64> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.value, self.unit)
    }
}

impl Mul<Unit> for f64 {
    type Output = Quantity<Self>;

    fn mul(self, unit: Unit) -> Self::Output {
        Quantity::new(self, unit)
    }
}

/// Scaling operation used to convert complete [`Quantity`] values.
pub trait ScaleValue: Sized {
    fn scaled(self, factor: f64) -> Self;
}

impl ScaleValue for f64 {
    fn scaled(self, factor: f64) -> Self {
        self * factor
    }
}

impl<T> ScaleValue for Option<T>
where
    T: ScaleValue,
{
    fn scaled(self, factor: f64) -> Self {
        self.map(|value| value.scaled(factor))
    }
}

impl<T, const N: usize> ScaleValue for [T; N]
where
    T: ScaleValue,
{
    fn scaled(self, factor: f64) -> Self {
        self.map(|value| value.scaled(factor))
    }
}

impl<T> ScaleValue for Vec<T>
where
    T: ScaleValue,
{
    fn scaled(self, factor: f64) -> Self {
        self.into_iter().map(|value| value.scaled(factor)).collect()
    }
}

pub const DIMENSIONLESS: Unit = Unit::named(Dimension::DIMENSIONLESS, 1.0, "1");

pub const METER: Unit = Unit::named(Dimension::LENGTH, 1.0, "m");
pub const NANOMETER: Unit = Unit::named(Dimension::LENGTH, 1.0e-9, "nm");
pub const ANGSTROM: Unit = Unit::named(Dimension::LENGTH, 1.0e-10, "A");
pub const BOHR: Unit = Unit::named(Dimension::LENGTH, 5.291_772_109_03e-11, "bohr");
pub const SQUARE_ANGSTROM: Unit =
    Unit::named(Dimension::new([2, 0, 0, 0, 0, 0, 0]), 1.0e-20, "A^2");

pub const KILOGRAM: Unit = Unit::named(Dimension::MASS, 1.0, "kg");
pub const DALTON: Unit = Unit::named(Dimension::MASS, 1.660_539_068_92e-27, "Da");

pub const SECOND: Unit = Unit::named(Dimension::TIME, 1.0, "s");
pub const PICOSECOND: Unit = Unit::named(Dimension::TIME, 1.0e-12, "ps");
pub const FEMTOSECOND: Unit = Unit::named(Dimension::TIME, 1.0e-15, "fs");

pub const KELVIN: Unit = Unit::named(Dimension::TEMPERATURE, 1.0, "K");
pub const MOLE: Unit = Unit::named(Dimension::AMOUNT, 1.0, "mol");

pub const COULOMB: Unit = Unit::named(Dimension::CHARGE, 1.0, "C");
pub const ELEMENTARY_CHARGE: Unit = Unit::named(Dimension::CHARGE, 1.602_176_634e-19, "e");

pub const RADIAN: Unit = Unit::named(Dimension::ANGLE, 1.0, "rad");
pub const DEGREE: Unit = Unit::named(Dimension::ANGLE, std::f64::consts::PI / 180.0, "deg");

pub const JOULE: Unit = Unit::named(Dimension::new([2, 1, -2, 0, 0, 0, 0]), 1.0, "J");
pub const KILOJOULE: Unit = Unit::named(JOULE.dimension, 1.0e3, "kJ");
pub const KILOCALORIE: Unit = Unit::named(JOULE.dimension, 4.184e3, "kcal");
pub const KILOJOULE_PER_MOLE: Unit =
    Unit::named(Dimension::new([2, 1, -2, 0, -1, 0, 0]), 1.0e3, "kJ/mol");
pub const KILOCALORIE_PER_MOLE: Unit =
    Unit::named(KILOJOULE_PER_MOLE.dimension, 4.184e3, "kcal/mol");
pub const KILOJOULE_PER_MOLE_PER_ANGSTROM: Unit =
    Unit::named(Dimension::new([1, 1, -2, 0, -1, 0, 0]), 1.0e13, "kJ/mol/A");
pub const KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM: Unit = Unit::named(
    Dimension::new([0, 1, -2, 0, -1, 0, 0]),
    1.0e23,
    "kJ/mol/A^2",
);

/// Preferred explicit units used by the fixed-topology modelling kernel.
pub const MODEL_LENGTH_UNIT: Unit = ANGSTROM;
pub const MODEL_ENERGY_UNIT: Unit = KILOJOULE_PER_MOLE;
pub const MODEL_GRADIENT_UNIT: Unit = KILOJOULE_PER_MOLE_PER_ANGSTROM;
pub const MODEL_FORCE_CONSTANT_UNIT: Unit = KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converts_scalar_and_collection_quantities() {
        let length = Quantity::new(1.0, NANOMETER);
        assert_eq!(length.value_in(ANGSTROM).unwrap(), 10.0);

        let positions = Quantity::new(vec![[0.1, 0.2, 0.3]], NANOMETER);
        assert_eq!(positions.value_in(ANGSTROM).unwrap(), vec![[1.0, 2.0, 3.0]]);
    }

    #[test]
    fn composes_and_checks_dimensions() {
        let gradient = KILOJOULE_PER_MOLE / ANGSTROM;
        assert!(gradient.is_compatible(MODEL_GRADIENT_UNIT));
        assert_eq!(
            gradient.conversion_factor_to(MODEL_GRADIENT_UNIT).unwrap(),
            1.0
        );
        assert!(ANGSTROM.conversion_factor_to(KELVIN).is_err());
    }

    #[test]
    fn scalar_arithmetic_preserves_or_composes_units() {
        let sum = Quantity::new(1.0, NANOMETER)
            .try_add(Quantity::new(5.0, ANGSTROM))
            .unwrap();
        assert_eq!(sum, Quantity::new(1.5, NANOMETER));

        let energy = Quantity::new(2.0, MODEL_FORCE_CONSTANT_UNIT)
            * Quantity::new(3.0, ANGSTROM)
            * Quantity::new(3.0, ANGSTROM);
        assert!(energy.unit().is_compatible(MODEL_ENERGY_UNIT));
        assert_eq!(energy.value_in(MODEL_ENERGY_UNIT).unwrap(), 18.0);
    }

    #[test]
    fn supports_checked_custom_units_and_explicit_equivalence() {
        let picometer = Unit::new(Dimension::LENGTH, 1.0e-12, Some("pm")).unwrap();
        assert_eq!(Quantity::new(100.0, picometer).value_in(ANGSTROM), Ok(1.0));
        assert!(Unit::new(Dimension::LENGTH, 0.0, None).is_err());
        assert!(Quantity::new(1.0, NANOMETER)
            .is_close(&Quantity::new(10.0, ANGSTROM), 1.0e-12, 0.0)
            .unwrap());
    }
}
