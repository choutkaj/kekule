# Physical Quantities and Units

## Summary

Represent scalar and collection-valued physical quantities as an explicit value
paired with a composable runtime unit.

## Behavior/API

- Exposes `units::Dimension`, `BaseDimension`, `Unit`, `Quantity<T>`,
  `ScaleValue`, and `UnitError`.
- Supports checked custom linear units plus predefined constants; symbols are
  static so `Unit` stays a small copyable value.
- Provides predefined molecular and SI-scale units for length, area, mass,
  time, temperature, amount, charge, angle, energy, molar energy, gradients,
  and force constants.
- Converts only between compatible dimensions and applies one scale factor to
  a complete scalar or collection value.
- Provides explicit exact-equivalence and tolerance-based scalar comparison;
  ordinary `PartialEq` compares the represented value and unit without hidden
  conversion.
- Composes units through multiplication, division, and integer powers.
- Keeps values in their declared unit until conversion is explicitly requested.

## Implementation Notes

- Dimensions use integer powers of seven independent base dimensions: length,
  mass, time, temperature, amount, charge, and angle.
- Unit scale factors refer to SI-scale base dimensions, but quantity values are
  not eagerly normalized to SI.
- `ScaleValue` supports recursive conversion of arrays, options, and vectors;
  geometry types implement it at their owning module boundary.
- The modelling layer declares its preferred length, molar-energy, gradient,
  and force-constant units as explicit public constants.

## Validation

- Focused unit tests cover scalar and collection conversion, incompatible
  dimensions, composite units, and scalar quantity arithmetic.
- This feature intentionally has no validation-harness requirement and remains
  no external parity result is recorded.

## Out Of Scope

- Affine or logarithmic units, arbitrary unit-string parsing, runtime unit or
  dimension registries, fractional powers, uncertainty propagation, automatic
  provenance propagation, and implicit conversion inside numerical kernels.

## Revision Notes

- v1: Add the runtime `Quantity<T>` and composable `Unit` foundation.
