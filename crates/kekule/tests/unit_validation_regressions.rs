use kekule::geometry::Point3;
use kekule::structure::{PositionError, Positions};
use kekule::units::{
    Dimension, Quantity, Unit, UnitError, ANGSTROM, DIMENSIONLESS, METER, NANOMETER,
};

fn unit(dimension: Dimension, scale: f64) -> Unit {
    Unit::new(dimension, scale, None).unwrap()
}

#[test]
fn dynamic_unit_composition_rejects_unrepresentable_scales() {
    let tiny_length = unit(Dimension::LENGTH, 1e-200);
    let tiny_factor = unit(Dimension::DIMENSIONLESS, 1e-200);
    let huge_factor = unit(Dimension::DIMENSIONLESS, 1e200);
    let huge_length = unit(Dimension::LENGTH, 1e200);
    for result in [
        tiny_length.try_mul(tiny_factor),
        tiny_length.try_div(huge_factor),
        tiny_length.try_powi(2),
        huge_length.try_mul(huge_factor),
        huge_length.try_div(tiny_factor),
        tiny_length.try_powi(-2),
    ] {
        assert!(matches!(result, Err(UnitError::InvalidScale(_))));
    }
    assert_eq!(tiny_length.try_powi(0), Ok(DIMENSIONLESS));
    assert_eq!(tiny_length.try_mul(huge_factor), Ok(METER));
    assert_eq!(tiny_length.try_div(tiny_factor), Ok(METER));

    let tiny_quantity = Quantity::new(1.0, tiny_length);
    assert!(matches!(
        tiny_quantity.try_mul(Quantity::new(1.0, tiny_factor)),
        Err(UnitError::InvalidScale(_))
    ));
    assert!(matches!(
        tiny_quantity.try_div(Quantity::new(1.0, huge_factor)),
        Err(UnitError::InvalidScale(_))
    ));
    assert_eq!(
        Quantity::new(2.0, NANOMETER).try_mul(Quantity::new(3.0, DIMENSIONLESS)),
        Ok(Quantity::new(6.0, NANOMETER))
    );
    assert_eq!(
        Quantity::new(6.0, NANOMETER).try_div(Quantity::new(3.0, DIMENSIONLESS)),
        Ok(Quantity::new(2.0, NANOMETER))
    );
}

#[test]
fn dimension_composition_never_wraps_in_any_build_profile() {
    let max = unit(Dimension::new([i32::MAX, 0, 0, 0, 0, 0, 0]), 1.0);
    let min = unit(Dimension::new([i32::MIN, 0, 0, 0, 0, 0, 0]), 1.0);
    for result in [
        max.try_mul(METER),
        min.try_div(METER),
        max.try_powi(2),
        min.try_powi(-1),
    ] {
        assert_eq!(result, Err(UnitError::DimensionOverflow));
    }
    assert_eq!(max.try_div(max), Ok(DIMENSIONLESS));
    assert_eq!(min.try_powi(0), Ok(DIMENSIONLESS));
}

#[test]
fn convenience_operators_fail_instead_of_constructing_invalid_units() {
    let tiny = unit(Dimension::LENGTH, 1e-200);
    let huge = unit(Dimension::LENGTH, 1e200);
    assert!(std::panic::catch_unwind(|| tiny * tiny).is_err());
    assert!(std::panic::catch_unwind(|| tiny / huge).is_err());
    assert!(std::panic::catch_unwind(|| tiny.powi(2)).is_err());
    assert!(
        std::panic::catch_unwind(|| { Quantity::new(1.0, tiny) * Quantity::new(1.0, tiny) })
            .is_err()
    );
    assert!(
        std::panic::catch_unwind(|| { Quantity::new(1.0, tiny) / Quantity::new(1.0, huge) })
            .is_err()
    );
    let max = unit(Dimension::new([i32::MAX, 0, 0, 0, 0, 0, 0]), 1.0);
    assert!(std::panic::catch_unwind(|| max * METER).is_err());
}

#[test]
fn conversions_reject_zero_and_infinite_intermediate_factors() {
    let small = unit(Dimension::LENGTH, 1e-300);
    let large = unit(Dimension::LENGTH, 1e100);
    for (from, to) in [(small, large), (large, small)] {
        assert_eq!(
            from.conversion_factor_to(to),
            Err(UnitError::UnrepresentableConversion { from, to })
        );
        assert_eq!(
            Quantity::new(1e300, from).value_in(to),
            Err(UnitError::UnrepresentableConversion { from, to })
        );
        assert!(matches!(
            Quantity::new(vec![1e300], from).into_unit(to),
            Err(UnitError::UnrepresentableConversion { .. })
        ));
    }
    // Positive subnormal ratios are representable and remain usable.
    let subnormal = unit(Dimension::LENGTH, f64::from_bits(1));
    assert_eq!(subnormal.conversion_factor_to(METER), Ok(f64::from_bits(1)));
    assert_eq!(small.conversion_factor_to(small), Ok(1.0));

    let too_large = unit(Dimension::LENGTH, 1e300);
    assert!(matches!(
        Positions::new(Quantity::new(vec![Point3::new(1.0, 0.0, 0.0)], too_large)),
        Err(PositionError::Unit(
            UnitError::UnrepresentableConversion { .. }
        ))
    ));
    let mut positions = Positions::zeros(1);
    let before = positions.clone();
    assert!(matches!(
        positions.set_position_at(0, Quantity::new(Point3::new(1.0, 0.0, 0.0), too_large)),
        Err(PositionError::Unit(
            UnitError::UnrepresentableConversion { .. }
        ))
    ));
    assert_eq!(positions, before);
}

#[test]
fn approximate_comparison_rejects_nonfinite_values_before_tolerance_arithmetic() {
    let finite = Quantity::new(1.0, DIMENSIONLESS);
    for tolerance in [-1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        assert_eq!(finite.is_close(&finite, tolerance, 0.0), Ok(false));
        assert_eq!(finite.is_close(&finite, 0.0, tolerance), Ok(false));
    }
    for invalid in [f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        let invalid = Quantity::new(invalid, DIMENSIONLESS);
        let finite = Quantity::new(1.0, DIMENSIONLESS);
        for tolerance in [0.0, 1e-6, 1.0, f64::MAX] {
            assert_eq!(finite.is_close(&invalid, tolerance, 0.0), Ok(false));
            assert_eq!(invalid.is_close(&finite, tolerance, 0.0), Ok(false));
            assert_eq!(invalid.is_close(&invalid, tolerance, 0.0), Ok(false));
        }
    }
    assert_eq!(
        Quantity::new(1.0, NANOMETER).is_close(&Quantity::new(f64::MAX, METER), 1e-6, 0.0),
        Ok(false)
    );
    assert!(matches!(
        Quantity::new(1.0, NANOMETER).is_close(
            &Quantity::new(1.0, unit(Dimension::LENGTH, 1e300)),
            1e-6,
            0.0
        ),
        Err(UnitError::UnrepresentableConversion { .. })
    ));
}

#[test]
fn approximate_comparison_handles_finite_extremes_without_overflow() {
    let positive = Quantity::new(f64::MAX, DIMENSIONLESS);
    let negative = Quantity::new(-f64::MAX, DIMENSIONLESS);
    assert_eq!(positive.is_close(&negative, 1.0, 0.0), Ok(false));
    assert_eq!(positive.is_close(&negative, 0.75, f64::MAX), Ok(false));
    assert_eq!(positive.is_close(&negative, 1.0, f64::MAX), Ok(true));
    assert_eq!(positive.is_close(&positive, 0.0, 0.0), Ok(true));
    assert_eq!(
        Quantity::new(0.0, DIMENSIONLESS).is_close(&Quantity::new(-0.0, DIMENSIONLESS), 0.0, 0.0),
        Ok(true)
    );
    let subnormal = f64::from_bits(1);
    assert_eq!(
        Quantity::new(subnormal, DIMENSIONLESS).is_close(
            &Quantity::new(0.0, DIMENSIONLESS),
            0.0,
            subnormal
        ),
        Ok(true)
    );
    assert_eq!(
        Quantity::new(1.0, NANOMETER).is_close(&Quantity::new(10.0, ANGSTROM), 1e-12, 0.0),
        Ok(true)
    );
}
