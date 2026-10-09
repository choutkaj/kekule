use kekule::geometry::Vector3;

#[test]
fn vector_norm_preserves_representable_extreme_magnitudes() {
    for scale in [1e-200, 1.0, 1e200] {
        let norm = Vector3::new(3.0 * scale, -4.0 * scale, 0.0).norm();
        assert!((norm / (5.0 * scale) - 1.0).abs() < 1e-15);
    }
    for magnitude in [f64::from_bits(1), f64::MIN_POSITIVE, f64::MAX] {
        assert_eq!(Vector3::new(0.0, 0.0, magnitude).norm(), magnitude);
    }
    assert_eq!(Vector3::zero().norm(), 0.0);
}
