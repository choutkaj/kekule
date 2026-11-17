use super::log_sum_exp;

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} differs from {expected}"
    );
}

#[test]
fn log_sum_exp_is_the_log_of_the_sum() {
    let values = [1.0_f64.ln(), 2.0_f64.ln(), 3.0_f64.ln()];
    assert_close(log_sum_exp(values.iter().copied()), 6.0_f64.ln());
    // exp(1000) overflows; the shifted sum does not.
    assert_close(
        log_sum_exp([1000.0, 1000.0].into_iter()),
        1000.0 + 2.0_f64.ln(),
    );
    assert_close(
        log_sum_exp([-1000.0, -1000.0].into_iter()),
        -1000.0 + 2.0_f64.ln(),
    );
    assert_eq!(log_sum_exp(std::iter::empty()), f64::NEG_INFINITY);
}
