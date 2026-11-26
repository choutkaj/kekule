//! Landscape estimators: exact constant groups, the cumulant series, and
//! values that no estimate represents.

use kekule::statistics::{CumulantOrder, Estimator};

const EVERY_ESTIMATOR: [Estimator; 4] = [
    Estimator::Exponential,
    Estimator::Cumulant(CumulantOrder::First),
    Estimator::Cumulant(CumulantOrder::Second),
    Estimator::Cumulant(CumulantOrder::Third),
];

fn assert_close(actual: f64, expected: f64, estimator: Estimator) {
    let tolerance = 8.0 * f64::EPSILON * actual.abs().max(expected.abs()).max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "{estimator:?}: {actual} != {expected}"
    );
}

#[test]
fn constant_groups_reweight_exactly_with_every_estimator() {
    for (value, count) in [(1.0, 2), (0.25, 4), (0.0, 1), (-3.5, 7), (812.5, 3)] {
        let reduced = vec![value; count];
        let exact = (count as f64).ln() + value;
        for estimator in EVERY_ESTIMATOR {
            assert_close(estimator.log_sum(&reduced).unwrap(), exact, estimator);
        }
    }
}

#[test]
fn estimators_follow_their_series() {
    // Two reduced values a ± d.
    let (a, d) = (1.2_f64, 0.4_f64);
    let reduced = [a - d, a + d];
    let cases = [
        (Estimator::Exponential, a + d.cosh().ln()),
        (Estimator::Cumulant(CumulantOrder::First), a),
        (Estimator::Cumulant(CumulantOrder::Second), a + d * d / 2.0),
        (Estimator::Cumulant(CumulantOrder::Third), a + d * d / 2.0),
    ];
    for (estimator, log_average) in cases {
        assert_close(
            estimator.log_sum(&reduced).unwrap(),
            2.0_f64.ln() + log_average,
            estimator,
        );
    }

    // Mean 1, variance (1 + 1 + 4) / 3 = 2, third central moment
    // (-1 - 1 + 8) / 3 = 2.
    let skewed = [0.0, 0.0, 3.0];
    let third = Estimator::Cumulant(CumulantOrder::Third);
    assert_close(
        third.log_sum(&skewed).unwrap(),
        3.0_f64.ln() + 1.0 + 1.0 + 2.0 / 6.0,
        third,
    );
}

#[test]
fn a_common_offset_shifts_every_estimate_by_itself() {
    let deviations = [-0.7, 0.1, 0.4, 1.3, -0.2];
    let offset = 1.0e6;
    let shifted = deviations.map(|value| value + offset);
    for estimator in EVERY_ESTIMATOR {
        let base = estimator.log_sum(&deviations).unwrap();
        let estimate = estimator.log_sum(&shifted).unwrap();
        // The shifted values themselves are rounded to the offset's spacing.
        assert!(
            (estimate - offset - base).abs() < 1e-9,
            "{estimator:?}: {estimate} - {offset} != {base}"
        );
    }
    // Neither sum overflows although exp(b) would.
    assert_close(
        Estimator::Exponential.log_sum(&[1000.0, 1000.0]).unwrap(),
        1000.0 + 2.0_f64.ln(),
        Estimator::Exponential,
    );
    // Groups of the largest finite value, of sizes that divide it inexactly
    // too.
    for count in [2, 3, 5, 7, 1000] {
        let largest = vec![f64::MAX; count];
        for estimator in EVERY_ESTIMATOR {
            assert_close(estimator.log_sum(&largest).unwrap(), f64::MAX, estimator);
        }
    }
    let negative = vec![-f64::MAX; 3];
    for estimator in EVERY_ESTIMATOR {
        assert_close(estimator.log_sum(&negative).unwrap(), -f64::MAX, estimator);
    }
}

#[test]
fn moments_are_finite_whenever_the_estimate_is() {
    let second = Estimator::Cumulant(CumulantOrder::Second);
    let third = Estimator::Cumulant(CumulantOrder::Third);
    // Each squared deviation is 1e308, so their sum overflows, but the
    // variance is 1e308 and the estimate ln n + 5e307.
    for count in [4, 7, 1000] {
        let reduced = (0..count)
            .map(|index| if index % 2 == 0 { -1.0e154 } else { 1.0e154 })
            .collect::<Vec<_>>();
        let n = count as f64;
        // An odd count leaves a mean of -1e154 / n.
        let mean = if count % 2 == 0 { 0.0 } else { -1.0e154 / n };
        let variance = 1.0e308 - mean * mean;
        let expected = n.ln() + mean + variance / 2.0;
        let estimate = second.log_sum(&reduced).unwrap();
        assert!(
            (estimate - expected).abs() <= 1e-12 * expected,
            "{count}: {estimate} != {expected}"
        );
    }
    // The cubes overflow, but the symmetric group's third central moment is
    // zero: ln 6 + 2e206 / 6 / 2.
    let symmetric = [-1.0e103, 0.0, 0.0, 0.0, 0.0, 1.0e103];
    let expected = 6.0_f64.ln() + 2.0e206 / 12.0;
    for estimator in [second, third] {
        let estimate = estimator.log_sum(&symmetric).unwrap();
        assert!(
            (estimate - expected).abs() <= 1e-12 * expected,
            "{estimator:?}: {estimate} != {expected}"
        );
    }
}

#[test]
fn empty_groups_and_values_without_an_estimate_have_none() {
    for estimator in EVERY_ESTIMATOR {
        assert_eq!(estimator.log_sum(&[]), None, "{estimator:?}");
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(estimator.log_sum(&[0.0, value]), None, "{estimator:?}");
        }
    }
    // The variance of values this far apart is not finite.
    let spread = [-1.0e200, 1.0e200];
    assert_eq!(
        Estimator::Cumulant(CumulantOrder::Second).log_sum(&spread),
        None
    );
    assert_eq!(
        Estimator::Cumulant(CumulantOrder::First).log_sum(&spread),
        Some(2.0_f64.ln())
    );
    assert_eq!(Estimator::Exponential.log_sum(&spread), Some(1.0e200));
}

#[test]
fn the_default_is_the_exact_exponential_average() {
    assert_eq!(Estimator::default(), Estimator::Exponential);
}

#[test]
fn the_second_cumulant_of_gaussian_values_is_their_exponential_average() {
    // Reduced values drawn from N(mu, sigma^2), for which
    // ln <exp b> = mu + sigma^2 / 2.
    let (mu, sigma) = (1.0, 0.5);
    let count = 4_000;
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut uniform = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1_u64 << 53) as f64
    };
    let values = (0..count)
        .map(|_| {
            let (radius, angle) = (1.0 - uniform(), uniform());
            mu + sigma * (-2.0 * radius.ln()).sqrt() * (std::f64::consts::TAU * angle).cos()
        })
        .collect::<Vec<_>>();
    let expected = (count as f64).ln() + mu + sigma * sigma / 2.0;
    for estimator in [
        Estimator::Exponential,
        Estimator::Cumulant(CumulantOrder::Second),
    ] {
        let estimate = estimator.log_sum(&values).unwrap();
        assert!(
            (estimate - expected).abs() < 0.03,
            "{estimator:?}: {estimate} differs from {expected}"
        );
    }
}
