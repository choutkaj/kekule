//! Estimators of a group's population from the reduced values of its
//! samples.

/// The highest cumulant that [`Estimator::Cumulant`] keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum CumulantOrder {
    /// The mean alone.
    First = 1,
    /// The mean and the variance; exact in expectation for Gaussian reduced
    /// values.
    Second = 2,
    /// The mean, the variance, and the third central moment.
    Third = 3,
}

/// How the population `Σ exp(b_i)` of a group is estimated from the reduced
/// values `b_i` of its samples.
///
/// A reduced value is a dimensionless log importance weight, such as a bias
/// potential over the thermal energy, `βΔV`. The exponential sum is exact but
/// dominated by the few largest values. The cumulant expansion trades bias for
/// lower variance (Miao et al., J. Chem. Theory Comput. 10, 2677, 2014); to the
/// second order, the default, it is exact in expectation for Gaussian reduced
/// values.
///
/// For a group of `n` samples:
///
/// | Estimator | `ln Σ exp(b_i)` |
/// | --- | --- |
/// | `Exponential` | `ln Σ exp(b_i)` |
/// | `Cumulant(First)` | `ln n + κ1` |
/// | `Cumulant(Second)` | `ln n + κ1 + κ2 / 2` |
/// | `Cumulant(Third)` | `ln n + κ1 + κ2 / 2 + κ3 / 6` |
///
/// `κ1` is the mean, `κ2` the variance, and `κ3` the third central moment of
/// the reduced values, each over `n`. Adding a constant to every reduced
/// value of a group adds it to the group's estimate, with every estimator.
///
/// An estimator gives the totals of groups only. How a group's population is
/// shared among its samples is not an estimate: each sample's share is
/// `exp(b_i) / Σ exp(b_j)`.
///
/// ```
/// use kekule::statistics::{CumulantOrder, Estimator};
///
/// let reduced = [0.5, 1.5];
/// let exact = Estimator::Exponential.log_sum(&reduced).unwrap();
/// assert!((exact - (0.5_f64.exp() + 1.5_f64.exp()).ln()).abs() < 1e-15);
/// // ln 2 + 1 + 0.25 / 2: the mean and half the variance.
/// assert_eq!(Estimator::default(), Estimator::Cumulant(CumulantOrder::Second));
/// let cumulant = Estimator::default().log_sum(&reduced).unwrap();
/// assert!((cumulant - (2.0_f64.ln() + 1.125)).abs() < 1e-15);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Estimator {
    /// The exponential sum itself.
    Exponential,
    /// The cumulant expansion of `ln⟨exp(b)⟩` to the given order.
    Cumulant(CumulantOrder),
}

impl Default for Estimator {
    /// The second-order cumulant expansion.
    fn default() -> Self {
        Self::Cumulant(CumulantOrder::Second)
    }
}

impl Estimator {
    /// The estimated `ln Σ exp(b_i)` over the reduced values of a group.
    ///
    /// `None` when the group is empty, a reduced value is not finite, or the
    /// estimate or one of its cumulant terms is not finite.
    pub fn log_sum(self, reduced: &[f64]) -> Option<f64> {
        if reduced.is_empty() || reduced.iter().any(|value| !value.is_finite()) {
            return None;
        }
        let estimate = match self {
            Self::Exponential => log_sum_exp(reduced.iter().copied()),
            Self::Cumulant(order) => (reduced.len() as f64).ln() + cumulant_sum(reduced, order),
        };
        estimate.is_finite().then_some(estimate)
    }
}

/// `κ1 + κ2 / 2 + κ3 / 6` up to `order`, with population moments.
///
/// The values are divided by a power of two near the largest magnitude, which
/// is exact and leaves every scaled value below 2, so no sum overflows. The
/// mean is refined by a second pass, so the central moments keep their digits
/// under a large common offset of the values. Each cumulant is scaled back
/// on its own, and a term that does not fit makes the sum not finite.
fn cumulant_sum(values: &[f64], order: CumulantOrder) -> f64 {
    let largest = values
        .iter()
        .fold(0.0_f64, |largest, value| largest.max(value.abs()));
    if largest == 0.0 {
        return 0.0;
    }
    if let [value] = values {
        return *value;
    }
    let scale = power_of_two_below(largest);
    let count = values.len() as f64;
    let scaled = || values.iter().map(|value| value / scale);
    let (lowest, highest) = scaled().fold(
        (f64::INFINITY, f64::NEG_INFINITY),
        |(lowest, highest), value| (lowest.min(value), highest.max(value)),
    );
    let estimate = scaled().map(|value| value / count).sum::<f64>();
    // The mean lies between the extremes; rounding must not push it past
    // them, where the largest finite values would overflow when scaled back.
    let mean = (estimate
        + scaled()
            .map(|value| (value - estimate) / count)
            .sum::<f64>())
    .clamp(lowest, highest);
    let mut sum = mean * scale;
    if order == CumulantOrder::First {
        return sum;
    }
    let (second, third) = scaled().fold((0.0, 0.0), |(second, third), value| {
        let deviation = value - mean;
        let share = deviation * deviation / count;
        (second + share, third + share * deviation)
    });
    sum += second / 2.0 * scale * scale;
    if order == CumulantOrder::Third {
        sum += third / 6.0 * scale * scale * scale;
    }
    sum
}

/// The largest power of two that is at most a positive, finite value, or the
/// smallest normal one for a subnormal value.
fn power_of_two_below(value: f64) -> f64 {
    const EXPONENT_SHIFT: u32 = 52;
    let exponent = value.to_bits() >> EXPONENT_SHIFT;
    f64::from_bits(exponent.max(1) << EXPONENT_SHIFT)
}

/// `ln Σ exp(x)` over finite values, shifted by their largest so that no
/// term overflows; negative infinity for none.
pub(crate) fn log_sum_exp(values: impl Iterator<Item = f64> + Clone) -> f64 {
    let maximum = values.clone().fold(f64::NEG_INFINITY, f64::max);
    if maximum == f64::NEG_INFINITY {
        return maximum;
    }
    maximum
        + values
            .map(|value| (value - maximum).exp())
            .sum::<f64>()
            .ln()
}

#[cfg(test)]
mod tests;
