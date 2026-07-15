use std::fmt;

use crate::core::Point3;
use crate::units::{Quantity, UnitError, MODEL_GRADIENT_UNIT, MODEL_LENGTH_UNIT};

use super::potential::{Potential, PotentialError, PotentialEvaluation, Vector3};
use super::{Model, PositionError};

#[derive(Debug, Clone, Copy, PartialEq)]
/// Controls normalized steepest descent with Armijo backtracking.
pub struct MinimizeOptions {
    /// Maximum number of accepted coordinate updates.
    pub max_iterations: usize,
    /// Convergence threshold for the maximum atom-gradient norm.
    pub gradient_tolerance: Quantity<f64>,
    /// Initial maximum atom displacement for each line search.
    pub initial_step: Quantity<f64>,
    /// Smallest line-search displacement to consider.
    pub minimum_step: Quantity<f64>,
    /// Multiplicative line-search step reduction, strictly between zero and one.
    pub backtracking_factor: f64,
    /// Armijo sufficient-decrease coefficient, strictly between zero and one.
    pub armijo_coefficient: f64,
    /// Maximum potential evaluations attempted by each line search.
    pub max_backtracks: usize,
}

impl Default for MinimizeOptions {
    fn default() -> Self {
        Self {
            max_iterations: 1_000,
            gradient_tolerance: Quantity::new(1.0e-4, MODEL_GRADIENT_UNIT),
            initial_step: Quantity::new(0.1, MODEL_LENGTH_UNIT),
            minimum_step: Quantity::new(1.0e-8, MODEL_LENGTH_UNIT),
            backtracking_factor: 0.5,
            armijo_coefficient: 1.0e-4,
            max_backtracks: 24,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Terminal state of a non-error minimization run.
pub enum MinimizationStatus {
    Converged,
    MaxIterations,
    LineSearchStalled,
}

#[derive(Debug, Clone, PartialEq)]
/// Minimized model and convergence diagnostics.
pub struct MinimizationResult {
    pub model: Model,
    pub initial_energy: Quantity<f64>,
    pub final_energy: Quantity<f64>,
    pub final_max_gradient: Quantity<f64>,
    pub iterations: usize,
    pub evaluations: usize,
    pub status: MinimizationStatus,
}

/// Minimize a cloned model while leaving the input unchanged.
///
/// Uses normalized steepest descent with an Armijo backtracking line search.
/// Trial evaluations with [`PotentialError::InvalidGeometry`] are rejected and
/// backtracked; other potential failures abort the minimization.
pub fn minimize(
    model: &Model,
    potential: &mut dyn Potential,
    options: MinimizeOptions,
) -> Result<MinimizationResult, MinimizationError> {
    let validated = validate_options(options)?;
    let mut working = model.clone();
    let mut evaluation = potential.evaluate(&working)?;
    let initial_energy = evaluation.energy();
    let mut evaluations = 1;
    let mut iterations = 0;

    loop {
        let max_gradient = maximum_gradient(evaluation.gradient().value());
        if max_gradient <= validated.gradient_tolerance {
            return Ok(result(
                working,
                initial_energy,
                &evaluation,
                max_gradient,
                iterations,
                evaluations,
                MinimizationStatus::Converged,
            ));
        }
        if iterations >= options.max_iterations {
            return Ok(result(
                working,
                initial_energy,
                &evaluation,
                max_gradient,
                iterations,
                evaluations,
                MinimizationStatus::MaxIterations,
            ));
        }

        let direction = evaluation
            .gradient()
            .value()
            .iter()
            .map(|gradient| {
                Vector3::new(
                    -gradient.x / max_gradient,
                    -gradient.y / max_gradient,
                    -gradient.z / max_gradient,
                )
            })
            .collect::<Vec<_>>();
        let directional_derivative = evaluation
            .gradient()
            .value()
            .iter()
            .zip(&direction)
            .map(|(gradient, direction)| gradient.dot(*direction))
            .sum::<f64>();
        let current_positions = working.positions_value().to_vec();
        let current_energy = evaluation.energy().into_value();
        let mut step = validated.initial_step;
        let mut accepted = None;

        for _ in 0..options.max_backtracks {
            if step < validated.minimum_step {
                break;
            }
            let trial_positions = displaced_positions(&current_positions, &direction, step);
            working.set_positions(Quantity::new(trial_positions, MODEL_LENGTH_UNIT))?;
            let trial_result = potential.evaluate(&working);
            evaluations += 1;
            let trial = match trial_result {
                Ok(trial) => trial,
                Err(error) if error.is_invalid_geometry() => {
                    working.set_positions(Quantity::new(&current_positions, MODEL_LENGTH_UNIT))?;
                    step *= options.backtracking_factor;
                    continue;
                }
                Err(error) => {
                    working.set_positions(Quantity::new(&current_positions, MODEL_LENGTH_UNIT))?;
                    return Err(MinimizationError::Potential(error));
                }
            };
            let armijo_limit =
                current_energy + options.armijo_coefficient * step * directional_derivative;
            if trial.energy().into_value() <= armijo_limit {
                accepted = Some(trial);
                break;
            }
            working.set_positions(Quantity::new(&current_positions, MODEL_LENGTH_UNIT))?;
            step *= options.backtracking_factor;
        }

        let Some(trial) = accepted else {
            working.set_positions(Quantity::new(&current_positions, MODEL_LENGTH_UNIT))?;
            return Ok(result(
                working,
                initial_energy,
                &evaluation,
                max_gradient,
                iterations,
                evaluations,
                MinimizationStatus::LineSearchStalled,
            ));
        };
        evaluation = trial;
        iterations += 1;
    }
}

struct ValidatedOptions {
    gradient_tolerance: f64,
    initial_step: f64,
    minimum_step: f64,
}

fn validate_options(options: MinimizeOptions) -> Result<ValidatedOptions, MinimizationError> {
    let gradient_tolerance = options.gradient_tolerance.value_in(MODEL_GRADIENT_UNIT)?;
    let initial_step = options.initial_step.value_in(MODEL_LENGTH_UNIT)?;
    let minimum_step = options.minimum_step.value_in(MODEL_LENGTH_UNIT)?;
    if !gradient_tolerance.is_finite() || gradient_tolerance <= 0.0 {
        return Err(MinimizationError::InvalidOptions(
            "gradient tolerance must be finite and positive",
        ));
    }
    if !initial_step.is_finite() || initial_step <= 0.0 {
        return Err(MinimizationError::InvalidOptions(
            "initial step must be finite and positive",
        ));
    }
    if !minimum_step.is_finite() || minimum_step <= 0.0 || minimum_step > initial_step {
        return Err(MinimizationError::InvalidOptions(
            "minimum step must be finite, positive, and no larger than the initial step",
        ));
    }
    if !options.backtracking_factor.is_finite()
        || options.backtracking_factor <= 0.0
        || options.backtracking_factor >= 1.0
    {
        return Err(MinimizationError::InvalidOptions(
            "backtracking factor must be finite and between zero and one",
        ));
    }
    if !options.armijo_coefficient.is_finite()
        || options.armijo_coefficient <= 0.0
        || options.armijo_coefficient >= 1.0
    {
        return Err(MinimizationError::InvalidOptions(
            "Armijo coefficient must be finite and between zero and one",
        ));
    }
    if options.max_backtracks == 0 {
        return Err(MinimizationError::InvalidOptions(
            "maximum backtracks must be at least one",
        ));
    }
    Ok(ValidatedOptions {
        gradient_tolerance,
        initial_step,
        minimum_step,
    })
}

fn displaced_positions(positions: &[Point3], direction: &[Vector3], step: f64) -> Vec<Point3> {
    positions
        .iter()
        .zip(direction)
        .map(|(position, direction)| {
            Point3::new(
                position.x + step * direction.x,
                position.y + step * direction.y,
                position.z + step * direction.z,
            )
        })
        .collect()
}

fn maximum_gradient(gradient: &[Vector3]) -> f64 {
    gradient
        .iter()
        .map(|gradient| gradient.norm())
        .fold(0.0, f64::max)
}

fn result(
    model: Model,
    initial_energy: Quantity<f64>,
    evaluation: &PotentialEvaluation,
    final_max_gradient: f64,
    iterations: usize,
    evaluations: usize,
    status: MinimizationStatus,
) -> MinimizationResult {
    MinimizationResult {
        model,
        initial_energy,
        final_energy: evaluation.energy(),
        final_max_gradient: Quantity::new(final_max_gradient, MODEL_GRADIENT_UNIT),
        iterations,
        evaluations,
        status,
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum MinimizationError {
    InvalidOptions(&'static str),
    Potential(PotentialError),
    Position(PositionError),
    Unit(UnitError),
}

impl fmt::Display for MinimizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => write!(f, "invalid minimization options: {message}"),
            Self::Potential(error) => write!(f, "potential evaluation failed: {error}"),
            Self::Position(error) => write!(f, "cannot update model positions: {error}"),
            Self::Unit(error) => write!(f, "invalid minimization quantity unit: {error}"),
        }
    }
}

impl std::error::Error for MinimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidOptions(_) => None,
            Self::Potential(error) => Some(error),
            Self::Position(error) => Some(error),
            Self::Unit(error) => Some(error),
        }
    }
}

impl From<PotentialError> for MinimizationError {
    fn from(error: PotentialError) -> Self {
        Self::Potential(error)
    }
}

impl From<PositionError> for MinimizationError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<UnitError> for MinimizationError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
