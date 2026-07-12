use std::fmt;

use crate::core::Point3;

use super::potential::{Potential, PotentialError, PotentialEvaluation, Vector3};
use super::{MolecularModel, PositionError};

#[derive(Debug, Clone, Copy, PartialEq)]
/// Controls normalized steepest descent with Armijo backtracking.
pub struct MinimizeOptions {
    /// Maximum number of accepted coordinate updates.
    pub max_iterations: usize,
    /// Convergence threshold for the maximum atom-gradient norm, in kJ/mol/angstrom.
    pub gradient_tolerance: f64,
    /// Initial maximum atom displacement for each line search, in angstroms.
    pub initial_step: f64,
    /// Smallest line-search displacement to consider, in angstroms.
    pub minimum_step: f64,
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
            gradient_tolerance: 1.0e-4,
            initial_step: 0.1,
            minimum_step: 1.0e-8,
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
    pub model: MolecularModel,
    pub initial_energy: f64,
    pub final_energy: f64,
    pub final_max_gradient: f64,
    pub iterations: usize,
    pub evaluations: usize,
    pub status: MinimizationStatus,
}

/// Minimize a cloned model while leaving the input unchanged.
///
/// Uses normalized steepest descent with an Armijo backtracking line search.
pub fn minimize(
    model: &MolecularModel,
    potential: &mut dyn Potential,
    options: MinimizeOptions,
) -> Result<MinimizationResult, MinimizationError> {
    validate_options(options)?;
    let mut working = model.clone();
    let mut evaluation = potential.evaluate(&working)?;
    let initial_energy = evaluation.energy();
    let mut evaluations = 1;
    let mut iterations = 0;

    loop {
        let max_gradient = maximum_gradient(evaluation.gradient());
        if max_gradient <= options.gradient_tolerance {
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
            .iter()
            .zip(&direction)
            .map(|(gradient, direction)| gradient.dot(*direction))
            .sum::<f64>();
        let current_positions = working.positions().to_vec();
        let current_energy = evaluation.energy();
        let mut step = options.initial_step;
        let mut accepted = None;

        for _ in 0..options.max_backtracks {
            if step < options.minimum_step {
                break;
            }
            let trial_positions = displaced_positions(&current_positions, &direction, step);
            working.set_positions(&trial_positions)?;
            let trial = match potential.evaluate(&working) {
                Ok(trial) => trial,
                Err(error) => {
                    working.set_positions(&current_positions)?;
                    return Err(MinimizationError::Potential(error));
                }
            };
            evaluations += 1;
            let armijo_limit =
                current_energy + options.armijo_coefficient * step * directional_derivative;
            if trial.energy() <= armijo_limit {
                accepted = Some(trial);
                break;
            }
            working.set_positions(&current_positions)?;
            step *= options.backtracking_factor;
        }

        let Some(trial) = accepted else {
            working.set_positions(&current_positions)?;
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

fn validate_options(options: MinimizeOptions) -> Result<(), MinimizationError> {
    if !options.gradient_tolerance.is_finite() || options.gradient_tolerance <= 0.0 {
        return Err(MinimizationError::InvalidOptions(
            "gradient tolerance must be finite and positive",
        ));
    }
    if !options.initial_step.is_finite() || options.initial_step <= 0.0 {
        return Err(MinimizationError::InvalidOptions(
            "initial step must be finite and positive",
        ));
    }
    if !options.minimum_step.is_finite()
        || options.minimum_step <= 0.0
        || options.minimum_step > options.initial_step
    {
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
    Ok(())
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
    model: MolecularModel,
    initial_energy: f64,
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
        final_max_gradient,
        iterations,
        evaluations,
        status,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MinimizationError {
    InvalidOptions(&'static str),
    Potential(PotentialError),
    Position(PositionError),
}

impl fmt::Display for MinimizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => write!(f, "invalid minimization options: {message}"),
            Self::Potential(error) => write!(f, "potential evaluation failed: {error}"),
            Self::Position(error) => write!(f, "cannot update model positions: {error}"),
        }
    }
}

impl std::error::Error for MinimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidOptions(_) => None,
            Self::Potential(error) => Some(error),
            Self::Position(error) => Some(error),
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
