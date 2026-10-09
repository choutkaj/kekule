//! Local geometry optimization over any [`Potential`].
//!
//! [`minimize`] runs limited-memory BFGS (L-BFGS) with a strong-Wolfe line
//! search on all Cartesian coordinates. It reads a borrowed [`ModelView`], so
//! models, ensemble members, and trajectory frames are minimized the same way,
//! and returns new [`Positions`] without changing its input. No constraints or
//! fixed atoms are applied. Periodic cells are passed through to the potential
//! unchanged.
//!
//! Each line search caps the largest atom displacement at
//! [`MinimizeOptions::max_step`]. A trial whose coordinates are singular for the
//! potential ([`EvaluationError::InvalidGeometry`]) is treated as an energy
//! increase and the step is shortened; every other evaluation failure aborts.
//! If a line search cannot reduce the energy, the curvature history is discarded
//! and a steepest-descent step is tried before reporting
//! [`MinimizationStatus::LineSearchFailed`].

use std::collections::VecDeque;
use std::fmt;
use std::ops::ControlFlow;
use std::sync::Arc;

use kekule::geometry::{Point3, Vector3};
use kekule::structure::{Model, ModelError, ModelView, PositionError, Positions};
use kekule::topology::Topology;
use kekule::units::{
    Quantity, UnitError, CANONICAL_GRADIENT_UNIT, CANONICAL_LENGTH_UNIT,
    KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER,
};

use crate::potential::{Energy, Evaluation, EvaluationError, Potential};

const ARMIJO: f64 = 1e-4;
const CURVATURE: f64 = 0.9;

/// Controls [`minimize`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MinimizeOptions {
    /// Maximum number of accepted coordinate updates.
    pub max_iterations: usize,
    /// Convergence threshold for the largest per-atom gradient norm.
    pub gradient_tolerance: Quantity<f64>,
    /// Largest displacement of any atom in one line search.
    pub max_step: Quantity<f64>,
    /// Number of curvature pairs retained by L-BFGS.
    pub memory: usize,
    /// Potential evaluations allowed in one line search.
    pub max_line_search_evaluations: usize,
}

impl Default for MinimizeOptions {
    fn default() -> Self {
        Self {
            max_iterations: 1_000,
            gradient_tolerance: Quantity::new(1.0, KILOJOULE_PER_MOLE_PER_NANOMETER),
            max_step: Quantity::new(0.03, NANOMETER),
            memory: 10,
            max_line_search_evaluations: 20,
        }
    }
}

/// Terminal state of a minimization that did not fail.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MinimizationStatus {
    /// The largest per-atom gradient norm reached the tolerance.
    Converged,
    MaxIterations,
    /// No step along a steepest-descent direction reduced the energy.
    LineSearchFailed,
    /// The observer requested an early stop.
    Stopped,
}

/// State after one accepted step, passed to a [`minimize_with_observer`] callback.
#[derive(Debug, Clone, Copy)]
pub struct MinimizationStep<'a> {
    /// Number of accepted steps so far, starting at one.
    pub iteration: usize,
    pub energy: &'a Energy,
    pub max_gradient: Quantity<f64>,
    pub positions: &'a Positions,
}

/// Minimized coordinates and convergence diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct Minimization {
    positions: Positions,
    initial_energy: Energy,
    evaluation: Evaluation,
    iterations: usize,
    evaluations: usize,
    status: MinimizationStatus,
}

impl Minimization {
    /// Final coordinates in dense topology atom order.
    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    pub fn into_positions(self) -> Positions {
        self.positions
    }

    pub fn initial_energy(&self) -> &Energy {
        &self.initial_energy
    }

    /// Energy and gradient at the final coordinates.
    pub fn final_evaluation(&self) -> &Evaluation {
        &self.evaluation
    }

    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// Potential evaluations performed, including rejected line-search trials.
    pub fn evaluations(&self) -> usize {
        self.evaluations
    }

    pub fn status(&self) -> MinimizationStatus {
        self.status
    }

    /// Copies `source` with the minimized coordinates.
    ///
    /// `source` must share the minimized topology snapshot. Its cell and
    /// realization properties are retained unchanged.
    pub fn to_model(&self, source: ModelView<'_>) -> Result<Model, MinimizationError> {
        if !self.evaluation.topology().shares_layout(source.topology()) {
            return Err(MinimizationError::IncompatibleTopology);
        }
        let mut model = source.to_model();
        model.set_positions(self.positions.values())?;
        Ok(model)
    }
}

/// Minimizes the energy of `model` under `potential`.
pub fn minimize<P: Potential + ?Sized>(
    potential: &P,
    model: ModelView<'_>,
    options: &MinimizeOptions,
) -> Result<Minimization, MinimizationError> {
    minimize_with_observer(potential, model, options, |_| ControlFlow::Continue(()))
}

/// Minimizes like [`minimize`], calling `observer` after every accepted step.
///
/// Returning [`ControlFlow::Break`] stops with [`MinimizationStatus::Stopped`]
/// at the step just observed. Recording the observed positions, for example
/// into a trajectory, is the caller's choice.
pub fn minimize_with_observer<P: Potential + ?Sized>(
    potential: &P,
    model: ModelView<'_>,
    options: &MinimizeOptions,
    mut observer: impl FnMut(&MinimizationStep<'_>) -> ControlFlow<()>,
) -> Result<Minimization, MinimizationError> {
    let limits = validate(options)?;
    let initial = potential.evaluate(model)?;
    let mut search = Search {
        potential,
        topology: model.shared_topology(),
        model,
        evaluations: 1,
    };
    let mut current = Point {
        x: flatten_points(model.positions().values().into_value()),
        gradient: flatten_vectors(initial.gradient().into_value()),
        positions: model.positions().clone(),
        evaluation: initial,
    };
    let initial_energy = current.evaluation.energy().clone();
    let mut history: VecDeque<Pair> = VecDeque::with_capacity(options.memory);
    let mut iterations = 0;

    let status = loop {
        let max_gradient = max_atom_norm(&current.gradient);
        if !max_gradient.is_finite() {
            return Err(MinimizationError::NumericalFailure(
                "maximum gradient norm is not representable",
            ));
        }
        if max_gradient <= limits.gradient_tolerance {
            break MinimizationStatus::Converged;
        }
        if iterations >= options.max_iterations {
            break MinimizationStatus::MaxIterations;
        }

        let mut direction = two_loop(&history, &current.gradient);
        let mut slope = dot(&current.gradient, &direction);
        if !(slope.is_finite() && slope < 0.0) {
            history.clear();
            direction = current.gradient.iter().map(|g| -g).collect();
            slope = dot(&current.gradient, &direction);
            if !(slope.is_finite() && slope < 0.0) {
                return Err(MinimizationError::NumericalFailure(
                    "steepest-descent slope is not finite and negative",
                ));
            }
        }
        let alpha_max = limits.max_step / max_atom_norm(&direction);
        if !(alpha_max.is_finite() && alpha_max > 0.0) {
            return Err(MinimizationError::NumericalFailure(
                "step length bound is not representable",
            ));
        }
        // Without curvature information the direction has no natural scale, so
        // the first trial uses the full displacement bound.
        let alpha = if history.is_empty() {
            alpha_max
        } else {
            alpha_max.min(1.0)
        };
        let Some(next) = search.line_search(
            &current,
            &direction,
            slope,
            alpha,
            alpha_max,
            options.max_line_search_evaluations,
        )?
        else {
            if history.is_empty() {
                break MinimizationStatus::LineSearchFailed;
            }
            history.clear();
            continue;
        };

        let s: Vec<f64> = next.x.iter().zip(&current.x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = next
            .gradient
            .iter()
            .zip(&current.gradient)
            .map(|(a, b)| a - b)
            .collect();
        let sy = dot(&s, &y);
        // Skip pairs without positive curvature to keep the inverse Hessian
        // approximation positive definite.
        if sy.is_finite() && sy > f64::EPSILON * norm(&s) * norm(&y) {
            if history.len() == options.memory {
                history.pop_front();
            }
            history.push_back(Pair {
                rho: 1.0 / sy,
                s,
                y,
            });
        }
        current = next;
        iterations += 1;
        let step = MinimizationStep {
            iteration: iterations,
            energy: current.evaluation.energy(),
            max_gradient: Quantity::new(max_atom_norm(&current.gradient), CANONICAL_GRADIENT_UNIT),
            positions: &current.positions,
        };
        if observer(&step).is_break() {
            break MinimizationStatus::Stopped;
        }
    };

    Ok(Minimization {
        positions: current.positions,
        initial_energy,
        evaluation: current.evaluation,
        iterations,
        evaluations: search.evaluations,
        status,
    })
}

struct Limits {
    gradient_tolerance: f64,
    max_step: f64,
}

fn validate(options: &MinimizeOptions) -> Result<Limits, MinimizationError> {
    let gradient_tolerance = options
        .gradient_tolerance
        .value_in(CANONICAL_GRADIENT_UNIT)?;
    let max_step = options.max_step.value_in(CANONICAL_LENGTH_UNIT)?;
    if !(gradient_tolerance.is_finite() && gradient_tolerance > 0.0) {
        return Err(MinimizationError::InvalidOptions(
            "gradient tolerance must be finite and positive",
        ));
    }
    if !(max_step.is_finite() && max_step > 0.0) {
        return Err(MinimizationError::InvalidOptions(
            "maximum step must be finite and positive",
        ));
    }
    if options.memory == 0 {
        return Err(MinimizationError::InvalidOptions(
            "L-BFGS memory must be at least one",
        ));
    }
    if options.max_line_search_evaluations == 0 {
        return Err(MinimizationError::InvalidOptions(
            "line searches need at least one evaluation",
        ));
    }
    Ok(Limits {
        gradient_tolerance,
        max_step,
    })
}

/// One accepted or trial configuration with flattened canonical coordinates.
struct Point {
    x: Vec<f64>,
    gradient: Vec<f64>,
    positions: Positions,
    evaluation: Evaluation,
}

impl Point {
    fn energy(&self) -> f64 {
        self.evaluation.energy().total().into_value()
    }
}

struct Pair {
    s: Vec<f64>,
    y: Vec<f64>,
    rho: f64,
}

/// The L-BFGS two-loop recursion: returns `-H g`.
fn two_loop(history: &VecDeque<Pair>, gradient: &[f64]) -> Vec<f64> {
    let mut q = gradient.to_vec();
    let mut alphas = Vec::with_capacity(history.len());
    for pair in history.iter().rev() {
        let alpha = pair.rho * dot(&pair.s, &q);
        axpy(&mut q, -alpha, &pair.y);
        alphas.push(alpha);
    }
    if let Some(last) = history.back() {
        let gamma = 1.0 / (last.rho * dot(&last.y, &last.y));
        q.iter_mut().for_each(|v| *v *= gamma);
    }
    for (pair, alpha) in history.iter().zip(alphas.into_iter().rev()) {
        let beta = pair.rho * dot(&pair.y, &q);
        axpy(&mut q, alpha - beta, &pair.s);
    }
    q.iter_mut().for_each(|v| *v = -*v);
    q
}

struct Search<'a, 'm, P: ?Sized> {
    potential: &'a P,
    topology: Arc<Topology>,
    model: ModelView<'m>,
    evaluations: usize,
}

/// A line-search sample: `phi = E(x + alpha d)` and `dphi = g . d`.
struct Sample {
    alpha: f64,
    phi: f64,
    dphi: f64,
    point: Option<Point>,
}

impl<P: Potential + ?Sized> Search<'_, '_, P> {
    /// Evaluates `x + alpha d`. Singular coordinates yield an infinite energy;
    /// coordinates that round back to `x` yield `None`.
    fn sample(
        &mut self,
        from: &Point,
        direction: &[f64],
        alpha: f64,
    ) -> Result<Option<Sample>, MinimizationError> {
        let x: Vec<f64> = from
            .x
            .iter()
            .zip(direction)
            .map(|(x, d)| x + alpha * d)
            .collect();
        if x == from.x {
            return Ok(None);
        }
        let positions = Positions::from_vec(Quantity::new(to_points(&x), CANONICAL_LENGTH_UNIT))?;
        let view = ModelView::new(
            &self.topology,
            &positions,
            self.model.cell(),
            self.model.properties(),
        )?;
        self.evaluations += 1;
        match self.potential.evaluate(view) {
            Ok(evaluation) => {
                let gradient = flatten_vectors(evaluation.gradient().into_value());
                let point = Point {
                    x,
                    gradient,
                    positions,
                    evaluation,
                };
                Ok(Some(Sample {
                    alpha,
                    phi: point.energy(),
                    dphi: dot(&point.gradient, direction),
                    point: Some(point),
                }))
            }
            Err(error) if error.is_invalid_geometry() => Ok(Some(Sample {
                alpha,
                phi: f64::INFINITY,
                dphi: f64::NAN,
                point: None,
            })),
            Err(error) => Err(error.into()),
        }
    }

    /// Strong-Wolfe line search (Nocedal and Wright, algorithms 3.5 and 3.6).
    ///
    /// Returns the accepted point, or the best sufficient-decrease point found
    /// when the evaluation budget runs out, or `None` without any decrease.
    fn line_search(
        &mut self,
        from: &Point,
        direction: &[f64],
        slope: f64,
        mut alpha: f64,
        alpha_max: f64,
        budget: usize,
    ) -> Result<Option<Point>, MinimizationError> {
        let phi0 = from.energy();
        let armijo = |s: &Sample| s.phi.is_finite() && s.phi <= phi0 + ARMIJO * s.alpha * slope;
        let curvature = |s: &Sample| s.dphi.abs() <= -CURVATURE * slope;
        let start = Sample {
            alpha: 0.0,
            phi: phi0,
            dphi: slope,
            point: None,
        };
        let mut used = 0;
        let mut previous = start;
        let (mut lo, mut hi) = loop {
            if used == budget {
                return Ok(previous.point);
            }
            used += 1;
            let Some(sample) = self.sample(from, direction, alpha)? else {
                return Ok(previous.point);
            };
            if !armijo(&sample) || (previous.alpha > 0.0 && sample.phi >= previous.phi) {
                break (previous, sample);
            }
            if curvature(&sample) {
                return Ok(sample.point);
            }
            if sample.dphi >= 0.0 {
                break (sample, previous);
            }
            if alpha >= alpha_max {
                // The displacement bound is reached with sufficient decrease.
                return Ok(sample.point);
            }
            alpha = (2.0 * alpha).min(alpha_max);
            previous = sample;
        };
        // `lo` always satisfies sufficient decrease and has the lowest energy.
        while used < budget {
            used += 1;
            let trial = interpolate(&lo, &hi);
            let Some(sample) = self.sample(from, direction, trial)? else {
                break;
            };
            if !armijo(&sample) || sample.phi >= lo.phi {
                hi = sample;
            } else {
                if curvature(&sample) {
                    return Ok(sample.point);
                }
                if sample.dphi * (hi.alpha - lo.alpha) >= 0.0 {
                    hi = lo;
                }
                lo = sample;
            }
        }
        Ok(lo.point)
    }
}

/// Safeguarded cubic interpolation of the minimum between two samples, with
/// bisection when the cubic is undefined or leaves the inner interval.
fn interpolate(lo: &Sample, hi: &Sample) -> f64 {
    let (a, b) = (lo.alpha, hi.alpha);
    let midpoint = 0.5 * (a + b);
    let (left, right) = (a.min(b), a.max(b));
    let margin = 0.1 * (right - left);
    if !(hi.phi.is_finite() && hi.dphi.is_finite()) {
        return midpoint;
    }
    let d1 = lo.dphi + hi.dphi - 3.0 * (lo.phi - hi.phi) / (a - b);
    let radicand = d1 * d1 - lo.dphi * hi.dphi;
    if radicand < 0.0 {
        return midpoint;
    }
    let d2 = (b - a).signum() * radicand.sqrt();
    let candidate = b - (b - a) * (hi.dphi + d2 - d1) / (hi.dphi - lo.dphi + 2.0 * d2);
    if candidate.is_finite() && candidate >= left + margin && candidate <= right - margin {
        candidate
    } else {
        midpoint
    }
}

fn flatten_points(points: &[Point3]) -> Vec<f64> {
    points.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

fn flatten_vectors(vectors: &[Vector3]) -> Vec<f64> {
    vectors.iter().flat_map(|v| [v.x, v.y, v.z]).collect()
}

fn to_points(x: &[f64]) -> Vec<Point3> {
    x.as_chunks::<3>()
        .0
        .iter()
        .map(|&[x, y, z]| Point3::new(x, y, z))
        .collect()
}

fn max_atom_norm(x: &[f64]) -> f64 {
    x.as_chunks::<3>()
        .0
        .iter()
        .map(|&[x, y, z]| Vector3::new(x, y, z).norm())
        .fold(0.0, f64::max)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

fn axpy(y: &mut [f64], a: f64, x: &[f64]) {
    for (y, x) in y.iter_mut().zip(x) {
        *y += a * x;
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum MinimizationError {
    InvalidOptions(&'static str),
    /// Finite evaluations cannot produce a representable descent step.
    NumericalFailure(&'static str),
    /// The initial evaluation failed, or a trial failed for a reason other
    /// than singular geometry.
    Evaluation(EvaluationError),
    /// [`Minimization::to_model`] received a view of a different topology.
    IncompatibleTopology,
    Model(Box<ModelError>),
    Position(PositionError),
    Unit(UnitError),
}

impl fmt::Display for MinimizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => write!(f, "invalid minimization options: {message}"),
            Self::NumericalFailure(message) => {
                write!(f, "minimization numerical failure: {message}")
            }
            Self::Evaluation(error) => write!(f, "potential evaluation failed: {error}"),
            Self::IncompatibleTopology => f.write_str(
                "model view belongs to a different topology layout than the minimization",
            ),
            Self::Model(error) => write!(f, "cannot build the minimized realization: {error}"),
            Self::Position(error) => write!(f, "cannot update positions: {error}"),
            Self::Unit(error) => write!(f, "invalid minimization quantity unit: {error}"),
        }
    }
}

impl std::error::Error for MinimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Evaluation(error) => Some(error),
            Self::Model(error) => Some(error.as_ref()),
            Self::Position(error) => Some(error),
            Self::Unit(error) => Some(error),
            Self::InvalidOptions(_) | Self::NumericalFailure(_) | Self::IncompatibleTopology => {
                None
            }
        }
    }
}

impl From<EvaluationError> for MinimizationError {
    fn from(error: EvaluationError) -> Self {
        Self::Evaluation(error)
    }
}

impl From<ModelError> for MinimizationError {
    fn from(error: ModelError) -> Self {
        Self::Model(Box::new(error))
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
