//! Local geometry optimization over any [`Potential`].
//!
//! [`minimize`] runs limited-memory BFGS (L-BFGS) with a strong-Wolfe line
//! search on all Cartesian coordinates. It reads a borrowed [`ModelView`], so
//! models, ensemble samples, and trajectory frames are minimized the same way,
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

use kekule::geometry::{Point3, Vector3};
use kekule::structure::{DenseArrayError, Model, ModelView, Positions, StructureError};
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
    max_iterations: usize,
    gradient_tolerance: Quantity<f64>,
    max_step: Quantity<f64>,
    memory: usize,
    max_line_search_evaluations: usize,
}

impl MinimizeOptions {
    /// Maximum number of accepted coordinate updates.
    pub const fn max_iterations(&self) -> usize {
        self.max_iterations
    }

    /// Sets [`Self::max_iterations`].
    #[must_use]
    pub const fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Convergence threshold for the largest per-atom gradient norm.
    pub const fn gradient_tolerance(&self) -> Quantity<f64> {
        self.gradient_tolerance
    }

    /// Sets [`Self::gradient_tolerance`].
    ///
    /// # Errors
    ///
    /// [`MinimizationError::Unit`] for an incompatible unit and
    /// [`MinimizationError::InvalidOptions`] unless the value is finite and
    /// positive.
    pub fn with_gradient_tolerance(
        mut self,
        gradient_tolerance: Quantity<f64>,
    ) -> Result<Self, MinimizationError> {
        let value = gradient_tolerance.value_in(CANONICAL_GRADIENT_UNIT)?;
        if !(value.is_finite() && value > 0.0) {
            return Err(MinimizationError::InvalidOptions(
                MinimizeOption::GradientTolerance,
            ));
        }
        self.gradient_tolerance = gradient_tolerance;
        Ok(self)
    }

    /// Largest displacement of any atom in one line search.
    pub const fn max_step(&self) -> Quantity<f64> {
        self.max_step
    }

    /// Sets [`Self::max_step`].
    ///
    /// # Errors
    ///
    /// [`MinimizationError::Unit`] for an incompatible unit and
    /// [`MinimizationError::InvalidOptions`] unless the value is finite and
    /// positive.
    pub fn with_max_step(mut self, max_step: Quantity<f64>) -> Result<Self, MinimizationError> {
        let value = max_step.value_in(CANONICAL_LENGTH_UNIT)?;
        if !(value.is_finite() && value > 0.0) {
            return Err(MinimizationError::InvalidOptions(MinimizeOption::MaxStep));
        }
        self.max_step = max_step;
        Ok(self)
    }

    /// Number of curvature pairs retained by L-BFGS.
    pub const fn memory(&self) -> usize {
        self.memory
    }

    /// Sets [`Self::memory`].
    ///
    /// # Errors
    ///
    /// [`MinimizationError::InvalidOptions`] when the value is zero.
    pub fn with_memory(mut self, memory: usize) -> Result<Self, MinimizationError> {
        if memory == 0 {
            return Err(MinimizationError::InvalidOptions(MinimizeOption::Memory));
        }
        self.memory = memory;
        Ok(self)
    }

    /// Potential evaluations allowed in one line search.
    pub const fn max_line_search_evaluations(&self) -> usize {
        self.max_line_search_evaluations
    }

    /// Sets [`Self::max_line_search_evaluations`].
    ///
    /// # Errors
    ///
    /// [`MinimizationError::InvalidOptions`] when the value is zero.
    pub fn with_max_line_search_evaluations(
        mut self,
        max_line_search_evaluations: usize,
    ) -> Result<Self, MinimizationError> {
        if max_line_search_evaluations == 0 {
            return Err(MinimizationError::InvalidOptions(
                MinimizeOption::MaxLineSearchEvaluations,
            ));
        }
        self.max_line_search_evaluations = max_line_search_evaluations;
        Ok(self)
    }
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
#[non_exhaustive]
pub enum MinimizationStatus {
    /// The largest per-atom gradient norm reached the tolerance.
    Converged,
    /// [`MinimizeOptions::max_iterations`] steps were accepted without
    /// convergence.
    MaxIterations,
    /// No step along a steepest-descent direction reduced the energy.
    LineSearchFailed,
    /// The observer requested an early stop.
    Stopped,
}

/// State after one accepted step, passed to a [`minimize_with_observer`] callback.
#[derive(Debug, Clone, Copy)]
pub struct MinimizationStep<'a> {
    pub(crate) iteration: usize,
    pub(crate) energy: &'a Energy,
    pub(crate) max_gradient: Quantity<f64>,
    pub(crate) positions: &'a Positions,
}

impl<'a> MinimizationStep<'a> {
    /// Number of accepted steps so far, starting at one.
    pub const fn iteration(&self) -> usize {
        self.iteration
    }

    /// Energy at the accepted coordinates.
    pub const fn energy(&self) -> &'a Energy {
        self.energy
    }

    /// Largest per-atom gradient norm at the accepted coordinates, in
    /// kJ/(mol nm).
    pub const fn max_gradient(&self) -> Quantity<f64> {
        self.max_gradient
    }

    /// The accepted coordinates in dense topology atom order.
    pub const fn positions(&self) -> &'a Positions {
        self.positions
    }
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

    /// Takes the final coordinates.
    pub fn into_positions(self) -> Positions {
        self.positions
    }

    /// Energy at the input coordinates.
    pub fn initial_energy(&self) -> &Energy {
        &self.initial_energy
    }

    /// Energy and gradient at the final coordinates.
    pub fn final_evaluation(&self) -> &Evaluation {
        &self.evaluation
    }

    /// Number of accepted steps.
    pub fn iterations(&self) -> usize {
        self.iterations
    }

    /// Potential evaluations performed, including rejected line-search trials.
    pub fn evaluations(&self) -> usize {
        self.evaluations
    }

    /// Why the minimization stopped.
    pub fn status(&self) -> MinimizationStatus {
        self.status
    }

    /// Copies `source` with the minimized coordinates.
    ///
    /// `source` must share the minimized topology snapshot. Its cell and
    /// structure properties are retained unchanged.
    ///
    /// # Errors
    ///
    /// [`MinimizationError::IncompatibleTopology`] if `source` does not share
    /// the minimized topology layout, and [`MinimizationError::Position`] if
    /// the coordinates cannot be stored.
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
///
/// A run that stops without converging is still `Ok`; check
/// [`Minimization::status`].
///
/// # Errors
///
/// - [`MinimizationError::InvalidOptions`] or [`MinimizationError::Unit`] if
///   `options` are unusable.
/// - [`MinimizationError::Evaluation`] if the initial evaluation fails, or a
///   trial fails for a reason other than singular geometry.
/// - [`MinimizationError::NumericalFailure`] if finite evaluations cannot
///   produce a representable descent step.
/// - [`MinimizationError::Position`] or [`MinimizationError::Structure`] if
///   trial coordinates cannot be stored.
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
///
/// # Errors
///
/// As [`minimize`].
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
        trial: model.to_model(),
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
                NumericalFailure::MaxGradientNorm,
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
                    NumericalFailure::SteepestDescentSlope,
                ));
            }
        }
        let alpha_max = limits.max_step / max_atom_norm(&direction);
        if !(alpha_max.is_finite() && alpha_max > 0.0) {
            return Err(MinimizationError::NumericalFailure(
                NumericalFailure::StepLengthBound,
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
            MinimizeOption::GradientTolerance,
        ));
    }
    if !(max_step.is_finite() && max_step > 0.0) {
        return Err(MinimizationError::InvalidOptions(MinimizeOption::MaxStep));
    }
    if options.memory == 0 {
        return Err(MinimizationError::InvalidOptions(MinimizeOption::Memory));
    }
    if options.max_line_search_evaluations == 0 {
        return Err(MinimizationError::InvalidOptions(
            MinimizeOption::MaxLineSearchEvaluations,
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

struct Search<'a, P: ?Sized> {
    potential: &'a P,
    /// Scratch model carrying the topology, cell, and annotations of the input.
    trial: Model,
    evaluations: usize,
}

/// A line-search sample: `phi = E(x + alpha d)` and `dphi = g . d`.
struct Sample {
    alpha: f64,
    phi: f64,
    dphi: f64,
    point: Option<Point>,
}

impl<P: Potential + ?Sized> Search<'_, P> {
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
        self.trial.set_positions(positions.values())?;
        let view = self.trial.as_model_view();
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

/// A [`MinimizeOptions`] field with an unusable value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum MinimizeOption {
    /// [`MinimizeOptions::gradient_tolerance`] is not finite and positive.
    GradientTolerance,
    /// [`MinimizeOptions::max_step`] is not finite and positive.
    MaxStep,
    /// [`MinimizeOptions::memory`] is zero.
    Memory,
    /// [`MinimizeOptions::max_line_search_evaluations`] is zero.
    MaxLineSearchEvaluations,
}

impl fmt::Display for MinimizeOption {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::GradientTolerance => "gradient tolerance must be finite and positive",
            Self::MaxStep => "maximum step must be finite and positive",
            Self::Memory => "L-BFGS memory must be at least one",
            Self::MaxLineSearchEvaluations => "line searches need at least one evaluation",
        })
    }
}

/// The quantity that stopped being representable during a minimization.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum NumericalFailure {
    /// The maximum per-atom gradient norm is not finite.
    MaxGradientNorm,
    /// Even the steepest-descent slope is not finite and negative.
    SteepestDescentSlope,
    /// The step-length bound is not finite and positive.
    StepLengthBound,
}

impl fmt::Display for NumericalFailure {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::MaxGradientNorm => "maximum gradient norm is not representable",
            Self::SteepestDescentSlope => "steepest-descent slope is not finite and negative",
            Self::StepLengthBound => "step length bound is not representable",
        })
    }
}

/// A minimization that could not run or could not continue.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum MinimizationError {
    /// An option has an unusable value.
    InvalidOptions(MinimizeOption),
    /// Finite evaluations cannot produce a representable descent step.
    NumericalFailure(NumericalFailure),
    /// The initial evaluation failed, or a trial failed for a reason other
    /// than singular geometry.
    Evaluation(EvaluationError),
    /// [`Minimization::to_model`] received a view of a different topology.
    IncompatibleTopology,
    /// The minimized structure could not be built; see the source.
    Structure(StructureError),
    /// The positions could not be updated; see the source.
    Position(DenseArrayError),
    /// An option could not be converted to the canonical unit; see the
    /// source.
    Unit(UnitError),
}

impl fmt::Display for MinimizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => write!(f, "invalid minimization options: {message}"),
            Self::NumericalFailure(message) => {
                write!(f, "minimization numerical failure: {message}")
            }
            Self::Evaluation(_) => f.write_str("potential evaluation failed"),
            Self::IncompatibleTopology => f.write_str(
                "model view belongs to a different topology layout than the minimization",
            ),
            Self::Structure(_) => f.write_str("cannot build the minimized structure"),
            Self::Position(_) => f.write_str("cannot update positions"),
            Self::Unit(_) => f.write_str("invalid minimization quantity unit"),
        }
    }
}

impl std::error::Error for MinimizationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Evaluation(error) => Some(error),
            Self::Structure(error) => Some(error),
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

impl From<StructureError> for MinimizationError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<DenseArrayError> for MinimizationError {
    fn from(error: DenseArrayError) -> Self {
        Self::Position(error)
    }
}

impl From<UnitError> for MinimizationError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
