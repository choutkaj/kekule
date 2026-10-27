//! Local geometry optimization over any [`Potential`].
//!
//! [`minimize`] runs limited-memory BFGS (L-BFGS) with a strong-Wolfe line
//! search on all Cartesian coordinates. It reads a borrowed [`ModelView`], so
//! models, ensemble samples, and trajectory frames are minimized the same way,
//! and returns new [`Positions`] without changing its input. Periodic cells are
//! passed through to the potential unchanged.
//!
//! [`minimize_scoped`] moves only the atoms of a [`MinimizationScope`]. The
//! other atoms keep their input coordinates bit for bit but still enter every
//! potential evaluation, so they shape the energy and the forces on the moving
//! atoms. A scope may add harmonic position restraints,
//! `k/2 |r - r_ref|^2` per restrained atom, to the minimized energy; restraints
//! use absolute Cartesian coordinates without periodic imaging. Convergence and
//! the step bound consider moving atoms only.
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
use kekule::structure::{DenseArrayError, Model, ModelView, Positions, StructureError};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::units::{
    Quantity, UnitError, CANONICAL_ENERGY_UNIT, CANONICAL_FORCE_CONSTANT_UNIT,
    CANONICAL_GRADIENT_UNIT, CANONICAL_LENGTH_UNIT, KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER,
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

/// The atoms a minimization moves and the position restraints it adds.
///
/// Atoms outside [`Self::moving_atoms`] are fixed: they keep their input
/// coordinates exactly and still contribute to the energy and to the forces on
/// the moving atoms. A scope binds the topology snapshot of its selection and
/// applies to any view sharing that layout.
///
/// ```
/// use kekule::geometry::Point3;
/// use kekule::smiles;
/// use kekule::topology::AtomSelection;
/// use kekule::units::{Quantity, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER, NANOMETER};
/// use kekule_potentials::MinimizationScope;
///
/// let topology = smiles::to_topology("CCO")?;
/// let oxygen = topology.atom_ids()[2];
/// let scope = MinimizationScope::moving(AtomSelection::from_atoms(topology, [oxygen])?)
///     .with_position_restraint(
///         oxygen,
///         Quantity::new(Point3::new(0.1, 0.0, 0.0), NANOMETER),
///         Quantity::new(1_000.0, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER),
///     )?;
/// assert_eq!(scope.moving_atoms().len(), 1);
/// assert_eq!(scope.position_restraints().len(), 1);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone)]
pub struct MinimizationScope {
    moving: AtomSelection,
    restraints: Vec<PositionRestraint>,
    /// Dense atom index of each restraint, parallel to `restraints`.
    restraint_indices: Vec<usize>,
}

impl MinimizationScope {
    /// Moves every atom of `topology`, without restraints, as [`minimize`] does.
    pub fn all(topology: impl Into<Arc<Topology>>) -> Self {
        Self::moving(AtomSelection::all(topology))
    }

    /// Moves the selected atoms and fixes every other atom.
    ///
    /// An empty selection fixes every atom, so the minimization converges
    /// without a step.
    pub fn moving(selection: AtomSelection) -> Self {
        Self {
            moving: selection,
            restraints: Vec::new(),
            restraint_indices: Vec::new(),
        }
    }

    /// Adds a harmonic restraint `k/2 |r - reference|^2` on `atom`.
    ///
    /// Restraints on the same atom add up. A restraint on a fixed atom adds a
    /// constant energy and no force.
    ///
    /// # Errors
    ///
    /// [`MinimizationError::Unit`] for an incompatible unit, and
    /// [`MinimizationError::InvalidRestraint`] if `atom` is not in the scope's
    /// topology, `reference` is not finite, or `force_constant` is not finite
    /// and non-negative.
    pub fn with_position_restraint(
        mut self,
        atom: InstanceAtomId,
        reference: Quantity<Point3>,
        force_constant: Quantity<f64>,
    ) -> Result<Self, MinimizationError> {
        let reference = reference.value_in(CANONICAL_LENGTH_UNIT)?;
        let force_constant = force_constant.value_in(CANONICAL_FORCE_CONSTANT_UNIT)?;
        let invalid = |problem| MinimizationError::InvalidRestraint { atom, problem };
        let index = self
            .moving
            .topology()
            .atom_index(atom)
            .ok_or(invalid(RestraintProblem::UnknownAtom))?;
        if ![reference.x, reference.y, reference.z]
            .iter()
            .all(|v| v.is_finite())
        {
            return Err(invalid(RestraintProblem::NonFiniteReference));
        }
        if !(force_constant.is_finite() && force_constant >= 0.0) {
            return Err(invalid(RestraintProblem::ForceConstant));
        }
        self.restraints.push(PositionRestraint {
            atom,
            reference,
            force_constant,
        });
        self.restraint_indices.push(index.index());
        Ok(self)
    }

    /// The atoms that move; every other atom is fixed.
    pub fn moving_atoms(&self) -> &AtomSelection {
        &self.moving
    }

    /// The position restraints in the order added.
    pub fn position_restraints(&self) -> &[PositionRestraint] {
        &self.restraints
    }
}

/// A harmonic restraint `k/2 |r - reference|^2` holding one atom near a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PositionRestraint {
    atom: InstanceAtomId,
    reference: Point3,
    force_constant: f64,
}

impl PositionRestraint {
    /// The restrained atom.
    pub const fn atom(&self) -> InstanceAtomId {
        self.atom
    }

    /// The point the atom is held near, in nm.
    pub const fn reference(&self) -> Quantity<Point3> {
        Quantity::new(self.reference, CANONICAL_LENGTH_UNIT)
    }

    /// The force constant `k`, in kJ/(mol nm^2).
    pub const fn force_constant(&self) -> Quantity<f64> {
        Quantity::new(self.force_constant, CANONICAL_FORCE_CONSTANT_UNIT)
    }
}

/// Terminal state of a minimization that did not fail.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MinimizationStatus {
    /// The largest per-atom gradient norm over the moving atoms reached the
    /// tolerance.
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
    pub(crate) restraint_energy: f64,
    pub(crate) max_gradient: Quantity<f64>,
    pub(crate) positions: &'a Positions,
}

impl<'a> MinimizationStep<'a> {
    /// Number of accepted steps so far, starting at one.
    pub const fn iteration(&self) -> usize {
        self.iteration
    }

    /// Potential energy at the accepted coordinates, without restraints.
    pub const fn energy(&self) -> &'a Energy {
        self.energy
    }

    /// Position-restraint energy at the accepted coordinates, in kJ/mol; zero
    /// without restraints. The minimized energy is the potential energy plus
    /// this.
    pub const fn restraint_energy(&self) -> Quantity<f64> {
        Quantity::new(self.restraint_energy, CANONICAL_ENERGY_UNIT)
    }

    /// Largest per-atom norm, over the moving atoms, of the minimized energy's
    /// gradient at the accepted coordinates, in kJ/(mol nm).
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
    restraint_energy: f64,
    max_gradient: f64,
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

    /// Potential energy at the input coordinates.
    pub fn initial_energy(&self) -> &Energy {
        &self.initial_energy
    }

    /// Potential energy and gradient at the final coordinates, for every atom
    /// and without restraints.
    pub fn final_evaluation(&self) -> &Evaluation {
        &self.evaluation
    }

    /// Position-restraint energy at the final coordinates, in kJ/mol; zero
    /// without restraints.
    pub fn restraint_energy(&self) -> Quantity<f64> {
        Quantity::new(self.restraint_energy, CANONICAL_ENERGY_UNIT)
    }

    /// Largest per-atom norm, over the moving atoms, of the minimized energy's
    /// gradient at the final coordinates, in kJ/(mol nm): the convergence
    /// measure. Without fixed atoms or restraints it equals
    /// [`Evaluation::max_gradient_norm`] of [`Self::final_evaluation`].
    pub fn max_gradient(&self) -> Quantity<f64> {
        Quantity::new(self.max_gradient, CANONICAL_GRADIENT_UNIT)
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
    observer: impl FnMut(&MinimizationStep<'_>) -> ControlFlow<()>,
) -> Result<Minimization, MinimizationError> {
    let scope = MinimizationScope::all(model.shared_topology());
    minimize_scoped_with_observer(potential, model, &scope, options, observer)
}

/// Minimizes the energy of `model` under `potential` plus the scope's position
/// restraints, moving only the scope's atoms.
///
/// Fixed atoms keep their input coordinates bit for bit and still contribute
/// to every evaluation. A run that stops without converging is still `Ok`;
/// check [`Minimization::status`].
///
/// # Errors
///
/// As [`minimize`], and [`MinimizationError::IncompatibleTopology`] if
/// `scope` does not share the layout of `model`.
pub fn minimize_scoped<P: Potential + ?Sized>(
    potential: &P,
    model: ModelView<'_>,
    scope: &MinimizationScope,
    options: &MinimizeOptions,
) -> Result<Minimization, MinimizationError> {
    minimize_scoped_with_observer(potential, model, scope, options, |_| {
        ControlFlow::Continue(())
    })
}

/// Minimizes like [`minimize_scoped`], calling `observer` after every accepted
/// step as [`minimize_with_observer`] does.
///
/// # Errors
///
/// As [`minimize_scoped`].
pub fn minimize_scoped_with_observer<P: Potential + ?Sized>(
    potential: &P,
    model: ModelView<'_>,
    scope: &MinimizationScope,
    options: &MinimizeOptions,
    mut observer: impl FnMut(&MinimizationStep<'_>) -> ControlFlow<()>,
) -> Result<Minimization, MinimizationError> {
    let limits = validate(options)?;
    if !scope.moving.topology().shares_layout(model.topology()) {
        return Err(MinimizationError::IncompatibleTopology);
    }
    let objective = Objective::new(scope, model.positions().values().into_value());
    let initial = potential.evaluate(model)?;
    let mut search = Search {
        potential,
        objective: &objective,
        trial: model.to_model(),
        evaluations: 1,
    };
    let mut current = objective.point(model.positions().clone(), initial)?;
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
            restraint_energy: current.restraint_energy,
            max_gradient: Quantity::new(max_atom_norm(&current.gradient), CANONICAL_GRADIENT_UNIT),
            positions: &current.positions,
        };
        if observer(&step).is_break() {
            break MinimizationStatus::Stopped;
        }
    };

    Ok(Minimization {
        max_gradient: max_atom_norm(&current.gradient),
        restraint_energy: current.restraint_energy,
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

/// The minimized function: the potential plus position restraints, over the
/// flattened canonical coordinates of the moving atoms only.
struct Objective {
    /// Input coordinates; fixed atoms are always copied from here.
    base: Vec<Point3>,
    /// Dense indices of the moving atoms, ascending.
    moving: Vec<usize>,
    restraints: Vec<ActiveRestraint>,
}

struct ActiveRestraint {
    index: usize,
    /// Position of the atom in `Objective::moving`, or `None` when it is fixed.
    slot: Option<usize>,
    reference: Point3,
    force_constant: f64,
}

impl Objective {
    fn new(scope: &MinimizationScope, base: &[Point3]) -> Self {
        let moving: Vec<usize> = scope.moving.indices().iter().map(|i| i.index()).collect();
        let restraints = scope
            .restraints
            .iter()
            .zip(&scope.restraint_indices)
            .map(|(restraint, &index)| ActiveRestraint {
                index,
                slot: moving.binary_search(&index).ok(),
                reference: restraint.reference,
                force_constant: restraint.force_constant,
            })
            .collect();
        Self {
            base: base.to_vec(),
            moving,
            restraints,
        }
    }

    /// Full coordinates with the moving atoms at `x`.
    fn positions(&self, x: &[f64]) -> Vec<Point3> {
        let mut points = self.base.clone();
        for (&index, &[x, y, z]) in self.moving.iter().zip(x.as_chunks::<3>().0) {
            points[index] = Point3::new(x, y, z);
        }
        points
    }

    /// Builds a point from a potential evaluation at `positions`.
    fn point(
        &self,
        positions: Positions,
        evaluation: Evaluation,
    ) -> Result<Point, EvaluationError> {
        let points = positions.values().into_value();
        let full = evaluation.gradient().into_value();
        let mut x = Vec::with_capacity(3 * self.moving.len());
        let mut gradient = Vec::with_capacity(3 * self.moving.len());
        for &index in &self.moving {
            let (p, g) = (points[index], full[index]);
            x.extend([p.x, p.y, p.z]);
            gradient.extend([g.x, g.y, g.z]);
        }
        let mut restraint_energy = 0.0;
        for restraint in &self.restraints {
            let d = points[restraint.index] - restraint.reference;
            restraint_energy += 0.5 * restraint.force_constant * d.norm_squared();
            if let Some(slot) = restraint.slot {
                let g = d * restraint.force_constant;
                gradient[3 * slot] += g.x;
                gradient[3 * slot + 1] += g.y;
                gradient[3 * slot + 2] += g.z;
            }
        }
        let energy = evaluation.energy().total().into_value() + restraint_energy;
        if !energy.is_finite() {
            return Err(EvaluationError::NonFiniteEnergy);
        }
        Ok(Point {
            x,
            gradient,
            energy,
            restraint_energy,
            positions,
            evaluation,
        })
    }
}

/// One accepted or trial configuration. `x` and `gradient` hold the
/// flattened canonical coordinates and minimized-energy gradient of the moving
/// atoms; `positions` and `evaluation` cover every atom.
struct Point {
    x: Vec<f64>,
    gradient: Vec<f64>,
    energy: f64,
    restraint_energy: f64,
    positions: Positions,
    evaluation: Evaluation,
}

impl Point {
    fn energy(&self) -> f64 {
        self.energy
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
    objective: &'a Objective,
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
        let positions = Positions::from_vec(Quantity::new(
            self.objective.positions(&x),
            CANONICAL_LENGTH_UNIT,
        ))?;
        self.trial.set_positions(positions.values())?;
        let view = self.trial.as_model_view();
        self.evaluations += 1;
        match self
            .potential
            .evaluate(view)
            .and_then(|evaluation| self.objective.point(positions, evaluation))
        {
            Ok(point) => Ok(Some(Sample {
                alpha,
                phi: point.energy(),
                dphi: dot(&point.gradient, direction),
                point: Some(point),
            })),
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

/// Why [`MinimizationScope::with_position_restraint`] rejected a restraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum RestraintProblem {
    /// The atom is not in the scope's topology.
    UnknownAtom,
    /// The reference position is not finite.
    NonFiniteReference,
    /// The force constant is not finite and non-negative.
    ForceConstant,
}

impl fmt::Display for RestraintProblem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::UnknownAtom => "atom is not in the topology",
            Self::NonFiniteReference => "reference position must be finite",
            Self::ForceConstant => "force constant must be finite and non-negative",
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
    /// A position restraint is unusable.
    InvalidRestraint {
        /// The restrained atom.
        atom: InstanceAtomId,
        /// What is wrong with the restraint.
        problem: RestraintProblem,
    },
    /// Finite evaluations cannot produce a representable descent step.
    NumericalFailure(NumericalFailure),
    /// The initial evaluation failed, or a trial failed for a reason other
    /// than singular geometry.
    Evaluation(EvaluationError),
    /// A [`MinimizationScope`] or the view given to [`Minimization::to_model`]
    /// belongs to a different topology layout.
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
            Self::InvalidRestraint { atom, problem } => {
                write!(f, "invalid position restraint on atom {atom}: {problem}")
            }
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
            Self::InvalidOptions(_)
            | Self::InvalidRestraint { .. }
            | Self::NumericalFailure(_)
            | Self::IncompatibleTopology => None,
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
