//! L-BFGS minimization: convergence, observers, limits, line-search failures,
//! options, fixed atoms and position restraints.

use std::ops::ControlFlow;
use std::sync::Arc;

use kekule::geometry::{Point3, Vector3};
use kekule::structure::{Model, ModelView, Positions};
use kekule::topology::{AtomSelection, Topology};
use kekule::units::{
    Quantity, ANGSTROM, CANONICAL_ENERGY_UNIT, CANONICAL_GRADIENT_UNIT,
    KILOJOULE_PER_MOLE_PER_NANOMETER, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER, NANOMETER,
};
use kekule_potentials::{
    minimize, minimize_scoped, minimize_scoped_with_observer, minimize_with_observer, Energy,
    Evaluation, EvaluationError, InteractionKind, MinimizationError, MinimizationScope,
    MinimizationStatus, MinimizeOption, MinimizeOptions, Potential, RestraintProblem,
    SingularGeometry,
};

/// `sum_i sum_axis k_axis * scale * (x - center_i)^2`, optionally singular
/// beyond a wall at `x > wall` and optionally with a constant energy.
struct Bowl {
    topology: Arc<Topology>,
    centers: Vec<Point3>,
    stiffness: [f64; 3],
    scale: f64,
    wall: f64,
    flat_energy: bool,
}

impl Bowl {
    fn new(model: &Model, centers: Vec<Point3>) -> Self {
        Self {
            topology: model.shared_topology(),
            centers,
            stiffness: [1.0, 100.0, 10_000.0],
            scale: 1.0,
            wall: f64::INFINITY,
            flat_energy: false,
        }
    }
}

impl Potential for Bowl {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn evaluate(&self, model: ModelView<'_>) -> Result<Evaluation, EvaluationError> {
        if !self.topology.shares_layout(model.topology()) {
            return Err(EvaluationError::IncompatibleTopology);
        }
        let x = model.positions().values().into_value();
        let mut energy = 0.0;
        let mut gradient = Vec::with_capacity(x.len());
        for (index, (p, c)) in x.iter().zip(&self.centers).enumerate() {
            if p.x > self.wall {
                return Err(EvaluationError::InvalidGeometry {
                    interaction: InteractionKind::Restraint,
                    atoms: vec![model.topology().atom_ids()[index]],
                    kind: SingularGeometry::CoincidentAtoms,
                });
            }
            let d = [p.x - c.x, p.y - c.y, p.z - c.z];
            let k = self.stiffness.map(|k| k * self.scale);
            energy += (0..3).map(|i| k[i] * d[i] * d[i]).sum::<f64>();
            gradient.push(Vector3::new(
                2.0 * k[0] * d[0],
                2.0 * k[1] * d[1],
                2.0 * k[2] * d[2],
            ));
        }
        if self.flat_energy {
            energy = 1.0;
        }
        Evaluation::new(
            model,
            Energy::new(Quantity::new(energy, CANONICAL_ENERGY_UNIT))?,
            Quantity::new(gradient, CANONICAL_GRADIENT_UNIT),
        )
    }
}

fn model(points: &[Point3]) -> Model {
    Model::new(
        kekule::smiles::to_topology("CCO").unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn start() -> Model {
    model(&[
        Point3::new(0.3, -0.2, 0.1),
        Point3::new(-0.1, 0.4, 0.05),
        Point3::new(0.2, 0.2, -0.3),
    ])
}

fn centers() -> Vec<Point3> {
    vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.1, 0.2, 0.0),
        Point3::new(0.4, -0.1, 0.2),
    ]
}

fn tight(tolerance: f64) -> MinimizeOptions {
    MinimizeOptions::default()
        .with_gradient_tolerance(Quantity::new(tolerance, KILOJOULE_PER_MOLE_PER_NANOMETER))
        .unwrap()
}

#[test]
fn lbfgs_converges_on_an_ill_conditioned_bowl() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let result = minimize(&bowl, model.as_model_view(), &tight(1e-6)).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);
    // Steepest descent needs on the order of the condition number (1e4) steps.
    assert!(
        result.iterations() < 100,
        "{} iterations",
        result.iterations()
    );
    for (p, c) in result
        .positions()
        .values()
        .into_value()
        .iter()
        .zip(centers())
    {
        assert!((p.x - c.x).abs() < 1e-6 && (p.y - c.y).abs() < 1e-7 && (p.z - c.z).abs() < 1e-9);
    }
    assert!(result.final_evaluation().max_gradient_norm().into_value() <= 1e-6);
    assert!(result.evaluations() > result.iterations());
}

#[test]
fn every_accepted_step_is_observed_and_observers_can_stop() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let mut energies = Vec::new();
    let result = minimize_with_observer(&bowl, model.as_model_view(), &tight(1e-6), |step| {
        energies.push(step.energy().total().into_value());
        assert_eq!(step.iteration(), energies.len());
        assert_eq!(step.positions().len(), 3);
        if step.iteration() == 3 {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        }
    })
    .unwrap();
    assert_eq!(result.status(), MinimizationStatus::Stopped);
    assert_eq!(result.iterations(), 3);
    assert!(energies.windows(2).all(|w| w[1] < w[0]));
    assert_eq!(
        result.final_evaluation().energy().total().into_value(),
        energies[2]
    );
}

#[test]
fn iteration_limit_is_reported() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let options = tight(1e-9).with_max_iterations(2);
    let result = minimize(&bowl, model.as_model_view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::MaxIterations);
    assert_eq!(result.iterations(), 2);
    let converged = minimize(
        &bowl,
        model.as_model_view(),
        &tight(1e9).with_max_iterations(0),
    )
    .unwrap();
    assert_eq!(converged.status(), MinimizationStatus::Converged);
    assert_eq!(converged.positions(), model.positions());
}

#[test]
fn singular_trials_are_backtracked_but_other_failures_abort() {
    let model = start();
    let mut bowl = Bowl::new(&model, centers());
    // The minimum at x = 0.4 lies just inside a wall the long trial steps cross.
    bowl.wall = 0.401;
    let options = tight(1e-6)
        .with_max_step(Quantity::new(1.0, NANOMETER))
        .unwrap();
    let result = minimize(&bowl, model.as_model_view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);

    bowl.wall = 0.2;
    let error = minimize(&bowl, model.as_model_view(), &options).unwrap_err();
    assert!(matches!(
        error,
        MinimizationError::Evaluation(EvaluationError::InvalidGeometry { .. })
    ));

    let unrelated = start();
    let error = minimize(
        &Bowl::new(&unrelated, centers()),
        model.as_model_view(),
        &options,
    )
    .unwrap_err();
    assert_eq!(
        error,
        MinimizationError::Evaluation(EvaluationError::IncompatibleTopology)
    );
}

#[test]
fn gradients_inconsistent_with_the_energy_stop_the_line_search() {
    let model = start();
    let mut bowl = Bowl::new(&model, centers());
    bowl.flat_energy = true;
    let result = minimize(&bowl, model.as_model_view(), &MinimizeOptions::default()).unwrap();
    assert_eq!(result.status(), MinimizationStatus::LineSearchFailed);
    assert_eq!(result.iterations(), 0);
    assert_eq!(result.positions(), model.positions());
}

#[test]
fn displacements_that_round_away_stop_the_line_search() {
    let far = Point3::new(1e15, 1e15, 1e15);
    let model = model(&[far; 3]);
    let mut bowl = Bowl::new(&model, vec![Point3::new(1e15 + 1.0, 1e15, 1e15); 3]);
    bowl.stiffness = [1.0; 3];
    let options = MinimizeOptions::default()
        .with_max_step(Quantity::new(1e-6, NANOMETER))
        .unwrap();
    let result = minimize(&bowl, model.as_model_view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::LineSearchFailed);
    assert_eq!(result.positions(), model.positions());
}

#[test]
fn progress_is_independent_of_the_gradient_scale() {
    let model = start();
    let ratio = |scale: f64| {
        let mut bowl = Bowl::new(&model, centers());
        bowl.scale = scale;
        let options = tight(1e-300).with_max_iterations(3);
        let result = minimize(&bowl, model.as_model_view(), &options).unwrap();
        assert_eq!(result.iterations(), 3, "scale {scale}");
        let initial = result.initial_energy().total().into_value();
        result.final_evaluation().energy().total().into_value() / initial
    };
    let reference = ratio(1.0);
    assert!(reference < 1.0);
    for scale in [1e-100, 1e100] {
        let observed = ratio(scale);
        assert!(
            (observed - reference).abs() < 1e-9,
            "scale {scale}: {observed} vs {reference}"
        );
    }
}

#[test]
fn invalid_options_are_rejected_when_set() {
    let tolerance = |value: f64| {
        MinimizeOptions::default()
            .with_gradient_tolerance(Quantity::new(value, KILOJOULE_PER_MOLE_PER_NANOMETER))
            .unwrap_err()
    };
    for (error, option) in [
        (tolerance(0.0), MinimizeOption::GradientTolerance),
        (tolerance(f64::NAN), MinimizeOption::GradientTolerance),
        (
            MinimizeOptions::default()
                .with_max_step(Quantity::new(0.0, NANOMETER))
                .unwrap_err(),
            MinimizeOption::MaxStep,
        ),
        (
            MinimizeOptions::default().with_memory(0).unwrap_err(),
            MinimizeOption::Memory,
        ),
        (
            MinimizeOptions::default()
                .with_max_line_search_evaluations(0)
                .unwrap_err(),
            MinimizeOption::MaxLineSearchEvaluations,
        ),
    ] {
        assert_eq!(error, MinimizationError::InvalidOptions(option));
    }
    assert!(matches!(
        MinimizeOptions::default().with_max_step(Quantity::new(1.0, CANONICAL_ENERGY_UNIT)),
        Err(MinimizationError::Unit(_))
    ));
}

/// `sum_i |x_{i+1} - x_i - offset|^2`: neighbours pull on each other, so a
/// moving atom between fixed ones settles at their midpoint.
struct Chain {
    topology: Arc<Topology>,
    offset: Vector3,
}

impl Potential for Chain {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn evaluate(&self, model: ModelView<'_>) -> Result<Evaluation, EvaluationError> {
        if !self.topology.shares_layout(model.topology()) {
            return Err(EvaluationError::IncompatibleTopology);
        }
        let x = model.positions().values().into_value();
        let mut energy = 0.0;
        let mut gradient = vec![Vector3::zero(); x.len()];
        for i in 1..x.len() {
            let d = x[i] - x[i - 1] - self.offset;
            energy += d.norm_squared();
            gradient[i] += 2.0 * d;
            gradient[i - 1] -= 2.0 * d;
        }
        Evaluation::new(
            model,
            Energy::new(Quantity::new(energy, CANONICAL_ENERGY_UNIT))?,
            Quantity::new(gradient, CANONICAL_GRADIENT_UNIT),
        )
    }
}

fn moving(model: &Model, indices: &[usize]) -> MinimizationScope {
    let atoms = indices.iter().map(|&i| model.topology().atom_ids()[i]);
    MinimizationScope::moving(AtomSelection::from_atoms(model.shared_topology(), atoms).unwrap())
}

#[test]
fn fixed_atoms_keep_their_bits_and_pull_on_the_moving_atoms() {
    // Signed zeros and subnormals expose any arithmetic on fixed coordinates.
    let model = model(&[
        Point3::new(-0.0, 5e-324, 0.1),
        Point3::new(0.7, -0.4, 0.3),
        Point3::new(0.2, -0.0, -0.3),
    ]);
    let chain = Chain {
        topology: model.shared_topology(),
        offset: Vector3::new(0.1, 0.0, 0.0),
    };
    let scope = moving(&model, &[1]);
    let result = minimize_scoped(&chain, model.as_model_view(), &scope, &tight(1e-8)).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);
    let input = model.positions().values().into_value();
    let output = result.positions().values().into_value();
    for i in [0, 2] {
        let bits = |p: &Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
        assert_eq!(bits(&output[i]), bits(&input[i]));
    }
    // The fixed neighbours set the minimum: their midpoint.
    let midpoint = input[0] + (input[2] - input[0]) / 2.0;
    assert!((output[1] - midpoint).norm() < 1e-9, "{:?}", output[1]);
    assert!(result.max_gradient().into_value() <= 1e-8);
    // Fixed atoms still feel forces; the convergence measure ignores them.
    assert!(result.final_evaluation().max_gradient_norm().into_value() > 1e-3);
    assert_eq!(result.restraint_energy().into_value(), 0.0);
    // Moving every atom reaches the unconstrained minimum instead.
    let free = minimize(&chain, model.as_model_view(), &tight(1e-8)).unwrap();
    let all = minimize_scoped(
        &chain,
        model.as_model_view(),
        &MinimizationScope::all(model.shared_topology()),
        &tight(1e-8),
    )
    .unwrap();
    assert_eq!(free, all);
    assert_eq!(
        free.max_gradient(),
        free.final_evaluation().max_gradient_norm()
    );
    assert!(free.final_evaluation().energy().total().into_value() < 1e-14);
}

#[test]
fn moving_atoms_converge_in_an_ill_conditioned_bowl_with_others_fixed() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let scope = moving(&model, &[0, 2]);
    let result = minimize_scoped(&bowl, model.as_model_view(), &scope, &tight(1e-6)).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);
    let output = result.positions().values().into_value();
    assert_eq!(output[1], model.positions().values().into_value()[1]);
    for i in [0, 2] {
        assert!((output[i] - centers()[i]).norm() < 1e-6);
    }
}

#[test]
fn a_restraint_holds_an_atom_near_its_reference() {
    let model = start();
    let mut bowl = Bowl::new(&model, centers());
    bowl.stiffness = [1.0, 1.0, 1.0];
    let atom = model.topology().atom_ids()[2];
    let reference = Point3::new(0.0, 0.0, 0.0);
    // `k/2 |r - r_ref|^2` against the bowl's `|r - c|^2` balances at
    // `(2c + k r_ref) / (2 + k)`.
    for k in [10.0, 1.0e6] {
        let scope = MinimizationScope::all(model.shared_topology())
            .with_position_restraint(
                atom,
                Quantity::new(reference, NANOMETER),
                Quantity::new(k, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER),
            )
            .unwrap();
        let result = minimize_scoped(&bowl, model.as_model_view(), &scope, &tight(1e-6)).unwrap();
        assert_eq!(result.status(), MinimizationStatus::Converged);
        let o = Point3::origin();
        let expected = o + (2.0 * (centers()[2] - o) + k * (reference - o)) / (2.0 + k);
        let position = result.positions().values().into_value()[2];
        assert!((position - expected).norm() < 1e-7, "k {k}: {position:?}");
        let restraint = 0.5 * k * (position - reference).norm_squared();
        assert!((result.restraint_energy().into_value() - restraint).abs() < 1e-12);
        // The potential's own evaluation excludes the restraint.
        let minimized = result.to_model(model.as_model_view()).unwrap();
        assert_eq!(
            &bowl.evaluate(minimized.as_model_view()).unwrap(),
            result.final_evaluation()
        );
    }
    // A stiff restraint holds the atom within 1e-5 nm.
    let stiff = MinimizationScope::all(model.shared_topology())
        .with_position_restraint(
            atom,
            Quantity::new(Point3::new(1.0, 1.0, 1.0), ANGSTROM),
            Quantity::new(1.0e7, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER),
        )
        .unwrap();
    let result = minimize_scoped(&bowl, model.as_model_view(), &stiff, &tight(1e-4)).unwrap();
    let position = result.positions().values().into_value()[2];
    assert!((position - Point3::new(0.1, 0.1, 0.1)).norm() < 1e-5);
    // A restraint on a fixed atom adds a constant energy and no force.
    let fixed = moving(&model, &[0, 1])
        .with_position_restraint(
            atom,
            Quantity::new(reference, NANOMETER),
            Quantity::new(2.0, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER),
        )
        .unwrap();
    let result = minimize_scoped(&bowl, model.as_model_view(), &fixed, &tight(1e-6)).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);
    let start = model.positions().values().into_value()[2];
    assert_eq!(result.positions().values().into_value()[2], start);
    let constant = (start - reference).norm_squared();
    assert!((result.restraint_energy().into_value() - constant).abs() < 1e-15);
}

#[test]
fn scoped_observers_see_restraint_energies_and_can_stop() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let atom = model.topology().atom_ids()[0];
    let scope = moving(&model, &[0, 1])
        .with_position_restraint(
            atom,
            Quantity::new(Point3::new(0.3, -0.2, 0.1), NANOMETER),
            Quantity::new(100.0, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER),
        )
        .unwrap();
    let mut seen = Vec::new();
    let result =
        minimize_scoped_with_observer(&bowl, model.as_model_view(), &scope, &tight(1e-6), |step| {
            seen.push((step.restraint_energy().into_value(), step.max_gradient()));
            assert_eq!(
                step.positions().values().into_value()[2],
                model.positions().values().into_value()[2]
            );
            if step.iteration() == 2 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        })
        .unwrap();
    assert_eq!(result.status(), MinimizationStatus::Stopped);
    assert_eq!(result.iterations(), 2);
    let (restraint, max_gradient) = seen[1];
    assert!(restraint > 0.0);
    assert_eq!(result.restraint_energy().into_value(), restraint);
    assert_eq!(result.max_gradient(), max_gradient);
}

#[test]
fn an_empty_scope_converges_without_moving() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let scope = MinimizationScope::moving(AtomSelection::empty(model.shared_topology()));
    let result = minimize_scoped(&bowl, model.as_model_view(), &scope, &tight(1e-6)).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);
    assert_eq!(result.iterations(), 0);
    assert_eq!(result.evaluations(), 1);
    assert_eq!(result.positions(), model.positions());
    assert_eq!(result.max_gradient().into_value(), 0.0);
}

#[test]
fn invalid_scopes_and_restraints_are_rejected() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let other = self::model(&centers());
    let foreign = MinimizationScope::all(other.shared_topology());
    assert_eq!(
        minimize_scoped(&bowl, model.as_model_view(), &foreign, &tight(1e-6)).unwrap_err(),
        MinimizationError::IncompatibleTopology
    );

    let atom = model.topology().atom_ids()[0];
    let scope = || MinimizationScope::all(model.shared_topology());
    let origin = Quantity::new(Point3::origin(), NANOMETER);
    let k = |value: f64| Quantity::new(value, KILOJOULE_PER_MOLE_PER_SQUARE_NANOMETER);
    let water = kekule::smiles::to_topology("O.O").unwrap();
    let unknown = water.atom_ids()[1];
    let problem = |error: MinimizationError| match error {
        MinimizationError::InvalidRestraint { atom, problem } => (atom, problem),
        other => panic!("unexpected {other:?}"),
    };
    assert_eq!(
        problem(
            scope()
                .with_position_restraint(unknown, origin, k(1.0))
                .unwrap_err()
        ),
        (unknown, RestraintProblem::UnknownAtom)
    );
    let nan = Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), NANOMETER);
    assert_eq!(
        problem(
            scope()
                .with_position_restraint(atom, nan, k(1.0))
                .unwrap_err()
        ),
        (atom, RestraintProblem::NonFiniteReference)
    );
    for value in [-1.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            problem(
                scope()
                    .with_position_restraint(atom, origin, k(value))
                    .unwrap_err()
            ),
            (atom, RestraintProblem::ForceConstant)
        );
    }
    assert!(matches!(
        scope().with_position_restraint(atom, origin, Quantity::new(1.0, NANOMETER)),
        Err(MinimizationError::Unit(_))
    ));
    let error = scope()
        .with_position_restraint(atom, origin, k(-1.0))
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        format!(
            "invalid position restraint on atom {atom}: force constant must be finite and non-negative"
        )
    );
    // A zero force constant is allowed and changes nothing.
    let zero = scope()
        .with_position_restraint(atom, origin, k(0.0))
        .unwrap();
    assert_eq!(zero.position_restraints()[0].force_constant(), k(0.0));
    assert_eq!(zero.position_restraints()[0].atom(), atom);
    assert_eq!(
        minimize_scoped(&bowl, model.as_model_view(), &zero, &tight(1e-6)).unwrap(),
        minimize(&bowl, model.as_model_view(), &tight(1e-6)).unwrap()
    );
}
