use std::ops::ControlFlow;
use std::sync::Arc;

use kekule::geometry::{Point3, Vector3};
use kekule::structure::{Model, ModelView, Positions};
use kekule::topology::Topology;
use kekule::units::{
    Quantity, CANONICAL_ENERGY_UNIT, CANONICAL_GRADIENT_UNIT, KILOJOULE_PER_MOLE_PER_NANOMETER,
    NANOMETER,
};
use kekule_potentials::{
    minimize, minimize_with_observer, Energy, Evaluation, EvaluationError, MinimizationError,
    MinimizationStatus, MinimizeOptions, Potential, SingularGeometry,
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
    fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }

    fn evaluate(&self, model: ModelView<'_>) -> Result<Evaluation, EvaluationError> {
        if !std::ptr::eq(self.topology.as_ref(), model.topology()) {
            return Err(EvaluationError::IncompatibleTopology);
        }
        let x = model.positions().values().into_value();
        let mut energy = 0.0;
        let mut gradient = Vec::with_capacity(x.len());
        for (index, (p, c)) in x.iter().zip(&self.centers).enumerate() {
            if p.x > self.wall {
                return Err(EvaluationError::InvalidGeometry {
                    interaction: "wall",
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
    MinimizeOptions {
        gradient_tolerance: Quantity::new(tolerance, KILOJOULE_PER_MOLE_PER_NANOMETER),
        ..MinimizeOptions::default()
    }
}

#[test]
fn lbfgs_converges_on_an_ill_conditioned_bowl() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    let result = minimize(&bowl, model.view(), &tight(1e-6)).unwrap();
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
    let result = minimize_with_observer(&bowl, model.view(), &tight(1e-6), |step| {
        energies.push(step.energy.total().into_value());
        assert_eq!(step.iteration, energies.len());
        assert_eq!(step.positions.len(), 3);
        if step.iteration == 3 {
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
    let options = MinimizeOptions {
        max_iterations: 2,
        ..tight(1e-9)
    };
    let result = minimize(&bowl, model.view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::MaxIterations);
    assert_eq!(result.iterations(), 2);
    let converged = minimize(
        &bowl,
        model.view(),
        &MinimizeOptions {
            max_iterations: 0,
            ..tight(1e9)
        },
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
    let options = MinimizeOptions {
        max_step: Quantity::new(1.0, NANOMETER),
        ..tight(1e-6)
    };
    let result = minimize(&bowl, model.view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::Converged);

    bowl.wall = 0.2;
    let error = minimize(&bowl, model.view(), &options).unwrap_err();
    assert!(matches!(
        error,
        MinimizationError::Evaluation(EvaluationError::InvalidGeometry { .. })
    ));

    let unrelated = start();
    let error = minimize(&Bowl::new(&unrelated, centers()), model.view(), &options).unwrap_err();
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
    let result = minimize(&bowl, model.view(), &MinimizeOptions::default()).unwrap();
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
    let options = MinimizeOptions {
        max_step: Quantity::new(1e-6, NANOMETER),
        ..MinimizeOptions::default()
    };
    let result = minimize(&bowl, model.view(), &options).unwrap();
    assert_eq!(result.status(), MinimizationStatus::LineSearchFailed);
    assert_eq!(result.positions(), model.positions());
}

#[test]
fn progress_is_independent_of_the_gradient_scale() {
    let model = start();
    let ratio = |scale: f64| {
        let mut bowl = Bowl::new(&model, centers());
        bowl.scale = scale;
        let options = MinimizeOptions {
            max_iterations: 3,
            ..tight(1e-300)
        };
        let result = minimize(&bowl, model.view(), &options).unwrap();
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
fn invalid_options_are_rejected_before_evaluation() {
    let model = start();
    let bowl = Bowl::new(&model, centers());
    for options in [
        tight(0.0),
        tight(f64::NAN),
        MinimizeOptions {
            max_step: Quantity::new(0.0, NANOMETER),
            ..MinimizeOptions::default()
        },
        MinimizeOptions {
            memory: 0,
            ..MinimizeOptions::default()
        },
        MinimizeOptions {
            max_line_search_evaluations: 0,
            ..MinimizeOptions::default()
        },
    ] {
        assert!(matches!(
            minimize(&bowl, model.view(), &options),
            Err(MinimizationError::InvalidOptions(_))
        ));
    }
    let wrong_unit = MinimizeOptions {
        max_step: Quantity::new(1.0, CANONICAL_ENERGY_UNIT),
        ..MinimizeOptions::default()
    };
    assert!(matches!(
        minimize(&bowl, model.view(), &wrong_unit),
        Err(MinimizationError::Unit(_))
    ));
}
