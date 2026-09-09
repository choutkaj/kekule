use kekule::geometry::{Point3, Vector3};
use kekule::modeling::potential::{
    HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError, PotentialEvaluation,
};
use kekule::modeling::{minimize, MinimizationError, MinimizationStatus, MinimizeOptions};
use kekule::structure::{Model, ModelView, Positions};
use kekule::units::{
    Quantity, CANONICAL_ENERGY_UNIT, CANONICAL_FORCE_CONSTANT_UNIT, CANONICAL_GRADIENT_UNIT,
    NANOMETER,
};

fn model(smiles: &str, positions: &[Point3]) -> Model {
    Model::new(
        kekule::smiles::to_topology(smiles).unwrap(),
        Positions::new(Quantity::new(positions, NANOMETER)).unwrap(),
    )
    .unwrap()
}

#[test]
fn vector_norm_preserves_representable_extreme_magnitudes() {
    for scale in [1e-200, 1.0, 1e200] {
        let norm = Vector3::new(3.0 * scale, -4.0 * scale, 0.0).norm();
        assert!((norm / (5.0 * scale) - 1.0).abs() < 1e-15);
    }
    for magnitude in [f64::from_bits(1), f64::MIN_POSITIVE, f64::MAX] {
        assert_eq!(Vector3::new(0.0, 0.0, magnitude).norm(), magnitude);
    }
    assert_eq!(Vector3::zero().norm(), 0.0);
}

#[test]
fn minimization_makes_progress_at_extreme_finite_gradient_scales() {
    let source = model("CC", &[Point3::origin(), Point3::new(2.0, 0.0, 0.0)]);
    let initial_positions = source.positions().clone();
    for (stiffness, tolerance) in [(1e-200, 1e-250), (1.0, 1e-4), (1e200, 1e-4)] {
        let topology = source.shared_topology();
        let mut potential = HarmonicBondPotential::new(
            &topology,
            [HarmonicBondParameter::new(
                topology.bond_ids()[0],
                Quantity::new(1.0, NANOMETER),
                Quantity::new(stiffness, CANONICAL_FORCE_CONSTANT_UNIT),
            )],
        )
        .unwrap();
        let result = minimize(
            &source,
            &mut potential,
            MinimizeOptions {
                max_iterations: 1,
                gradient_tolerance: Quantity::new(tolerance, CANONICAL_GRADIENT_UNIT),
                ..MinimizeOptions::default()
            },
        )
        .unwrap();

        assert_eq!(result.status, MinimizationStatus::MaxIterations);
        assert_eq!(result.iterations, 1);
        assert!(result.final_energy < result.initial_energy);
        assert!((result.final_max_gradient.into_value() / stiffness - 0.98).abs() < 1e-14);
        assert_ne!(result.model.positions(), source.positions());
        assert_eq!(source.positions(), &initial_positions);
    }
}

struct LinearPotential(Vector3);

impl Potential for LinearPotential {
    fn evaluate(&mut self, model: ModelView<'_>) -> Result<PotentialEvaluation, PotentialError> {
        let energy = model
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| self.0.dot(Vector3::new(point.x, point.y, point.z)))
            .sum();
        PotentialEvaluation::new(
            model,
            Quantity::new(energy, CANONICAL_ENERGY_UNIT),
            Quantity::new(vec![self.0; model.atom_count()], CANONICAL_GRADIENT_UNIT),
        )
    }
}

#[test]
fn minimization_rejects_unrepresentable_gradient_norm_and_derivative() {
    for (source, slope, diagnostic) in [
        (
            model("C", &[Point3::origin()]),
            Vector3::new(f64::MAX, f64::MAX, 0.0),
            "gradient norm",
        ),
        (
            model("CC", &[Point3::origin(), Point3::origin()]),
            Vector3::new(1e308, 0.0, 0.0),
            "directional derivative",
        ),
    ] {
        let original = source.clone();
        let error = minimize(
            &source,
            &mut LinearPotential(slope),
            MinimizeOptions {
                max_iterations: 1,
                ..MinimizeOptions::default()
            },
        )
        .unwrap_err();
        assert!(matches!(error, MinimizationError::NumericalFailure(_)));
        assert!(error.to_string().contains(diagnostic), "{error}");
        assert_eq!(source, original);
    }
}

#[test]
fn minimization_stalls_when_displacements_round_to_unchanged_positions() {
    let source = model("C", &[Point3::new(1e100, 0.0, 0.0)]);
    let result = minimize(
        &source,
        &mut LinearPotential(Vector3::new(1.0, 0.0, 0.0)),
        MinimizeOptions {
            max_iterations: 1,
            ..MinimizeOptions::default()
        },
    )
    .unwrap();
    assert_eq!(result.status, MinimizationStatus::LineSearchStalled);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.evaluations, 1);
    assert_eq!(result.model, source);
}
