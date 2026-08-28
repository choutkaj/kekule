use super::*;
use crate::core::{Atom, AtomId, BondId, BondOrder, Element, Molecule};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::modeling::potential::{
    HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError, PotentialEvaluation,
    PotentialGeometryError,
};
use crate::structure::{Ensemble, Model, ModelBuildError, ModelView, PositionError, Positions};
use crate::topology::AtomSiteMetadata;
use crate::topology::{InstanceAtomId, InstanceBondId, MoleculeInstanceId, TopologyBuildError};
use crate::units::{
    Quantity, ANGSTROM, CANONICAL_ENERGY_UNIT, CANONICAL_GRADIENT_UNIT, CANONICAL_LENGTH_UNIT,
    KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM, NANOMETER,
};

fn two_atom_small(distance: f64) -> (Molecule, Positions, AtomId, AtomId, BondId) {
    let mut graph = crate::core::MoleculeEditor::new();
    let carbon = Element::from_symbol("C").unwrap();
    let a = graph
        .add_atom(Atom::new(carbon))
        .expect("atom identifier capacity");
    let tombstone = graph
        .add_atom(Atom::new(carbon))
        .expect("atom identifier capacity");
    graph.delete_atom(tombstone).unwrap();
    let b = graph
        .add_atom(Atom::new(carbon))
        .expect("atom identifier capacity");
    let bond = graph.add_bond(a, b, BondOrder::Single).unwrap();
    let positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(distance, 0.0, 0.0)],
        ANGSTROM,
    ))
    .unwrap();
    (
        graph.finish().expect("connected test molecule"),
        positions,
        a,
        b,
        bond,
    )
}

fn one_atom_macro() -> (Molecule, Positions, AtomId) {
    let mut graph = crate::core::MoleculeEditor::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol("N").unwrap()))
        .expect("atom identifier capacity");
    let positions = Positions::new(Quantity::new([Point3::new(2.0, 0.0, 0.0)], ANGSTROM)).unwrap();
    (graph.finish().unwrap(), positions, atom)
}

#[test]
fn model_preserves_local_ids_and_dense_round_trips() {
    let (small, positions, a, b, _) = two_atom_small(1.5);
    let mut builder = Model::builder();
    let instance = builder.add_molecule(&small, &positions).unwrap();
    let model = builder.build().unwrap();
    let qa = InstanceAtomId::new(instance, a);
    let qb = InstanceAtomId::new(instance, b);
    assert_eq!(model.topology().atom_ids(), &[qa, qb]);
    assert_eq!(
        model
            .topology()
            .atom_id(model.topology().atom_index(qb).unwrap()),
        Some(qb)
    );
    assert_eq!(
        model.position(qb).unwrap(),
        Quantity::new(Point3::new(0.15, 0.0, 0.0), CANONICAL_LENGTH_UNIT)
    );
    assert!(model
        .topology()
        .atom(InstanceAtomId::new(instance, AtomId::new(1)))
        .is_err());
}

#[test]
fn positions_convert_source_units_before_model_construction() {
    let mut graph = crate::core::MoleculeEditor::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let source = [Point3::new(0.15, 0.0, 0.0)];
    let positions = Positions::new(Quantity::new(source, NANOMETER)).unwrap();
    let small = graph.finish().unwrap();

    let model = Model::from_molecule(&small, &positions).unwrap();
    let qualified = InstanceAtomId::new(MoleculeInstanceId::new(0), atom);
    assert_eq!(
        model.position(qualified).unwrap().unit(),
        CANONICAL_LENGTH_UNIT
    );
    assert_eq!(model.position(qualified).unwrap().x, 0.15);
    assert_eq!(positions.values().unit(), CANONICAL_LENGTH_UNIT);
    assert_eq!(source[0].x, 0.15);
}

#[test]
fn topology_allocation_is_shared_only_by_model_clones() {
    let (small, positions, _, b, _) = two_atom_small(1.5);
    let model = Model::from_molecule(&small, &positions).unwrap();
    let mut cloned = model.clone();
    cloned
        .set_position(
            InstanceAtomId::new(MoleculeInstanceId::new(0), b),
            Quantity::new(Point3::new(2.0, 0.0, 0.0), ANGSTROM),
        )
        .unwrap();
    let rebuilt = Model::from_molecule(&small, &positions).unwrap();

    assert!(std::sync::Arc::ptr_eq(
        &model.shared_topology(),
        &cloned.shared_topology()
    ));
    assert!(!std::sync::Arc::ptr_eq(
        &model.shared_topology(),
        &rebuilt.shared_topology()
    ));
    assert!(model.topology().same_layout(rebuilt.topology()));
    assert_ne!(model, cloned);
    assert_ne!(model, rebuilt);
}

#[test]
fn mixed_instances_and_hierarchy_use_qualified_ids() {
    let (small, small_positions, _, _, _) = two_atom_small(1.0);
    let (macromolecule, macro_positions, atom) = one_atom_macro();
    let mut builder = Model::builder();
    let small_id = builder.add_molecule(&small, &small_positions).unwrap();
    let macro_id = builder
        .add_molecule(&macromolecule, &macro_positions)
        .unwrap();
    let chain = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_chain("A", None)
        .unwrap();
    let residue = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    let site = builder
        .topology_builder_mut()
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(macro_id, atom),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    let model = builder.build().unwrap();
    assert_ne!(small_id, macro_id);
    assert!(std::ptr::eq(
        model.topology().molecule(small_id).unwrap().molecule(),
        model
            .topology()
            .definition_for_instance(small_id)
            .unwrap()
            .molecule()
    ));
    let hierarchy = model.topology();
    assert_eq!(
        hierarchy.atom_for_site(site).unwrap(),
        InstanceAtomId::new(macro_id, atom)
    );
}

#[test]
fn repeated_molecules_get_distinct_instance_ids() {
    let (small, positions, atom, _, _) = two_atom_small(1.0);
    let mut builder = Model::builder();
    let first = builder.add_molecule(&small, &positions).unwrap();
    let second = builder.add_molecule(&small, &positions).unwrap();
    let model = builder.build().unwrap();
    assert_ne!(first, second);
    assert_ne!(
        InstanceAtomId::new(first, atom),
        InstanceAtomId::new(second, atom)
    );
    assert_eq!(model.topology().instance_count(), 2);
}

#[test]
fn construction_copies_positions_and_preserves_sources() {
    let (small, positions, _a, _, _) = two_atom_small(1.0);
    let source = small.clone();
    let mut model = Model::from_molecule(&small, &positions).unwrap();
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(0));
    model
        .set_position(atom, Quantity::new(Point3::new(3.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    assert_eq!(small, source);
    assert_eq!(
        positions.position_at(0).unwrap(),
        Quantity::new(Point3::new(0.0, 0.0, 0.0), CANONICAL_LENGTH_UNIT)
    );
    assert_eq!(model.topology().definition_count(), 1);
}

#[test]
fn construction_rejects_empty_and_mismatched_positions_transactionally() {
    assert!(matches!(
        Model::builder().build(),
        Err(ModelBuildError::Topology(
            TopologyBuildError::NoMoleculeInstances
        ))
    ));
    let (small, _positions, _, _, _) = two_atom_small(1.0);
    let wrong_count = Positions::new(Quantity::new([Point3::origin()], ANGSTROM)).unwrap();
    let mut builder = Model::builder();
    assert!(matches!(
        builder.add_molecule(&small, &wrong_count),
        Err(ModelBuildError::InstancePositionCountMismatch {
            expected: 2,
            actual: 1
        })
    ));
    assert!(matches!(
        builder.build(),
        Err(ModelBuildError::Topology(
            TopologyBuildError::NoMoleculeInstances
        ))
    ));
}

#[test]
fn position_updates_are_complete_finite_and_transactional() {
    let (small, positions, _a, _, _) = two_atom_small(1.0);
    let mut model = Model::from_molecule(&small, &positions).unwrap();
    let original = model.positions().values().value().to_vec();
    assert!(matches!(
        model.set_positions(Quantity::new(&[Point3::default()], ANGSTROM)),
        Err(PositionError::PositionCountMismatch { .. })
    ));
    assert_eq!(model.positions().values().to_value(), original.as_slice());
    let mut invalid = original.clone();
    invalid[0] = Point3::new(f64::INFINITY, 0.0, 0.0);
    assert!(matches!(
        model.set_positions(Quantity::new(&invalid, ANGSTROM)),
        Err(PositionError::NonFinitePosition { index: 0 })
    ));
    assert_eq!(model.positions().values().to_value(), original.as_slice());
}

#[test]
fn harmonic_potential_and_minimization_use_instance_qualified_topology() {
    let (small, positions, _, _, bond) = two_atom_small(2.0);
    let model = Model::from_molecule(&small, &positions).unwrap();
    let qualified = InstanceBondId::new(MoleculeInstanceId::new(0), bond);
    let mut potential = HarmonicBondPotential::new(
        &model.shared_topology(),
        [HarmonicBondParameter::new(
            qualified,
            Quantity::new(1.0, ANGSTROM),
            Quantity::new(100.0, KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM),
        )],
    )
    .unwrap();
    let initial = potential.evaluate(model.view()).unwrap();
    assert!((initial.energy().to_value() - 50.0).abs() < 1.0e-10);
    let result = minimize(&model, &mut potential, MinimizeOptions::default()).unwrap();
    assert!(result.final_energy < result.initial_energy);
    assert!((model.positions().values().value()[1].x - 0.2).abs() < 1.0e-15);

    let rebuilt = Model::from_molecule(&small, &positions).unwrap();
    assert_eq!(
        potential.evaluate(rebuilt.view()),
        Err(PotentialError::IncompatibleTopology)
    );

    let mut coincident = model.clone();
    let instance = MoleculeInstanceId::new(0);
    coincident
        .set_position(
            InstanceAtomId::new(instance, AtomId::new(2)),
            Quantity::new(
                coincident.positions().values().value()[0],
                CANONICAL_LENGTH_UNIT,
            ),
        )
        .unwrap();
    assert_eq!(
        potential.evaluate(coincident.view()),
        Err(PotentialError::InvalidGeometry {
            interaction: "harmonic bond",
            atoms: vec![
                InstanceAtomId::new(instance, AtomId::new(0)),
                InstanceAtomId::new(instance, AtomId::new(2)),
            ],
            kind: PotentialGeometryError::CoincidentAtoms,
        })
    );
}

#[test]
fn harmonic_potential_rejects_periodic_model_and_ensemble_state() {
    let (small, positions, _, _, bond) = two_atom_small(9.8);
    let mut model = Model::from_molecule(&small, &positions).unwrap();
    model
        .set_positions(Quantity::new(
            [Point3::new(0.1, 0.0, 0.0), Point3::new(9.9, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap();
    let qualified = InstanceBondId::new(MoleculeInstanceId::new(0), bond);
    let mut potential = HarmonicBondPotential::new(
        &model.shared_topology(),
        [HarmonicBondParameter::new(
            qualified,
            Quantity::new(1.0, ANGSTROM),
            Quantity::new(100.0, KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM),
        )],
    )
    .unwrap();

    assert!(potential.evaluate(model.view()).is_ok());
    let nonperiodic_ensemble = Ensemble::from_models(&[model.clone()]).unwrap();
    assert!(potential
        .evaluate(nonperiodic_ensemble.member(0).unwrap().as_model())
        .is_ok());
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    let mut periodic_model = model.clone();
    periodic_model.set_cell(Some(cell));
    assert_eq!(
        potential.evaluate(periodic_model.view()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );

    let periodic_ensemble = Ensemble::from_models(&[periodic_model.clone()]).unwrap();
    assert_eq!(
        potential.evaluate(periodic_ensemble.member(0).unwrap().as_model()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );
    let mut independent = Model::from_molecule(&small, &positions).unwrap();
    independent.set_cell(Some(cell));
    assert_eq!(
        potential.evaluate(independent.view()),
        Err(PotentialError::IncompatibleTopology)
    );
}

struct RecoverableGeometryPotential;

impl Potential for RecoverableGeometryPotential {
    fn evaluate(&mut self, model: ModelView<'_>) -> Result<PotentialEvaluation, PotentialError> {
        let coordinate = model.positions().values().value()[1].x;
        if coordinate <= 0.025 {
            return Err(PotentialError::invalid_geometry(
                "test coordinate",
                [model.topology().atom_ids()[1]],
                PotentialGeometryError::CoincidentAtoms,
            ));
        }
        PotentialEvaluation::new(
            model,
            Quantity::new(0.5 * coordinate * coordinate, CANONICAL_ENERGY_UNIT),
            Quantity::new(
                vec![Vector3::zero(), Vector3::new(coordinate, 0.0, 0.0)],
                CANONICAL_GRADIENT_UNIT,
            ),
        )
    }
}

struct BackendFailurePotential {
    calls: usize,
}

impl Potential for BackendFailurePotential {
    fn evaluate(&mut self, model: ModelView<'_>) -> Result<PotentialEvaluation, PotentialError> {
        self.calls += 1;
        if self.calls > 1 {
            return Err(PotentialError::backend("test backend", "evaluation failed"));
        }
        PotentialEvaluation::new(
            model,
            Quantity::new(0.5, CANONICAL_ENERGY_UNIT),
            Quantity::new(
                vec![Vector3::zero(), Vector3::new(1.0, 0.0, 0.0)],
                CANONICAL_GRADIENT_UNIT,
            ),
        )
    }
}

#[test]
fn minimization_backtracks_invalid_geometry_but_propagates_backend_failures() {
    let (small, positions, _, _, _) = two_atom_small(1.0);
    let model = Model::from_molecule(&small, &positions).unwrap();
    let options = MinimizeOptions {
        max_iterations: 1,
        initial_step: Quantity::new(1.0, ANGSTROM),
        ..MinimizeOptions::default()
    };

    let result = minimize(&model, &mut RecoverableGeometryPotential, options).unwrap();
    assert_eq!(result.status, MinimizationStatus::MaxIterations);
    assert_eq!(result.iterations, 1);
    assert_eq!(result.evaluations, 3);
    assert!((result.model.positions().values().value()[1].x - 0.05).abs() < 1.0e-15);
    assert!((model.positions().values().value()[1].x - 0.1).abs() < 1.0e-15);

    let error = minimize(&model, &mut BackendFailurePotential { calls: 0 }, options).unwrap_err();
    assert!(matches!(
        error,
        MinimizationError::Potential(PotentialError::Backend {
            backend: "test backend",
            ..
        })
    ));
}
