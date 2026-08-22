use super::*;
use crate::bio::{Hierarchy, SmcraAtomSiteId, SmcraAtomSiteMetadata};
use crate::core::{Atom, AtomId, BondId, BondOrder, Conformer, Element, Molecule};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::modeling::potential::{
    HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError, PotentialEvaluation,
    PotentialGeometryError,
};
use crate::structure::{Ensemble, Model, ModelBuildError, ModelView, PositionError};
use crate::topology::{
    InstanceAtomId, InstanceAtomSiteId, InstanceBondId, MoleculeInstanceId, TopologyBuildError,
};
use crate::units::{
    Quantity, ANGSTROM, MODEL_ENERGY_UNIT, MODEL_FORCE_CONSTANT_UNIT, MODEL_GRADIENT_UNIT,
    NANOMETER,
};

fn two_atom_small(distance: f64) -> (Molecule, Conformer, AtomId, AtomId, BondId) {
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
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            a,
            crate::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            b,
            crate::units::Quantity::new(Point3::new(distance, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    (
        graph.finish().expect("connected test molecule"),
        conformer,
        a,
        b,
        bond,
    )
}

fn one_atom_macro() -> (Molecule, Conformer, AtomId, SmcraAtomSiteId) {
    let mut graph = crate::core::MoleculeEditor::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol("N").unwrap()))
        .expect("atom identifier capacity");
    let mut conformer = Conformer::new(crate::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            atom,
            crate::units::Quantity::new(Point3::new(2.0, 0.0, 0.0), crate::units::ANGSTROM),
        )
        .unwrap();
    let mut hierarchy = Hierarchy::new();
    let chain = hierarchy.add_chain("A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    let site = hierarchy
        .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
        .unwrap();
    *graph.hierarchy_mut() = hierarchy;
    (graph.finish().unwrap(), conformer, atom, site)
}

#[test]
fn model_preserves_local_ids_and_dense_round_trips() {
    let (small, conformer, a, b, _) = two_atom_small(1.5);
    let mut builder = Model::builder();
    let instance = builder.add_molecule(&small, &conformer).unwrap();
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
        Quantity::new(Point3::new(1.5, 0.0, 0.0), ANGSTROM)
    );
    assert!(model
        .topology()
        .atom(InstanceAtomId::new(instance, AtomId::new(1)))
        .is_err());
}

#[test]
fn model_converts_source_conformer_units_once_without_mutating_the_source() {
    let mut graph = crate::core::MoleculeEditor::new();
    let atom = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let mut conformer = Conformer::new(NANOMETER).unwrap();
    conformer
        .set_position(atom, Quantity::new(Point3::new(0.15, 0.0, 0.0), NANOMETER))
        .unwrap();
    let small = graph.finish().unwrap();

    let model = Model::from_molecule(&small, &conformer).unwrap();
    let qualified = InstanceAtomId::new(MoleculeInstanceId::new(0), atom);
    assert_eq!(model.position(qualified).unwrap().unit(), ANGSTROM);
    assert_eq!(model.position(qualified).unwrap().x, 1.5);
    assert_eq!(conformer.unit(), NANOMETER);
    assert_eq!(conformer.position(atom).unwrap().x, 0.15);
}

#[test]
fn topology_allocation_is_shared_only_by_model_clones() {
    let (small, conformer, _, b, _) = two_atom_small(1.5);
    let model = Model::from_molecule(&small, &conformer).unwrap();
    let mut cloned = model.clone();
    cloned
        .set_position(
            InstanceAtomId::new(MoleculeInstanceId::new(0), b),
            Quantity::new(Point3::new(2.0, 0.0, 0.0), ANGSTROM),
        )
        .unwrap();
    let rebuilt = Model::from_molecule(&small, &conformer).unwrap();

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
    let (small, small_conformer, _, _, _) = two_atom_small(1.0);
    let (macromolecule, macro_conformer, atom, site) = one_atom_macro();
    let mut builder = Model::builder();
    let small_id = builder.add_molecule(&small, &small_conformer).unwrap();
    let macro_id = builder
        .add_molecule(&macromolecule, &macro_conformer)
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
    let hierarchy = model.topology().hierarchy(macro_id).unwrap().unwrap();
    assert_eq!(
        hierarchy
            .atom_for_site(InstanceAtomSiteId::new(macro_id, site))
            .unwrap(),
        InstanceAtomId::new(macro_id, atom)
    );
}

#[test]
fn repeated_molecules_get_distinct_instance_ids() {
    let (small, conformer, atom, _, _) = two_atom_small(1.0);
    let mut builder = Model::builder();
    let first = builder.add_molecule(&small, &conformer).unwrap();
    let second = builder.add_molecule(&small, &conformer).unwrap();
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
    let (small, conformer, a, _, _) = two_atom_small(1.0);
    let source = small.clone();
    let mut model = Model::from_molecule(&small, &conformer).unwrap();
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(0), a);
    model
        .set_position(atom, Quantity::new(Point3::new(3.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    assert_eq!(small, source);
    assert_eq!(
        conformer.position(a),
        Some(Quantity::new(Point3::new(0.0, 0.0, 0.0), ANGSTROM))
    );
    assert_eq!(model.topology().definition_count(), 1);
}

#[test]
fn construction_rejects_empty_missing_and_nonfinite_inputs_transactionally() {
    assert!(matches!(
        Model::builder().build(),
        Err(ModelBuildError::Topology(
            TopologyBuildError::NoMoleculeInstances
        ))
    ));
    let (small, mut conformer, a, _, _) = two_atom_small(1.0);
    conformer
        .set_position(a, Quantity::new(Point3::new(f64::NAN, 0.0, 0.0), ANGSTROM))
        .unwrap();
    let mut builder = Model::builder();
    assert!(
        matches!(builder.add_molecule(&small, &conformer), Err(ModelBuildError::NonFinitePosition { atom }) if atom == a)
    );
    assert!(matches!(
        builder.build(),
        Err(ModelBuildError::Topology(
            TopologyBuildError::NoMoleculeInstances
        ))
    ));
}

#[test]
fn position_updates_are_complete_finite_and_transactional() {
    let (small, conformer, a, _, _) = two_atom_small(1.0);
    let mut model = Model::from_molecule(&small, &conformer).unwrap();
    let original = model.positions().values().value().to_vec();
    assert!(matches!(
        model.set_positions(Quantity::new(&[Point3::default()], ANGSTROM)),
        Err(PositionError::PositionCountMismatch { .. })
    ));
    assert_eq!(model.positions().values().to_value(), original.as_slice());
    let mut invalid = original.clone();
    invalid[0] = Point3::new(f64::INFINITY, 0.0, 0.0);
    assert!(
        matches!(model.set_positions(Quantity::new(&invalid, ANGSTROM)), Err(PositionError::NonFinitePosition { atom }) if atom.atom() == a)
    );
    assert_eq!(model.positions().values().to_value(), original.as_slice());
}

#[test]
fn harmonic_potential_and_minimization_use_instance_qualified_topology() {
    let (small, conformer, _, _, bond) = two_atom_small(2.0);
    let model = Model::from_molecule(&small, &conformer).unwrap();
    let qualified = InstanceBondId::new(MoleculeInstanceId::new(0), bond);
    let mut potential = HarmonicBondPotential::new(
        &model.shared_topology(),
        [HarmonicBondParameter::new(
            qualified,
            Quantity::new(1.0, ANGSTROM),
            Quantity::new(100.0, MODEL_FORCE_CONSTANT_UNIT),
        )],
    )
    .unwrap();
    let initial = potential.evaluate(model.view()).unwrap();
    assert!((initial.energy().to_value() - 50.0).abs() < 1.0e-10);
    let result = minimize(&model, &mut potential, MinimizeOptions::default()).unwrap();
    assert!(result.final_energy < result.initial_energy);
    assert_eq!(
        model.positions().values().value()[1],
        Point3::new(2.0, 0.0, 0.0)
    );

    let rebuilt = Model::from_molecule(&small, &conformer).unwrap();
    assert_eq!(
        potential.evaluate(rebuilt.view()),
        Err(PotentialError::IncompatibleTopology)
    );

    let mut coincident = model.clone();
    let instance = MoleculeInstanceId::new(0);
    coincident
        .set_position(
            InstanceAtomId::new(instance, AtomId::new(2)),
            Quantity::new(coincident.positions().values().value()[0], ANGSTROM),
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
    let (small, conformer, _, _, bond) = two_atom_small(9.8);
    let mut model = Model::from_molecule(&small, &conformer).unwrap();
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
            Quantity::new(100.0, MODEL_FORCE_CONSTANT_UNIT),
        )],
    )
    .unwrap();

    assert!(potential.evaluate(model.view()).is_ok());
    let nonperiodic_ensemble = Ensemble::from_models(&[model.clone()]).unwrap();
    assert!(potential
        .evaluate(nonperiodic_ensemble.views().next().unwrap())
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
        potential.evaluate(periodic_ensemble.views().next().unwrap()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );
    let mut independent = Model::from_molecule(&small, &conformer).unwrap();
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
        if coordinate <= 0.25 {
            return Err(PotentialError::invalid_geometry(
                "test coordinate",
                [model.topology().atom_ids()[1]],
                PotentialGeometryError::CoincidentAtoms,
            ));
        }
        PotentialEvaluation::new(
            model,
            Quantity::new(0.5 * coordinate * coordinate, MODEL_ENERGY_UNIT),
            Quantity::new(
                vec![Vector3::zero(), Vector3::new(coordinate, 0.0, 0.0)],
                MODEL_GRADIENT_UNIT,
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
            Quantity::new(0.5, MODEL_ENERGY_UNIT),
            Quantity::new(
                vec![Vector3::zero(), Vector3::new(1.0, 0.0, 0.0)],
                MODEL_GRADIENT_UNIT,
            ),
        )
    }
}

#[test]
fn minimization_backtracks_invalid_geometry_but_propagates_backend_failures() {
    let (small, conformer, _, _, _) = two_atom_small(1.0);
    let model = Model::from_molecule(&small, &conformer).unwrap();
    let options = MinimizeOptions {
        max_iterations: 1,
        initial_step: Quantity::new(1.0, ANGSTROM),
        ..MinimizeOptions::default()
    };

    let result = minimize(&model, &mut RecoverableGeometryPotential, options).unwrap();
    assert_eq!(result.status, MinimizationStatus::MaxIterations);
    assert_eq!(result.iterations, 1);
    assert_eq!(result.evaluations, 3);
    assert_eq!(result.model.positions().values().value()[1].x, 0.5);
    assert_eq!(model.positions().values().value()[1].x, 1.0);

    let error = minimize(&model, &mut BackendFailurePotential { calls: 0 }, options).unwrap_err();
    assert!(matches!(
        error,
        MinimizationError::Potential(PotentialError::Backend {
            backend: "test backend",
            ..
        })
    ));
}
