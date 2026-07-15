use super::*;
use crate::bio::{AtomSiteId, AtomSiteMetadata, BioHierarchy, MacroMolecule};
use crate::core::{
    Atom, AtomId, BondId, BondOrder, Conformer, ConformerId, Element, Molecule, Point3,
};
use crate::modeling::potential::{
    HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError, PotentialEvaluation,
    PotentialGeometryError, Vector3,
};
use crate::small::SmallMolecule;

fn two_atom_small(distance: f64) -> (SmallMolecule, ConformerId, AtomId, AtomId, BondId) {
    let mut graph = Molecule::new();
    let carbon = Element::from_symbol("C").unwrap();
    let a = graph.add_atom(Atom::new(carbon));
    let tombstone = graph.add_atom(Atom::new(carbon));
    graph.delete_atom(tombstone).unwrap();
    let b = graph.add_atom(Atom::new(carbon));
    let bond = graph.add_bond(a, b, BondOrder::Single).unwrap();
    let mut conformer = Conformer::new();
    conformer.set_position(a, Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(b, Point3::new(distance, 0.0, 0.0));
    let conformer = graph.add_conformer(conformer).expect("valid conformer");
    (SmallMolecule::from_graph(graph), conformer, a, b, bond)
}

fn one_atom_macro() -> (MacroMolecule, ConformerId, AtomId, AtomSiteId) {
    let mut graph = Molecule::new();
    let atom = graph.add_atom(Atom::new(Element::from_symbol("N").unwrap()));
    let mut conformer = Conformer::new();
    conformer.set_position(atom, Point3::new(2.0, 0.0, 0.0));
    let conformer = graph.add_conformer(conformer).expect("valid conformer");
    let mut hierarchy = BioHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    let site = hierarchy
        .add_atom_site(residue, atom, AtomSiteMetadata::default())
        .unwrap();
    (
        MacroMolecule::from_parts(graph, hierarchy),
        conformer,
        atom,
        site,
    )
}

#[test]
fn model_preserves_local_ids_and_dense_round_trips() {
    let (small, conformer, a, b, _) = two_atom_small(1.5);
    let mut builder = MolecularModel::builder();
    let instance = builder.add_small_molecule(&small, conformer).unwrap();
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
    assert_eq!(model.position(qb).unwrap(), Point3::new(1.5, 0.0, 0.0));
    assert!(model
        .topology()
        .atom(InstanceAtomId::new(instance, AtomId::new(1)))
        .is_err());
}

#[test]
fn model_definition_identity_is_shared_only_by_clones() {
    let (small, conformer, _, b, _) = two_atom_small(1.5);
    let model = MolecularModel::from_small_molecule(&small, conformer).unwrap();
    let mut cloned = model.clone();
    cloned
        .set_position(
            InstanceAtomId::new(MoleculeInstanceId::new(0), b),
            Point3::new(2.0, 0.0, 0.0),
        )
        .unwrap();
    let rebuilt = MolecularModel::from_small_molecule(&small, conformer).unwrap();

    assert_eq!(model.definition_key(), cloned.definition_key());
    assert_ne!(model.definition_key(), rebuilt.definition_key());
    assert_ne!(model, cloned);
    assert_eq!(model, rebuilt);
}

#[test]
fn mixed_instances_and_hierarchy_use_qualified_ids() {
    let (small, small_conformer, _, _, _) = two_atom_small(1.0);
    let (macromolecule, macro_conformer, atom, site) = one_atom_macro();
    let mut metadata = MoleculeInstanceMetadata::default();
    metadata.insert_role(MoleculeRole::Ligand);
    let mut builder = MolecularModel::builder();
    let small_id = builder
        .add_small_molecule_with_metadata(&small, small_conformer, metadata)
        .unwrap();
    let macro_id = builder
        .add_macro_molecule(&macromolecule, macro_conformer)
        .unwrap();
    let model = builder.build().unwrap();
    assert_ne!(small_id, macro_id);
    assert!(model
        .topology()
        .molecule(small_id)
        .unwrap()
        .has_role(MoleculeRole::Ligand));
    let hierarchy = model
        .topology()
        .molecule(macro_id)
        .unwrap()
        .bio_hierarchy()
        .unwrap();
    assert_eq!(
        hierarchy.atom_for_site(site).unwrap(),
        InstanceAtomId::new(macro_id, atom)
    );
}

#[test]
fn repeated_molecules_get_distinct_instance_ids() {
    let (small, conformer, atom, _, _) = two_atom_small(1.0);
    let mut builder = MolecularModel::builder();
    let first = builder.add_small_molecule(&small, conformer).unwrap();
    let second = builder.add_small_molecule(&small, conformer).unwrap();
    let model = builder.build().unwrap();
    assert_ne!(first, second);
    assert_ne!(
        InstanceAtomId::new(first, atom),
        InstanceAtomId::new(second, atom)
    );
    assert_eq!(model.topology().molecule_count(), 2);
}

#[test]
fn construction_copies_positions_and_preserves_sources() {
    let (small, conformer, a, _, _) = two_atom_small(1.0);
    let source = small.clone();
    let mut model = MolecularModel::from_small_molecule(&small, conformer).unwrap();
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(0), a);
    model
        .set_position(atom, Point3::new(3.0, 0.0, 0.0))
        .unwrap();
    assert_eq!(small, source);
    assert_eq!(
        small.graph().conformer(conformer).unwrap().position(a),
        Some(Point3::new(0.0, 0.0, 0.0))
    );
    assert_eq!(
        model
            .topology()
            .molecule(MoleculeInstanceId::new(0))
            .unwrap()
            .graph()
            .conformers()
            .count(),
        0
    );
}

#[test]
fn construction_rejects_empty_missing_and_nonfinite_inputs_transactionally() {
    assert_eq!(
        MolecularModel::builder().build(),
        Err(ModelBuildError::EmptyModel)
    );
    let empty = SmallMolecule::new();
    let mut builder = MolecularModel::builder();
    assert!(matches!(
        builder.add_small_molecule(&empty, ConformerId::new(0)),
        Err(ModelBuildError::EmptyMolecule)
    ));
    assert_eq!(builder.build(), Err(ModelBuildError::EmptyModel));

    let (mut small, conformer, a, _, _) = two_atom_small(1.0);
    small
        .graph_mut_raw()
        .conformer_mut(conformer)
        .unwrap()
        .set_position(a, Point3::new(f64::NAN, 0.0, 0.0));
    let mut builder = MolecularModel::builder();
    assert!(
        matches!(builder.add_small_molecule(&small, conformer), Err(ModelBuildError::NonFinitePosition { atom }) if atom == a)
    );
    assert_eq!(builder.build(), Err(ModelBuildError::EmptyModel));
}

#[test]
fn position_updates_are_complete_finite_and_transactional() {
    let (small, conformer, a, _, _) = two_atom_small(1.0);
    let mut model = MolecularModel::from_small_molecule(&small, conformer).unwrap();
    let original = model.positions().to_vec();
    assert!(matches!(
        model.set_positions(&[Point3::default()]),
        Err(PositionError::PositionCountMismatch { .. })
    ));
    assert_eq!(model.positions(), original);
    let mut invalid = original.clone();
    invalid[0] = Point3::new(f64::INFINITY, 0.0, 0.0);
    assert!(
        matches!(model.set_positions(&invalid), Err(PositionError::NonFinitePosition { atom }) if atom.atom() == a)
    );
    assert_eq!(model.positions(), original);
}

#[test]
fn harmonic_potential_and_minimization_use_instance_qualified_topology() {
    let (small, conformer, _, _, bond) = two_atom_small(2.0);
    let model = MolecularModel::from_small_molecule(&small, conformer).unwrap();
    let qualified = InstanceBondId::new(MoleculeInstanceId::new(0), bond);
    let mut potential =
        HarmonicBondPotential::new(&model, [HarmonicBondParameter::new(qualified, 1.0, 100.0)])
            .unwrap();
    let initial = potential.evaluate(&model).unwrap();
    assert!((initial.energy() - 50.0).abs() < 1.0e-10);
    let result = minimize(&model, &mut potential, MinimizeOptions::default()).unwrap();
    assert!(result.final_energy < result.initial_energy);
    assert_eq!(model.positions()[1], Point3::new(2.0, 0.0, 0.0));

    let rebuilt = MolecularModel::from_small_molecule(&small, conformer).unwrap();
    assert_eq!(
        potential.evaluate(&rebuilt),
        Err(PotentialError::IncompatibleModel)
    );

    let mut coincident = model.clone();
    let instance = MoleculeInstanceId::new(0);
    coincident
        .set_position(
            InstanceAtomId::new(instance, AtomId::new(2)),
            coincident.positions()[0],
        )
        .unwrap();
    assert_eq!(
        potential.evaluate(&coincident),
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

struct RecoverableGeometryPotential;

impl Potential for RecoverableGeometryPotential {
    fn evaluate(&mut self, model: &MolecularModel) -> Result<PotentialEvaluation, PotentialError> {
        let coordinate = model.positions()[1].x;
        if coordinate <= 0.25 {
            return Err(PotentialError::invalid_geometry(
                "test coordinate",
                [model.topology().atom_ids()[1]],
                PotentialGeometryError::CoincidentAtoms,
            ));
        }
        PotentialEvaluation::new(
            model,
            0.5 * coordinate * coordinate,
            vec![Vector3::zero(), Vector3::new(coordinate, 0.0, 0.0)],
        )
    }
}

struct BackendFailurePotential {
    calls: usize,
}

impl Potential for BackendFailurePotential {
    fn evaluate(&mut self, model: &MolecularModel) -> Result<PotentialEvaluation, PotentialError> {
        self.calls += 1;
        if self.calls > 1 {
            return Err(PotentialError::backend("test backend", "evaluation failed"));
        }
        PotentialEvaluation::new(
            model,
            0.5,
            vec![Vector3::zero(), Vector3::new(1.0, 0.0, 0.0)],
        )
    }
}

#[test]
fn minimization_backtracks_invalid_geometry_but_propagates_backend_failures() {
    let (small, conformer, _, _, _) = two_atom_small(1.0);
    let model = MolecularModel::from_small_molecule(&small, conformer).unwrap();
    let options = MinimizeOptions {
        max_iterations: 1,
        initial_step: 1.0,
        ..MinimizeOptions::default()
    };

    let result = minimize(&model, &mut RecoverableGeometryPotential, options).unwrap();
    assert_eq!(result.status, MinimizationStatus::MaxIterations);
    assert_eq!(result.iterations, 1);
    assert_eq!(result.evaluations, 3);
    assert_eq!(result.model.positions()[1].x, 0.5);
    assert_eq!(model.positions()[1].x, 1.0);

    let error = minimize(&model, &mut BackendFailurePotential { calls: 0 }, options).unwrap_err();
    assert!(matches!(
        error,
        MinimizationError::Potential(PotentialError::Backend {
            backend: "test backend",
            ..
        })
    ));
}
