use super::potential::*;
use super::*;
use crate::core::*;
use crate::small::SmallMolecule;

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).expect("test element"))
}

fn diatomic(distance: f64) -> (SmallMolecule, ConformerId, BondId) {
    let mut graph = Molecule::new();
    let a = graph.add_atom(atom("C"));
    let b = graph.add_atom(atom("O"));
    let bond = graph.add_bond(a, b, BondOrder::Single).expect("bond");
    let mut conformer = Conformer::new();
    conformer.set_position(a, Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(b, Point3::new(distance, 0.0, 0.0));
    let conformer = graph.add_conformer(conformer);
    (SmallMolecule::from_graph(graph), conformer, bond)
}

fn diatomic_at(a: Point3, b: Point3) -> (SmallMolecule, ConformerId, BondId) {
    let mut graph = Molecule::new();
    let a_id = graph.add_atom(atom("C"));
    let b_id = graph.add_atom(atom("O"));
    let bond = graph.add_bond(a_id, b_id, BondOrder::Single).expect("bond");
    let mut conformer = Conformer::new();
    conformer.set_position(a_id, a);
    conformer.set_position(b_id, b);
    let conformer = graph.add_conformer(conformer);
    (SmallMolecule::from_graph(graph), conformer, bond)
}

#[test]
fn model_from_conformer_copies_one_complete_coordinate_set() {
    let (mut source, conformer, _) = diatomic(1.5);
    let mut second = Conformer::new();
    second.set_position(AtomId::new(0), Point3::new(10.0, 0.0, 0.0));
    second.set_position(AtomId::new(1), Point3::new(12.0, 0.0, 0.0));
    source.graph_mut().add_conformer(second);

    let model = MolecularModel::from_conformer(&source, conformer).expect("model");

    assert_eq!(model.atom_count(), 2);
    assert_eq!(model.positions()[1], Point3::new(1.5, 0.0, 0.0));
    assert_eq!(model.topology().conformers().count(), 0);
    assert_eq!(source.graph().conformers().count(), 2);
    assert_eq!(model.components().count(), 1);
}

#[test]
fn builder_remaps_tombstones_and_multiple_components() {
    let mut graph = Molecule::new();
    let first = graph.add_atom(atom("C"));
    let deleted = graph.add_atom(atom("N"));
    let second = graph.add_atom(atom("O"));
    graph.delete_atom(deleted).expect("delete tombstone");
    let source_bond = graph
        .add_bond(first, second, BondOrder::Double)
        .expect("bond across source IDs");
    let mut conformer = Conformer::new();
    conformer.set_position(first, Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(second, Point3::new(1.2, 0.0, 0.0));
    let conformer = graph.add_conformer(conformer);
    let source = SmallMolecule::from_graph(graph);
    let (other, other_conformer, _) = diatomic(2.0);

    let mut builder = MolecularModel::builder();
    let first_mapping = builder
        .add_component(&source, conformer)
        .expect("first component");
    let second_mapping = builder
        .add_component(&other, other_conformer)
        .expect("second component");
    let model = builder.build().expect("model");

    assert_eq!(first_mapping.atom(first), Some(AtomId::new(0)));
    assert_eq!(first_mapping.atom(second), Some(AtomId::new(1)));
    assert_eq!(first_mapping.bond(source_bond), Some(BondId::new(0)));
    assert_eq!(second_mapping.component(), ComponentId::new(1));
    assert_eq!(model.atom_count(), 4);
    assert_eq!(model.components().count(), 2);
    assert_eq!(
        model
            .component_for_atom(AtomId::new(3))
            .expect("membership")
            .id(),
        ComponentId::new(1)
    );
}

#[test]
fn model_copy_preserves_payload_properties_and_stored_stereo() {
    let mut graph = Molecule::new();
    graph
        .props_mut()
        .insert("name".into(), PropValue::String("chiral".into()));
    let center = graph.add_atom(atom("C"));
    let fluorine = graph.add_atom(atom("F"));
    let chlorine = graph.add_atom(atom("Cl"));
    let bromine = graph.add_atom(atom("Br"));
    let f_bond = graph
        .add_bond(center, fluorine, BondOrder::Single)
        .expect("C-F");
    graph
        .add_bond(center, chlorine, BondOrder::Single)
        .expect("C-Cl");
    graph
        .add_bond(center, bromine, BondOrder::Single)
        .expect("C-Br");
    graph
        .atom_mut(center)
        .expect("center")
        .props
        .insert("atom-note".into(), PropValue::String("preserve".into()));
    graph
        .bond_mut(f_bond)
        .expect("bond")
        .props
        .insert("bond-note".into(), PropValue::Bool(true));
    let element = graph
        .add_stereo_element(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers: vec![
                    StereoCarrier::Atom(fluorine),
                    StereoCarrier::Atom(chlorine),
                    StereoCarrier::Atom(bromine),
                    StereoCarrier::ImplicitHydrogen,
                ],
                orientation: TetrahedralOrientation::Clockwise,
            }),
            StereoSource::User,
        ))
        .expect("stereo element");
    graph
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: vec![element],
        })
        .expect("stereo group");
    graph
        .set_stereo_bond_mark(StereoBondMark {
            bond: f_bond,
            kind: StereoBondMarkKind::WedgeUp,
            source: StereoSource::User,
        })
        .expect("bond mark");
    graph
        .stereo_element_mut(element)
        .expect("descriptor")
        .descriptor = Some(StereoDescriptor::R);
    let mut conformer = Conformer::new();
    for (index, atom) in [center, fluorine, chlorine, bromine]
        .into_iter()
        .enumerate()
    {
        conformer.set_position(atom, Point3::new(index as f64, 0.0, 0.0));
    }
    let conformer = graph.add_conformer(conformer);
    let source = SmallMolecule::from_graph(graph);

    let model = MolecularModel::from_conformer(&source, conformer).expect("model");
    let component = model.component(ComponentId::new(0)).expect("component");
    let (_, copied_stereo) = model
        .topology()
        .stereo_elements()
        .next()
        .expect("copied stereo");

    assert_eq!(
        component.props().get("name"),
        Some(&PropValue::String("chiral".into()))
    );
    assert_eq!(
        model
            .topology()
            .atom(AtomId::new(0))
            .unwrap()
            .props
            .get("atom-note"),
        Some(&PropValue::String("preserve".into()))
    );
    assert_eq!(
        model
            .topology()
            .bond(BondId::new(0))
            .unwrap()
            .props
            .get("bond-note"),
        Some(&PropValue::Bool(true))
    );
    assert_eq!(model.topology().stereo_groups().count(), 1);
    assert_eq!(model.topology().stereo_bond_marks().count(), 1);
    assert_eq!(copied_stereo.descriptor, None);
}

#[test]
fn model_construction_errors_are_transactional() {
    let (valid, valid_conformer, _) = diatomic(1.0);
    let mut incomplete_graph = Molecule::new();
    let a = incomplete_graph.add_atom(atom("C"));
    incomplete_graph.add_atom(atom("O"));
    let mut incomplete = Conformer::new();
    incomplete.set_position(a, Point3::new(0.0, 0.0, 0.0));
    let incomplete_id = incomplete_graph.add_conformer(incomplete);
    let incomplete = SmallMolecule::from_graph(incomplete_graph);

    let mut builder = MolecularModel::builder();
    builder
        .add_component(&valid, valid_conformer)
        .expect("valid component");
    assert_eq!(
        builder.add_component(&incomplete, incomplete_id),
        Err(ModelBuildError::MissingPosition {
            atom: AtomId::new(1)
        })
    );
    let model = builder.build().expect("staged valid component survives");
    assert_eq!(model.components().count(), 1);

    assert_eq!(
        MolecularModel::from_conformer(&valid, ConformerId::new(99)),
        Err(ModelBuildError::InvalidConformerId(ConformerId::new(99)))
    );
    assert_eq!(
        MolecularModel::builder().build(),
        Err(ModelBuildError::EmptyModel)
    );

    let mut empty_graph = Molecule::new();
    let empty_conformer = empty_graph.add_conformer(Conformer::new());
    assert_eq!(
        MolecularModel::from_conformer(&SmallMolecule::from_graph(empty_graph), empty_conformer),
        Err(ModelBuildError::EmptyComponent)
    );
}

#[test]
fn model_rejects_non_finite_and_non_transactional_position_updates() {
    let (mut source, conformer, _) = diatomic(1.0);
    source
        .graph_mut()
        .conformer_mut(conformer)
        .expect("conformer")
        .set_position(AtomId::new(1), Point3::new(f64::NAN, 0.0, 0.0));
    assert_eq!(
        MolecularModel::from_conformer(&source, conformer),
        Err(ModelBuildError::NonFinitePosition {
            atom: AtomId::new(1)
        })
    );

    let (source, conformer, _) = diatomic(1.0);
    let mut model = MolecularModel::from_conformer(&source, conformer).expect("model");
    let original = model.positions().to_vec();
    let invalid = [original[0], Point3::new(f64::INFINITY, 0.0, 0.0)];
    assert_eq!(
        model.set_positions(&invalid),
        Err(PositionError::NonFinitePosition {
            atom: AtomId::new(1)
        })
    );
    assert_eq!(model.positions(), original);
    assert_eq!(
        model.set_positions(&original[..1]),
        Err(PositionError::PositionCountMismatch {
            expected: 2,
            actual: 1
        })
    );
}

#[test]
fn harmonic_bond_energy_and_gradient_match_finite_difference() {
    let (source, conformer, source_bond) =
        diatomic_at(Point3::new(0.2, -0.4, 0.7), Point3::new(1.1, 1.6, 2.2));
    let mut builder = MolecularModel::builder();
    let mapping = builder
        .add_component(&source, conformer)
        .expect("component");
    let model = builder.build().expect("model");
    let bond = mapping.bond(source_bond).expect("mapped bond");
    let parameter = HarmonicBondParameter::new(bond, 1.5, 12.0);
    let mut potential = HarmonicBondPotential::new(&model, [parameter]).expect("potential");
    let evaluation = potential.evaluate(&model).expect("evaluation");

    let epsilon = 1.0e-6;
    for atom_index in 0..2 {
        for axis in 0..3 {
            let mut plus = model.clone();
            let mut minus = model.clone();
            let mut plus_point = plus.positions()[atom_index];
            let mut minus_point = minus.positions()[atom_index];
            coordinate_mut(&mut plus_point, axis, epsilon);
            coordinate_mut(&mut minus_point, axis, -epsilon);
            plus.set_position(AtomId::new(atom_index as u32), plus_point)
                .unwrap();
            minus
                .set_position(AtomId::new(atom_index as u32), minus_point)
                .unwrap();
            let finite_difference = (potential.evaluate(&plus).unwrap().energy()
                - potential.evaluate(&minus).unwrap().energy())
                / (2.0 * epsilon);
            let analytic = vector_coordinate(evaluation.gradient()[atom_index], axis);
            assert!((analytic - finite_difference).abs() < 1.0e-6);
        }
    }
}

fn coordinate_mut(point: &mut Point3, axis: usize, delta: f64) {
    match axis {
        0 => point.x += delta,
        1 => point.y += delta,
        2 => point.z += delta,
        _ => unreachable!(),
    }
}

fn vector_coordinate(vector: Vector3, axis: usize) -> f64 {
    match axis {
        0 => vector.x,
        1 => vector.y,
        2 => vector.z,
        _ => unreachable!(),
    }
}

#[test]
fn harmonic_potential_rejects_invalid_terms_and_evaluations() {
    let (source, conformer, source_bond) = diatomic(1.0);
    let mut builder = MolecularModel::builder();
    let mapping = builder.add_component(&source, conformer).unwrap();
    let model = builder.build().unwrap();
    let bond = mapping.bond(source_bond).unwrap();

    assert_eq!(
        HarmonicBondPotential::new(
            &model,
            [HarmonicBondParameter::new(BondId::new(99), 1.0, 1.0)]
        ),
        Err(PotentialError::InvalidBondId(BondId::new(99)))
    );
    assert!(matches!(
        HarmonicBondPotential::new(
            &model,
            [
                HarmonicBondParameter::new(bond, 1.0, 1.0),
                HarmonicBondParameter::new(bond, 1.0, 1.0)
            ]
        ),
        Err(PotentialError::DuplicateBondParameter(_))
    ));
    assert!(matches!(
        HarmonicBondPotential::new(&model, [HarmonicBondParameter::new(bond, 0.0, 1.0)]),
        Err(PotentialError::InvalidBondParameter { .. })
    ));
    assert_eq!(
        PotentialEvaluation::new(&model, 0.0, vec![Vector3::zero()]),
        Err(PotentialError::GradientLengthMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        PotentialEvaluation::new(
            &model,
            0.0,
            vec![Vector3::zero(), Vector3::new(f64::NAN, 0.0, 0.0)]
        ),
        Err(PotentialError::NonFiniteGradient {
            atom: AtomId::new(1)
        })
    );
    assert_eq!(
        PotentialEvaluation::new(
            &model,
            f64::INFINITY,
            vec![Vector3::zero(), Vector3::zero()]
        ),
        Err(PotentialError::NonFiniteEnergy)
    );

    let mut mismatched_graph = Molecule::new();
    let a = mismatched_graph.add_atom(atom("C"));
    let b = mismatched_graph.add_atom(atom("O"));
    mismatched_graph
        .add_bond(a, b, BondOrder::Double)
        .expect("different bond order");
    let mut mismatched_conformer = Conformer::new();
    mismatched_conformer.set_position(a, Point3::new(0.0, 0.0, 0.0));
    mismatched_conformer.set_position(b, Point3::new(1.0, 0.0, 0.0));
    let mismatched_conformer = mismatched_graph.add_conformer(mismatched_conformer);
    let mismatched = MolecularModel::from_conformer(
        &SmallMolecule::from_graph(mismatched_graph),
        mismatched_conformer,
    )
    .unwrap();
    let mut bound =
        HarmonicBondPotential::new(&model, [HarmonicBondParameter::new(bond, 1.0, 1.0)]).unwrap();
    assert_eq!(
        bound.evaluate(&mismatched),
        Err(PotentialError::ModelTopologyMismatch(bond))
    );

    let (coincident, conformer, source_bond) = diatomic(0.0);
    let mut builder = MolecularModel::builder();
    let mapping = builder.add_component(&coincident, conformer).unwrap();
    let coincident = builder.build().unwrap();
    let bond = mapping.bond(source_bond).unwrap();
    let mut potential =
        HarmonicBondPotential::new(&coincident, [HarmonicBondParameter::new(bond, 1.0, 1.0)])
            .unwrap();
    assert_eq!(
        potential.evaluate(&coincident),
        Err(PotentialError::CoincidentBondAtoms(bond))
    );
}

#[test]
fn minimization_converges_without_mutating_input() {
    let (source, conformer, source_bond) = diatomic(2.0);
    let mut builder = MolecularModel::builder();
    let mapping = builder.add_component(&source, conformer).unwrap();
    let model = builder.build().unwrap();
    let original = model.clone();
    let bond = mapping.bond(source_bond).unwrap();
    let mut potential =
        HarmonicBondPotential::new(&model, [HarmonicBondParameter::new(bond, 1.0, 100.0)]).unwrap();

    let result = minimize(&model, &mut potential, MinimizeOptions::default()).unwrap();
    let distance = (result.model.positions()[1].x - result.model.positions()[0].x).abs();

    assert_eq!(result.status, MinimizationStatus::Converged);
    assert!(result.final_energy < result.initial_energy);
    assert!((distance - 1.0).abs() < 1.0e-6);
    assert!(result.final_max_gradient <= 1.0e-4);
    assert_eq!(model, original);
}

#[test]
fn minimization_handles_multiple_components_and_limits() {
    let (first, first_conformer, first_bond) = diatomic(2.0);
    let (second, second_conformer, second_bond) = diatomic(3.0);
    let mut builder = MolecularModel::builder();
    let first_mapping = builder.add_component(&first, first_conformer).unwrap();
    let second_mapping = builder.add_component(&second, second_conformer).unwrap();
    let model = builder.build().unwrap();
    let mut potential = HarmonicBondPotential::new(
        &model,
        [
            HarmonicBondParameter::new(first_mapping.bond(first_bond).unwrap(), 1.0, 50.0),
            HarmonicBondParameter::new(second_mapping.bond(second_bond).unwrap(), 1.5, 50.0),
        ],
    )
    .unwrap();
    let result = minimize(&model, &mut potential, MinimizeOptions::default()).unwrap();
    assert_eq!(result.status, MinimizationStatus::Converged);
    assert!((distance(&result.model, 0, 1) - 1.0).abs() < 1.0e-6);
    assert!((distance(&result.model, 2, 3) - 1.5).abs() < 1.0e-6);

    let mut potential = HarmonicBondPotential::new(
        &model,
        [HarmonicBondParameter::new(
            first_mapping.bond(first_bond).unwrap(),
            1.0,
            50.0,
        )],
    )
    .unwrap();
    let result = minimize(
        &model,
        &mut potential,
        MinimizeOptions {
            max_iterations: 0,
            ..MinimizeOptions::default()
        },
    )
    .unwrap();
    assert_eq!(result.status, MinimizationStatus::MaxIterations);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.evaluations, 1);
}

fn distance(model: &MolecularModel, a: usize, b: usize) -> f64 {
    let a = model.positions()[a];
    let b = model.positions()[b];
    ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
}

struct ConstantEnergyPotential;

impl Potential for ConstantEnergyPotential {
    fn evaluate(
        &mut self,
        model: &MolecularModel,
    ) -> std::result::Result<PotentialEvaluation, PotentialError> {
        PotentialEvaluation::new(
            model,
            1.0,
            vec![Vector3::new(1.0, 0.0, 0.0); model.atom_count()],
        )
    }
}

#[test]
fn minimization_reports_line_search_stall_and_invalid_options() {
    let (source, conformer, _) = diatomic(1.0);
    let model = MolecularModel::from_conformer(&source, conformer).unwrap();
    let mut potential = ConstantEnergyPotential;
    let result = minimize(&model, &mut potential, MinimizeOptions::default()).unwrap();
    assert_eq!(result.status, MinimizationStatus::LineSearchStalled);
    assert_eq!(result.model, model);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.evaluations, 25);

    let error = minimize(
        &model,
        &mut potential,
        MinimizeOptions {
            gradient_tolerance: 0.0,
            ..MinimizeOptions::default()
        },
    )
    .expect_err("invalid tolerance");
    assert!(matches!(error, MinimizationError::InvalidOptions(_)));
}
