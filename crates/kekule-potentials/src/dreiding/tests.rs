use kekule::bio::{Hierarchy, SmcraAtomSiteMetadata};
use kekule::core::{Atom, AtomId, BondOrder, Conformer, Element, HydrogenDeclaration, Molecule};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::modeling::potential::{Potential, PotentialError};
use kekule::structure::{Ensemble, Model};
use kekule::topology::{InstanceAtomId, MoleculeInstanceId};
use kekule_traj::{FrameBuffer, TrajectoryFrame};

use super::{DreidingPotential, DreidingPrepareError, DreidingPrepareOptions, QeqGrouping};

fn explicit_atom(symbol: &str) -> Atom {
    let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
    atom.hydrogens = HydrogenDeclaration::Fixed(0);
    atom
}

fn molecule(
    elements: &[&str],
    bonds: &[(usize, usize, BondOrder)],
    positions: &[Point3],
) -> (Molecule, Conformer) {
    let mut graph = kekule::core::MoleculeEditor::new();
    let atoms = elements
        .iter()
        .map(|symbol| {
            graph
                .add_atom(explicit_atom(symbol))
                .expect("atom identifier capacity")
        })
        .collect::<Vec<_>>();
    for &(a, b, order) in bonds {
        graph.add_bond(atoms[a], atoms[b], order).unwrap();
    }
    let mut conformer = Conformer::new(kekule::units::ANGSTROM).unwrap();
    for (&atom, &position) in atoms.iter().zip(positions) {
        conformer
            .set_position(
                atom,
                kekule::units::Quantity::new(position, kekule::units::ANGSTROM),
            )
            .unwrap();
    }
    let graph = graph.finish().expect("connected test molecule");
    (graph, conformer)
}

fn water(offset: f64) -> (Molecule, Conformer) {
    molecule(
        &["O", "H", "H"],
        &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        &[
            Point3::new(offset, 0.0, 0.0),
            Point3::new(offset + 0.9575, 0.0, 0.0),
            Point3::new(offset - 0.2399, 0.9272, 0.0),
        ],
    )
}

#[test]
fn preparation_and_evaluation_are_finite() {
    let (water, conformer) = water(0.0);
    let model = Model::from_molecule(&water, &conformer).unwrap();
    let mut potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();
    let evaluation = potential.evaluate(model.view()).unwrap();
    let oxygen = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(0));
    assert!(evaluation.energy().is_finite());
    assert_eq!(evaluation.gradient().len(), 3);
    assert!(potential.atom_type(oxygen).is_some());
    assert!(potential.partial_charge(oxygen).unwrap().is_finite());
}

#[test]
fn qeq_is_prepared_per_molecule_instance() {
    let (first, first_conf) = water(0.0);
    let (second, second_conf) = water(5.0);
    let mut builder = Model::builder();
    let first_id = builder.add_molecule(&first, &first_conf).unwrap();
    let second_id = builder.add_molecule(&second, &second_conf).unwrap();
    let model = builder.build().unwrap();
    let potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();
    assert_eq!(potential.qeq_grouping(), QeqGrouping::MoleculeInstances);
    for instance in [first_id, second_id] {
        let total = (0..3)
            .map(|atom| {
                potential
                    .partial_charge(InstanceAtomId::new(instance, AtomId::new(atom)))
                    .unwrap()
                    .to_value()
            })
            .sum::<f64>();
        assert!(total.abs() < 1.0e-8);
    }
    assert_eq!(
        potential.nonbonded.len(),
        9,
        "two waters have nine inter-instance pairs and no intramolecular nonbonded pairs"
    );
}

#[test]
fn prepared_potential_evaluates_models_ensembles_and_frames_sharing_topology() {
    let (water, conformer) = water(0.0);
    let model = Model::from_molecule(&water, &conformer).unwrap();
    let mut displaced = model.clone();
    let hydrogen = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(1));
    displaced
        .set_position(
            hydrogen,
            kekule::units::Quantity::new(Point3::new(1.05, 0.0, 0.0), kekule::units::ANGSTROM),
        )
        .unwrap();
    let ensemble = Ensemble::from_models(&[model.clone(), displaced.clone()]).unwrap();
    let frame = TrajectoryFrame::new(
        displaced.positions().clone(),
        displaced.topology().bond_count(),
    );
    let mut potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();

    let energies = ensemble
        .views()
        .map(|view| potential.evaluate(view).unwrap().energy().to_value())
        .collect::<Vec<_>>();
    assert_eq!(energies.len(), 2);
    assert!(energies.iter().all(|energy| energy.is_finite()));
    let topology = model.shared_topology();
    let frame_view = frame.view(&topology).unwrap();
    assert!(potential
        .evaluate(frame_view.model_view())
        .unwrap()
        .energy()
        .is_finite());
}

#[test]
fn periodic_state_is_rejected_during_preparation_and_across_structural_views() {
    let (molecule, conformer) = molecule(
        &["C", "C"],
        &[(0, 1, BondOrder::Single)],
        &[Point3::new(0.1, 0.0, 0.0), Point3::new(9.9, 0.0, 0.0)],
    );
    let model = Model::from_molecule(&molecule, &conformer).unwrap();
    let cell = PeriodicCell::orthorhombic(
        kekule::units::Quantity::new(Vector3::new(10.0, 10.0, 10.0), kekule::units::ANGSTROM),
        [true; 3],
    )
    .unwrap();
    let mut periodic_model = model.clone();
    periodic_model.set_cell(Some(cell));

    assert!(matches!(
        DreidingPotential::prepare(
            &periodic_model.shared_topology(),
            periodic_model.view(),
            DreidingPrepareOptions::default(),
        ),
        Err(DreidingPrepareError::UnsupportedPeriodicCell)
    ));

    let mut potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();
    assert!(potential.evaluate(model.view()).is_ok());
    let mut nonperiodic_buffer = FrameBuffer::new(model.shared_topology());
    nonperiodic_buffer
        .set_positions(model.positions().values())
        .unwrap();
    assert!(potential.evaluate(nonperiodic_buffer.model_view()).is_ok());
    assert_eq!(
        potential.evaluate(periodic_model.view()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );

    let periodic_ensemble = Ensemble::from_models(&[periodic_model.clone()]).unwrap();
    assert_eq!(
        potential.evaluate(periodic_ensemble.views().next().unwrap()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );
    let mut periodic_frame = TrajectoryFrame::new(
        periodic_model.positions().clone(),
        periodic_model.topology().bond_count(),
    );
    periodic_frame.set_cell(periodic_model.cell().copied());
    assert_eq!(
        potential.evaluate(
            periodic_frame
                .view(&model.shared_topology())
                .unwrap()
                .model_view(),
        ),
        Err(PotentialError::UnsupportedPeriodicCell)
    );
    let mut periodic_buffer = FrameBuffer::new(model.shared_topology());
    periodic_buffer
        .set_positions(model.positions().values())
        .unwrap();
    periodic_buffer.set_cell(Some(cell));
    assert_eq!(
        potential.evaluate(periodic_buffer.model_view()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );

    let mut independent = Model::from_molecule(&molecule, &conformer).unwrap();
    independent.set_cell(Some(cell));
    assert_eq!(
        potential.evaluate(independent.view()),
        Err(PotentialError::IncompatibleTopology)
    );
}

#[test]
fn qeq_grouping_policy_is_explicit() {
    let (water, conformer) = water(0.0);
    let model = Model::from_molecule(&water, &conformer).unwrap();
    for grouping in [QeqGrouping::WholeTopology, QeqGrouping::MoleculeInstances] {
        let potential = DreidingPotential::prepare(
            &model.shared_topology(),
            model.view(),
            DreidingPrepareOptions {
                qeq_grouping: grouping,
            },
        )
        .unwrap();
        assert_eq!(potential.qeq_grouping(), grouping);
    }
}

#[test]
fn preparation_maps_tombstoned_local_ids_to_dense_adjacency() {
    let mut graph = kekule::core::MoleculeEditor::new();
    let oxygen = graph
        .add_atom(explicit_atom("O"))
        .expect("atom identifier capacity");
    let tombstone = graph
        .add_atom(explicit_atom("H"))
        .expect("atom identifier capacity");
    let first_hydrogen = graph
        .add_atom(explicit_atom("H"))
        .expect("atom identifier capacity");
    let second_hydrogen = graph
        .add_atom(explicit_atom("H"))
        .expect("atom identifier capacity");
    graph.delete_atom(tombstone).unwrap();
    graph
        .add_bond(oxygen, first_hydrogen, BondOrder::Single)
        .unwrap();
    graph
        .add_bond(oxygen, second_hydrogen, BondOrder::Single)
        .unwrap();
    let mut conformer = Conformer::new(kekule::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            oxygen,
            kekule::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), kekule::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            first_hydrogen,
            kekule::units::Quantity::new(Point3::new(0.9575, 0.0, 0.0), kekule::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            second_hydrogen,
            kekule::units::Quantity::new(
                Point3::new(-0.2399, 0.9272, 0.0),
                kekule::units::ANGSTROM,
            ),
        )
        .unwrap();
    let molecule = graph.finish().expect("connected water");
    let model = Model::from_molecule(&molecule, &conformer).unwrap();

    let potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();
    assert!(potential.nonbonded.is_empty());
}

#[test]
fn hierarchy_bearing_molecules_are_supported() {
    let (small, conformer) = water(0.0);
    let mut hierarchy = Hierarchy::new();
    let chain = hierarchy.add_chain("A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "HOH", None, None, None)
        .unwrap();
    for atom in small.atom_ids() {
        hierarchy
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .unwrap();
    }
    let mut editor = small.edit();
    *editor.hierarchy_mut() = hierarchy;
    let macromolecule = editor.finish().unwrap();
    let model = Model::from_molecule(&macromolecule, &conformer).unwrap();
    let mut potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();
    assert!(potential
        .evaluate(model.view())
        .unwrap()
        .energy()
        .is_finite());
}

#[test]
fn unresolved_or_counted_hydrogens_are_rejected_with_qualified_ids() {
    let mut atom = Atom::new(Element::from_symbol("C").unwrap());
    let mut graph = kekule::core::MoleculeEditor::new();
    let id = graph
        .add_atom(atom.clone())
        .expect("atom identifier capacity");
    let mut conformer = Conformer::new(kekule::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            id,
            kekule::units::Quantity::new(Point3::default(), kekule::units::ANGSTROM),
        )
        .unwrap();
    let molecule = graph.finish().expect("single atom molecule");
    let model = Model::from_molecule(&molecule, &conformer).unwrap();
    assert!(matches!(
        DreidingPotential::prepare(
            &model.shared_topology(),
            model.view(),
            DreidingPrepareOptions::default(),
        ),
        Err(DreidingPrepareError::UnresolvedImplicitHydrogens { atom })
            if atom == InstanceAtomId::new(MoleculeInstanceId::new(0), id)
    ));

    atom.hydrogens = HydrogenDeclaration::Fixed(1);
    let mut graph = kekule::core::MoleculeEditor::new();
    let id = graph.add_atom(atom).expect("atom identifier capacity");
    let mut conformer = Conformer::new(kekule::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            id,
            kekule::units::Quantity::new(Point3::default(), kekule::units::ANGSTROM),
        )
        .unwrap();
    let molecule = graph.finish().expect("single atom molecule");
    let model = Model::from_molecule(&molecule, &conformer).unwrap();
    assert!(matches!(
        DreidingPotential::prepare(
            &model.shared_topology(),
            model.view(),
            DreidingPrepareOptions::default(),
        ),
        Err(DreidingPrepareError::CountedHydrogens { .. })
    ));
}

#[test]
fn prepared_potential_uses_exact_shared_topology() {
    let (combined, combined_conf) = molecule(
        &["C", "C"],
        &[(0, 1, BondOrder::Single)],
        &[Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
    );
    let combined_model = Model::from_molecule(&combined, &combined_conf).unwrap();
    let mut potential = DreidingPotential::prepare(
        &combined_model.shared_topology(),
        combined_model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();

    let (one, one_conf) = molecule(&["C"], &[], &[Point3::new(0.0, 0.0, 0.0)]);
    let mut builder = Model::builder();
    builder.add_molecule(&one, &one_conf).unwrap();
    builder.add_molecule(&one, &one_conf).unwrap();
    let split_model = builder.build().unwrap();
    assert_eq!(
        potential.evaluate(split_model.view()),
        Err(PotentialError::IncompatibleTopology)
    );

    let mut singular = combined_model.clone();
    singular
        .set_position(
            InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(1)),
            kekule::units::Quantity::new(
                singular.positions().values().value()[0],
                kekule::units::ANGSTROM,
            ),
        )
        .unwrap();
    assert!(matches!(
        potential.evaluate(singular.view()),
        Err(PotentialError::InvalidGeometry { .. })
    ));
}
