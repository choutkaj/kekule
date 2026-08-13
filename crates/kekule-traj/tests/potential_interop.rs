use kekule::core::{Atom, BondOrder, Conformer, Element, Molecule};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::modeling::potential::{
    HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError,
};
use kekule::small::SmallMolecule;
use kekule::structure::Model;
use kekule::topology::{InstanceBondId, MoleculeInstanceId};
use kekule::units::{Quantity, ANGSTROM, MODEL_FORCE_CONSTANT_UNIT};
use kekule_traj::{FrameBuffer, TrajectoryFrame};

fn bonded_model() -> (Model, InstanceBondId) {
    let mut molecule = Molecule::builder();
    let carbon = Element::from_symbol("C").unwrap();
    let first = molecule.add_atom(Atom::new(carbon)).unwrap();
    let second = molecule.add_atom(Atom::new(carbon)).unwrap();
    let bond = molecule.add_bond(first, second, BondOrder::Single).unwrap();
    let mut conformer = Conformer::new(ANGSTROM).unwrap();
    conformer
        .set_position(first, Quantity::new(Point3::new(0.0, 0.0, 0.0), ANGSTROM))
        .unwrap();
    conformer
        .set_position(second, Quantity::new(Point3::new(1.1, 0.0, 0.0), ANGSTROM))
        .unwrap();
    let mut molecule = molecule.build().unwrap();
    let conformer = molecule.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(molecule);
    let model = Model::from_small_molecule(&molecule, conformer).unwrap();
    (model, InstanceBondId::new(MoleculeInstanceId::new(0), bond))
}

#[test]
fn kekule_potentials_consume_trajectory_views_without_coordinate_copies() {
    let (model, bond) = bonded_model();
    let topology = model.shared_topology();
    let mut potential = HarmonicBondPotential::new(
        &topology,
        [HarmonicBondParameter::new(
            bond,
            Quantity::new(1.0, ANGSTROM),
            Quantity::new(100.0, MODEL_FORCE_CONSTANT_UNIT),
        )],
    )
    .unwrap();

    let frame = TrajectoryFrame::new(model.positions().clone());
    assert!(potential
        .evaluate(frame.view(&topology).unwrap().model_view())
        .is_ok());

    let mut buffer = FrameBuffer::new(topology);
    buffer.set_positions(model.positions().values()).unwrap();
    assert!(potential.evaluate(buffer.model_view()).is_ok());

    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    buffer.set_cell(Some(cell));
    assert_eq!(
        potential.evaluate(buffer.model_view()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );
}
