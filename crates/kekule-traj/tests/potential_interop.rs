use kekule::core::{Atom, BondOrder, Element};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::modeling::potential::{
    HarmonicBondParameter, HarmonicBondPotential, Potential, PotentialError,
};
use kekule::structure::{Model, Positions};
use kekule::topology::{InstanceBondId, MoleculeInstanceId};
use kekule::units::{Quantity, ANGSTROM, KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM};
use kekule_traj::{FrameBuffer, TrajectoryFrame};

fn bonded_model() -> (Model, InstanceBondId) {
    let mut molecule = kekule::core::MoleculeEditor::new();
    let carbon = Element::from_symbol("C").unwrap();
    let first = molecule.add_atom(Atom::new(carbon)).unwrap();
    let second = molecule.add_atom(Atom::new(carbon)).unwrap();
    let bond = molecule.add_bond(first, second, BondOrder::Single).unwrap();
    let positions = Positions::new(Quantity::new(
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.1, 0.0, 0.0)],
        ANGSTROM,
    ))
    .unwrap();
    let molecule = molecule.finish().unwrap();
    let model = Model::from_molecule(&molecule, &positions).unwrap();
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
            Quantity::new(100.0, KILOJOULE_PER_MOLE_PER_SQUARE_ANGSTROM),
        )],
    )
    .unwrap();

    let frame = TrajectoryFrame::new(model.positions().clone(), model.topology().bond_count());
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
