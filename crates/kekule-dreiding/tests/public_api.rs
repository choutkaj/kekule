use kekule::core::{Atom, AtomId, BondOrder, Conformer, Element, Molecule};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::modeling::potential::{Potential, PotentialError};
use kekule::small::SmallMolecule;
use kekule::structure::Model;
use kekule::topology::{InstanceAtomId, MoleculeInstanceId};
use kekule_dreiding::{DreidingPotential, DreidingPrepareError, DreidingPrepareOptions};

#[test]
fn downstream_preparation_and_evaluation() {
    let mut graph = Molecule::new();
    let mut explicit_atom = |symbol: &str| {
        let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
        atom.no_implicit_hydrogens = true;
        graph.add_atom(atom).expect("atom identifier capacity")
    };
    let oxygen = explicit_atom("O");
    let first_hydrogen = explicit_atom("H");
    let second_hydrogen = explicit_atom("H");
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
    let conformer = graph.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(graph);
    let model = Model::from_small_molecule(&molecule, conformer).unwrap();
    let independently_built = Model::from_small_molecule(&molecule, conformer).unwrap();
    let mut periodic = model.clone();
    periodic.set_cell(Some(
        PeriodicCell::orthorhombic(
            kekule::units::Quantity::new(Vector3::new(10.0, 10.0, 10.0), kekule::units::ANGSTROM),
            [true; 3],
        )
        .unwrap(),
    ));
    assert!(matches!(
        DreidingPotential::prepare(
            periodic.topology(),
            periodic.view(),
            DreidingPrepareOptions::default(),
        ),
        Err(DreidingPrepareError::UnsupportedPeriodicCell)
    ));
    let mut potential = DreidingPotential::prepare(
        model.topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )
    .unwrap();
    let evaluation = potential.evaluate(model.view()).unwrap();
    let oxygen = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(0));
    assert!(evaluation.energy().is_finite());
    assert_eq!(evaluation.gradient().len(), model.atom_count());
    assert!(potential.atom_type(oxygen).is_some());
    assert!(potential.partial_charge(oxygen).unwrap().is_finite());
    assert_eq!(
        potential.evaluate(periodic.view()),
        Err(PotentialError::UnsupportedPeriodicCell)
    );
    assert_eq!(
        potential.evaluate(independently_built.view()),
        Err(PotentialError::IncompatibleTopology)
    );
}
