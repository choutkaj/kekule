use molecules::core::{Atom, AtomId, BondOrder, Conformer, Element, Molecule, Point3};
use molecules::modeling::potential::{Potential, PotentialError};
use molecules::modeling::{InstanceAtomId, Model, MoleculeInstanceId};
use molecules::small::SmallMolecule;
use molecules_dreiding::DreidingPotential;

#[test]
fn downstream_preparation_and_evaluation() {
    let mut graph = Molecule::new();
    let mut explicit_atom = |symbol: &str| {
        let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
        atom.no_implicit_hydrogens = true;
        graph.add_atom(atom)
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
    let mut conformer = Conformer::new(molecules::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            oxygen,
            molecules::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), molecules::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            first_hydrogen,
            molecules::units::Quantity::new(
                Point3::new(0.9575, 0.0, 0.0),
                molecules::units::ANGSTROM,
            ),
        )
        .unwrap();
    conformer
        .set_position(
            second_hydrogen,
            molecules::units::Quantity::new(
                Point3::new(-0.2399, 0.9272, 0.0),
                molecules::units::ANGSTROM,
            ),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).unwrap();
    let molecule = SmallMolecule::from_graph(graph);
    let model = Model::from_small_molecule(&molecule, conformer).unwrap();
    let independently_built = Model::from_small_molecule(&molecule, conformer).unwrap();
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    let evaluation = potential.evaluate(&model).unwrap();
    let oxygen = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(0));
    assert!(evaluation.energy().is_finite());
    assert_eq!(evaluation.gradient().len(), model.atom_count());
    assert!(potential.atom_type(oxygen).is_some());
    assert!(potential.partial_charge(oxygen).unwrap().is_finite());
    assert_eq!(
        potential.evaluate(&independently_built),
        Err(PotentialError::IncompatibleModel)
    );
}
