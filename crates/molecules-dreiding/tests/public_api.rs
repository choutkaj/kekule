use molecules::core::{Atom, AtomId, BondOrder, Conformer, Element, Molecule, Point3};
use molecules::modeling::MolecularModel;
use molecules::modeling::potential::Potential;
use molecules::small::SmallMolecule;
use molecules_dreiding::DreidingPotential;

#[test]
fn downstream_preparation_and_evaluation() {
    let mut graph = Molecule::new();
    let oxygen = graph.add_atom(Atom::new(Element::from_symbol("O").unwrap()));
    let first_hydrogen = graph.add_atom(Atom::new(Element::from_symbol("H").unwrap()));
    let second_hydrogen = graph.add_atom(Atom::new(Element::from_symbol("H").unwrap()));
    graph
        .add_bond(oxygen, first_hydrogen, BondOrder::Single)
        .unwrap();
    graph
        .add_bond(oxygen, second_hydrogen, BondOrder::Single)
        .unwrap();
    let mut conformer = Conformer::new();
    conformer.set_position(oxygen, Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(first_hydrogen, Point3::new(0.9575, 0.0, 0.0));
    conformer.set_position(second_hydrogen, Point3::new(-0.2399, 0.9272, 0.0));
    let conformer = graph.add_conformer(conformer);
    let model =
        MolecularModel::from_conformer(&SmallMolecule::from_graph(graph), conformer).unwrap();

    let mut potential = DreidingPotential::prepare(&model).unwrap();
    let evaluation = potential.evaluate(&model).unwrap();

    assert!(evaluation.energy().is_finite());
    assert_eq!(evaluation.gradient().len(), model.atom_count());
    assert!(potential.atom_type(AtomId::new(0)).is_some());
    assert!(
        potential
            .partial_charge(AtomId::new(0))
            .unwrap()
            .is_finite()
    );
}
