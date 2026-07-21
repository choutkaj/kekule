use molecular::bio::{MacroMolecule, SmcraAtomSiteMetadata, SmcraHierarchy};
use molecular::core::{Atom, AtomId, BondOrder, Conformer, Element, Molecule, Point3};
use molecular::modeling::potential::{Potential, PotentialError};
use molecular::modeling::{InstanceAtomId, Model, MoleculeInstanceId};
use molecular::small::SmallMolecule;

use crate::{DreidingPotential, DreidingPrepareError};

fn explicit_atom(symbol: &str) -> Atom {
    let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
    atom.no_implicit_hydrogens = true;
    atom
}

fn molecule(
    elements: &[&str],
    bonds: &[(usize, usize, BondOrder)],
    positions: &[Point3],
) -> (SmallMolecule, molecular::core::ConformerId) {
    let mut graph = Molecule::new();
    let atoms = elements
        .iter()
        .map(|symbol| graph.add_atom(explicit_atom(symbol)))
        .collect::<Vec<_>>();
    for &(a, b, order) in bonds {
        graph.add_bond(atoms[a], atoms[b], order).unwrap();
    }
    let mut conformer = Conformer::new(molecular::units::ANGSTROM).unwrap();
    for (&atom, &position) in atoms.iter().zip(positions) {
        conformer
            .set_position(
                atom,
                molecular::units::Quantity::new(position, molecular::units::ANGSTROM),
            )
            .unwrap();
    }
    let conformer = graph.add_conformer(conformer).expect("valid conformer");
    (SmallMolecule::from_graph(graph), conformer)
}

fn water(offset: f64) -> (SmallMolecule, molecular::core::ConformerId) {
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
    let model = Model::from_small_molecule(&water, conformer).unwrap();
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    let evaluation = potential.evaluate(&model).unwrap();
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
    let first_id = builder.add_small_molecule(&first, first_conf).unwrap();
    let second_id = builder.add_small_molecule(&second, second_conf).unwrap();
    let model = builder.build().unwrap();
    let potential = DreidingPotential::prepare(&model).unwrap();
    for instance in [first_id, second_id] {
        let total = (0..3)
            .map(|atom| {
                potential
                    .partial_charge(InstanceAtomId::new(instance, AtomId::new(atom)))
                    .unwrap()
                    .into_value()
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
fn preparation_maps_tombstoned_local_ids_to_dense_adjacency() {
    let mut graph = Molecule::new();
    let oxygen = graph.add_atom(explicit_atom("O"));
    let tombstone = graph.add_atom(explicit_atom("H"));
    let first_hydrogen = graph.add_atom(explicit_atom("H"));
    let second_hydrogen = graph.add_atom(explicit_atom("H"));
    graph.delete_atom(tombstone).unwrap();
    graph
        .add_bond(oxygen, first_hydrogen, BondOrder::Single)
        .unwrap();
    graph
        .add_bond(oxygen, second_hydrogen, BondOrder::Single)
        .unwrap();
    let mut conformer = Conformer::new(molecular::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            oxygen,
            molecular::units::Quantity::new(Point3::new(0.0, 0.0, 0.0), molecular::units::ANGSTROM),
        )
        .unwrap();
    conformer
        .set_position(
            first_hydrogen,
            molecular::units::Quantity::new(
                Point3::new(0.9575, 0.0, 0.0),
                molecular::units::ANGSTROM,
            ),
        )
        .unwrap();
    conformer
        .set_position(
            second_hydrogen,
            molecular::units::Quantity::new(
                Point3::new(-0.2399, 0.9272, 0.0),
                molecular::units::ANGSTROM,
            ),
        )
        .unwrap();
    let conformer = graph.add_conformer(conformer).expect("valid conformer");
    let model = Model::from_small_molecule(&SmallMolecule::from_graph(graph), conformer).unwrap();

    let potential = DreidingPotential::prepare(&model).unwrap();
    assert!(potential.nonbonded.is_empty());
}

#[test]
fn eligible_macro_molecules_are_supported() {
    let (small, conformer) = water(0.0);
    let mut hierarchy = SmcraHierarchy::new();
    let model = hierarchy.add_model("1");
    let chain = hierarchy.add_chain(model, "A", None).unwrap();
    let residue = hierarchy
        .add_residue(chain, "HOH", None, None, None)
        .unwrap();
    for atom in small.graph().atom_ids() {
        hierarchy
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .unwrap();
    }
    let macromolecule = MacroMolecule::try_from_parts(small.graph().clone(), hierarchy).unwrap();
    let model = Model::from_macro_molecule(&macromolecule, conformer).unwrap();
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    assert!(potential.evaluate(&model).unwrap().energy().is_finite());
}

#[test]
fn unresolved_or_counted_hydrogens_are_rejected_with_qualified_ids() {
    let mut atom = Atom::new(Element::from_symbol("C").unwrap());
    let mut graph = Molecule::new();
    let id = graph.add_atom(atom.clone());
    let mut conformer = Conformer::new(molecular::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            id,
            molecular::units::Quantity::new(Point3::default(), molecular::units::ANGSTROM),
        )
        .unwrap();
    let conformer_id = graph.add_conformer(conformer).expect("valid conformer");
    let model =
        Model::from_small_molecule(&SmallMolecule::from_graph(graph), conformer_id).unwrap();
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::UnresolvedImplicitHydrogens { atom })
            if atom == InstanceAtomId::new(MoleculeInstanceId::new(0), id)
    ));

    atom.no_implicit_hydrogens = true;
    atom.explicit_hydrogens = 1;
    let mut graph = Molecule::new();
    let id = graph.add_atom(atom);
    let mut conformer = Conformer::new(molecular::units::ANGSTROM).unwrap();
    conformer
        .set_position(
            id,
            molecular::units::Quantity::new(Point3::default(), molecular::units::ANGSTROM),
        )
        .unwrap();
    let conformer_id = graph.add_conformer(conformer).expect("valid conformer");
    let model =
        Model::from_small_molecule(&SmallMolecule::from_graph(graph), conformer_id).unwrap();
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::CountedHydrogens { .. })
    ));
}

#[test]
fn prepared_potential_uses_model_definition_identity() {
    let (combined, combined_conf) = molecule(
        &["C", "C"],
        &[],
        &[Point3::new(0.0, 0.0, 0.0), Point3::new(4.0, 0.0, 0.0)],
    );
    let combined_model = Model::from_small_molecule(&combined, combined_conf).unwrap();
    let mut potential = DreidingPotential::prepare(&combined_model).unwrap();

    let (one, one_conf) = molecule(&["C"], &[], &[Point3::new(0.0, 0.0, 0.0)]);
    let mut builder = Model::builder();
    builder.add_small_molecule(&one, one_conf).unwrap();
    builder.add_small_molecule(&one, one_conf).unwrap();
    let split_model = builder.build().unwrap();
    assert_eq!(
        potential.evaluate(&split_model),
        Err(PotentialError::IncompatibleModel)
    );

    let mut singular = combined_model.clone();
    singular
        .set_position(
            InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(1)),
            molecular::units::Quantity::new(singular.positions()[0], molecular::units::ANGSTROM),
        )
        .unwrap();
    assert!(matches!(
        potential.evaluate(&singular),
        Err(PotentialError::InvalidGeometry { .. })
    ));
}
