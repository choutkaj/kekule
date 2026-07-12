use dreid_forge::{
    Atom as ForgeAtom, Bond as ForgeBond, BondOrder as ForgeBondOrder, BondPotential, ChargeMethod,
    ForgeConfig, QeqConfig, System, forge,
};
use molecules::core::{Atom, AtomId, AtomRadical, BondOrder, Conformer, Element, Molecule, Point3};
use molecules::modeling::potential::{Potential, Vector3};
use molecules::modeling::{MinimizationStatus, MinimizeOptions, MolecularModel, minimize};
use molecules::small::SmallMolecule;

use crate::{DreidingPotential, DreidingPrepareError};

fn atom(symbol: &str) -> Atom {
    let mut atom = Atom::new(Element::from_symbol(symbol).unwrap());
    atom.implicit_hydrogens = Some(0);
    atom
}

fn component(
    elements: &[&str],
    bonds: &[(usize, usize, BondOrder)],
    positions: &[Point3],
) -> (SmallMolecule, molecules::core::ConformerId) {
    assert_eq!(elements.len(), positions.len());
    let mut graph = Molecule::new();
    let atoms = elements
        .iter()
        .map(|symbol| graph.add_atom(atom(symbol)))
        .collect::<Vec<_>>();
    for &(first, second, order) in bonds {
        graph.add_bond(atoms[first], atoms[second], order).unwrap();
    }
    let mut conformer = Conformer::new();
    for (&atom, &position) in atoms.iter().zip(positions) {
        conformer.set_position(atom, position);
    }
    let conformer = graph.add_conformer(conformer);
    (SmallMolecule::from_graph(graph), conformer)
}

fn water(offset: [f64; 3]) -> (SmallMolecule, molecules::core::ConformerId) {
    let translate = |point: [f64; 3]| {
        Point3::new(
            point[0] + offset[0],
            point[1] + offset[1],
            point[2] + offset[2],
        )
    };
    component(
        &["O", "H", "H"],
        &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        &[
            translate([0.0, 0.0, 0.0]),
            translate([0.9575, 0.0, 0.0]),
            translate([-0.2399, 0.9272, 0.0]),
        ],
    )
}

fn ethanol() -> (SmallMolecule, molecules::core::ConformerId) {
    component(
        &["C", "C", "O", "H", "H", "H", "H", "H", "H"],
        &[
            (0, 1, BondOrder::Single),
            (1, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
            (0, 4, BondOrder::Single),
            (0, 5, BondOrder::Single),
            (1, 6, BondOrder::Single),
            (1, 7, BondOrder::Single),
            (2, 8, BondOrder::Single),
        ],
        &[
            Point3::new(-1.270, 0.248, 0.011),
            Point3::new(0.139, -0.308, -0.027),
            Point3::new(1.036, 0.789, 0.083),
            Point3::new(-1.317, 0.885, 0.883),
            Point3::new(-1.317, 0.885, -0.883),
            Point3::new(-2.030, -0.533, 0.113),
            Point3::new(0.358, -0.920, 0.876),
            Point3::new(0.358, -0.920, -0.876),
            Point3::new(1.939, 0.473, 0.191),
        ],
    )
}

fn formaldehyde() -> (SmallMolecule, molecules::core::ConformerId) {
    component(
        &["C", "O", "H", "H"],
        &[
            (0, 1, BondOrder::Double),
            (0, 2, BondOrder::Single),
            (0, 3, BondOrder::Single),
        ],
        &[
            Point3::new(0.0, 0.0, 0.04),
            Point3::new(1.22, 0.02, -0.01),
            Point3::new(-0.58, 0.92, 0.03),
            Point3::new(-0.62, -0.89, -0.02),
        ],
    )
}

fn localized_benzene() -> (SmallMolecule, molecules::core::ConformerId) {
    let mut graph = Molecule::new();
    let carbons = (0..6)
        .map(|_| {
            let id = graph.add_atom(atom("C"));
            graph.atom_mut(id).unwrap().aromatic = true;
            id
        })
        .collect::<Vec<_>>();
    let hydrogens = (0..6)
        .map(|_| graph.add_atom(atom("H")))
        .collect::<Vec<_>>();

    for index in 0..6 {
        let order = if index % 2 == 0 {
            BondOrder::Double
        } else {
            BondOrder::Single
        };
        let bond = graph
            .add_bond(carbons[index], carbons[(index + 1) % 6], order)
            .unwrap();
        graph.bond_mut(bond).unwrap().aromatic = true;
    }
    for index in 0..6 {
        graph
            .add_bond(carbons[index], hydrogens[index], BondOrder::Single)
            .unwrap();
    }

    let carbon_radius = 1.397;
    let hydrogen_radius = 2.48;
    let mut conformer = Conformer::new();
    for index in 0..6 {
        let angle = index as f64 * std::f64::consts::TAU / 6.0;
        conformer.set_position(
            carbons[index],
            Point3::new(
                carbon_radius * angle.cos(),
                carbon_radius * angle.sin(),
                0.0,
            ),
        );
        conformer.set_position(
            hydrogens[index],
            Point3::new(
                hydrogen_radius * angle.cos(),
                hydrogen_radius * angle.sin(),
                0.0,
            ),
        );
    }
    let conformer = graph.add_conformer(conformer);
    (SmallMolecule::from_graph(graph), conformer)
}

fn model_of(source: &SmallMolecule, conformer: molecules::core::ConformerId) -> MolecularModel {
    MolecularModel::from_conformer(source, conformer).unwrap()
}

#[test]
fn preparation_matches_direct_forge_for_single_component() {
    let (source, conformer) = water([0.0, 0.0, 0.0]);
    let model = model_of(&source, conformer);
    let original = model.clone();
    let prepared = DreidingPotential::prepare(&model).unwrap();

    let mut system = System::new();
    system
        .atoms
        .push(ForgeAtom::new("O".parse().unwrap(), [0.0, 0.0, 0.0]));
    system
        .atoms
        .push(ForgeAtom::new("H".parse().unwrap(), [0.9575, 0.0, 0.0]));
    system
        .atoms
        .push(ForgeAtom::new("H".parse().unwrap(), [-0.2399, 0.9272, 0.0]));
    system
        .bonds
        .push(ForgeBond::new(0, 1, ForgeBondOrder::Single));
    system
        .bonds
        .push(ForgeBond::new(0, 2, ForgeBondOrder::Single));
    let forged = forge(
        &system,
        &ForgeConfig {
            charge_method: ChargeMethod::Qeq(QeqConfig::default()),
            ..ForgeConfig::default()
        },
    )
    .unwrap();

    for index in 0..3 {
        let expected_type = &forged.atom_types[forged.atom_properties[index].type_idx];
        assert_eq!(
            prepared.atom_type(AtomId::new(index as u32)),
            Some(expected_type.as_str())
        );
        assert!(
            (prepared.partial_charge(AtomId::new(index as u32)).unwrap()
                - forged.atom_properties[index].charge)
                .abs()
                < 1.0e-12
        );
    }
    let BondPotential::Harmonic { k_half, .. } = forged.potentials.bonds[0] else {
        panic!("fixed configuration must produce harmonic bonds");
    };
    assert!((prepared.bonds[0].k_half - k_half * crate::prepare::KCAL_TO_KJ).abs() < 1.0e-12);
    assert_eq!(prepared.nonbonded.len(), 0, "water has only 1-2/1-3 pairs");
    assert_eq!(model, original);
}

#[test]
fn preparation_accepts_localized_aromatic_bonds() {
    let (source, conformer) = localized_benzene();
    let model = model_of(&source, conformer);
    let original = model.clone();

    let prepared = DreidingPotential::prepare(&model).unwrap();

    for index in 0..6 {
        assert_eq!(prepared.atom_type(AtomId::new(index)), Some("C_R"));
    }
    let localized_orders = model
        .topology()
        .bonds()
        .take(6)
        .map(|(_, bond)| (bond.order, bond.aromatic))
        .collect::<Vec<_>>();
    assert_eq!(
        localized_orders,
        vec![
            (BondOrder::Double, true),
            (BondOrder::Single, true),
            (BondOrder::Double, true),
            (BondOrder::Single, true),
            (BondOrder::Double, true),
            (BondOrder::Single, true),
        ]
    );
    assert_eq!(model, original);
}

#[test]
fn component_qeq_charges_are_isolated_and_cross_pairs_are_present() {
    let (first, first_conformer) = water([0.0, 0.0, 0.0]);
    let (second, second_conformer) = water([2.8, 0.1, 0.2]);
    let mut builder = MolecularModel::builder();
    builder.add_component(&first, first_conformer).unwrap();
    builder.add_component(&second, second_conformer).unwrap();
    let model = builder.build().unwrap();
    let prepared = DreidingPotential::prepare(&model).unwrap();

    for atoms in [[0_u32, 1, 2], [3, 4, 5]] {
        let total = atoms
            .into_iter()
            .map(|atom| prepared.partial_charge(AtomId::new(atom)).unwrap())
            .sum::<f64>();
        assert!(total.abs() < 1.0e-9);
    }
    assert_eq!(prepared.nonbonded.len(), 9);
    assert!(!prepared.hydrogen_bonds.is_empty());
}

#[test]
fn component_qeq_respects_nonzero_formal_charge_sums() {
    let (mut sodium, sodium_conformer) = component(&["Na"], &[], &[Point3::new(0.0, 0.0, 0.0)]);
    sodium
        .graph_mut()
        .atom_mut(AtomId::new(0))
        .unwrap()
        .formal_charge = 1;
    let (mut chloride, chloride_conformer) = component(&["Cl"], &[], &[Point3::new(3.0, 0.0, 0.0)]);
    chloride
        .graph_mut()
        .atom_mut(AtomId::new(0))
        .unwrap()
        .formal_charge = -1;
    let mut builder = MolecularModel::builder();
    builder.add_component(&sodium, sodium_conformer).unwrap();
    builder
        .add_component(&chloride, chloride_conformer)
        .unwrap();
    let model = builder.build().unwrap();
    let prepared = DreidingPotential::prepare(&model).unwrap();

    assert!((prepared.partial_charge(AtomId::new(0)).unwrap() - 1.0).abs() < 1.0e-12);
    assert!((prepared.partial_charge(AtomId::new(1)).unwrap() + 1.0).abs() < 1.0e-12);
    assert_eq!(prepared.nonbonded.len(), 1);
}

#[test]
fn complete_ethanol_gradient_matches_central_difference() {
    let (source, conformer) = ethanol();
    let model = model_of(&source, conformer);
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    assert!(!potential.bonds.is_empty());
    assert!(!potential.angles.is_empty());
    assert!(!potential.torsions.is_empty());
    assert!(!potential.nonbonded.is_empty());
    assert_gradient_matches(&model, &mut potential, 3.0e-3);
    assert_zero_total_gradient(potential.evaluate(&model).unwrap().gradient(), 1.0e-7);
}

#[test]
fn inversion_and_hydrogen_bond_gradients_match_central_difference() {
    let (source, conformer) = formaldehyde();
    let model = model_of(&source, conformer);
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    assert!(!potential.inversions.is_empty());
    assert_gradient_matches(&model, &mut potential, 3.0e-3);

    let (first, first_conformer) = water([0.0, 0.0, 0.0]);
    let (second, second_conformer) = water([2.75, 0.2, 0.15]);
    let mut builder = MolecularModel::builder();
    builder.add_component(&first, first_conformer).unwrap();
    builder.add_component(&second, second_conformer).unwrap();
    let model = builder.build().unwrap();
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    assert!(!potential.hydrogen_bonds.is_empty());
    assert_gradient_matches(&model, &mut potential, 5.0e-3);
}

#[test]
fn preparation_rejects_unrepresented_hydrogens_radicals_and_bond_orders() {
    let (mut source, conformer) = water([0.0, 0.0, 0.0]);
    source
        .graph_mut()
        .atom_mut(AtomId::new(0))
        .unwrap()
        .implicit_hydrogens = None;
    let model = model_of(&source, conformer);
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::UnresolvedImplicitHydrogens { .. })
    ));

    let (mut source, conformer) = water([0.0, 0.0, 0.0]);
    source
        .graph_mut()
        .atom_mut(AtomId::new(0))
        .unwrap()
        .implicit_hydrogens = Some(1);
    let model = model_of(&source, conformer);
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::CountedHydrogens { .. })
    ));

    let (mut source, conformer) = water([0.0, 0.0, 0.0]);
    source.graph_mut().atom_mut(AtomId::new(0)).unwrap().radical = Some(AtomRadical::Doublet);
    let model = model_of(&source, conformer);
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::RadicalAtom { .. })
    ));

    let (source, conformer) = component(
        &["C", "C"],
        &[(0, 1, BondOrder::Dative)],
        &[Point3::new(0.0, 0.0, 0.0), Point3::new(1.4, 0.0, 0.0)],
    );
    let model = model_of(&source, conformer);
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::UnsupportedBondOrder { .. })
    ));

    let (mut source, conformer) = component(
        &["C", "C"],
        &[(0, 1, BondOrder::Aromatic)],
        &[Point3::new(0.0, 0.0, 0.0), Point3::new(1.4, 0.0, 0.0)],
    );
    source
        .graph_mut()
        .bond_mut(molecules::core::BondId::new(0))
        .unwrap()
        .aromatic = false;
    let model = model_of(&source, conformer);
    assert!(matches!(
        DreidingPotential::prepare(&model),
        Err(DreidingPrepareError::InconsistentAromaticBond { .. })
    ));
}

#[test]
fn no_implicit_hydrogens_assertion_resolves_an_unset_count() {
    let (mut source, conformer) = water([0.0, 0.0, 0.0]);
    {
        let mut oxygen = source.graph_mut().atom_mut(AtomId::new(0)).unwrap();
        oxygen.implicit_hydrogens = None;
        oxygen.no_implicit_hydrogens = true;
    }
    let model = model_of(&source, conformer);

    DreidingPotential::prepare(&model).unwrap();
}

#[test]
fn evaluation_rejects_incompatible_topology_and_singular_geometry() {
    let (source, conformer) = water([0.0, 0.0, 0.0]);
    let model = model_of(&source, conformer);
    let mut potential = DreidingPotential::prepare(&model).unwrap();

    let (other, other_conformer) = formaldehyde();
    let other = model_of(&other, other_conformer);
    assert!(matches!(
        potential.evaluate(&other),
        Err(molecules::modeling::potential::PotentialError::IncompatibleModel)
    ));

    let mut singular = model.clone();
    singular
        .set_position(AtomId::new(1), singular.positions()[0])
        .unwrap();
    assert!(matches!(
        potential.evaluate(&singular),
        Err(molecules::modeling::potential::PotentialError::InvalidGeometry { .. })
    ));
}

#[test]
fn minimization_decreases_energy_without_mutating_input() {
    let (source, conformer) = component(
        &["O", "H", "H"],
        &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.25, 0.0, 0.0),
            Point3::new(-0.15, 1.18, 0.12),
        ],
    );
    let model = model_of(&source, conformer);
    let original = model.clone();
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    let result = minimize(
        &model,
        &mut potential,
        MinimizeOptions {
            max_iterations: 400,
            ..MinimizeOptions::default()
        },
    )
    .unwrap();

    assert!(result.final_energy < result.initial_energy);
    assert!(matches!(
        result.status,
        MinimizationStatus::Converged | MinimizationStatus::MaxIterations
    ));
    assert_eq!(model, original);
}

#[test]
fn multi_component_minimization_uses_fixed_component_charges() {
    let (first, first_conformer) = water([0.0, 0.0, 0.0]);
    let (second, second_conformer) = water([3.4, 0.3, 0.2]);
    let mut builder = MolecularModel::builder();
    builder.add_component(&first, first_conformer).unwrap();
    builder.add_component(&second, second_conformer).unwrap();
    let model = builder.build().unwrap();
    let original = model.clone();
    let mut potential = DreidingPotential::prepare(&model).unwrap();
    let charges = (0..model.atom_count())
        .map(|index| potential.partial_charge(AtomId::new(index as u32)).unwrap())
        .collect::<Vec<_>>();
    let result = minimize(
        &model,
        &mut potential,
        MinimizeOptions {
            max_iterations: 80,
            ..MinimizeOptions::default()
        },
    )
    .unwrap();

    assert!(result.final_energy < result.initial_energy);
    assert_eq!(model, original);
    for (index, charge) in charges.into_iter().enumerate() {
        assert_eq!(
            potential.partial_charge(AtomId::new(index as u32)),
            Some(charge)
        );
    }
}

fn assert_gradient_matches(
    model: &MolecularModel,
    potential: &mut DreidingPotential,
    tolerance: f64,
) {
    let analytic = potential.evaluate(model).unwrap();
    let epsilon = 1.0e-6;
    for atom in 0..model.atom_count() {
        for axis in 0..3 {
            let mut plus = model.clone();
            let mut minus = model.clone();
            let mut plus_point = plus.positions()[atom];
            let mut minus_point = minus.positions()[atom];
            coordinate_mut(&mut plus_point, axis, epsilon);
            coordinate_mut(&mut minus_point, axis, -epsilon);
            plus.set_position(AtomId::new(atom as u32), plus_point)
                .unwrap();
            minus
                .set_position(AtomId::new(atom as u32), minus_point)
                .unwrap();
            let numerical = (potential.evaluate(&plus).unwrap().energy()
                - potential.evaluate(&minus).unwrap().energy())
                / (2.0 * epsilon);
            let expected = coordinate(analytic.gradient()[atom], axis);
            assert!(
                (expected - numerical).abs() < tolerance,
                "atom {atom} axis {axis}: analytic={expected}, numerical={numerical}"
            );
        }
    }
}

fn assert_zero_total_gradient(gradient: &[Vector3], tolerance: f64) {
    let total = gradient.iter().fold(Vector3::zero(), |mut total, value| {
        total.x += value.x;
        total.y += value.y;
        total.z += value.z;
        total
    });
    assert!(total.norm() < tolerance, "net gradient is {total:?}");
}

fn coordinate_mut(point: &mut Point3, axis: usize, delta: f64) {
    match axis {
        0 => point.x += delta,
        1 => point.y += delta,
        2 => point.z += delta,
        _ => unreachable!(),
    }
}

fn coordinate(vector: Vector3, axis: usize) -> f64 {
    match axis {
        0 => vector.x,
        1 => vector.y,
        2 => vector.z,
        _ => unreachable!(),
    }
}
