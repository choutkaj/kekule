use kekule::{
    alignment::{
        kabsch, kabsch_with_options, AlignmentError, AlignmentWeighting, KabschOptions,
        PeriodicAlignmentPolicy, RigidAlignment,
    },
    core::{Atom, Element, Molecule},
    geometry::Point3,
    small::SmallMolecule,
    structure::{Configuration, Model, Positions},
    topology::{AtomSelection, MoleculeInstanceMetadata, TopologyBuilder},
    units::{Quantity, ANGSTROM, MODEL_LENGTH_UNIT},
};

#[test]
fn focused_alignment_facade_is_downstream_usable() -> Result<(), Box<dyn std::error::Error>> {
    let mut graph = Molecule::new();
    for _ in 0..4 {
        graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()))?;
    }
    let molecule = SmallMolecule::from_graph(graph);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule)?;
    builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
    let topology = builder.build()?;
    let moving_points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 2.0, 0.0),
        Point3::new(0.0, 0.0, 3.0),
    ];
    let reference_points =
        moving_points.map(|point| Point3::new(-point.y + 4.0, point.x - 2.0, point.z + 0.5));
    let moving = Model::new(
        topology.clone(),
        Configuration::new(Positions::new(
            &topology,
            Quantity::new(moving_points, ANGSTROM),
        )?),
    )?;
    let reference = Model::new(
        topology.clone(),
        Configuration::new(Positions::new(
            &topology,
            Quantity::new(reference_points, ANGSTROM),
        )?),
    )?;
    let selection = AtomSelection::from_atoms(&topology, topology.atom_ids().iter().copied())?;

    let result: RigidAlignment = kabsch(moving.view(), reference.view(), &selection)?;
    assert_eq!(result.selected_atom_count(), 4);
    assert_eq!(result.rmsd().unit(), MODEL_LENGTH_UNIT);
    for (moving, reference) in moving_points.into_iter().zip(reference_points) {
        let aligned = result.transform().transform_point(moving);
        assert!((aligned.x - reference.x).abs() < 1.0e-12);
        assert!((aligned.y - reference.y).abs() < 1.0e-12);
        assert!((aligned.z - reference.z).abs() < 1.0e-12);
    }

    let weights = [1.0; 4];
    let weighted = kabsch_with_options(
        moving.view(),
        reference.view(),
        &selection,
        KabschOptions {
            weighting: AlignmentWeighting::Explicit(&weights),
            periodic_policy: PeriodicAlignmentPolicy::RejectPeriodic,
        },
    )?;
    assert!(weighted.rmsd().into_value() < 1.0e-12);

    let error = kabsch_with_options(
        moving.view(),
        reference.view(),
        &selection,
        KabschOptions {
            weighting: AlignmentWeighting::Explicit(&weights[..3]),
            ..KabschOptions::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        error,
        AlignmentError::WeightCountMismatch {
            expected: 4,
            actual: 3
        }
    ));
    Ok(())
}
