use kekule::alignment::PeriodicAlignmentPolicy;
use kekule::core::{Atom, BondOrder, Element, Molecule};
use kekule::geometry::Point3;
use kekule::small::SmallMolecule;
use kekule::structure::Positions;
use kekule::topology::{AtomSelection, MoleculeInstanceMetadata, Topology, TopologyBuilder};
use kekule::units::{Quantity, ANGSTROM};
use kekule_traj::analysis::{
    AlignedRmsdOptions, PeriodicRmsdPolicy, RmsdOptions, SuperpositionOptions,
};
use kekule_traj::{Trajectory, TrajectoryFrame};

fn topology() -> Arc<Topology> {
    let mut graph = Molecule::builder();
    let mut previous = None;
    for _ in 0..3 {
        let atom = graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(previous) = previous {
            graph.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    let molecule = SmallMolecule::from_graph(graph.build().unwrap());
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule).unwrap();
    builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    Arc::new(builder.build().unwrap())
}

fn frame(topology: &Arc<Topology>, points: [Point3; 3]) -> TrajectoryFrame {
    TrajectoryFrame::new(Positions::new(topology, Quantity::new(points, ANGSTROM)).unwrap())
}

#[test]
fn downstream_code_can_split_or_fuse_superposition_and_rmsd() {
    let topology = topology();
    let reference = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let moving = reference.map(|point| Point3::new(point.x + 4.0, point.y - 2.0, point.z + 1.0));
    let trajectory = Trajectory::from_frames(
        Arc::clone(&topology),
        [frame(&topology, reference), frame(&topology, moving)],
    )
    .unwrap();
    let selection =
        AtomSelection::from_atoms(&topology, topology.atom_ids().iter().copied()).unwrap();

    let direct = trajectory
        .rmsd_to_frame_with_options(
            0,
            &selection,
            RmsdOptions {
                periodic_policy: PeriodicRmsdPolicy::RejectPeriodic,
                ..RmsdOptions::default()
            },
        )
        .unwrap();
    assert!(direct.value()[1] > 4.0);

    let fused = trajectory
        .aligned_rmsd_to_frame_with_options(
            0,
            &selection,
            &selection,
            AlignedRmsdOptions {
                superposition: SuperpositionOptions {
                    periodic_policy: PeriodicAlignmentPolicy::RejectPeriodic,
                    ..SuperpositionOptions::default()
                },
                ..AlignedRmsdOptions::default()
            },
        )
        .unwrap();
    assert!(fused.value()[1] < 1.0e-12);

    let mut split = trajectory;
    let report = split.superpose_to_frame(0, &selection).unwrap();
    assert_eq!(report.alignments().len(), split.len());
    let measured = split.rmsd_to_frame(0, &selection).unwrap();
    assert!(measured.value()[1] < 1.0e-12);
}
use std::sync::Arc;
