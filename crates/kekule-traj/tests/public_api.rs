use std::fs;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use kekule::core::{Atom, Element, Molecule};
use kekule::small::SmallMolecule;
use kekule::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use kekule_traj::io::{
    open_indexed_trajectory, open_trajectory, FieldAvailability, RandomAccessCapability,
    TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryTopologyBinding,
};
use kekule_traj::{
    AtomOrderAssertion, FrameBuffer, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryReader,
};

fn topology() -> Topology {
    let mut graph = Molecule::builder();
    graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let molecule = SmallMolecule::from_graph(graph.build().unwrap());
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule).unwrap();
    builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    builder.build().unwrap()
}

fn binding(topology: &Topology) -> TrajectoryTopologyBinding {
    TrajectoryTopologyBinding::new(
        topology.clone(),
        AtomOrderAssertion::assert_file_uses_topology_order(topology),
    )
    .unwrap()
}

fn temporary_xyz() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "kekule-trajectory-public-api-{}-{nonce}.xyz",
        std::process::id()
    ))
}

#[test]
fn format_agnostic_public_api_opens_sequential_and_indexed_readers() {
    let path = temporary_xyz();
    fs::write(&path, b"1\npublic API\nC 1.0 2.0 3.0\n").unwrap();
    let topology = topology();
    let options = TrajectoryOpenOptions::default()
        .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Xyz));

    let (mut sequential, report) =
        open_trajectory(&path, binding(&topology), options.clone()).unwrap();
    assert_eq!(report.selected_format(), TrajectoryFormat::Xyz);
    assert_eq!(
        sequential.metadata().fields().positions,
        FieldAvailability::Required
    );
    assert_eq!(
        sequential.metadata().random_access(),
        RandomAccessCapability::SequentialOnly
    );
    let mut destination = FrameBuffer::new(topology.clone());
    assert!(sequential.read_next(&mut destination).unwrap());
    assert!(!sequential.read_next(&mut destination).unwrap());

    let (mut indexed, _) = open_indexed_trajectory(&path, binding(&topology), options).unwrap();
    assert_eq!(indexed.frame_count(), Some(1));
    assert_eq!(
        indexed.metadata().random_access(),
        RandomAccessCapability::Indexed
    );
    indexed.read_frame(0, &mut destination).unwrap();
    assert_eq!(
        destination.configuration().positions().values().value()[0].x,
        1.0
    );

    fs::remove_file(path).unwrap();
}
