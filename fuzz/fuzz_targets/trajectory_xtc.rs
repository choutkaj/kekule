#![no_main]

use std::io::Cursor;
use std::sync::OnceLock;

use libfuzzer_sys::fuzz_target;
use molecular::core::{Atom, Element, Molecule};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use molecular::trajectory::{AtomOrderAssertion, FrameBuffer, TrajectoryReader};
use molecular_trajectory_io::xtc::{XtcReadOptions, XtcReader};
use molecular_trajectory_io::{TrajectoryIoLimits, TrajectoryTopologyBinding};

fn topology() -> &'static Topology {
    static TOPOLOGY: OnceLock<Topology> = OnceLock::new();
    TOPOLOGY.get_or_init(|| {
        let mut graph = Molecule::new();
        for _ in 0..12 {
            graph
                .add_atom(Atom::new(Element::from_symbol("C").expect("element")))
                .expect("atom");
        }
        let molecule = SmallMolecule::from_graph(graph);
        let mut builder = TopologyBuilder::new();
        let definition = builder
            .add_small_molecule_definition(&molecule)
            .expect("definition");
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .expect("instance");
        builder.build().expect("topology")
    })
}

fuzz_target!(|data: &[u8]| {
    let topology = topology();
    let binding = TrajectoryTopologyBinding::new(
        topology.clone(),
        AtomOrderAssertion::assert_file_uses_topology_order(topology),
    )
    .expect("binding");
    let limits = TrajectoryIoLimits {
        max_atoms: 16,
        max_frames: 8,
        max_frame_bytes: 8192,
        max_record_bytes: 8192,
        max_scratch_bytes: 8192,
        max_index_entries: 8,
        max_index_bytes: 64,
        ..TrajectoryIoLimits::default()
    };
    if let Ok(mut reader) = XtcReader::new(
        Cursor::new(data),
        binding,
        XtcReadOptions::default(),
        limits,
        "fuzz-input.xtc",
    ) {
        let mut buffer = FrameBuffer::new(topology.clone());
        for _ in 0..8 {
            match reader.read_next(&mut buffer) {
                Ok(true) => {}
                Ok(false) | Err(_) => break,
            }
        }
    }
});
