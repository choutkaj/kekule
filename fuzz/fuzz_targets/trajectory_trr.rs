#![no_main]

use std::io::Cursor;
use std::sync::OnceLock;

use libfuzzer_sys::fuzz_target;
use molecular::core::{Atom, Element, Molecule};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use molecular::trajectory::{AtomOrderAssertion, FrameBuffer, TrajectoryReader};
use molecular_trajectory_io::trr::{TrrReadOptions, TrrReader};
use molecular_trajectory_io::{TrajectoryIoLimits, TrajectoryTopologyBinding};

fn topology() -> &'static Topology {
    static TOPOLOGY: OnceLock<Topology> = OnceLock::new();
    TOPOLOGY.get_or_init(|| {
        let mut graph = Molecule::new();
        for symbol in ["C", "H", "O"] {
            graph
                .add_atom(Atom::new(Element::from_symbol(symbol).expect("element")))
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
        max_atoms: 8,
        max_frames: 8,
        max_frame_bytes: 4096,
        max_record_bytes: 4096,
        max_scratch_bytes: 4096,
        max_index_entries: 8,
        max_index_bytes: 64,
        max_text_line_bytes: 256,
        ..TrajectoryIoLimits::default()
    };
    if let Ok(mut reader) = TrrReader::new(
        Cursor::new(data),
        binding,
        TrrReadOptions::default(),
        limits,
        "fuzz-input.trr",
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
