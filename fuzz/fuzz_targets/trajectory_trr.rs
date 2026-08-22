#![no_main]

use std::io::Cursor;
use std::sync::{Arc, OnceLock};

use libfuzzer_sys::fuzz_target;
use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::topology::{Topology, TopologyBuilder};
use kekule_traj::{AtomOrderAssertion, FrameBuffer, TrajectoryReader};
use kekule_traj::io::trr::{TrrReadOptions, TrrReader};
use kekule_traj::io::{TrajectoryIoLimits, TrajectoryTopologyBinding};

fn topology() -> &'static Arc<Topology> {
    static TOPOLOGY: OnceLock<Arc<Topology>> = OnceLock::new();
    TOPOLOGY.get_or_init(|| {
        let mut graph = MoleculeEditor::new();
        let mut previous = None;
        for symbol in ["C", "H", "O"] {
            let atom = graph
                .add_atom(Atom::new(Element::from_symbol(symbol).expect("element")))
                .expect("atom");
            if let Some(previous) = previous {
                graph
                    .add_bond(previous, atom, BondOrder::Single)
                    .expect("bond");
            }
            previous = Some(atom);
        }
        let molecule = graph.finish().expect("molecule");
        let mut builder = TopologyBuilder::new();
        let definition = builder
            .add_molecule_definition(&molecule)
            .expect("definition");
        builder
            .add_instance(definition)
            .expect("instance");
        Arc::new(builder.build().expect("topology"))
    })
}

fuzz_target!(|data: &[u8]| {
    let topology = topology();
    let binding = TrajectoryTopologyBinding::new(
        Arc::clone(topology),
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
        let mut buffer = FrameBuffer::new(Arc::clone(topology));
        for _ in 0..8 {
            match reader.read_next(&mut buffer) {
                Ok(true) => {}
                Ok(false) | Err(_) => break,
            }
        }
    }
});
