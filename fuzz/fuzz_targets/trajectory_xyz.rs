#![no_main]

use std::io::{BufReader, Cursor};
use std::sync::{Arc, OnceLock};

use libfuzzer_sys::fuzz_target;
use kekule::core::{Atom, BondOrder, Element, Molecule};
use kekule::small::SmallMolecule;
use kekule::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use kekule_traj::{AtomOrderAssertion, FrameBuffer, TrajectoryReader};
use kekule_traj::io::xyz::{XyzReadOptions, XyzReader};
use kekule_traj::io::{TrajectoryIoLimits, TrajectoryTopologyBinding};

fn topology() -> &'static Arc<Topology> {
    static TOPOLOGY: OnceLock<Arc<Topology>> = OnceLock::new();
    TOPOLOGY.get_or_init(|| {
        let mut graph = Molecule::builder();
        let mut previous = None;
        for symbol in ["C", "H"] {
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
        let molecule = SmallMolecule::from_graph(graph.build().expect("molecule"));
        let mut builder = TopologyBuilder::new();
        let definition = builder
            .add_small_molecule_definition(&molecule)
            .expect("definition");
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
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
    let mut limits = TrajectoryIoLimits::default();
    limits.max_atoms = 8;
    limits.max_frames = 8;
    limits.max_frame_bytes = 4096;
    limits.max_scratch_bytes = 4096;
    limits.max_text_line_bytes = 512;
    limits.max_comment_bytes = 512;

    if let Ok(mut reader) = XyzReader::new(
        BufReader::new(Cursor::new(data)),
        binding,
        XyzReadOptions::default(),
        limits,
        "fuzz-input.xyz",
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
