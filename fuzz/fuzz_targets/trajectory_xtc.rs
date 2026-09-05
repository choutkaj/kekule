#![no_main]

use std::io::Cursor;
use std::sync::{Arc, OnceLock};

use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::topology::{Topology, TopologyBuilder};
use kekule_traj::io::xtc::{XtcReadOptions, XtcReader};
use kekule_traj::io::TrajectoryIoLimits;
use kekule_traj::{FrameBuffer, TrajectoryReader};
use libfuzzer_sys::fuzz_target;

fn topology() -> &'static Arc<Topology> {
    static TOPOLOGY: OnceLock<Arc<Topology>> = OnceLock::new();
    TOPOLOGY.get_or_init(|| {
        let mut graph = MoleculeEditor::new();
        let mut previous = None;
        for _ in 0..12 {
            let atom = graph
                .add_atom(Atom::new(Element::from_symbol("C").expect("element")))
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
        builder.add_instance(definition).expect("instance");
        Arc::new(builder.build().expect("topology"))
    })
}

fuzz_target!(|data: &[u8]| {
    let topology = topology();
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
        Arc::clone(topology),
        XtcReadOptions::default()
            .with_limits(limits)
            .with_source_label("fuzz-input.xtc"),
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
