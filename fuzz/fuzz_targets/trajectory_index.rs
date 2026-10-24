#![no_main]

//! Decodes untrusted frame indexes, and opens real files of every trajectory
//! format through checksummed indexes with fuzzed offsets.

use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::structure::TrajectoryFrame;
use kekule::topology::{Topology, TopologyBuilder};
use kekule_traj::io::{
    open_indexed_trajectory, open_indexed_trajectory_from_index, TrajectoryFrameIndex,
    TrajectoryIoLimits,
};
use kekule_traj::{SeekableTrajectoryReader, TrajectoryReader};
use libfuzzer_sys::fuzz_target;

/// One indexed file: its path, topology, canonical index, and every frame.
struct Case {
    path: PathBuf,
    topology: Arc<Topology>,
    index: TrajectoryFrameIndex,
    frames: Vec<TrajectoryFrame>,
}

fn topology(symbols: &[&str]) -> Arc<Topology> {
    let mut graph = MoleculeEditor::new();
    let mut previous = None;
    for symbol in symbols {
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
        .add_molecule_definition(molecule)
        .expect("definition");
    builder.add_instance(definition).expect("instance");
    Arc::new(builder.build().expect("topology"))
}

fn cases() -> &'static [Case] {
    static CASES: OnceLock<Vec<Case>> = OnceLock::new();
    CASES.get_or_init(|| {
        let fixtures: [(&str, &[u8], Arc<Topology>); 4] = [
            (
                "xyz",
                include_bytes!("../../crates/kekule-traj/tests/fixtures/ase-3.26.0-water.xyz"),
                topology(&["O", "H", "H"]),
            ),
            (
                "dcd",
                include_bytes!(
                    "../../crates/kekule-traj/tests/fixtures/mdanalysis-2.9.0-three-atoms.dcd"
                ),
                topology(&["C", "H", "O"]),
            ),
            (
                "trr",
                include_bytes!(
                    "../../crates/kekule-traj/tests/fixtures/mdanalysis-2.9.0-three-atoms.trr"
                ),
                topology(&["C", "H", "O"]),
            ),
            (
                "xtc",
                include_bytes!(
                    "../../crates/kekule-traj/tests/fixtures/mdanalysis-2.9.0-twelve-atoms.xtc"
                ),
                topology(&["C"; 12]),
            ),
        ];
        fixtures
            .into_iter()
            .map(|(extension, bytes, topology)| {
                let path = std::env::temp_dir().join(format!(
                    "kekule-fuzz-trajectory-index-{}.{extension}",
                    std::process::id()
                ));
                std::fs::write(&path, bytes).expect("fixture copy");
                let mut reader =
                    open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
                        .expect("fixture indexes");
                let mut buffer = reader.frame_buffer();
                let mut frames = Vec::new();
                while reader.read_next(&mut buffer).expect("fixture frame") {
                    frames.push(buffer.frame_view().payload().clone());
                }
                Case {
                    path,
                    topology,
                    index: reader.to_frame_index(),
                    frames,
                }
            })
            .collect()
    })
}

/// The documented FNV-1a 64 checksum of an encoded index.
fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

fuzz_target!(|data: &[u8]| {
    let limits = TrajectoryIoLimits::default()
        .with_max_frames(64)
        .with_max_index_entries(64)
        .with_max_index_bytes(512);

    // Any accepted encoding is canonical.
    if let Ok(index) = TrajectoryFrameIndex::from_bytes(data, &limits) {
        assert_eq!(index.to_bytes(), data);
    }

    // Fuzzed, strictly increasing offsets in an otherwise valid index of a
    // real file, each with the step and time of the file's frame at that
    // offset or else of another frame: every frame read through it must be
    // a frame of the file.
    let Some((&selector, deltas)) = data.split_first() else {
        return;
    };
    let case = &cases()[usize::from(selector) % cases().len()];
    let size = case.index.fingerprint().byte_count();
    let mut starts = Vec::new();
    let mut next = 0_u64;
    for delta in deltas.chunks_exact(2).take(16) {
        next += u64::from(u16::from_le_bytes([delta[0], delta[1]]));
        if next >= size {
            break;
        }
        starts.push(next);
        next += 1;
    }
    let offsets = case.index.frame_offsets();
    let frames = starts
        .iter()
        .enumerate()
        .map(|(forged, start)| {
            offsets
                .iter()
                .position(|offset| offset == start)
                .unwrap_or(forged % offsets.len())
        })
        .collect::<Vec<_>>();
    let mut forged = case.index.to_bytes()[..48].to_vec();
    if starts.is_empty() {
        // An empty index records no steps or times.
        forged[13] &= 1;
    }
    forged.extend_from_slice(&(starts.len() as u64).to_le_bytes());
    for start in &starts {
        forged.extend_from_slice(&start.to_le_bytes());
    }
    if let Some(steps) = case.index.steps() {
        for &frame in &frames {
            forged.extend_from_slice(&steps[frame].to_le_bytes());
        }
    }
    if let Some(times) = case.index.times() {
        for &frame in &frames {
            forged.extend_from_slice(&times.value()[frame].to_le_bytes());
        }
    }
    let checksum = fnv1a(&forged);
    forged.extend_from_slice(&checksum.to_le_bytes());
    let index = TrajectoryFrameIndex::from_bytes(&forged, &limits).expect("forged index decodes");
    let Ok(mut reader) = open_indexed_trajectory_from_index(
        &case.path,
        Arc::clone(&case.topology),
        index,
        Default::default(),
    ) else {
        return;
    };
    let mut buffer = reader.frame_buffer();
    for frame in 0..reader.frame_count().unwrap_or(0) {
        if reader.read_frame(frame, &mut buffer).is_ok() {
            assert!(case.frames.contains(buffer.frame_view().payload()));
        }
    }
});
