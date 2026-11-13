use std::sync::Arc;

use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, Positions, StructureSet, Velocities};
use kekule::topology::TopologyBuilder;
use kekule::traj::periodic::{self, MoleculeImager, PeriodicError};
use kekule::traj::{FrameBuffer, FrameView};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::{self, dynamics, time_ordered, topology, Frame, Recorded};

fn cubic() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

fn frame(points: &[[f64; 3]], cell: Option<PeriodicCell>, index: usize) -> Frame {
    let points: Vec<_> = points
        .iter()
        .map(|p| Point3::new(p[0], p[1], p[2]))
        .collect();
    let count = points.len();
    let mut frame = Frame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap());
    frame.set_cell(cell);
    frame
        .set_time(Some(Quantity::new(index as f64, PICOSECOND)))
        .unwrap();
    frame.set_step(Some(index as u64 * 10));
    frame
        .set_velocities(Some(Velocities::zeros(count)))
        .unwrap();
    frame.set_forces(Some(Forces::zeros(count))).unwrap();
    frame
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("label").unwrap(),
            PropertyValue::Int(index as i64),
        )
        .unwrap();
    frame
}

fn points(buffer: &FrameBuffer) -> Vec<Point3> {
    buffer.positions().values().value().to_vec()
}

fn close(actual: &[Point3], expected: &[[f64; 3]]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        let e = Point3::new(e[0], e[1], e[2]);
        assert!((*a - e).norm() < 1.0e-12, "{actual:?} != {expected:?}");
    }
}

/// Images `set[frame]` onto `set[reference]` through a buffer, whose
/// positions and cell are the next ones `recorded`.
fn imaged(
    set: &StructureSet,
    frame: usize,
    reference: usize,
    recorded: &mut Recorded,
) -> Result<FrameBuffer, PeriodicError> {
    let imager = MoleculeImager::new(set.shared_topology());
    let mut buffer = FrameBuffer::new(set.shared_topology());
    buffer
        .copy_from(FrameView::from_geometry(set.get(frame).unwrap()).unwrap())
        .unwrap();
    imager.image_onto(frame, &mut buffer, &set.get(reference).unwrap())?;
    recorded.check(&(buffer.positions().clone(), buffer.cell().copied()));
    Ok(buffer)
}

#[test]
fn a_molecule_crossing_the_boundary_is_imaged_onto_references_before_and_after_it() {
    let mut recorded = Recorded::new(&[
        0x6fbf_099b_3ffa_e369, // frame 1 onto frame 0
        0x2d6f_a386_f5b4_b55c, // frame 0 onto frame 1
        0x0569_a700_e1df_54bc, // in memory
        0x6f0f_4ecd_deef_5786, // streamed, frame 0
        0x6fbf_099b_3ffa_e369, // streamed, frame 1
        0xf3e7_1cdd_903b_e37b, // frame 1 onto itself
    ]);
    // An ethane-like C–C molecule and a lone O atom. Between frames 0 and 1 the
    // C–C molecule crosses x = 1 and is wrapped back; the O atom stays put.
    let mut builder = TopologyBuilder::new();
    for smiles in ["CC", "O"] {
        let molecule = kekule::smiles::to_molecules(smiles).unwrap().pop().unwrap();
        builder.add_molecule(molecule).unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    let before = frame(
        &[[0.9, 0.5, 0.5], [1.05, 0.5, 0.5], [0.5, 0.5, 0.5]],
        Some(cubic()),
        0,
    );
    let after = frame(
        &[[0.0, 0.5, 0.5], [0.15, 0.5, 0.5], [0.5, 0.5, 0.5]],
        Some(cubic()),
        1,
    );
    let set = time_ordered(topology, [before, after]);

    // Reference before the frame: the molecule continues past x = 1.
    let forward = imaged(&set, 1, 0, &mut recorded).unwrap();
    close(
        &points(&forward),
        &[[1.0, 0.5, 0.5], [1.15, 0.5, 0.5], [0.5, 0.5, 0.5]],
    );
    // Reference after the frame: the molecule continues below x = 0.
    let backward = imaged(&set, 0, 1, &mut recorded).unwrap();
    close(
        &points(&backward),
        &[[-0.1, 0.5, 0.5], [0.05, 0.5, 0.5], [0.5, 0.5, 0.5]],
    );

    // Everything but positions is kept.
    for (buffer, index) in [(&forward, 1), (&backward, 0)] {
        let original = set.get(index).unwrap();
        assert_eq!(buffer.cell(), original.cell());
        assert_eq!(Some(buffer.frame_view().place()), original.place());
        assert_eq!(buffer.velocities(), original.velocities());
        assert_eq!(buffer.forces(), original.forces());
        assert_eq!(buffer.properties(), original.properties());
    }

    // The in-memory form matches the streaming form frame by frame.
    let reference = set.get(0).unwrap().to_model();
    let mut loaded = set.clone();
    periodic::image_onto(&mut loaded, &reference).unwrap();
    recorded.check(&dynamics(&loaded));
    for index in 0..loaded.len() {
        let streamed = imaged(&set, index, 0, &mut recorded).unwrap();
        assert_eq!(
            Frame::of(&streamed),
            Frame::of_set_sample(loaded.get(index).unwrap())
        );
    }
    // Imaging is idempotent: a frame imaged onto itself is unchanged.
    let itself = imaged(&set, 1, 1, &mut recorded).unwrap();
    assert_eq!(
        Frame::of(&itself),
        Frame::of_set_sample(set.get(1).unwrap())
    );
}

#[test]
fn whole_molecules_move_together_even_when_their_atoms_straddle_half_a_cell() {
    let mut recorded = Recorded::new(&[
        0xccac_9a7b_dc92_2df7, // frame 1 onto frame 0
    ]);
    let topology = topology(&["C", "C"], &[(0, 1)]);
    let reference = frame(&[[0.0, 0.5, 0.5], [0.15, 0.5, 0.5]], Some(cubic()), 0);
    // Atom 0 is 0.44 and atom 1 is 0.55 from its reference position: per-atom
    // images would move atom 1 alone and split the molecule.
    let moved = frame(&[[0.44, 0.5, 0.5], [0.7, 0.5, 0.5]], Some(cubic()), 1);
    let set = time_ordered(topology, [reference, moved]);
    close(
        &points(&imaged(&set, 1, 0, &mut recorded).unwrap()),
        &[[0.44, 0.5, 0.5], [0.7, 0.5, 0.5]],
    );
}

#[test]
fn triclinic_and_partially_periodic_cells_use_exact_lattice_images() {
    let mut recorded = Recorded::new(&[
        0x54b8_165a_0012_597d, // fully periodic
        0xe1ed_8b43_bc2b_7010, // nonperiodic c
    ]);
    let a = Vector3::new(1.0, 0.0, 0.0);
    let b = Vector3::new(0.5, 0.9, 0.0);
    let c = Vector3::new(0.2, 0.1, 1.0);
    let cell = |axes| PeriodicCell::new(Quantity::new([a, b, c], NANOMETER), axes).unwrap();
    let topology = topology(&["O"], &[]);
    let near = Point3::new(0.15, 0.12, 0.5);
    let shifted = near + b - a + c;
    let at = |p: Point3| [p.x, p.y, p.z];

    for (axes, expected) in [
        // Every lattice translation is removed.
        ([true; 3], near),
        // Along the nonperiodic c axis no image exists, so c remains.
        ([true, true, false], near + c),
    ] {
        let set = time_ordered(
            topology.clone(),
            [
                frame(&[at(Point3::new(0.1, 0.1, 0.5))], Some(cell(axes)), 0),
                frame(&[at(shifted)], Some(cell(axes)), 1),
            ],
        );
        close(
            &points(&imaged(&set, 1, 0, &mut recorded).unwrap()),
            &[at(expected)],
        );
    }
}

#[test]
fn imaging_onto_a_reference_requires_the_frame_cell_and_is_transactional() {
    let mut recorded = Recorded::new(&[
        0x2efe_b0c7_7a55_6272, // frame 0 onto frame 1
    ]);
    let topology = topology(&["C", "C"], &[(0, 1)]);
    let with_cell = frame(&[[0.0, 0.5, 0.5], [0.15, 0.5, 0.5]], Some(cubic()), 0);
    let without_cell = frame(&[[0.9, 0.5, 0.5], [1.05, 0.5, 0.5]], None, 1);
    let set = time_ordered(topology.clone(), [with_cell, without_cell]);
    let imager = MoleculeImager::new(topology.clone());

    // A frame without a cell has no lattice to image by.
    let mut buffer = FrameBuffer::new(topology.clone());
    buffer
        .copy_from(FrameView::from_geometry(set.get(1).unwrap()).unwrap())
        .unwrap();
    let snapshot = format!("{buffer:?}");
    assert_eq!(
        imager.image_onto(1, &mut buffer, &set.get(0).unwrap()),
        Err(PeriodicError::MissingCell { frame: 1 })
    );
    assert_eq!(format!("{buffer:?}"), snapshot);
    let mut loaded = set.clone();
    let reference = set.get(0).unwrap().to_model();
    assert_eq!(
        periodic::image_onto(&mut loaded, &reference),
        Err(PeriodicError::MissingCell { frame: 1 })
    );
    assert_eq!(format!("{loaded:?}"), format!("{set:?}"));

    // A reference without a cell is fine: only the frame's cell is used.
    close(
        &points(&imaged(&set, 0, 1, &mut recorded).unwrap()),
        &[[1.0, 0.5, 0.5], [1.15, 0.5, 0.5]],
    );

    // Foreign layouts fail without changing the buffer.
    let foreign = time_ordered(
        support::topology(&["C", "C"], &[(0, 1)]),
        [frame(&[[0.0; 3], [0.15, 0.0, 0.0]], Some(cubic()), 0)],
    );
    let mut buffer = FrameBuffer::new(topology);
    buffer
        .copy_from(FrameView::from_geometry(set.get(0).unwrap()).unwrap())
        .unwrap();
    let snapshot = format!("{buffer:?}");
    assert_eq!(
        imager.image_onto(0, &mut buffer, &foreign.get(0).unwrap()),
        Err(PeriodicError::ReferenceTopologyMismatch)
    );
    let mut foreign_buffer = FrameBuffer::new(foreign.shared_topology());
    assert_eq!(
        imager.image_onto(4, &mut foreign_buffer, &set.get(0).unwrap()),
        Err(PeriodicError::TopologyMismatch { frame: 4 })
    );
    assert_eq!(format!("{buffer:?}"), snapshot);
}
