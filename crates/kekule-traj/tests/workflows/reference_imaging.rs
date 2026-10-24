use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, Positions, Trajectory, TrajectoryFrame, Velocities};
use kekule::topology::TopologyBuilder;
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::periodic::{self, MoleculeImager, PeriodicError};
use kekule_traj::FrameBuffer;

use crate::support;
use support::topology;

fn cubic() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

fn frame(points: &[[f64; 3]], cell: Option<PeriodicCell>, index: usize) -> TrajectoryFrame {
    let points: Vec<_> = points
        .iter()
        .map(|p| Point3::new(p[0], p[1], p[2]))
        .collect();
    let count = points.len();
    let mut frame = TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap());
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

/// Images `trajectory[frame]` onto `trajectory[reference]` through a buffer.
fn imaged(
    trajectory: &Trajectory,
    frame: usize,
    reference: usize,
) -> Result<FrameBuffer, PeriodicError> {
    let imager = MoleculeImager::new(trajectory.shared_topology());
    let mut buffer = FrameBuffer::new(trajectory.shared_topology());
    buffer.copy_from(trajectory.get(frame).unwrap()).unwrap();
    imager.image_onto(frame, &mut buffer, &trajectory.get(reference).unwrap())?;
    Ok(buffer)
}

#[test]
fn a_molecule_crossing_the_boundary_is_imaged_onto_references_before_and_after_it() {
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
    let trajectory = Trajectory::from_items(topology, [before, after]).unwrap();

    // Reference before the frame: the molecule continues past x = 1.
    let forward = imaged(&trajectory, 1, 0).unwrap();
    close(
        &points(&forward),
        &[[1.0, 0.5, 0.5], [1.15, 0.5, 0.5], [0.5, 0.5, 0.5]],
    );
    // Reference after the frame: the molecule continues below x = 0.
    let backward = imaged(&trajectory, 0, 1).unwrap();
    close(
        &points(&backward),
        &[[-0.1, 0.5, 0.5], [0.05, 0.5, 0.5], [0.5, 0.5, 0.5]],
    );

    // Everything but positions is kept.
    for (buffer, index) in [(&forward, 1), (&backward, 0)] {
        let original = trajectory.get(index).unwrap();
        assert_eq!(buffer.cell(), original.cell());
        assert_eq!(buffer.time(), original.time());
        assert_eq!(buffer.step(), original.step());
        assert_eq!(buffer.velocities(), original.velocities());
        assert_eq!(buffer.forces(), original.forces());
        assert_eq!(buffer.properties(), original.properties());
    }

    // The in-memory form matches the streaming form frame by frame.
    let reference = trajectory.get(0).unwrap().to_model();
    let mut loaded = trajectory.clone();
    periodic::image_onto(&mut loaded, &reference).unwrap();
    for index in 0..loaded.len() {
        let streamed = imaged(&trajectory, index, 0).unwrap();
        assert_eq!(
            streamed.frame_view().payload(),
            loaded.get(index).unwrap().payload()
        );
    }
    // Imaging is idempotent: a frame imaged onto itself is unchanged.
    let itself = imaged(&trajectory, 1, 1).unwrap();
    assert_eq!(
        itself.frame_view().payload(),
        trajectory.get(1).unwrap().payload()
    );
}

#[test]
fn whole_molecules_move_together_even_when_their_atoms_straddle_half_a_cell() {
    let topology = topology(&["C", "C"], &[(0, 1)]);
    let reference = frame(&[[0.0, 0.5, 0.5], [0.15, 0.5, 0.5]], Some(cubic()), 0);
    // Atom 0 is 0.44 and atom 1 is 0.55 from its reference position: per-atom
    // images would move atom 1 alone and split the molecule.
    let moved = frame(&[[0.44, 0.5, 0.5], [0.7, 0.5, 0.5]], Some(cubic()), 1);
    let trajectory = Trajectory::from_items(topology, [reference, moved]).unwrap();
    close(
        &points(&imaged(&trajectory, 1, 0).unwrap()),
        &[[0.44, 0.5, 0.5], [0.7, 0.5, 0.5]],
    );
}

#[test]
fn triclinic_and_partially_periodic_cells_use_exact_lattice_images() {
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
        let trajectory = Trajectory::from_items(
            topology.clone(),
            [
                frame(&[at(Point3::new(0.1, 0.1, 0.5))], Some(cell(axes)), 0),
                frame(&[at(shifted)], Some(cell(axes)), 1),
            ],
        )
        .unwrap();
        close(
            &points(&imaged(&trajectory, 1, 0).unwrap()),
            &[at(expected)],
        );
    }
}

#[test]
fn imaging_onto_a_reference_requires_the_frame_cell_and_is_transactional() {
    let topology = topology(&["C", "C"], &[(0, 1)]);
    let with_cell = frame(&[[0.0, 0.5, 0.5], [0.15, 0.5, 0.5]], Some(cubic()), 0);
    let without_cell = frame(&[[0.9, 0.5, 0.5], [1.05, 0.5, 0.5]], None, 1);
    let trajectory = Trajectory::from_items(topology.clone(), [with_cell, without_cell]).unwrap();
    let imager = MoleculeImager::new(topology.clone());

    // A frame without a cell has no lattice to image by.
    let mut buffer = FrameBuffer::new(topology.clone());
    buffer.copy_from(trajectory.get(1).unwrap()).unwrap();
    let snapshot = format!("{buffer:?}");
    assert_eq!(
        imager.image_onto(1, &mut buffer, &trajectory.get(0).unwrap()),
        Err(PeriodicError::MissingCell { frame: 1 })
    );
    assert_eq!(format!("{buffer:?}"), snapshot);
    let mut loaded = trajectory.clone();
    let reference = trajectory.get(0).unwrap().to_model();
    assert_eq!(
        periodic::image_onto(&mut loaded, &reference),
        Err(PeriodicError::MissingCell { frame: 1 })
    );
    assert_eq!(format!("{loaded:?}"), format!("{trajectory:?}"));

    // A reference without a cell is fine: only the frame's cell is used.
    close(
        &points(&imaged(&trajectory, 0, 1).unwrap()),
        &[[1.0, 0.5, 0.5], [1.15, 0.5, 0.5]],
    );

    // Foreign layouts fail without changing the buffer.
    let foreign = Trajectory::from_items(
        support::topology(&["C", "C"], &[(0, 1)]),
        [frame(&[[0.0; 3], [0.15, 0.0, 0.0]], Some(cubic()), 0)],
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    buffer.copy_from(trajectory.get(0).unwrap()).unwrap();
    let snapshot = format!("{buffer:?}");
    assert_eq!(
        imager.image_onto(0, &mut buffer, &foreign.get(0).unwrap()),
        Err(PeriodicError::ReferenceTopologyMismatch)
    );
    let mut foreign_buffer = FrameBuffer::new(foreign.shared_topology());
    assert_eq!(
        imager.image_onto(4, &mut foreign_buffer, &trajectory.get(0).unwrap()),
        Err(PeriodicError::TopologyMismatch { frame: 4 })
    );
    assert_eq!(format!("{buffer:?}"), snapshot);
}
