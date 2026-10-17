//! The in-memory test reader and the in-memory stream contracts.

use crate::TrajectoryIoErrorKind;
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Forces, Positions, TrajectoryFrame, Velocities};
use kekule::topology::TopologyBuilder;
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use super::*;

/// Sequential and seekable reader over an in-memory trajectory.
struct MemoryTrajectoryReader<'a> {
    trajectory: &'a Trajectory,
    cursor: usize,
}

impl fmt::Debug for MemoryTrajectoryReader<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MemoryTrajectoryReader")
            .field("frame_count", &self.trajectory.len())
            .field("cursor", &self.cursor)
            .finish_non_exhaustive()
    }
}

impl<'a> MemoryTrajectoryReader<'a> {
    fn new(trajectory: &'a Trajectory) -> Self {
        Self {
            trajectory,
            cursor: 0,
        }
    }
}

impl TrajectoryReader for MemoryTrajectoryReader<'_> {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        let Some(frame) = self.trajectory.get(self.cursor) else {
            return Ok(false);
        };
        destination.copy_from(frame)?;
        self.cursor += 1;
        Ok(true)
    }
}

impl SeekableTrajectoryReader for MemoryTrajectoryReader<'_> {
    fn frame_count(&self) -> Option<usize> {
        Some(self.trajectory.len())
    }

    fn read_frame(
        &mut self,
        index: usize,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        let frame = self
            .trajectory
            .get(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        destination.copy_from(frame)?;
        Ok(())
    }
}

fn linear_carbon_topology(atom_count: usize) -> Arc<Topology> {
    let mut graph = MoleculeEditor::new();
    let atoms = (0..atom_count)
        .map(|_| {
            graph
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    for pair in atoms.windows(2) {
        graph.add_bond(pair[0], pair[1], BondOrder::Single).unwrap();
    }
    let mut builder = TopologyBuilder::new();
    let definition = builder
        .add_molecule_definition(graph.finish().unwrap())
        .unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

fn annotated() -> Trajectory {
    let mut trajectory = Trajectory::new(linear_carbon_topology(3));
    trajectory
        .properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(7))
        .unwrap();
    for index in 0..4 {
        let mut frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(
                [
                    Point3::new(index as f64, 0.0, 0.0),
                    Point3::new(0.1, 0.0, 0.0),
                    Point3::new(0.0, 0.1, 0.0),
                ],
                NANOMETER,
            ))
            .unwrap(),
        );
        frame
            .set_time(Some(Quantity::new(index as f64 * 0.5, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index * 10));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
        frame.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        let mut properties = frame.properties_mut();
        properties
            .owner_mut()
            .insert(
                PropertyKey::new("frame").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        properties
            .atoms_mut()
            .insert(
                PropertyKey::new("atom").unwrap(),
                PropertyColumn::Int(vec![Some(1), None, Some(3)]),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
        trajectory
            .get_mut(index as usize)
            .unwrap()
            .properties_mut()
            .bonds_mut()
            .insert(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
    }
    trajectory
}

#[test]
fn memory_random_reads_preserve_sequential_cursor_and_error_destinations() {
    let trajectory = annotated();
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let mut destination = reader.frame_buffer();
    reader.read_frame(3, &mut destination).unwrap();
    assert_eq!(
        (*destination).clone(),
        trajectory.get(3).unwrap().payload().clone()
    );
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(
        (*destination).clone(),
        trajectory.get(0).unwrap().payload().clone()
    );

    reader.read_frame(2, &mut destination).unwrap();
    let before = (*destination).clone();
    assert_eq!(
        reader.read_frame(99, &mut destination),
        Err(TrajectoryIoError::frame_index_out_of_range(99))
    );
    assert_eq!((*destination).clone(), before);
    let unrelated = annotated();
    let mut unrelated_buffer = MemoryTrajectoryReader::new(&unrelated).frame_buffer();
    let unrelated_before = (*unrelated_buffer).clone();
    assert_eq!(
        reader.read_frame(0, &mut unrelated_buffer),
        Err(TrajectoryIoError::new(
            TrajectoryIoErrorKind::TopologyMismatch
        ))
    );
    assert_eq!((*unrelated_buffer).clone(), unrelated_before);
    for index in 1..4 {
        assert!(reader.read_next(&mut destination).unwrap());
        assert_eq!(
            (*destination).clone(),
            trajectory.get(index).unwrap().payload().clone()
        );
    }
    reader.read_frame(0, &mut destination).unwrap();
    assert!(!reader.read_next(&mut destination).unwrap());
}

#[test]
fn memory_streams_debug_their_progress_without_their_frames() {
    let trajectory = annotated();
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let mut destination = reader.frame_buffer();
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(
        format!("{reader:?}"),
        "MemoryTrajectoryReader { frame_count: 4, cursor: 1, .. }"
    );
    let writer = MemoryTrajectoryWriter::new(trajectory.shared_topology());
    assert_eq!(
        format!("{writer:?}"),
        "MemoryTrajectoryWriter { frame_count: 0, .. }"
    );
}
