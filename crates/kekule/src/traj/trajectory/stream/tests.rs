//! The in-memory test reader and the in-memory stream contracts.

use crate::math::{PeriodicCell, Point3, Vector3};
use crate::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use crate::properties::{PropertyColumn, PropertyKey, PropertyValue};
use crate::structure::{Forces, Geometry, Place, Positions, Velocities};
use crate::topology::TopologyBuilder;
use crate::traj::TrajectoryIoErrorKind;
use crate::units::{Quantity, NANOMETER, PICOSECOND};

use super::*;

/// Sequential and seekable reader over the geometries of a time-ordered set.
struct MemoryTrajectoryReader<'a> {
    set: &'a StructureSet,
    cursor: usize,
}

impl fmt::Debug for MemoryTrajectoryReader<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MemoryTrajectoryReader")
            .field("frame_count", &self.set.len())
            .field("cursor", &self.cursor)
            .finish_non_exhaustive()
    }
}

impl<'a> MemoryTrajectoryReader<'a> {
    fn new(set: &'a StructureSet) -> Self {
        Self { set, cursor: 0 }
    }
}

impl TrajectoryReader for MemoryTrajectoryReader<'_> {
    fn topology(&self) -> &Topology {
        self.set.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.set.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        let Some(geometry) = self.set.get(self.cursor) else {
            return Ok(false);
        };
        destination.copy_from(FrameView::from_geometry(geometry)?)?;
        self.cursor += 1;
        Ok(true)
    }
}

impl SeekableTrajectoryReader for MemoryTrajectoryReader<'_> {
    fn frame_count(&self) -> Option<usize> {
        Some(self.set.len())
    }

    fn read_frame(
        &mut self,
        index: usize,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        let geometry = self
            .set
            .get(index)
            .ok_or(TrajectoryIoError::frame_index_out_of_range(index))?;
        destination.copy_from(FrameView::from_geometry(geometry)?)?;
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

/// Four frames carrying every kind of frame state, at times `0.5 i` ps and
/// steps `10 i`.
fn annotated() -> StructureSet {
    let mut set = StructureSet::new(linear_carbon_topology(3));
    let times = Quantity::new(Vec::new(), PICOSECOND);
    set.set_order(Some(Order::Time(
        TimeAxis::new(Some(times), Some(Vec::new())).unwrap(),
    )))
    .unwrap();
    set.properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(7))
        .unwrap();
    for index in 0..4 {
        let mut geometry = Geometry::new(
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
        geometry.set_velocities(Some(Velocities::zeros(3))).unwrap();
        geometry.set_forces(Some(Forces::zeros(3))).unwrap();
        geometry.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        let mut properties = geometry.properties_mut();
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
        let place = Place::Time {
            time: Some(Quantity::new(index as f64 * 0.5, PICOSECOND)),
            step: Some(index * 10),
        };
        set.push(geometry, Some(place)).unwrap();
        set.get_mut(index as usize)
            .unwrap()
            .properties_mut()
            .bonds_mut()
            .insert(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
    }
    set
}

/// The buffered frame's geometry and place.
fn buffered(buffer: &FrameBuffer) -> (Geometry, Place) {
    let frame = buffer.frame_view();
    (frame.payload().clone(), frame.place())
}

/// Geometry `index` of `set` with its place.
fn stored(set: &StructureSet, index: usize) -> (Geometry, Place) {
    let geometry = set.get(index).unwrap();
    (geometry.payload().clone(), geometry.place().unwrap())
}

#[test]
fn memory_random_reads_preserve_sequential_cursor_and_error_destinations() {
    let set = annotated();
    let mut reader = MemoryTrajectoryReader::new(&set);
    let mut destination = reader.frame_buffer();
    reader.read_frame(3, &mut destination).unwrap();
    assert_eq!(buffered(&destination), stored(&set, 3));
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(buffered(&destination), stored(&set, 0));

    reader.read_frame(2, &mut destination).unwrap();
    let before = buffered(&destination);
    assert_eq!(
        reader.read_frame(99, &mut destination),
        Err(TrajectoryIoError::frame_index_out_of_range(99))
    );
    assert_eq!(buffered(&destination), before);
    let unrelated = annotated();
    let mut unrelated_buffer = MemoryTrajectoryReader::new(&unrelated).frame_buffer();
    let unrelated_before = buffered(&unrelated_buffer);
    assert_eq!(
        reader.read_frame(0, &mut unrelated_buffer),
        Err(TrajectoryIoError::new(
            TrajectoryIoErrorKind::TopologyMismatch
        ))
    );
    assert_eq!(buffered(&unrelated_buffer), unrelated_before);
    for index in 1..4 {
        assert!(reader.read_next(&mut destination).unwrap());
        assert_eq!(buffered(&destination), stored(&set, index));
    }
    reader.read_frame(0, &mut destination).unwrap();
    assert!(!reader.read_next(&mut destination).unwrap());
}

#[test]
fn memory_streams_debug_their_progress_without_their_frames() {
    let set = annotated();
    let mut reader = MemoryTrajectoryReader::new(&set);
    let mut destination = reader.frame_buffer();
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(
        format!("{reader:?}"),
        "MemoryTrajectoryReader { frame_count: 4, cursor: 1, .. }"
    );
    let writer = MemoryTrajectoryWriter::new(set.shared_topology());
    assert_eq!(
        format!("{writer:?}"),
        "MemoryTrajectoryWriter { frame_count: 0, .. }"
    );
}

#[test]
fn geometries_off_a_time_order_are_not_frames() {
    let set = annotated();
    let mut reader = MemoryTrajectoryReader::new(&set);
    let mut destination = reader.frame_buffer();
    assert!(reader.read_next(&mut destination).unwrap());
    let before = buffered(&destination);
    let mut stepped = annotated();
    stepped
        .set_order(Some(Order::Step(vec![0, 10, 20, 30])))
        .unwrap();
    let error = MemoryTrajectoryReader::new(&stepped)
        .read_frame(0, &mut destination)
        .unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::UnsupportedField);
    assert_eq!(buffered(&destination), before);
}

/// The geometry payloads of `set` in order.
fn payloads(set: &StructureSet) -> Vec<Geometry> {
    set.iter()
        .map(|geometry| geometry.payload().clone())
        .collect()
}

#[test]
fn sampled_sets_are_the_same_in_memory_and_streamed() {
    use crate::topology::AtomSelection;

    let set = annotated();
    let mut reader = MemoryTrajectoryReader::new(&set);
    let mut buffer = reader.frame_buffer();
    let mut writer = MemoryTrajectoryWriter::new(reader.shared_topology());
    while reader.read_next(&mut buffer).unwrap() {
        writer.write_frame(buffer.frame_view()).unwrap();
    }
    let streamed = writer.into_structure_set();
    assert_eq!(streamed.order(), set.order());
    assert_eq!(payloads(&streamed), payloads(&set));
    // A reader has no set properties to keep.
    assert!(streamed.properties().is_empty());
    let ids = set.topology().atom_ids();
    let atoms = AtomSelection::from_atoms(set.shared_topology(), [ids[0], ids[1]]).unwrap();

    // Every atom, requested out of order and with a duplicate.
    let loaded = StructureSet::from_reader(&set, [3, 1, 3]).unwrap();
    let frame_key = PropertyKey::new("frame").unwrap();
    // Frame `i` has its first atom at x = i nm and its frame property `i`.
    assert_eq!(
        loaded
            .iter()
            .map(|geometry| (
                geometry.positions().values().value()[0].x,
                geometry.properties().owner().get(&frame_key).cloned()
            ))
            .collect::<Vec<_>>(),
        [3, 1, 3].map(|frame| (f64::from(frame), Some(PropertyValue::Int(frame.into()))))
    );
    assert_eq!(loaded.len(), 3);
    assert!(loaded.topology().shares_layout(set.topology()));
    assert_eq!(loaded.get(0).unwrap().place(), set.get(3).unwrap().place());
    let positions = |geometry: usize| loaded.get(geometry).unwrap().positions().clone();
    assert_eq!(positions(0), *set.get(3).unwrap().positions());
    assert_eq!(positions(1), *set.get(1).unwrap().positions());
    assert_eq!(positions(2), positions(0));
    let first = loaded.get(0).unwrap();
    let frame = set.get(3).unwrap();
    assert_eq!(first.cell(), frame.cell());
    assert_eq!(first.properties(), frame.properties());
    let from_stream = StructureSet::from_reader(&streamed, [3, 1, 3]).unwrap();
    assert_eq!(payloads(&from_stream), payloads(&loaded));

    // An atom subset with a stride.
    let loaded = StructureSet::from_reader(&set, (1..4).step_by(2))
        .unwrap()
        .subset(&atoms)
        .unwrap();
    // Frames 1 and 3, with the first two atoms of each, in their own box, and
    // without frame properties, which belong to the whole geometry.
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    for (geometry, frame) in loaded.iter().zip([1.0, 3.0]) {
        assert_eq!(
            geometry.positions().values().value().to_vec(),
            [Point3::new(frame, 0.0, 0.0), Point3::new(0.1, 0.0, 0.0)]
        );
        assert_eq!(geometry.cell(), Some(&cell));
        assert!(geometry.properties().owner().is_empty());
        assert_eq!(
            geometry.velocities().map(|velocities| velocities.len()),
            Some(2)
        );
        assert_eq!(geometry.forces().map(|forces| forces.len()), Some(2));
    }
    // The first two carbons of the chain and the bond between them.
    assert_eq!(loaded.topology().atom_count(), 2);
    assert_eq!(loaded.topology().bond_count(), 1);
    let view = loaded.get(0).unwrap().as_structure_view();
    let carbons = view.atoms().collect::<Vec<_>>();
    assert!(carbons.iter().all(|atom| atom.element.symbol() == "C"));
    assert_eq!(
        carbons[0]
            .neighbors()
            .map(|atom| atom.id())
            .collect::<Vec<_>>(),
        [carbons[1].id()]
    );
    assert!(loaded.properties().is_empty());
    let geometry = loaded.get(1).unwrap();
    let frame = set.get(3).unwrap();
    assert_eq!(geometry.place(), frame.place());
    assert_eq!(
        geometry.positions().values().value().to_vec(),
        frame.positions().values().value()[..2].to_vec()
    );
    assert_eq!(geometry.cell(), frame.cell());
    assert!(geometry.properties().owner().is_empty());
    let atoms_table = geometry.properties().atoms();
    let atom_key = PropertyKey::new("atom").unwrap();
    let atom_values = [0, 1].map(|index| {
        atoms_table
            .value(&atom_key, crate::topology::TopologyAtomIndex::new(index))
            .unwrap()
    });
    assert_eq!(atom_values, [Some(PropertyValue::Int(1)), None]);
    let bond_key = PropertyKey::new("bond").unwrap();
    let bond = geometry
        .properties()
        .bonds()
        .value(&bond_key, crate::topology::TopologyBondIndex::new(0))
        .unwrap();
    assert_eq!(bond, Some(PropertyValue::Int(4)));
    let from_stream = StructureSet::from_reader(&streamed, (1..4).step_by(2))
        .unwrap()
        .subset(&atoms)
        .unwrap();
    assert_eq!(payloads(&from_stream), payloads(&loaded));
    assert_eq!(from_stream.topology().atom_count(), 2);
    assert_eq!(from_stream.topology().bond_count(), 1);

    // No frames make an empty set.
    assert!(StructureSet::from_reader(&set, []).unwrap().is_empty());
    assert_eq!(
        StructureSet::from_reader(&set, [1, 5, 4]).unwrap_err(),
        crate::structure::StructureSetReadError::IndexOutOfRange { index: 5, len: 4 }
    );
}
