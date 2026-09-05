use super::*;
use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{Properties, PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::Positions;
use kekule::topology::{AtomSelection, Topology, TopologyBuilder};
use kekule::units::{
    Quantity, ANGSTROM, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT, DIMENSIONLESS, KELVIN,
    KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND, SQUARE_ANGSTROM, SQUARE_NANOMETER,
};
use std::sync::Arc;

fn make_topology(with_bond: bool) -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    let first = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    if with_bond {
        let second = editor
            .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
            .unwrap();
        editor.add_bond(first, second, BondOrder::Single).unwrap();
    }
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

fn positions(values: &[Point3]) -> Positions {
    Positions::new(Quantity::new(values, NANOMETER)).unwrap()
}

fn key(value: &str) -> PropertyKey {
    PropertyKey::new(value).unwrap()
}

fn real(value: f64) -> PropertyValue {
    PropertyValue::Real {
        value,
        unit: DIMENSIONLESS,
    }
}

#[test]
fn vector_arrays_are_topology_free_unit_aware_and_equal_by_values() {
    let vectors = [Vector3::new(1.0, 2.0, 3.0)];
    let velocities = Velocities::new(Quantity::new(vectors, ANGSTROM / PICOSECOND)).unwrap();
    let same = Velocities::new(Quantity::new(vectors, CANONICAL_VELOCITY_UNIT)).unwrap();
    let converted_velocity = velocities.values().value()[0];
    assert!((converted_velocity.x - 0.1).abs() < 1.0e-12);
    assert!((converted_velocity.y - 0.2).abs() < 1.0e-12);
    assert!((converted_velocity.z - 0.3).abs() < 1.0e-12);
    assert_ne!(velocities, same);
    assert_eq!(velocities.len(), 1);
    assert!(!velocities.is_empty());

    let forces = Forces::new(Quantity::new(vectors, KILOJOULE_PER_MOLE / ANGSTROM)).unwrap();
    let converted_force = forces.values().value()[0];
    assert!((converted_force.x - 10.0).abs() < 1.0e-12);
    assert!((converted_force.y - 20.0).abs() < 1.0e-12);
    assert!((converted_force.z - 30.0).abs() < 1.0e-12);
    assert!(matches!(
        Velocities::new(Quantity::new(
            [Vector3::new(f64::NAN, 0.0, 0.0)],
            CANONICAL_VELOCITY_UNIT
        )),
        Err(FrameError::NonFiniteVector { index: 0 })
    ));
    assert!(matches!(
        Forces::new(Quantity::new(vectors, KELVIN)),
        Err(FrameError::Unit(UnitError::IncompatibleUnits { .. }))
    ));
}

#[test]
fn frames_validate_all_dense_dimensions_at_the_owner_boundary() {
    let topology = make_topology(true);
    let valid_positions = positions(&[Point3::origin(), Point3::origin()]);
    let mut frame = TrajectoryFrame::new(valid_positions, topology.bond_count());
    frame
        .set_velocities(Some(Velocities::zeros(topology.atom_count())))
        .unwrap();
    frame
        .set_forces(Some(Forces::zeros(topology.atom_count())))
        .unwrap();
    frame.validate(&topology).unwrap();

    assert!(matches!(
        TrajectoryFrame::new(positions(&[Point3::origin()]), topology.bond_count())
            .validate(&topology),
        Err(FrameError::AtomCountMismatch {
            expected: 2,
            actual: 1
        })
    ));
    assert!(matches!(
        frame.set_velocities(Some(Velocities::zeros(1))),
        Err(FrameError::AtomCountMismatch {
            expected: 2,
            actual: 1
        })
    ));
    assert!(matches!(
        frame.set_forces(Some(Forces::zeros(1))),
        Err(FrameError::AtomCountMismatch {
            expected: 2,
            actual: 1
        })
    ));
    assert!(matches!(
        frame.set_properties(Properties::realization(1, topology.bond_count())),
        Err(FrameError::AtomCountMismatch { .. })
    ));
    assert!(matches!(
        frame.set_properties(Properties::realization(topology.atom_count(), 0)),
        Err(FrameError::BondCountMismatch { .. })
    ));
}

#[test]
fn frame_view_borrows_model_state_and_preserves_time_and_step() {
    let topology = make_topology(false);
    let atom = topology.atom_ids()[0];
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let score = key("score");
    let mut frame = TrajectoryFrame::new(
        positions(&[Point3::new(3.0, 0.0, 0.0)]),
        topology.bond_count(),
    );
    frame.set_cell(Some(cell));
    frame
        .insert_property(key("method"), PropertyValue::String("md".into()))
        .unwrap();
    frame
        .set_atom_property(0, score.clone(), Some(real(0.8)))
        .unwrap();
    frame
        .set_velocities(Some(Velocities::zeros(topology.atom_count())))
        .unwrap();
    frame
        .set_forces(Some(Forces::zeros(topology.atom_count())))
        .unwrap();
    frame
        .set_time(Some(Quantity::new(2.5, PICOSECOND)))
        .unwrap();
    frame.set_step(Some(7));

    let view = frame.view(&topology).unwrap();
    assert_eq!(view.as_model().position(atom).unwrap().value().x, 3.0);
    assert_eq!(
        view.atom_property(0, &key("score")).unwrap(),
        Some(PropertyValue::Real {
            value: 0.8,
            unit: DIMENSIONLESS,
        })
    );
    assert_eq!(view.time(), Some(Quantity::new(2.5, PICOSECOND)));
    assert_eq!(view.step(), Some(7));
    assert!(view.velocities().is_some());
    assert!(view.forces().is_some());
    assert_eq!(
        view.positions().values().value().as_ptr(),
        frame.positions().values().value().as_ptr()
    );
    let model = view.to_model();
    assert!(Arc::ptr_eq(&model.shared_topology(), &topology));
    assert_eq!(model.positions(), frame.positions());
    assert_eq!(model.cell(), frame.cell());
    assert_eq!(model.properties(), frame.properties());

    frame
        .set_atom_property(0, score.clone(), Some(real(0.2)))
        .unwrap();
    assert_eq!(
        model.atom_properties().value(&score, 0).unwrap(),
        Some(real(0.8))
    );
}

#[test]
fn canonical_trajectory_constructors_accept_owned_and_shared_topology() {
    let owned_new = Trajectory::new(Arc::try_unwrap(make_topology(false)).unwrap());
    let shared = owned_new.shared_topology();
    let shared_new = Trajectory::new(Arc::clone(&shared));
    assert!(Arc::ptr_eq(&shared_new.shared_topology(), &shared));

    let owned_topology = Arc::try_unwrap(make_topology(false)).unwrap();
    let owned_frames = Trajectory::from_frames(
        owned_topology,
        [TrajectoryFrame::new(
            positions(&[Point3::new(1.0, 0.0, 0.0)]),
            0,
        )],
    )
    .unwrap();
    let shared_frames = Trajectory::from_frames(
        owned_frames.shared_topology(),
        [TrajectoryFrame::new(
            positions(&[Point3::new(2.0, 0.0, 0.0)]),
            0,
        )],
    )
    .unwrap();
    assert!(owned_frames
        .topology()
        .same_layout(shared_frames.topology()));
}

#[test]
fn trajectory_frame_views_are_topology_bound_and_stably_ordered() {
    let topology = make_topology(false);
    let mut first = TrajectoryFrame::new(positions(&[Point3::new(1.0, 0.0, 0.0)]), 0);
    first.set_step(Some(10));
    let mut second = TrajectoryFrame::new(positions(&[Point3::new(2.0, 0.0, 0.0)]), 0);
    second.set_step(Some(20));
    let trajectory = Trajectory::from_frames(Arc::clone(&topology), [first, second]).unwrap();

    assert_eq!(trajectory.frame(0).unwrap().step(), Some(10));
    assert_eq!(trajectory.frame(1).unwrap().step(), Some(20));
    assert!(trajectory.frame(2).is_none());
    assert_eq!(
        trajectory
            .frames()
            .map(|frame| frame.step())
            .collect::<Vec<_>>(),
        [Some(10), Some(20)]
    );
    assert!(Arc::ptr_eq(
        &trajectory.frame(0).unwrap().shared_topology(),
        &topology
    ));
    let model = trajectory.frame(0).unwrap().to_model();
    assert!(Arc::ptr_eq(&model.shared_topology(), &topology));
    assert_eq!(model.positions().values().value()[0].x, 1.0);
}

#[test]
fn trajectory_rejects_bad_frame_dimensions_and_keeps_one_topology() {
    let topology = make_topology(true);
    let mut trajectory = Trajectory::new(Arc::clone(&topology));
    let mut valid = TrajectoryFrame::new(
        positions(&[Point3::origin(), Point3::origin()]),
        topology.bond_count(),
    );
    valid.set_step(Some(1));
    trajectory.push(valid).unwrap();
    assert!(Arc::ptr_eq(&trajectory.shared_topology(), &topology));
    assert_eq!(trajectory.frames().next().unwrap().step(), Some(1));

    let wrong = TrajectoryFrame::new(positions(&[Point3::origin()]), topology.bond_count());
    assert!(matches!(
        trajectory.push(wrong),
        Err(TrajectoryError::Frame(error))
            if matches!(*error, FrameError::AtomCountMismatch { .. })
    ));
}

#[test]
fn frame_and_buffer_scope_canonical_atom_properties_to_semantic_apis() {
    let topology = make_topology(true);
    let mut frame = TrajectoryFrame::new(
        positions(&[Point3::origin(), Point3::origin()]),
        topology.bond_count(),
    );
    assert!(matches!(
        frame.set_atom_property(0, key("occupancy"), Some(PropertyValue::Int(1))),
        Err(FrameError::Property(PropertyError::ReservedKey(_)))
    ));
    frame.set_occupancy_at(0, Some(0.8)).unwrap();
    frame
        .set_b_factor_at(0, Some(Quantity::new(25.0, SQUARE_ANGSTROM)))
        .unwrap();
    assert_eq!(frame.occupancy_at(0).unwrap(), Some(0.8));
    let b_factor = frame.b_factor_at(0).unwrap().unwrap();
    assert_eq!(b_factor.unit(), SQUARE_NANOMETER);
    assert!((*b_factor.value() - 0.25).abs() < 1.0e-12);
    assert!(frame.set_occupancy_at(0, Some(f64::NAN)).is_err());
    assert!(frame
        .set_b_factor_at(0, Some(Quantity::new(1.0, KELVIN)))
        .is_err());

    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(matches!(
        buffer.insert_atom_property_column(
            key("b_factor"),
            PropertyColumn::String(vec![Some("bad".into()), None]),
        ),
        Err(FrameError::Property(PropertyError::ReservedKey(_)))
    ));
    buffer.set_properties(frame.properties().clone()).unwrap();
    assert_eq!(buffer.occupancy_at(0).unwrap(), Some(0.8));
    assert_eq!(buffer.b_factor_at(0).unwrap(), Some(b_factor));
    buffer
        .set_atom_property(0, key("buffer_score"), Some(real(0.4)))
        .unwrap();
    assert_eq!(
        buffer.atom_property(0, &key("buffer_score")).unwrap(),
        Some(real(0.4))
    );
    assert!(matches!(
        buffer.remove_atom_property_column(&key("occupancy")),
        Err(FrameError::Property(PropertyError::ReservedKey(_)))
    ));
}

#[test]
fn trajectory_slice_transfers_every_per_atom_frame_field() {
    let topology = make_topology(true);
    let atoms = topology.atom_ids();
    let selection = AtomSelection::from_atoms(&topology, [atoms[1]]).unwrap();
    let mut frame = TrajectoryFrame::new(
        positions(&[Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)]),
        topology.bond_count(),
    );
    frame
        .set_atom_property(0, key("score"), Some(real(0.25)))
        .unwrap();
    frame
        .set_atom_property(1, key("score"), Some(real(0.75)))
        .unwrap();
    frame
        .set_bond_property(0, key("score"), Some(real(9.0)))
        .unwrap();
    frame
        .insert_property(key("frame_energy"), real(12.0))
        .unwrap();
    frame
        .set_velocities(Some(
            Velocities::new(Quantity::new(
                [Vector3::new(3.0, 0.0, 0.0), Vector3::new(4.0, 0.0, 0.0)],
                CANONICAL_VELOCITY_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
    frame
        .set_forces(Some(
            Forces::new(Quantity::new(
                [Vector3::new(5.0, 0.0, 0.0), Vector3::new(6.0, 0.0, 0.0)],
                CANONICAL_FORCE_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
    frame
        .set_time(Some(Quantity::new(2.0, PICOSECOND)))
        .unwrap();
    frame.set_step(Some(8));
    let mut trajectory = Trajectory::from_frames(Arc::clone(&topology), [frame]).unwrap();
    trajectory
        .insert_property(
            key("collection_source"),
            PropertyValue::String("test".into()),
        )
        .unwrap();

    let sliced = trajectory.slice(&selection).unwrap();
    assert!(sliced.properties().owner_is_empty());
    assert_eq!(sliced.topology().atom_count(), 1);
    assert_eq!(sliced.topology().bond_count(), 0);
    let frame = sliced.frame(0).unwrap();
    assert!(frame.properties().owner_is_empty());
    assert_eq!(
        frame.positions().values().value(),
        &[Point3::new(2.0, 0.0, 0.0)]
    );
    assert_eq!(
        frame.atom_properties().value(&key("score"), 0).unwrap(),
        Some(real(0.75))
    );
    assert!(!frame.bond_properties().has_data());
    assert_eq!(frame.velocities().unwrap().value()[0].x, 4.0);
    assert_eq!(frame.forces().unwrap().value()[0].x, 6.0);
    assert_eq!(frame.time(), Some(Quantity::new(2.0, PICOSECOND)));
    assert_eq!(frame.step(), Some(8));
}

#[test]
fn frame_buffer_publication_is_transactional_reuses_arrays_and_clears_optionals() {
    let topology = make_topology(true);
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    let position_ptr = buffer.positions().values().value().as_ptr();
    let velocity_ptr = buffer.velocities.values().value().as_ptr();
    let force_ptr = buffer.forces.values().value().as_ptr();
    let points = [Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)];
    let vectors = [Vector3::new(1.0, 0.0, 0.0); 2];
    let mut properties = Properties::realization(topology.atom_count(), topology.bond_count());
    properties
        .set_realization_atom_value(key("atom_score"), 0, Some(real(0.75)))
        .unwrap();
    properties
        .set_realization_bond_value(key("bond_score"), 0, Some(real(2.0)))
        .unwrap();
    properties
        .insert(key("codec_value"), PropertyValue::Int(7))
        .unwrap();

    buffer
        .replace_from_data(
            FrameBufferData::new(Quantity::new(points.as_slice(), ANGSTROM))
                .with_velocities(Quantity::new(vectors.as_slice(), CANONICAL_VELOCITY_UNIT))
                .with_forces(Quantity::new(vectors.as_slice(), CANONICAL_FORCE_UNIT))
                .with_time(Quantity::new(1.0, PICOSECOND))
                .with_step(4)
                .with_properties(&properties),
        )
        .unwrap();
    assert!(buffer.frame_view().velocities().is_some());
    assert!(buffer.frame_view().forces().is_some());
    assert_eq!(buffer.positions().values().value().as_ptr(), position_ptr);
    assert_eq!(buffer.velocities.values().value().as_ptr(), velocity_ptr);
    assert_eq!(buffer.forces.values().value().as_ptr(), force_ptr);
    assert!(buffer.atom_properties().has_data());
    assert!(buffer.bond_properties().has_data());
    assert_eq!(
        buffer.properties().get(&key("codec_value")),
        Some(&PropertyValue::Int(7))
    );

    let before = buffer.frame_view().positions().values().value().to_vec();
    assert!(matches!(
        buffer.replace_from_data(FrameBufferData::new(Quantity::new(
            [Point3::new(f64::NAN, 0.0, 0.0)].as_slice(),
            ANGSTROM
        ))),
        Err(FrameError::Position(_))
    ));
    assert_eq!(*buffer.positions().values().value(), before.as_slice());

    buffer
        .replace_from_data(FrameBufferData::new(Quantity::new(
            points.as_slice(),
            ANGSTROM,
        )))
        .unwrap();
    assert!(buffer.frame_view().velocities().is_none());
    assert!(buffer.frame_view().forces().is_none());
    assert!(buffer.frame_view().time().is_none());
    assert!(buffer.frame_view().step().is_none());
    assert!(buffer.properties().is_empty());
}

#[test]
fn frame_buffer_copy_requires_the_exact_frame_topology() {
    let topology = make_topology(false);
    let independent = make_topology(false);
    let frame = TrajectoryFrame::new(positions(&[Point3::origin()]), topology.bond_count());
    let view = frame.view(&topology).unwrap();
    let mut correct = FrameBuffer::new(Arc::clone(&topology));
    correct.copy_from(view).unwrap();

    let mut wrong = FrameBuffer::new(independent);
    assert_eq!(wrong.copy_from(view), Err(FrameError::TopologyMismatch));
}

#[test]
fn memory_reader_and_writer_round_trip_validated_frames() {
    let topology = make_topology(true);
    let mut trajectory = Trajectory::new(Arc::clone(&topology));
    let mut frame = TrajectoryFrame::new(
        positions(&[Point3::new(5.0, 0.0, 0.0), Point3::new(6.0, 0.0, 0.0)]),
        topology.bond_count(),
    );
    frame.set_step(Some(9));
    frame
        .set_atom_property(0, key("atom_score"), Some(real(0.6)))
        .unwrap();
    frame
        .set_bond_property(0, key("bond_score"), Some(real(4.0)))
        .unwrap();
    assert_eq!(
        frame.atom_property(0, &key("atom_score")).unwrap(),
        Some(real(0.6))
    );
    assert_eq!(
        frame.bond_property(0, &key("bond_score")).unwrap(),
        Some(real(4.0))
    );
    trajectory.push(frame).unwrap();

    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(buffer.positions().values().value()[0].x, 5.0);
    assert_eq!(buffer.frame_view().step(), Some(9));
    assert_eq!(
        buffer.atom_property(0, &key("atom_score")).unwrap(),
        Some(real(0.6))
    );
    assert_eq!(
        buffer.bond_property(0, &key("bond_score")).unwrap(),
        Some(real(4.0))
    );
    assert!(!reader.read_next(&mut buffer).unwrap());

    let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
    writer.write_frame(buffer.frame_view()).unwrap();
    let written = writer.to_trajectory();
    assert_eq!(written.len(), 1);
    assert_eq!(written.frames().next().unwrap().step(), Some(9));
    assert_eq!(
        written
            .frame(0)
            .unwrap()
            .atom_properties()
            .value(&key("atom_score"), 0)
            .unwrap(),
        Some(real(0.6))
    );
    assert_eq!(
        written
            .frame(0)
            .unwrap()
            .bond_properties()
            .value(&key("bond_score"), 0)
            .unwrap(),
        Some(real(4.0))
    );

    let independent = make_topology(true);
    let other = TrajectoryFrame::new(
        positions(&[Point3::origin(), Point3::origin()]),
        independent.bond_count(),
    );
    let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
    assert!(matches!(
        writer.write_frame(other.view(&independent).unwrap()),
        Err(TrajectoryError::TopologyMismatch)
    ));
}

#[test]
fn coordinate_reader_uses_topology_order_and_rejects_mismatched_buffers() {
    let topology = make_topology(false);
    let mut reader = CoordinateFrameReader::new(
        Arc::clone(&topology),
        [Quantity::new(vec![Point3::new(6.0, 0.0, 0.0)], ANGSTROM)],
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!((buffer.positions().values().value()[0].x - 0.6).abs() < 1.0e-15);

    let mut reader = CoordinateFrameReader::new(
        Arc::clone(&topology),
        [Quantity::new(vec![Point3::origin()], ANGSTROM)],
    )
    .unwrap();
    let mut wrong_buffer = FrameBuffer::new(make_topology(false));
    assert!(matches!(
        reader.read_next(&mut wrong_buffer),
        Err(TrajectoryError::TopologyMismatch)
    ));

    assert!(matches!(
        CoordinateFrameReader::new(
            Arc::clone(&topology),
            [Quantity::new(vec![Point3::origin(); 2], ANGSTROM)]
        ),
        Err(TrajectoryError::Frame(error))
            if matches!(*error, FrameError::AtomCountMismatch { expected: 1, actual: 2 })
    ));
}
