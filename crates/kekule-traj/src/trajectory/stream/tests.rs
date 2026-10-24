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

fn sample_payloads(
    ensemble: &kekule::structure::Ensemble,
) -> Vec<kekule::structure::EnsembleSample> {
    ensemble
        .iter()
        .map(|sample| sample.payload().clone())
        .collect()
}

fn read_both(
    trajectory: &Trajectory,
    sampling: &kekule::structure::TrajectorySampling,
) -> [Result<kekule::structure::Ensemble, TrajectoryIoError>; 2] {
    [
        crate::read_ensemble(
            &mut MemoryTrajectoryReader::new(trajectory),
            sampling,
            Default::default(),
        ),
        crate::read_ensemble_seekable(
            &mut MemoryTrajectoryReader::new(trajectory),
            sampling,
            Default::default(),
        ),
    ]
}

#[test]
fn sampled_ensembles_are_the_same_in_memory_and_streamed() {
    use kekule::structure::TrajectorySampling;
    use kekule::topology::AtomSelection;

    let trajectory = annotated();
    let ids = trajectory.topology().atom_ids();
    let atoms = AtomSelection::from_atoms(trajectory.shared_topology(), [ids[0], ids[1]]).unwrap();

    // Every atom, requested out of order and with a duplicate.
    let sampling = TrajectorySampling::new([3, 1, 3]);
    let loaded = trajectory.to_ensemble(&sampling).unwrap();
    assert_eq!(loaded.len(), 3);
    assert_eq!(loaded.state_count(), 1);
    assert!(loaded.topology().shares_layout(trajectory.topology()));
    assert_eq!(loaded.properties(), trajectory.properties());
    let positions = |sample: usize| loaded.get(sample).unwrap().positions().clone();
    assert_eq!(positions(0), *trajectory.get(3).unwrap().positions());
    assert_eq!(positions(1), *trajectory.get(1).unwrap().positions());
    assert_eq!(positions(2), positions(0));
    let first = loaded.get(0).unwrap();
    let frame = trajectory.get(3).unwrap();
    assert_eq!(first.weight(), 1.0);
    assert_eq!(first.cell(), frame.cell());
    assert_eq!(first.properties(), frame.properties());
    for streamed in read_both(&trajectory, &sampling) {
        let streamed = streamed.unwrap();
        assert_eq!(sample_payloads(&streamed), sample_payloads(&loaded));
        assert_eq!(streamed.partition(), loaded.partition());
        assert!(streamed.topology().shares_layout(trajectory.topology()));
        // A reader has no collection properties to keep.
        assert!(streamed.properties().is_empty());
    }

    // An atom subset with a stride and explicit weights.
    let sampling = TrajectorySampling::new((1..4).step_by(2))
        .with_weights([0.25, 2.0])
        .unwrap()
        .with_atoms(atoms);
    let loaded = trajectory.to_ensemble(&sampling).unwrap();
    assert_eq!(loaded.topology().atom_count(), 2);
    assert_eq!(loaded.topology().bond_count(), 1);
    assert!(loaded.properties().is_empty());
    let weights = loaded
        .iter()
        .map(|sample| sample.weight())
        .collect::<Vec<_>>();
    assert_eq!(weights, [0.25, 2.0]);
    let sample = loaded.get(1).unwrap();
    let frame = trajectory.get(3).unwrap();
    assert_eq!(
        sample.positions().values().value().to_vec(),
        frame.positions().values().value()[..2].to_vec()
    );
    assert_eq!(sample.cell(), frame.cell());
    assert!(sample.properties().owner().is_empty());
    let atoms = sample.properties().atoms();
    let atom_key = PropertyKey::new("atom").unwrap();
    let atom_values = [0, 1].map(|index| {
        atoms
            .value(&atom_key, kekule::topology::TopologyAtomIndex::new(index))
            .unwrap()
    });
    assert_eq!(atom_values, [Some(PropertyValue::Int(1)), None]);
    let bond_key = PropertyKey::new("bond").unwrap();
    let bond = sample
        .properties()
        .bonds()
        .value(&bond_key, kekule::topology::TopologyBondIndex::new(0))
        .unwrap();
    assert_eq!(bond, Some(PropertyValue::Int(4)));
    for streamed in read_both(&trajectory, &sampling) {
        let streamed = streamed.unwrap();
        assert_eq!(sample_payloads(&streamed), sample_payloads(&loaded));
        assert_eq!(streamed.topology().atom_count(), 2);
        assert_eq!(streamed.topology().bond_count(), 1);
    }

    // No frames make an empty ensemble.
    let sampling = TrajectorySampling::new([]);
    assert!(trajectory.to_ensemble(&sampling).unwrap().is_empty());
    for streamed in read_both(&trajectory, &sampling) {
        assert!(streamed.unwrap().is_empty());
    }
}

#[test]
fn ensemble_reads_stop_after_the_last_requested_frame_and_keep_the_seek_cursor() {
    let trajectory = annotated();
    let sampling = kekule::structure::TrajectorySampling::new([2, 0]);
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    crate::read_ensemble(&mut reader, &sampling, Default::default()).unwrap();
    assert_eq!(reader.cursor, 3);
    crate::read_ensemble_seekable(&mut reader, &sampling, Default::default()).unwrap();
    assert_eq!(reader.cursor, 3);
}

#[test]
fn ensemble_sampling_rejects_invalid_weights_frames_and_atoms() {
    use kekule::structure::{TrajectoryError, TrajectorySampler, TrajectorySampling};
    use kekule::topology::AtomSelection;

    assert_eq!(
        TrajectorySampling::new([0, 1])
            .with_weights([1.0])
            .unwrap_err(),
        TrajectoryError::WeightCountMismatch {
            expected: 2,
            actual: 1
        }
    );
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            TrajectorySampling::new([0, 1])
                .with_weights([1.0, weight])
                .unwrap_err(),
            TrajectoryError::InvalidWeight { sample: 1 }
        );
    }

    let trajectory = annotated();
    let sampling = TrajectorySampling::new([1, 5, 4]);
    assert_eq!(
        trajectory.to_ensemble(&sampling).unwrap_err(),
        TrajectoryError::IndexOutOfRange { index: 5, len: 4 }
    );
    for streamed in read_both(&trajectory, &sampling) {
        let error = streamed.unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::FrameIndexOutOfRange);
        assert_eq!(error.frame(), Some(4));
    }

    let empty = TrajectorySampling::new([0])
        .with_atoms(AtomSelection::from_atoms(trajectory.shared_topology(), []).unwrap());
    assert!(matches!(
        trajectory.to_ensemble(&empty),
        Err(TrajectoryError::Subset(_))
    ));
    for streamed in read_both(&trajectory, &empty) {
        assert_eq!(
            streamed.unwrap_err().kind(),
            TrajectoryIoErrorKind::InvalidTrajectory
        );
    }

    let mut sampler = TrajectorySampler::new(
        trajectory.shared_topology(),
        &TrajectorySampling::new([0, 2]),
    )
    .unwrap();
    assert_eq!(sampler.frames().collect::<Vec<_>>(), [0, 2]);
    assert!(!sampler.add(1, trajectory.get(1).unwrap()).unwrap());
    let other = Trajectory::from_items(
        linear_carbon_topology(3),
        [TrajectoryFrame::new(Positions::zeros(3))],
    )
    .unwrap();
    assert_eq!(
        sampler.add(0, other.get(0).unwrap()),
        Err(TrajectoryError::TopologyMismatch)
    );
    assert!(sampler.add(2, trajectory.get(2).unwrap()).unwrap());
    assert_eq!(
        sampler.finish().unwrap_err(),
        TrajectoryError::FrameNotSampled { frame: 0 }
    );
}

#[test]
fn ensemble_reads_report_progress_and_stop_on_cancellation() {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Mutex;

    use kekule::diagnostics::Cancellation;

    let trajectory = annotated();
    let sampling = kekule::structure::TrajectorySampling::new([2, 0, 2]);
    let reports = Arc::new(Mutex::new(Vec::new()));
    let options = |reports: &Arc<Mutex<Vec<(usize, usize)>>>| {
        let reports = Arc::clone(reports);
        crate::EnsembleReadOptions::default().with_progress(move |progress| {
            reports
                .lock()
                .unwrap()
                .push((progress.frame_count(), progress.total_frame_count()));
        })
    };

    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    crate::read_ensemble(&mut reader, &sampling, options(&reports)).unwrap();
    assert_eq!(
        std::mem::take(&mut *reports.lock().unwrap()),
        [(0, 3), (1, 3), (2, 3), (3, 3)]
    );
    crate::read_ensemble_seekable(&mut reader, &sampling, options(&reports)).unwrap();
    assert_eq!(*reports.lock().unwrap(), [(0, 2), (1, 2), (2, 2)]);

    // Cancel once the first frame has been read.
    let read = Arc::new(AtomicUsize::new(0));
    let progress = Arc::clone(&read);
    let check = Arc::clone(&read);
    let options = crate::EnsembleReadOptions::default()
        .with_progress(move |step| progress.store(step.frame_count(), Ordering::Relaxed))
        .with_cancellation(Cancellation::new(move || {
            check.load(Ordering::Relaxed) >= 1
        }));
    assert!(options.cancellation().is_some());
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let error = crate::read_ensemble(&mut reader, &sampling, options.clone()).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::Cancelled);
    assert_eq!(reader.cursor, 1);
    let error = crate::read_ensemble_seekable(&mut reader, &sampling, options).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::Cancelled);
}
