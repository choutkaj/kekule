use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::topology::Topology;
use kekule::units::{Quantity, ANGSTROM};
use kekule_traj::{
    AtomOrderAssertion, FrameBuffer, FrameBufferData, SeekableTrajectoryReader, TrajectoryError,
    TrajectoryFrameView, TrajectoryReader, TrajectoryWriter,
};

mod support;
use support::topology as build_topology;

fn topology() -> Arc<Topology> {
    build_topology(&["C"], &[])
}

struct CompanionSequential {
    topology: Arc<Topology>,
}

impl TrajectoryReader for CompanionSequential {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, _destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        Ok(false)
    }
}

struct CompanionIndexed {
    topology: Arc<Topology>,
}

impl TrajectoryReader for CompanionIndexed {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, _destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        Ok(false)
    }
}

impl SeekableTrajectoryReader for CompanionIndexed {
    fn frame_count(&self) -> Option<u64> {
        Some(0)
    }

    fn read_frame(
        &mut self,
        index: u64,
        _destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        Err(TrajectoryError::FrameIndexOutOfRange(index))
    }
}

struct CompanionWriter {
    topology: Arc<Topology>,
}

impl TrajectoryWriter for CompanionWriter {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !std::ptr::eq(self.topology.as_ref(), frame.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        Ok(())
    }
}

fn assert_sequential<T: TrajectoryReader>(_reader: &mut T) {}
fn assert_seekable<T: SeekableTrajectoryReader>(_reader: &mut T) {}
fn assert_writer<T: TrajectoryWriter>(_writer: &mut T) {}

#[test]
fn companion_crate_can_publish_frames_and_implement_all_streaming_traits() {
    let topology = topology();
    let assertion = AtomOrderAssertion::assert_file_uses_topology_order(&topology);
    assert!(assertion.is_compatible(&topology));

    let points = [Point3::new(1.0, 2.0, 3.0)];
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    buffer
        .replace_from_data(FrameBufferData::new(
            &topology,
            Quantity::new(&points[..], ANGSTROM),
        ))
        .unwrap();

    let mut sequential = CompanionSequential {
        topology: Arc::clone(&topology),
    };
    let mut indexed = CompanionIndexed {
        topology: Arc::clone(&topology),
    };
    let mut writer = CompanionWriter { topology };
    assert_sequential(&mut sequential);
    assert_seekable(&mut indexed);
    assert_writer(&mut writer);
    writer.write_frame(buffer.frame_view()).unwrap();
}
