//! The failure contract every codec writer shares: a rejected frame writes
//! nothing, and a write cut off partway refuses every further frame and the
//! finish, so a corrupt stream is never extended or reported as complete.

use std::io::{Cursor, ErrorKind, Seek, Write};
use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::topology::Topology;
use kekule::units::{Quantity, ANGSTROM, NANOMETER, PICOSECOND};
use kekule_traj::io::dcd::{DcdWriteOptions, DcdWriter};
use kekule_traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions, TrrWriter};
use kekule_traj::io::xtc::{XtcWriteOptions, XtcWriter};
use kekule_traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule_traj::{
    FrameBuffer, TrajectoryCodecErrorKind, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorContext, TrajectoryIoOperation, TrajectoryWriter,
};

use crate::support;
use support::{codec_error, ByteBudgetSink};

const CODECS: [TrajectoryFormat; 4] = [
    TrajectoryFormat::Xyz,
    TrajectoryFormat::Dcd,
    TrajectoryFormat::Trr,
    TrajectoryFormat::Xtc,
];

/// A codec writer behind one interface, so each test runs over every codec.
trait CodecWriter<W>: TrajectoryWriter {
    fn sink(&self) -> &W;
    fn finish_stream(self: Box<Self>) -> Result<W, TrajectoryIoError>;
}

macro_rules! codec_writer {
    ($writer:ident, $($bound:path),+) => {
        impl<W: $($bound +)+> CodecWriter<W> for $writer<W> {
            fn sink(&self) -> &W {
                self.writer()
            }

            fn finish_stream(self: Box<Self>) -> Result<W, TrajectoryIoError> {
                (*self).finish()
            }
        }
    };
}

codec_writer!(XyzWriter, Write);
codec_writer!(DcdWriter, Write, Seek);
codec_writer!(TrrWriter, Write);
codec_writer!(XtcWriter, Write);

fn label(format: TrajectoryFormat) -> String {
    format!("contract.{}", format.name().to_ascii_lowercase())
}

fn topology() -> Arc<Topology> {
    support::topology(&["C", "H", "O"], &[(0, 1), (0, 2)])
}

fn open<W: Write + Seek + 'static>(
    format: TrajectoryFormat,
    sink: W,
    topology: &Arc<Topology>,
) -> Box<dyn CodecWriter<W>> {
    let topology = Arc::clone(topology);
    let label = label(format);
    match format {
        TrajectoryFormat::Xyz => {
            Box::new(XyzWriter::new(sink, topology, XyzWriteOptions::default(), label).unwrap())
        }
        TrajectoryFormat::Dcd => Box::new(
            DcdWriter::new(
                sink,
                topology,
                DcdWriteOptions::default().with_cells(true),
                label,
            )
            .unwrap(),
        ),
        TrajectoryFormat::Trr => Box::new(
            TrrWriter::new(
                sink,
                topology,
                TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
                label,
            )
            .unwrap(),
        ),
        TrajectoryFormat::Xtc => {
            Box::new(XtcWriter::new(sink, topology, XtcWriteOptions::default(), label).unwrap())
        }
        _ => unreachable!("every supported codec is listed"),
    }
}

/// Frame `index` carrying exactly the state each codec requires: DCD the
/// step sequence and a canonical cell, TRR step and time, XTC step, time
/// and a cell; XYZ positions only.
fn valid_frame(format: TrajectoryFormat, topology: &Arc<Topology>, index: u64) -> FrameBuffer {
    let shift = index as f64;
    let mut frame = FrameBuffer::new(Arc::clone(topology));
    let mut frame_state = frame.frame_mut();
    frame_state
        .set_positions(Quantity::new(
            [
                Point3::new(shift, 1.0, 2.0),
                Point3::new(3.0, 4.0 + shift, 5.0),
                Point3::new(6.0, 7.0, 8.0 + shift),
            ],
            ANGSTROM,
        ))
        .unwrap();
    match format {
        TrajectoryFormat::Xyz => {}
        TrajectoryFormat::Dcd => {
            frame_state.set_step(Some(index));
            frame_state.set_cell(Some(cell(
                [[10.0, 0.0, 0.0], [2.0, 11.0, 0.0], [1.0, 3.0, 12.0]],
                ANGSTROM,
            )));
        }
        TrajectoryFormat::Trr | TrajectoryFormat::Xtc => {
            frame_state.set_step(Some(index));
            frame_state
                .set_time(Some(Quantity::new(0.5 * shift, PICOSECOND)))
                .unwrap();
            if format == TrajectoryFormat::Xtc {
                frame_state.set_cell(Some(cell(
                    [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]],
                    NANOMETER,
                )));
            }
        }
        _ => unreachable!("every supported codec is listed"),
    }
    frame
}

fn cell(vectors: [[f64; 3]; 3], unit: kekule::units::Unit) -> PeriodicCell {
    PeriodicCell::new(
        Quantity::new(vectors.map(|[x, y, z]| Vector3::new(x, y, z)), unit),
        [true; 3],
    )
    .unwrap()
}

/// A frame each codec rejects at its last check before writing, with the
/// exact rejection: XYZ cannot carry time, DCD cannot store a cell whose `a`
/// axis leaves x, TRR needs a time, and XTC a cell.
fn rejected_frame(
    format: TrajectoryFormat,
    topology: &Arc<Topology>,
) -> (FrameBuffer, TrajectoryIoError) {
    let mut frame = valid_frame(format, topology, 1);
    let mut frame_state = frame.frame_mut();
    let (kind, detail, frame_index) = match format {
        TrajectoryFormat::Xyz => {
            frame_state
                .set_time(Some(Quantity::new(1.0, PICOSECOND)))
                .unwrap();
            (
                TrajectoryCodecErrorKind::UnsupportedField,
                "XYZ cannot preserve time",
                Some(1),
            )
        }
        TrajectoryFormat::Dcd => {
            frame_state.set_cell(Some(cell(
                [[0.0, 10.0, 0.0], [-11.0, 2.0, 0.0], [3.0, 1.0, 12.0]],
                ANGSTROM,
            )));
            (
                TrajectoryCodecErrorKind::UnsupportedVariant,
                "DCD cells require a along positive x, b in the xy plane with positive y, and c with positive z",
                None,
            )
        }
        TrajectoryFormat::Trr => {
            frame_state.set_time(None).unwrap();
            (
                TrajectoryCodecErrorKind::InconsistentMetadata,
                "TRR requires explicit frame time",
                None,
            )
        }
        TrajectoryFormat::Xtc => {
            frame_state.set_cell(None);
            (
                TrajectoryCodecErrorKind::InconsistentMetadata,
                "XTC writer requires a periodic cell; XtcCellPolicy::ZeroMatrixAsAbsent writes a zero box instead",
                None,
            )
        }
        _ => unreachable!("every supported codec is listed"),
    };
    let error = codec_error(
        kind,
        TrajectoryIoOperation::WriteFrame,
        format,
        &label(format),
        frame_index,
        detail,
    );
    (frame, error)
}

fn partial_stream_error(
    format: TrajectoryFormat,
    operation: TrajectoryIoOperation,
) -> TrajectoryIoError {
    codec_error(
        TrajectoryCodecErrorKind::InvalidFrame,
        operation,
        format,
        &label(format),
        None,
        &format!("an earlier {format} write failed partway; the stream may be incomplete"),
    )
}

#[test]
fn every_codec_writer_refuses_to_finish_without_frames() {
    let topology = topology();
    for format in CODECS {
        let writer = open(format, Cursor::new(Vec::new()), &topology);
        assert_eq!(
            writer.finish_stream().map(Cursor::into_inner),
            Err(codec_error(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                format,
                &label(format),
                None,
                &format!("{format} production trajectories must contain at least one frame"),
            )),
            "{format}"
        );
    }
}

#[test]
fn every_codec_writer_rejects_a_frame_without_writing_any_of_it() {
    let topology = topology();
    for format in CODECS {
        let mut reference = open(format, Cursor::new(Vec::new()), &topology);
        let mut subject = open(format, Cursor::new(Vec::new()), &topology);
        let (rejected, expected) = rejected_frame(format, &topology);
        for index in [0, 1] {
            let frame = valid_frame(format, &topology, index);
            reference.write_frame(frame.frame_view()).unwrap();
            subject.write_frame(frame.frame_view()).unwrap();
            if index == 0 {
                assert_eq!(
                    subject.write_frame(rejected.frame_view()),
                    Err(expected.clone()),
                    "{format}"
                );
            }
        }
        let reference = reference.finish_stream().unwrap().into_inner();
        let subject = subject.finish_stream().unwrap().into_inner();
        assert_eq!(subject, reference, "{format}");
    }
}

#[test]
fn every_codec_writer_refuses_to_continue_after_a_write_cut_off_partway() {
    let topology = topology();
    for format in CODECS {
        let mut probe = open(format, Cursor::new(Vec::new()), &topology);
        probe
            .write_frame(valid_frame(format, &topology, 0).frame_view())
            .unwrap();
        let first_frame_end = probe.sink().get_ref().len();
        // Three bytes into the second frame, the sink fails.
        let budget = first_frame_end + 3;
        let mut writer = open(format, ByteBudgetSink::new(budget), &topology);
        writer
            .write_frame(valid_frame(format, &topology, 0).frame_view())
            .unwrap();
        let io_failure: TrajectoryIoError = TrajectoryIoErrorContext::new(
            TrajectoryIoOperation::WriteFrame,
            ErrorKind::Other,
            ByteBudgetSink::INJECTED_FAILURE,
        )
        .with_format(format)
        .with_source_label(label(format))
        .into();
        assert_eq!(
            writer.write_frame(valid_frame(format, &topology, 1).frame_view()),
            Err(io_failure),
            "{format}"
        );
        assert_eq!(writer.sink().bytes().len(), budget, "{format}");
        assert_eq!(
            writer.write_frame(valid_frame(format, &topology, 1).frame_view()),
            Err(partial_stream_error(
                format,
                TrajectoryIoOperation::WriteFrame
            )),
            "{format}"
        );
        assert_eq!(writer.sink().bytes().len(), budget, "{format}");
        assert_eq!(
            writer.finish_stream().map(|sink| sink.bytes().len()),
            Err(partial_stream_error(format, TrajectoryIoOperation::Finish)),
            "{format}"
        );
    }
}
