#![allow(dead_code)]

use std::cell::RefCell;
use std::io::{self, BufRead, Cursor, Read, Seek, SeekFrom, Write};
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element};
use kekule::structure::TrajectoryFrame;
use kekule::topology::{Topology, TopologyBuilder};
use kekule_traj::{
    FrameBuffer, SeekableTrajectoryReader, TrajectoryCodecErrorContext, TrajectoryCodecErrorKind,
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoOperation, TrajectoryReader,
};

/// Reject another topology context without changing its buffer or consuming input.
/// Call before reading and checking the first expected frame in the caller.
pub(crate) fn assert_rejects_unrelated_buffer(reader: &mut impl TrajectoryReader) {
    let mut wrong = FrameBuffer::new(linear_carbon_topology(reader.topology().atom_count()));
    wrong.frame_mut().set_step(Some(99));
    let before = buffer_snapshot(&wrong);
    assert!(matches!(
        reader.read_next(&mut wrong),
        Err(TrajectoryIoError::TopologyMismatch)
    ));
    assert_eq!(buffer_snapshot(&wrong), before);
    assert!(Arc::ptr_eq(
        &reader.shared_topology(),
        &reader.frame_buffer().shared_topology()
    ));
}

pub(crate) fn topology(symbols: &[&str], bonds: &[(usize, usize)]) -> Arc<Topology> {
    let mut graph = kekule::core::MoleculeEditor::new();
    let atoms = symbols
        .iter()
        .map(|symbol| {
            graph
                .add_atom(Atom::new(Element::from_symbol(symbol).unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    for &(source, target) in bonds {
        graph
            .add_bond(atoms[source], atoms[target], BondOrder::Single)
            .unwrap();
    }

    let molecule = graph.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

pub(crate) fn linear_carbon_topology(atom_count: usize) -> Arc<Topology> {
    let symbols = vec!["C"; atom_count];
    let bonds = (1..atom_count)
        .map(|index| (index - 1, index))
        .collect::<Vec<_>>();
    topology(&symbols, &bonds)
}

pub(crate) fn codec_kind(error: &TrajectoryIoError) -> Option<TrajectoryCodecErrorKind> {
    match error {
        TrajectoryIoError::Codec(context) => Some(context.kind()),
        _ => None,
    }
}

/// The exact codec error a reader or writer reports, for whole-value assertions.
pub(crate) fn codec_error(
    kind: TrajectoryCodecErrorKind,
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
    frame: Option<u64>,
    detail: &str,
) -> TrajectoryIoError {
    let mut context = TrajectoryCodecErrorContext::new(kind, operation, Some(format))
        .with_source_label(source_label)
        .with_detail(detail);
    if let Some(frame) = frame {
        context = context.with_frame(frame);
    }
    context.into()
}

pub(crate) fn x_coordinates(buffer: &FrameBuffer) -> Vec<f64> {
    buffer
        .positions()
        .values()
        .value()
        .iter()
        .map(|point| point.x)
        .collect()
}

pub(crate) fn buffer_snapshot(buffer: &FrameBuffer) -> String {
    format!("{buffer:#?}")
}

/// Bytes ahead of a trajectory, so a reader opens partway into its stream.
pub(crate) const LEADING_BYTES: &[u8] = b"bytes before the trajectory";

/// A stream holding [`LEADING_BYTES`] then `trajectory`, positioned at its first frame.
pub(crate) fn cursor_after_leading_bytes(trajectory: &[u8]) -> Cursor<Vec<u8>> {
    let mut bytes = LEADING_BYTES.to_vec();
    bytes.extend_from_slice(trajectory);
    let mut cursor = Cursor::new(bytes);
    cursor.set_position(LEADING_BYTES.len() as u64);
    cursor
}

/// Assert that an indexed reader holds exactly `expected`, maps index `i` to
/// frame `i` before and after a sequential read, and resumes sequential reads
/// where they stopped after random reads.
pub(crate) fn assert_indexes_every_frame_from_its_start(
    reader: &mut impl SeekableTrajectoryReader,
    expected: &[TrajectoryFrame],
) {
    let mut destination = FrameBuffer::new(reader.shared_topology());
    assert_eq!(reader.frame_count(), Some(expected.len() as u64));
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(destination.frame_view().payload(), &expected[0]);
    for index in (0..expected.len()).rev() {
        reader.read_frame(index as u64, &mut destination).unwrap();
        assert_eq!(destination.frame_view().payload(), &expected[index]);
    }
    for frame in &expected[1..] {
        assert!(reader.read_next(&mut destination).unwrap());
        assert_eq!(destination.frame_view().payload(), frame);
    }
    assert!(!reader.read_next(&mut destination).unwrap());
}

/// Every frame a sequential reader decodes, in order.
pub(crate) fn read_all(reader: &mut impl TrajectoryReader) -> Vec<TrajectoryFrame> {
    let mut destination = FrameBuffer::new(reader.shared_topology());
    let mut frames = Vec::new();
    while reader.read_next(&mut destination).unwrap() {
        frames.push(destination.frame_view().payload().clone());
    }
    frames
}

/// A seekable test stream whose backing bytes can change between reader calls.
#[derive(Clone)]
pub(crate) struct SharedCursor(Rc<RefCell<Cursor<Vec<u8>>>>);

impl SharedCursor {
    pub(crate) fn new(bytes: Vec<u8>) -> Self {
        Self(Rc::new(RefCell::new(Cursor::new(bytes))))
    }

    pub(crate) fn overwrite(&self, offset: usize, bytes: &[u8]) {
        self.0.borrow_mut().get_mut()[offset..offset + bytes.len()].copy_from_slice(bytes);
    }
}

impl Read for SharedCursor {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        self.0.borrow_mut().read(buffer)
    }
}

impl Seek for SharedCursor {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.0.borrow_mut().seek(position)
    }
}

#[derive(Clone)]
pub(crate) struct RestoreSeekControl {
    position: Arc<AtomicU64>,
    target: Arc<AtomicU64>,
    armed: Arc<AtomicBool>,
}

impl RestoreSeekControl {
    pub(crate) fn arm_at_current_position(&self) {
        self.target
            .store(self.position.load(Ordering::SeqCst), Ordering::SeqCst);
        self.armed.store(true, Ordering::SeqCst);
    }
}

pub(crate) struct RestoreSeekFailure {
    inner: Cursor<Vec<u8>>,
    control: RestoreSeekControl,
}

impl RestoreSeekFailure {
    pub(crate) fn new(bytes: Vec<u8>) -> (Self, RestoreSeekControl) {
        let control = RestoreSeekControl {
            position: Arc::new(AtomicU64::new(0)),
            target: Arc::new(AtomicU64::new(u64::MAX)),
            armed: Arc::new(AtomicBool::new(false)),
        };
        (
            Self {
                inner: Cursor::new(bytes),
                control: control.clone(),
            },
            control,
        )
    }

    fn update_position(&self) {
        self.control
            .position
            .store(self.inner.position(), Ordering::SeqCst);
    }
}

impl Read for RestoreSeekFailure {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.update_position();
        Ok(read)
    }
}

impl Seek for RestoreSeekFailure {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        if let SeekFrom::Start(target) = position {
            if self.control.armed.load(Ordering::SeqCst)
                && target == self.control.target.load(Ordering::SeqCst)
            {
                self.control.armed.store(false, Ordering::SeqCst);
                return Err(io::Error::other("injected restoration seek failure"));
            }
        }
        let position = self.inner.seek(position)?;
        self.update_position();
        Ok(position)
    }
}

#[derive(Clone)]
pub(crate) struct GuardControl {
    violation: Arc<AtomicBool>,
    probed_bytes: Arc<AtomicU64>,
}

impl GuardControl {
    pub(crate) fn violated(&self) -> bool {
        self.violation.load(Ordering::SeqCst)
    }

    pub(crate) fn probed_bytes(&self) -> u64 {
        self.probed_bytes.load(Ordering::SeqCst)
    }
}

pub(crate) struct GuardedCursor {
    inner: Cursor<Vec<u8>>,
    guard_offset: u64,
    control: GuardControl,
}

impl GuardedCursor {
    pub(crate) fn new(bytes: Vec<u8>, guard_offset: u64) -> (Self, GuardControl) {
        let control = GuardControl {
            violation: Arc::new(AtomicBool::new(false)),
            probed_bytes: Arc::new(AtomicU64::new(0)),
        };
        (
            Self {
                inner: Cursor::new(bytes),
                guard_offset,
                control: control.clone(),
            },
            control,
        )
    }
}

impl Read for GuardedCursor {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        if self.inner.position() >= self.guard_offset {
            if buffer.len() != 1 {
                self.control.violation.store(true, Ordering::SeqCst);
                return Err(io::Error::other(
                    "attempted to decode beyond the guarded frame boundary",
                ));
            }
            self.control.probed_bytes.fetch_add(1, Ordering::SeqCst);
        }
        self.inner.read(buffer)
    }
}

impl BufRead for GuardedCursor {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.inner.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        if amount > 0 && self.inner.position() >= self.guard_offset {
            self.control.violation.store(true, Ordering::SeqCst);
        }
        self.inner.consume(amount);
    }
}

impl Seek for GuardedCursor {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.inner.seek(position)
    }
}

#[derive(Clone)]
pub(crate) struct NoBackwardSeekControl {
    armed: Arc<AtomicBool>,
}

impl NoBackwardSeekControl {
    pub(crate) fn arm(&self) {
        self.armed.store(true, Ordering::SeqCst);
    }
}

pub(crate) struct NoBackwardSeekCursor {
    inner: Cursor<Vec<u8>>,
    control: NoBackwardSeekControl,
}

impl NoBackwardSeekCursor {
    pub(crate) fn new(bytes: Vec<u8>) -> (Self, NoBackwardSeekControl) {
        let control = NoBackwardSeekControl {
            armed: Arc::new(AtomicBool::new(false)),
        };
        (
            Self {
                inner: Cursor::new(bytes),
                control: control.clone(),
            },
            control,
        )
    }
}

impl Read for NoBackwardSeekCursor {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        self.inner.read(buffer)
    }
}

impl Seek for NoBackwardSeekCursor {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        if let SeekFrom::Start(target) = position {
            if self.control.armed.load(Ordering::SeqCst) && target < self.inner.position() {
                return Err(io::Error::other(
                    "ordinary frame decoding attempted a backward EOF-probe seek",
                ));
            }
        }
        self.inner.seek(position)
    }
}

/// A seekable sink that accepts `budget` bytes in total and then fails every
/// write, so a writer can be cut off partway through a frame.
pub(crate) struct ByteBudgetSink {
    inner: Cursor<Vec<u8>>,
    remaining: usize,
}

impl ByteBudgetSink {
    pub(crate) const INJECTED_FAILURE: &'static str = "injected write failure";

    pub(crate) fn new(budget: usize) -> Self {
        Self {
            inner: Cursor::new(Vec::new()),
            remaining: budget,
        }
    }

    pub(crate) fn bytes(&self) -> &[u8] {
        self.inner.get_ref()
    }
}

impl Write for ByteBudgetSink {
    fn write(&mut self, buffer: &[u8]) -> io::Result<usize> {
        if self.remaining == 0 && !buffer.is_empty() {
            return Err(io::Error::other(Self::INJECTED_FAILURE));
        }
        let accepted = buffer.len().min(self.remaining);
        let written = self.inner.write(&buffer[..accepted])?;
        self.remaining -= written;
        Ok(written)
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Seek for ByteBudgetSink {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.inner.seek(position)
    }
}
