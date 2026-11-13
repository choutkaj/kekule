#![allow(dead_code)]

use std::cell::RefCell;
use std::fmt::Debug;
use std::io::{self, BufRead, Cursor, Read, Seek, SeekFrom, Write};
use std::rc::Rc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;

use kekule::math::PeriodicCell;
use kekule::molecule::{Atom, BondOrder, Element};
use kekule::structure::{
    Forces, Geometry, GeometryError, GeometryView, Order, Place, Positions, StructureSet,
    StructureSetRead, StructureSetReadError, TimeAxis, Velocities,
};
use kekule::topology::{Topology, TopologyBuilder};
use kekule::traj::io::dcd::DcdEndian;
use kekule::traj::{
    FrameBuffer, FrameView, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader,
};
use kekule::units::{Quantity, CANONICAL_TIME_UNIT};

/// A frame a test expects or a reader decodes: a geometry with its time and
/// MD step, compared as a whole.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Frame {
    pub(crate) geometry: Geometry,
    pub(crate) time: Option<Quantity<f64>>,
    pub(crate) step: Option<u64>,
}

impl Frame {
    /// A frame of `positions` without a cell, time, or step.
    pub(crate) fn new(positions: Positions) -> Self {
        Self::of_sample(Geometry::new(positions))
    }

    /// A frame of `geometry` without a time or step.
    pub(crate) fn of_sample(geometry: Geometry) -> Self {
        Self {
            geometry,
            time: None,
            step: None,
        }
    }

    /// The frame a buffer holds.
    pub(crate) fn of(buffer: &FrameBuffer) -> Self {
        Self::of_view(buffer.frame_view())
    }

    /// The frame a frame view shows.
    pub(crate) fn of_view(frame: FrameView<'_>) -> Self {
        Self {
            geometry: frame.payload().clone(),
            time: frame.time(),
            step: frame.step(),
        }
    }

    /// A geometry of a set as a frame: the time and step of a time order.
    pub(crate) fn of_set_sample(view: GeometryView<'_>) -> Self {
        let (time, step) = match view.place() {
            Some(Place::Time { time, step }) => (time, step),
            _ => (None, None),
        };
        Self {
            geometry: view.payload().clone(),
            time,
            step,
        }
    }

    #[must_use]
    pub(crate) fn with_cell(mut self, cell: kekule::math::PeriodicCell) -> Self {
        self.geometry = self.geometry.with_cell(cell);
        self
    }

    /// Sets the time in the canonical time unit, as frames hold it; a time
    /// that is not a finite time fails, as for a frame buffer.
    pub(crate) fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), GeometryError> {
        self.time = time
            .map(|time| {
                let time = time.into_unit(CANONICAL_TIME_UNIT)?;
                if time.value().is_finite() {
                    Ok(time)
                } else {
                    Err(GeometryError::NonFiniteTime)
                }
            })
            .transpose()?;
        Ok(())
    }

    pub(crate) fn set_step(&mut self, step: Option<u64>) {
        self.step = step;
    }

    pub(crate) const fn time(&self) -> Option<Quantity<f64>> {
        self.time
    }

    pub(crate) const fn step(&self) -> Option<u64> {
        self.step
    }

    /// The frame with the annotation rows of `topology`, as a set or a
    /// buffer binds a geometry.
    #[must_use]
    pub(crate) fn bound(self, topology: &Arc<Topology>) -> Self {
        let mut geometry = self.geometry;
        geometry.bind(topology).unwrap();
        Self { geometry, ..self }
    }

    /// The bound frame as a writer takes it.
    pub(crate) fn view<'a>(&'a self, topology: &'a Arc<Topology>) -> FrameView<'a> {
        FrameView::new(topology, &self.geometry, self.time, self.step).unwrap()
    }
}

impl std::ops::Deref for Frame {
    type Target = Geometry;

    fn deref(&self) -> &Geometry {
        &self.geometry
    }
}

impl std::ops::DerefMut for Frame {
    fn deref_mut(&mut self) -> &mut Geometry {
        &mut self.geometry
    }
}

/// Frames in a time order whose columns the first frame's time and step
/// give, as a reader of them returns them.
pub(crate) fn time_ordered(
    topology: impl Into<Arc<Topology>>,
    frames: impl IntoIterator<Item = Frame>,
) -> StructureSet {
    let mut set = StructureSet::new(topology);
    for (index, frame) in frames.into_iter().enumerate() {
        if index == 0 {
            let axis = TimeAxis::new(
                frame
                    .time
                    .map(|_| Quantity::new(Vec::new(), CANONICAL_TIME_UNIT)),
                frame.step.map(|_| Vec::new()),
            )
            .unwrap();
            set.set_order(Some(Order::Time(axis))).unwrap();
        }
        let place = Place::Time {
            time: frame.time,
            step: frame.step,
        };
        set.push(frame.geometry, Some(place)).unwrap();
    }
    set
}

/// The frames of a set, in set order.
pub(crate) fn frames_of(set: &StructureSet) -> Vec<Frame> {
    set.iter().map(Frame::of_set_sample).collect()
}

/// FNV-1a 64 over `bytes`.
pub(crate) fn fingerprint_bytes(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// FNV-1a 64 over the `Debug` text of `value`. `Debug` writes every `f64`
/// as the shortest text that reads back to the same bits.
pub(crate) fn fingerprint(value: &impl Debug) -> u64 {
    fingerprint_bytes(format!("{value:?}").as_bytes())
}

/// Checks a sequence of results against their recorded fingerprints, in
/// order, and that every recorded fingerprint was checked.
pub(crate) struct Recorded {
    expected: &'static [u64],
    next: usize,
}

impl Recorded {
    pub(crate) const fn new(expected: &'static [u64]) -> Self {
        Self { expected, next: 0 }
    }

    pub(crate) fn check(&mut self, value: &impl Debug) {
        assert_eq!(
            Some(&fingerprint(value)),
            self.expected.get(self.next),
            "recorded result {}",
            self.next
        );
        self.next += 1;
    }
}

impl Drop for Recorded {
    fn drop(&mut self) {
        if !std::thread::panicking() {
            assert_eq!(self.next, self.expected.len(), "unchecked recorded results");
        }
    }
}

/// What a set holds apart from its topology: its properties, order,
/// distribution, and geometries.
pub(crate) fn contents(set: &StructureSet) -> impl Debug + '_ {
    (
        set.properties(),
        set.order(),
        set.distribution(),
        set.iter()
            .map(|geometry| geometry.payload())
            .collect::<Vec<_>>(),
    )
}

/// The positions, cell, velocities, and forces of a geometry.
pub(crate) type Dynamics = (
    Positions,
    Option<PeriodicCell>,
    Option<Velocities>,
    Option<Forces>,
);

/// The dynamics of every geometry: the state periodic operations and fits
/// change.
pub(crate) fn dynamics(set: &StructureSet) -> Vec<Dynamics> {
    set.iter()
        .map(|geometry| {
            (
                geometry.positions().clone(),
                geometry.cell().copied(),
                geometry.velocities().cloned(),
                geometry.forces().cloned(),
            )
        })
        .collect()
}

/// Every geometry `reader` reads, in order.
pub(crate) fn reads(
    reader: &(impl StructureSetRead + ?Sized),
) -> Vec<Result<Geometry, StructureSetReadError>> {
    (0..reader.len()).map(|index| reader.read(index)).collect()
}

/// Reject another topology context without changing its buffer or consuming input.
/// Call before reading and checking the first expected frame in the caller.
pub(crate) fn assert_rejects_unrelated_buffer(reader: &mut impl TrajectoryReader) {
    let mut wrong = FrameBuffer::new(linear_carbon_topology(reader.topology().atom_count()));
    wrong.set_step(Some(99));
    let before = buffer_snapshot(&wrong);
    assert!(matches!(
        reader.read_next(&mut wrong),
        Err(error) if error.kind() == TrajectoryIoErrorKind::TopologyMismatch
    ));
    assert_eq!(buffer_snapshot(&wrong), before);
    assert!(Arc::ptr_eq(
        &reader.shared_topology(),
        &reader.frame_buffer().shared_topology()
    ));
}

pub(crate) fn topology(symbols: &[&str], bonds: &[(usize, usize)]) -> Arc<Topology> {
    let mut graph = kekule::molecule::MoleculeEditor::new();
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

pub(crate) fn codec_kind(error: &TrajectoryIoError) -> Option<TrajectoryIoErrorKind> {
    Some(error.kind())
}

/// The exact codec error a reader or writer reports, for whole-value assertions.
pub(crate) fn codec_error(
    kind: TrajectoryIoErrorKind,
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> TrajectoryIoError {
    let mut error = TrajectoryIoError::codec(kind, operation, Some(format))
        .with_source_label(source_label)
        .with_detail(detail);
    if let Some(frame) = frame {
        error = error.with_frame(frame);
    }
    error
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
    expected: &[Frame],
) {
    let mut destination = FrameBuffer::new(reader.shared_topology());
    assert_eq!(reader.frame_count(), Some(expected.len()));
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(Frame::of(&destination), expected[0]);
    for index in (0..expected.len()).rev() {
        reader.read_frame(index, &mut destination).unwrap();
        assert_eq!(Frame::of(&destination), expected[index]);
    }
    for frame in &expected[1..] {
        assert!(reader.read_next(&mut destination).unwrap());
        assert_eq!(&Frame::of(&destination), frame);
    }
    assert!(!reader.read_next(&mut destination).unwrap());
}

/// Every frame a sequential reader decodes, in order.
pub(crate) fn read_all(reader: &mut impl TrajectoryReader) -> Vec<Frame> {
    let mut destination = FrameBuffer::new(reader.shared_topology());
    let mut frames = Vec::new();
    while reader.read_next(&mut destination).unwrap() {
        frames.push(Frame::of(&destination));
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

/// A DCD of three atoms, atoms 1 and 3 fixed, in two frames: frame 0
/// holds every atom at x = y = z = 0, 1, and 2 Å, and frame 1 holds the
/// free atom alone at (10, 11, 12) Å.
pub(crate) fn fixed_atom_fixture(endian: DcdEndian) -> Vec<u8> {
    let i32_bytes = |value: i32| match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
    };
    let f32_bytes = |value: f32| match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
    };
    let mut bytes = Vec::new();
    let mut record = |payload: &[u8]| {
        bytes.extend(i32_bytes(payload.len() as i32));
        bytes.extend(payload);
        bytes.extend(i32_bytes(payload.len() as i32));
    };

    let mut header = [0_u8; 84];
    header[..4].copy_from_slice(b"CORD");
    let mut controls = [0_i32; 20];
    controls[0] = 2;
    controls[2] = 1;
    controls[8] = 2;
    controls[19] = 24;
    for (chunk, value) in header[4..].as_chunks_mut::<4>().0.iter_mut().zip(controls) {
        chunk.copy_from_slice(&i32_bytes(value));
    }
    header[40..44].copy_from_slice(&f32_bytes(1.0));
    record(&header);
    let mut title = [0_u8; 84];
    title[..4].copy_from_slice(&i32_bytes(1));
    title[4..9].copy_from_slice(b"fixed");
    record(&title);
    record(&i32_bytes(3));
    record(&i32_bytes(2));
    for values in [
        vec![0.0_f32, 1.0, 2.0],
        vec![0.0, 1.0, 2.0],
        vec![0.0, 1.0, 2.0],
        vec![10.0],
        vec![11.0],
        vec![12.0],
    ] {
        let payload = values.into_iter().flat_map(f32_bytes).collect::<Vec<_>>();
        record(&payload);
    }
    bytes
}
