#![allow(dead_code)]

use std::io::{self, BufRead, Cursor, Read, Seek, SeekFrom};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::Arc;

use molecular::trajectory::FrameBuffer;

pub fn buffer_snapshot(buffer: &FrameBuffer) -> String {
    format!("{buffer:#?}")
}

#[derive(Clone)]
pub struct RestoreSeekControl {
    position: Arc<AtomicU64>,
    target: Arc<AtomicU64>,
    armed: Arc<AtomicBool>,
}

impl RestoreSeekControl {
    pub fn arm_at_current_position(&self) {
        self.target
            .store(self.position.load(Ordering::SeqCst), Ordering::SeqCst);
        self.armed.store(true, Ordering::SeqCst);
    }
}

pub struct RestoreSeekFailure {
    inner: Cursor<Vec<u8>>,
    control: RestoreSeekControl,
}

impl RestoreSeekFailure {
    pub fn new(bytes: Vec<u8>) -> (Self, RestoreSeekControl) {
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
pub struct GuardControl {
    violation: Arc<AtomicBool>,
    probed_bytes: Arc<AtomicU64>,
}

impl GuardControl {
    pub fn violated(&self) -> bool {
        self.violation.load(Ordering::SeqCst)
    }

    pub fn probed_bytes(&self) -> u64 {
        self.probed_bytes.load(Ordering::SeqCst)
    }
}

pub struct GuardedCursor {
    inner: Cursor<Vec<u8>>,
    guard_offset: u64,
    control: GuardControl,
}

impl GuardedCursor {
    pub fn new(bytes: Vec<u8>, guard_offset: u64) -> (Self, GuardControl) {
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
pub struct NoBackwardSeekControl {
    armed: Arc<AtomicBool>,
}

impl NoBackwardSeekControl {
    pub fn arm(&self) {
        self.armed.store(true, Ordering::SeqCst);
    }
}

pub struct NoBackwardSeekCursor {
    inner: Cursor<Vec<u8>>,
    control: NoBackwardSeekControl,
}

impl NoBackwardSeekCursor {
    pub fn new(bytes: Vec<u8>) -> (Self, NoBackwardSeekControl) {
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
