//! One open file read by several readers at once: [`SharedFile`], through
//! the [`FileSlot`] they share.

use std::fs::{File, Metadata};
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::{Arc, PoisonError, RwLock, RwLockReadGuard, RwLockWriteGuard};

/// The one handle of a file that readers share, and the path it was opened
/// from. A host can release the handle and reopen the file, at the same
/// path or another, for every reader at once.
///
/// A read of a frame, and a refresh, [enters](Self::enter) the slot: it
/// holds the path shared while it reads. Releasing and reopening hold it
/// exclusively, so they wait for the reads running and no read sees the
/// file close halfway. The handle has a lock of its own, taken shared by
/// each positional read, so that a refresh can hand a new handle to the
/// readers while they read.
#[derive(Debug)]
pub(crate) struct FileSlot {
    path: RwLock<PathBuf>,
    file: RwLock<Option<File>>,
}

fn shared<T>(lock: &RwLock<T>) -> RwLockReadGuard<'_, T> {
    lock.read().unwrap_or_else(PoisonError::into_inner)
}

fn exclusive<T>(lock: &RwLock<T>) -> RwLockWriteGuard<'_, T> {
    lock.write().unwrap_or_else(PoisonError::into_inner)
}

impl FileSlot {
    /// The path of the open file while it is open, held until the guard
    /// drops; `None` once released. A thread that holds the guard must not
    /// enter again, since a waiting release would block it.
    pub(crate) fn enter(&self) -> Option<RwLockReadGuard<'_, PathBuf>> {
        let path = shared(&self.path);
        shared(&self.file).is_some().then_some(path)
    }

    /// The path the file was opened from, or last reopened at.
    pub(crate) fn path(&self) -> PathBuf {
        shared(&self.path).clone()
    }

    /// Closes the handle, once the reads running finish. Releasing twice
    /// does nothing.
    pub(crate) fn release(&self) {
        let _path = exclusive(&self.path);
        drop(exclusive(&self.file).take());
    }

    /// Takes the file `open` opens at `path` as the handle, and `path` as
    /// the file's path, once the reads running finish. When `open` fails,
    /// nothing changes.
    pub(crate) fn reopen<E>(
        &self,
        path: &Path,
        open: impl FnOnce(&Path) -> Result<File, E>,
    ) -> Result<(), E> {
        let mut place = exclusive(&self.path);
        let file = open(path)?;
        *exclusive(&self.file) = Some(file);
        path.clone_into(&mut place);
        Ok(())
    }

    /// Takes `file` as the handle; for a refresh that entered the slot.
    fn install(&self, file: File) {
        *exclusive(&self.file) = Some(file);
    }

    fn take(&self) -> Option<File> {
        exclusive(&self.file).take()
    }

    fn read_at(&self, buffer: &mut [u8], position: u64) -> io::Result<usize> {
        match &*shared(&self.file) {
            Some(file) => read_at(file, buffer, position),
            None => Err(released()),
        }
    }

    fn metadata(&self) -> io::Result<Metadata> {
        match &*shared(&self.file) {
            Some(file) => file.metadata(),
            None => Err(released()),
        }
    }
}

/// The failure of a read of a released file. Reads enter the slot first,
/// and fail as released there, so none reads a released file.
fn released() -> io::Error {
    io::Error::other("The file is closed.")
}

/// A reading position in an open file that other positions share.
///
/// Every read is positional, so the positions never move each other and
/// readers on several threads read the same file object at once. Readers
/// that share its slot read the same handle, whatever the path names later,
/// and follow the slot when a host releases or reopens it.
#[derive(Debug)]
pub(crate) struct SharedFile {
    slot: Arc<FileSlot>,
    position: u64,
}

impl SharedFile {
    /// The start of `file`, opened from `path`, in a slot of its own.
    pub(crate) fn new(file: File, path: PathBuf) -> Self {
        Self {
            slot: Arc::new(FileSlot {
                path: RwLock::new(path),
                file: RwLock::new(Some(file)),
            }),
            position: 0,
        }
    }

    /// Another position in the same file, at `position`.
    pub(crate) fn at(&self, position: u64) -> Self {
        Self {
            slot: Arc::clone(&self.slot),
            position,
        }
    }

    /// The slot of the file.
    pub(crate) const fn slot(&self) -> &Arc<FileSlot> {
        &self.slot
    }

    /// Hands this position's file to `slot` as its handle, and reads
    /// through `slot` from now on.
    pub(crate) fn hand_over(&mut self, slot: &Arc<FileSlot>) {
        if let Some(file) = self.slot.take() {
            slot.install(file);
        }
        self.slot = Arc::clone(slot);
    }

    /// Where the next read starts.
    pub(crate) const fn position(&self) -> u64 {
        self.position
    }

    /// Whether both read through the same slot.
    pub(crate) fn same_file(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.slot, &other.slot)
    }

    /// The file's metadata now.
    pub(crate) fn metadata(&self) -> io::Result<Metadata> {
        self.slot.metadata()
    }

    /// The file's length now.
    pub(crate) fn len(&self) -> io::Result<u64> {
        self.metadata().map(|metadata| metadata.len())
    }
}

#[cfg(test)]
thread_local! {
    /// The reads of shared files on this thread, which tests count.
    pub(crate) static READS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
    /// Runs before each read of a shared file on this thread, for tests
    /// that pause a read.
    pub(crate) static BEFORE_READ: std::cell::RefCell<Option<Box<dyn FnMut()>>> =
        const { std::cell::RefCell::new(None) };
}

impl Read for SharedFile {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        #[cfg(test)]
        READS.with(|reads| reads.set(reads.get() + 1));
        #[cfg(test)]
        BEFORE_READ.with(|hook| {
            if let Some(hook) = hook.borrow_mut().as_mut() {
                hook();
            }
        });
        let read = self.slot.read_at(buffer, self.position)?;
        self.position += read as u64;
        Ok(read)
    }
}

impl Seek for SharedFile {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        let (base, offset) = match position {
            SeekFrom::Start(position) => {
                self.position = position;
                return Ok(position);
            }
            SeekFrom::Current(offset) => (self.position, offset),
            SeekFrom::End(offset) => (self.len()?, offset),
        };
        self.position = base.checked_add_signed(offset).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "The position is outside the file.",
            )
        })?;
        Ok(self.position)
    }

    fn stream_position(&mut self) -> io::Result<u64> {
        Ok(self.position)
    }
}

#[cfg(unix)]
fn read_at(file: &File, buffer: &mut [u8], position: u64) -> io::Result<usize> {
    std::os::unix::fs::FileExt::read_at(file, buffer, position)
}

#[cfg(windows)]
fn read_at(file: &File, buffer: &mut [u8], position: u64) -> io::Result<usize> {
    // This also moves the handle's own position, which no reader uses.
    std::os::windows::fs::FileExt::seek_read(file, buffer, position)
}

#[cfg(not(any(unix, windows)))]
fn read_at(mut file: &File, buffer: &mut [u8], position: u64) -> io::Result<usize> {
    static SEEK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _seek = SEEK
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    file.seek(SeekFrom::Start(position))?;
    file.read(buffer)
}
