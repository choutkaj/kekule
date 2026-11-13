//! One open file read by several readers at once: [`SharedFile`].

use std::fs::{File, Metadata};
use std::io::{self, Read, Seek, SeekFrom};
use std::sync::Arc;

/// A reading position in an open file that other positions share.
///
/// Every read is positional, so the positions never move each other and
/// readers on several threads read the same file object at once. Readers
/// that share it read the same bytes, whatever the path names later.
#[derive(Debug)]
pub(crate) struct SharedFile {
    file: Arc<File>,
    position: u64,
}

impl SharedFile {
    /// The start of `file`.
    pub(crate) fn new(file: File) -> Self {
        Self {
            file: Arc::new(file),
            position: 0,
        }
    }

    /// Another position in the same file, at `position`.
    pub(crate) fn at(&self, position: u64) -> Self {
        Self {
            file: Arc::clone(&self.file),
            position,
        }
    }

    /// Where the next read starts.
    pub(crate) const fn position(&self) -> u64 {
        self.position
    }

    /// Whether both read the same file object.
    pub(crate) fn same_file(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.file, &other.file)
    }

    /// The file's metadata now.
    pub(crate) fn metadata(&self) -> io::Result<Metadata> {
        self.file.metadata()
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
}

impl Read for SharedFile {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        #[cfg(test)]
        READS.with(|reads| reads.set(reads.get() + 1));
        let read = read_at(&self.file, buffer, self.position)?;
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
