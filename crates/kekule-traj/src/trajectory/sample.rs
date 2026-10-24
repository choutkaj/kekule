use std::fmt;
use std::sync::Arc;

use kekule::diagnostics::Cancellation;
use kekule::structure::{Ensemble, TrajectorySampler, TrajectorySampling};

use super::{
    SeekableTrajectoryReader, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
    TrajectoryReader,
};

/// Cancellation and progress reporting for [`read_ensemble`] and
/// [`read_ensemble_seekable`]. Neither is set by default.
#[derive(Debug, Clone, Default)]
pub struct EnsembleReadOptions {
    cancellation: Option<Cancellation>,
    progress: Option<ProgressCallback>,
}

impl EnsembleReadOptions {
    /// Sets [`Self::cancellation`].
    #[must_use]
    pub fn with_cancellation(mut self, cancellation: Cancellation) -> Self {
        self.cancellation = Some(cancellation);
        self
    }

    /// Sets the callback called with an [`EnsembleReadProgress`] before the
    /// first frame and after every frame read, on the reading thread. It
    /// should be cheap, such as storing atomics a progress bar reads. Clones
    /// of these options share the callback.
    #[must_use]
    pub fn with_progress(
        mut self,
        progress: impl Fn(EnsembleReadProgress) + Send + Sync + 'static,
    ) -> Self {
        self.progress = Some(ProgressCallback(Arc::new(progress)));
        self
    }

    /// The cancellation polled before the first frame and after every frame
    /// read; once requested, reading fails with
    /// [`TrajectoryIoErrorKind::Cancelled`] and returns no ensemble.
    pub const fn cancellation(&self) -> Option<&Cancellation> {
        self.cancellation.as_ref()
    }

    /// Reports `frame_count` of `total_frame_count` frames read, then polls
    /// the cancellation.
    fn checkpoint(
        &self,
        frame_count: usize,
        total_frame_count: usize,
    ) -> Result<(), TrajectoryIoError> {
        if let Some(progress) = &self.progress {
            (progress.0)(EnsembleReadProgress {
                frame_count,
                total_frame_count,
            });
        }
        if let Some(cancellation) = &self.cancellation {
            cancellation.check().map_err(|cancelled| {
                TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::Cancelled,
                    TrajectoryIoOperation::ReadFrame,
                    None,
                )
                .with_source(cancelled)
            })?;
        }
        Ok(())
    }
}

/// Progress of reading frames into an ensemble.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnsembleReadProgress {
    frame_count: usize,
    total_frame_count: usize,
}

impl EnsembleReadProgress {
    /// Frames read so far.
    pub const fn frame_count(&self) -> usize {
        self.frame_count
    }

    /// Frames the operation reads in all, which [`Self::frame_count`]
    /// reaches when it succeeds: through the last requested frame for
    /// [`read_ensemble`], and the distinct requested frames for
    /// [`read_ensemble_seekable`].
    pub const fn total_frame_count(&self) -> usize {
        self.total_frame_count
    }
}

/// A host's observer of reading progress; clones share one callback.
#[derive(Clone)]
struct ProgressCallback(Arc<dyn Fn(EnsembleReadProgress) + Send + Sync>);

impl fmt::Debug for ProgressCallback {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("ProgressCallback")
    }
}

// As for `Cancellation`: this crate only calls the callback, so a panic inside
// it cannot leave crate state broken, and options holding one stay unwind safe.
impl std::panic::UnwindSafe for ProgressCallback {}
impl std::panic::RefUnwindSafe for ProgressCallback {}

/// Reads frames sequentially into an ensemble of one state, as
/// [`kekule::structure::Trajectory::to_ensemble`] makes from a loaded
/// trajectory.
///
/// Frame 0 is the next frame the reader publishes. Reading stops after the
/// last requested frame, so frames after it are never decoded; frames before
/// it are decoded and skipped unless requested. Memory beyond the ensemble is
/// one [`FrameBuffer`](crate::FrameBuffer). For a reader with random access,
/// [`read_ensemble_seekable`] decodes only the requested frames.
///
/// ```
/// use std::sync::Arc;
/// use kekule::structure::{Positions, Trajectory, TrajectoryFrame, TrajectorySampling};
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::io::{open_trajectory, write_trajectory};
/// use kekule_traj::read_ensemble;
///
/// let topology = Arc::new(Topology::from_molecule(smiles::to_molecules("[He]")?.remove(0))?);
/// let frames = (0..10).map(|_| TrajectoryFrame::new(Positions::zeros(1)));
/// let trajectory = Trajectory::from_items(Arc::clone(&topology), frames)?;
/// let path = std::env::temp_dir().join(format!("kekule-traj-sample-{}.xyz", std::process::id()));
/// # let _ = std::fs::remove_file(&path);
/// write_trajectory(&path, &trajectory, Default::default())?;
///
/// // Leave out two frames of equilibration, then take every third frame.
/// let sampling = TrajectorySampling::new((2..10).step_by(3));
/// let mut reader = open_trajectory(&path, topology, Default::default())?;
/// let ensemble = read_ensemble(&mut reader, &sampling, Default::default());
/// std::fs::remove_file(&path)?;
/// assert_eq!(ensemble?.len(), 3);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::FrameIndexOutOfRange`] naming the
/// smallest requested frame when the stream ends before it,
/// [`TrajectoryIoErrorKind::Cancelled`] when `options` cancels,
/// [`TrajectoryIoErrorKind::InvalidTrajectory`] when the atom subset cannot be
/// published, and the reader's errors. No partial ensemble is returned.
pub fn read_ensemble<R: TrajectoryReader + ?Sized>(
    reader: &mut R,
    sampling: &TrajectorySampling,
    options: EnsembleReadOptions,
) -> Result<Ensemble, TrajectoryIoError> {
    let mut sampler = TrajectorySampler::new(reader.shared_topology(), sampling)?;
    let total = sampling.frames().iter().max().map_or(0, |last| last + 1);
    let mut buffer = reader.frame_buffer();
    options.checkpoint(0, total)?;
    for index in 0..total {
        if !reader.read_next(&mut buffer)? {
            let missing = sampler
                .frames()
                .find(|&frame| frame >= index)
                .expect("the last requested frame is at or after an unread frame");
            return Err(TrajectoryIoError::frame_index_out_of_range(missing));
        }
        sampler.add(index, buffer.frame_view())?;
        options.checkpoint(index + 1, total)?;
    }
    Ok(sampler.finish()?)
}

/// Reads only the requested frames, each once and in ascending order, into an
/// ensemble of one state, as [`read_ensemble`] does sequentially.
///
/// Frames are numbered as [`SeekableTrajectoryReader::read_frame`] numbers
/// them, and the sequential cursor is kept. Memory beyond the ensemble is one
/// [`FrameBuffer`](crate::FrameBuffer).
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::FrameIndexOutOfRange`] naming a requested
/// frame the reader does not have, before reading any frame when the reader
/// knows its frame count; otherwise as [`read_ensemble`].
pub fn read_ensemble_seekable<R: SeekableTrajectoryReader + ?Sized>(
    reader: &mut R,
    sampling: &TrajectorySampling,
    options: EnsembleReadOptions,
) -> Result<Ensemble, TrajectoryIoError> {
    let mut sampler = TrajectorySampler::new(reader.shared_topology(), sampling)?;
    let frames = sampler.frames().collect::<Vec<_>>();
    if let (Some(count), Some(&last)) = (reader.frame_count(), frames.last()) {
        if last >= count {
            let missing = frames
                .iter()
                .copied()
                .find(|&frame| frame >= count)
                .expect("the last requested frame is out of range");
            return Err(TrajectoryIoError::frame_index_out_of_range(missing));
        }
    }
    let mut buffer = reader.frame_buffer();
    options.checkpoint(0, frames.len())?;
    for (read, &index) in frames.iter().enumerate() {
        reader.read_frame(index, &mut buffer)?;
        sampler.add(index, buffer.frame_view())?;
        options.checkpoint(read + 1, frames.len())?;
    }
    Ok(sampler.finish()?)
}
