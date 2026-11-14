use crate::diagnostics::Cancellation;

use super::frame_index::IndexProgressCallback;
use super::{
    dcd, trr, xtc, xyz, TrajectoryFormatHint, TrajectoryIndexProgress, TrajectoryIoLimits,
};

/// Format-agnostic configuration for opening readers or loading entire trajectories.
///
/// Path readers use these top-level limits and the path as their diagnostic
/// source label, overriding limits and source labels in format-specific options.
#[derive(Debug, Clone)]
pub struct TrajectoryOpenOptions {
    pub(super) format_hint: TrajectoryFormatHint,
    pub(super) limits: TrajectoryIoLimits,
    pub(super) xyz: xyz::XyzReadOptions,
    pub(super) dcd: dcd::DcdReadOptions,
    pub(super) trr: trr::TrrReadOptions,
    pub(super) xtc: xtc::XtcReadOptions,
    pub(super) index_cancellation: Option<Cancellation>,
    pub(super) index_progress: Option<IndexProgressCallback>,
    pub(super) growing: bool,
}

impl Default for TrajectoryOpenOptions {
    fn default() -> Self {
        Self {
            format_hint: TrajectoryFormatHint::Auto,
            limits: TrajectoryIoLimits::default(),
            xyz: xyz::XyzReadOptions::default(),
            dcd: dcd::DcdReadOptions::default(),
            trr: trr::TrrReadOptions::default(),
            xtc: xtc::XtcReadOptions::default(),
            index_cancellation: None,
            index_progress: None,
            growing: false,
        }
    }
}

impl TrajectoryOpenOptions {
    /// Sets [`Self::format_hint`].
    #[must_use]
    pub fn with_format_hint(mut self, format_hint: TrajectoryFormatHint) -> Self {
        self.format_hint = format_hint;
        self
    }

    /// Sets [`Self::limits`].
    #[must_use]
    pub fn with_limits(mut self, limits: TrajectoryIoLimits) -> Self {
        self.limits = limits;
        self
    }

    /// Sets [`Self::xyz_options`].
    #[must_use]
    pub fn with_xyz_options(mut self, options: xyz::XyzReadOptions) -> Self {
        self.xyz = options;
        self
    }

    /// Sets [`Self::dcd_options`].
    #[must_use]
    pub fn with_dcd_options(mut self, options: dcd::DcdReadOptions) -> Self {
        self.dcd = options;
        self
    }

    /// Sets [`Self::trr_options`].
    #[must_use]
    pub fn with_trr_options(mut self, options: trr::TrrReadOptions) -> Self {
        self.trr = options;
        self
    }

    /// Sets [`Self::xtc_options`].
    #[must_use]
    pub fn with_xtc_options(mut self, options: xtc::XtcReadOptions) -> Self {
        self.xtc = options;
        self
    }

    /// Sets [`Self::index_cancellation`].
    #[must_use]
    pub fn with_index_cancellation(mut self, cancellation: Cancellation) -> Self {
        self.index_cancellation = Some(cancellation);
        self
    }

    /// Sets the callback that [`open_indexed_trajectory`] calls with a
    /// [`TrajectoryIndexProgress`] before the first frame and after every
    /// verified frame while it builds its index, on the opening thread. It
    /// should be cheap, such as storing atomics a progress bar reads. Clones
    /// of these options share the callback. None by default.
    ///
    /// [`open_indexed_trajectory`]: super::open_indexed_trajectory
    #[must_use]
    pub fn with_index_progress(
        mut self,
        progress: impl Fn(TrajectoryIndexProgress) + Send + Sync + 'static,
    ) -> Self {
        self.index_progress = Some(IndexProgressCallback::new(progress));
        self
    }

    /// The cancellation that [`open_indexed_trajectory`] polls before the
    /// first frame and after every verified frame while it builds its index;
    /// once requested, opening fails with
    /// [`TrajectoryIoErrorKind::Cancelled`] and returns no reader. Other
    /// operations do not poll it. None by default.
    ///
    /// [`open_indexed_trajectory`]: super::open_indexed_trajectory
    /// [`TrajectoryIoErrorKind::Cancelled`]: crate::traj::TrajectoryIoErrorKind::Cancelled
    pub const fn index_cancellation(&self) -> Option<&Cancellation> {
        self.index_cancellation.as_ref()
    }

    /// Sets [`Self::growing`].
    #[must_use]
    pub const fn with_growing(mut self, growing: bool) -> Self {
        self.growing = growing;
        self
    }

    /// Whether the file may still be growing: a writer, such as a running
    /// simulation in another process, is appending frames to it. False by
    /// default.
    ///
    /// Indexed readers of a growing file take a partially written final
    /// frame as not yet written rather than as corruption: they index the
    /// complete frames before it and report the bytes after them as
    /// [`pending`](super::IndexedFileTrajectoryReader::pending_byte_count),
    /// and [`refresh`](super::IndexedFileTrajectoryReader::refresh) indexes
    /// it once it is complete. Corrupt data is still an error wherever it
    /// lies. DCD readers count the frames in the file rather than check
    /// them against the header's `NSET`, which writers update only as they
    /// flush or when they finish. XYZ lines count only once their line
    /// break is written, so a growing XYZ file's last line must end in one.
    /// A file that holds no complete frame yet, or for DCD no complete
    /// header, fails to open with [`TrajectoryIoErrorKind::Incomplete`].
    ///
    /// Without it, a partial final frame is an error, as for a complete
    /// file. Sequential readers and [`read_structure_set`](super::read_structure_set)
    /// ignore it.
    ///
    /// [`TrajectoryIoErrorKind::Incomplete`]: crate::traj::TrajectoryIoErrorKind::Incomplete
    pub const fn growing(&self) -> bool {
        self.growing
    }

    /// How the format is selected; [`TrajectoryFormatHint::Auto`] by default.
    pub const fn format_hint(&self) -> TrajectoryFormatHint {
        self.format_hint
    }

    /// Limits for path readers, which override those of the format options.
    pub const fn limits(&self) -> &TrajectoryIoLimits {
        &self.limits
    }

    /// Options for the XYZ reader; path readers override its limits
    /// and source label.
    pub const fn xyz_options(&self) -> &xyz::XyzReadOptions {
        &self.xyz
    }

    /// Options for the DCD reader; path readers override its limits
    /// and source label.
    pub const fn dcd_options(&self) -> &dcd::DcdReadOptions {
        &self.dcd
    }

    /// Options for the TRR reader; path readers override its limits
    /// and source label.
    pub const fn trr_options(&self) -> &trr::TrrReadOptions {
        &self.trr
    }

    /// Options for the XTC reader; path readers override its limits
    /// and source label.
    pub const fn xtc_options(&self) -> &xtc::XtcReadOptions {
        &self.xtc
    }
}
