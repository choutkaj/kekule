use std::path::Path;

use crate::traj::TrajectoryFormat;

use super::{dcd, trr, xtc, xyz, TrajectoryFormatHint, TrajectoryIoLimits};

/// Format-agnostic configuration for opening readers or loading entire trajectories.
///
/// The codec readers, such as [`xyz::XyzReader`], take these options too, or
/// their format's options alone through [`From`], and use their format's
/// options, the limits, and the source label; they ignore the format hint
/// and [`Self::growing`].
#[derive(Debug, Clone)]
pub struct TrajectoryOpenOptions {
    pub(super) format_hint: TrajectoryFormatHint,
    pub(super) limits: TrajectoryIoLimits,
    pub(super) source_label: Option<String>,
    pub(super) xyz: xyz::XyzReadOptions,
    pub(super) dcd: dcd::DcdReadOptions,
    pub(super) trr: trr::TrrReadOptions,
    pub(super) xtc: xtc::XtcReadOptions,
    pub(super) growing: bool,
}

impl Default for TrajectoryOpenOptions {
    fn default() -> Self {
        Self {
            format_hint: TrajectoryFormatHint::Auto,
            limits: TrajectoryIoLimits::default(),
            source_label: None,
            xyz: xyz::XyzReadOptions::default(),
            dcd: dcd::DcdReadOptions::default(),
            trr: trr::TrrReadOptions::default(),
            xtc: xtc::XtcReadOptions::default(),
            growing: false,
        }
    }
}

impl From<xyz::XyzReadOptions> for TrajectoryOpenOptions {
    fn from(options: xyz::XyzReadOptions) -> Self {
        Self::default().with_xyz_options(options)
    }
}

impl From<dcd::DcdReadOptions> for TrajectoryOpenOptions {
    fn from(options: dcd::DcdReadOptions) -> Self {
        Self::default().with_dcd_options(options)
    }
}

impl From<trr::TrrReadOptions> for TrajectoryOpenOptions {
    fn from(options: trr::TrrReadOptions) -> Self {
        Self::default().with_trr_options(options)
    }
}

impl From<xtc::XtcReadOptions> for TrajectoryOpenOptions {
    fn from(options: xtc::XtcReadOptions) -> Self {
        Self::default().with_xtc_options(options)
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

    /// Sets [`Self::source_label`].
    #[must_use]
    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.source_label = Some(source_label.into());
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

    /// Limits applied before allocation, scanning, and seeking.
    pub const fn limits(&self) -> &TrajectoryIoLimits {
        &self.limits
    }

    /// The name of the source in errors. Unset by default: path readers then
    /// name the path, and codec readers the format's stream, such as
    /// `xyz stream`.
    pub fn source_label(&self) -> Option<&str> {
        self.source_label.as_deref()
    }

    /// Options for the XYZ reader.
    pub const fn xyz_options(&self) -> &xyz::XyzReadOptions {
        &self.xyz
    }

    /// Options for the DCD reader.
    pub const fn dcd_options(&self) -> &dcd::DcdReadOptions {
        &self.dcd
    }

    /// Options for the TRR reader.
    pub const fn trr_options(&self) -> &trr::TrrReadOptions {
        &self.trr
    }

    /// Options for the XTC reader.
    pub const fn xtc_options(&self) -> &xtc::XtcReadOptions {
        &self.xtc
    }

    /// The label of a path reader's errors about `path`.
    pub(super) fn path_label(&self, path: &Path) -> String {
        self.source_label
            .clone()
            .unwrap_or_else(|| path.display().to_string())
    }

    /// The source a codec reader of `format` reads from.
    pub(super) fn read_source(&self, format: TrajectoryFormat) -> ReadSource {
        ReadSource {
            label: self
                .source_label
                .clone()
                .unwrap_or_else(|| stream_label(format)),
            limits: self.limits.clone(),
        }
    }
}

/// The label of a format's stream when the caller sets none.
pub(super) fn stream_label(format: TrajectoryFormat) -> String {
    format!("{} stream", format.name().to_ascii_lowercase())
}

/// What a codec reader names in errors and the limits it reads within.
#[derive(Debug, Clone)]
pub(crate) struct ReadSource {
    pub(crate) label: String,
    pub(crate) limits: TrajectoryIoLimits,
}
