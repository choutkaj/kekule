use super::*;
use crate::io::xyz::write::{finite_in_unit, topology_omissions, write_frame};
use crate::traj::trajectory::reason::{Part, Reason};

/// Strict multi-frame XYZ writer over a caller-supplied stream.
pub struct XyzWriter<W> {
    pub(super) writer: W,
    pub(super) topology: Arc<Topology>,
    pub(super) options: XyzWriteOptions,
    pub(super) source_label: String,
    pub(super) frame_count: usize,
    /// Set while a frame's bytes are being written and left set if writing
    /// fails, so a stream that may end in a partial frame is never extended
    /// or reported as finished.
    pub(super) partial: bool,
    pub(super) report: WriteReport,
}

impl<W> fmt::Debug for XyzWriter<W> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("XyzWriter")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source_label)
            .field("frame_count", &self.frame_count)
            .field("partial", &self.partial)
            .finish_non_exhaustive()
    }
}

impl<W: Write> XyzWriter<W> {
    /// A writer of `topology` frames to `writer`. Nothing is written until the first frame.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::LimitExceeded`] for a topology without
    /// atoms.
    pub fn new(
        writer: W,
        topology: impl Into<Arc<Topology>>,
        options: impl Into<TrajectoryWriteOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source_label = options.stream_label(FORMAT);
        Self::open(writer, topology.into(), options.xyz, source_label)
    }

    /// Opens the writer like [`Self::new`] with the format options and
    /// label taken apart.
    pub(crate) fn open(
        writer: W,
        topology: Arc<Topology>,
        options: XyzWriteOptions,
        source_label: String,
    ) -> Result<Self, TrajectoryIoError> {
        if topology.atom_count() == 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                Reason::TopologyHasNoAtoms,
            ));
        }
        Ok(Self {
            writer,
            options,
            source_label,
            frame_count: 0,
            partial: false,
            report: topology_omissions(&topology),
            topology,
        })
    }

    /// What the stream written so far leaves out: from construction, the
    /// topology content [`kekule::xyz::write`](crate::xyz::write) reports,
    /// and then the frame state omitted under [`XyzWriteOptions::cell`],
    /// [`XyzWriteOptions::geometry_metadata`], and
    /// [`XyzWriteOptions::properties`], counted over every frame written.
    /// The writer writes the same text whatever it reports.
    pub fn report(&self) -> &WriteReport {
        &self.report
    }

    /// The topology every written frame must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The underlying stream.
    pub fn writer(&self) -> &W {
        &self.writer
    }

    /// Takes the stream as holding `frame_count` frames already, which the
    /// next frame follows, as when appending to a file.
    pub(crate) const fn resumed(mut self, frame_count: usize) -> Self {
        self.frame_count = frame_count;
        self
    }

    pub(crate) fn validate_finish(&self) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Xyz,
            &self.source_label,
        )?;
        require_nonempty_writer(self.frame_count, TrajectoryFormat::Xyz, &self.source_label)
    }

    /// Flushes and returns the completed nonempty XYZ stream.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidFrame`] when no frame was
    /// written or an earlier write left a partial frame, and
    /// [`TrajectoryIoErrorKind::Io`] when flushing fails.
    pub fn finish(mut self) -> Result<W, TrajectoryIoError> {
        self.validate_finish()?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::Finish,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })?;
        Ok(self.writer)
    }
}

impl<W: Write> TrajectoryWriter for XyzWriter<W> {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn write_frame(&mut self, frame: FrameView<'_>) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            &self.source_label,
        )?;
        if !self.topology().shares_layout(frame.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let omitted = check_frame_state(
            TrajectoryFormat::Xyz,
            &self.source_label,
            self.frame_count,
            frame,
            self.options.state_policies(FIELDS.carries()),
        )?;
        let factor = self.options.output_factor();
        self.check_positions(frame.positions(), factor)?;
        self.write_positions(frame.positions(), factor)?;
        self.record(omitted);
        Ok(())
    }

    /// Flushes the stream, which then ends after the last whole frame.
    fn flush(&mut self) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            &self.source_label,
        )?;
        self.writer.flush().map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })
    }
}

impl<W: Write> XyzWriter<W> {
    /// Adds `omissions` to [`Self::report`].
    pub(super) fn record(&mut self, omissions: impl IntoIterator<Item = WriteOmission>) {
        record(&mut self.report, omissions);
    }

    /// Fails unless every coordinate stays finite after conversion.
    pub(super) fn check_positions(
        &self,
        positions: &Positions,
        factor: f64,
    ) -> Result<(), TrajectoryIoError> {
        if finite_in_unit(positions, factor) {
            return Ok(());
        }
        Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::WriteFrame,
            Some(FORMAT),
            &self.source_label,
            Reason::NotFiniteInOutputUnit,
        ))
    }

    /// Writes one frame of checked positions.
    pub(super) fn write_positions(
        &mut self,
        positions: &Positions,
        factor: f64,
    ) -> Result<(), TrajectoryIoError> {
        self.partial = true;
        write_frame(
            &mut self.writer,
            &self.topology,
            positions,
            factor,
            &self.options,
        )
        .map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })?;
        self.partial = false;
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                Reason::AddressRange(Part::FrameCount),
            )
        })?;
        Ok(())
    }
}
