use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use kekule::topology::Topology;

use crate::{
    FrameBuffer, MemoryTrajectoryWriter, SeekableTrajectoryReader, TrajectoryFormat,
    TrajectoryIoError, TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use kekule::structure::Trajectory;

use super::{
    dcd, detect, io_context, trr, xtc, xyz, FileTrajectoryMetadata, TrajectoryOpenOptions,
    TrajectoryOpenReport,
};

enum SequentialReaderInner {
    Xyz(xyz::XyzReader<BufReader<File>>),
    Dcd(dcd::DcdReader<BufReader<File>>),
    Trr(Box<trr::TrrReader<BufReader<File>>>),
    Xtc(xtc::XtcReader<BufReader<File>>),
}

/// Format-agnostic path-backed sequential reader retaining one file handle.
pub struct SequentialFileTrajectoryReader {
    inner: SequentialReaderInner,
    metadata: FileTrajectoryMetadata,
    open_report: TrajectoryOpenReport,
}

impl fmt::Debug for SequentialFileTrajectoryReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SequentialFileTrajectoryReader")
            .field("metadata", &self.metadata)
            .field("open_report", &self.open_report)
            .finish_non_exhaustive()
    }
}

impl SequentialFileTrajectoryReader {
    /// Format detection and non-fatal facts recorded when this reader was opened.
    pub fn open_report(&self) -> &TrajectoryOpenReport {
        &self.open_report
    }

    /// Returns metadata verified through the most recent successful read.
    ///
    /// In particular, mixed-width TRR input is reported as mixed as soon as
    /// the second scalar width has been observed.
    pub fn metadata(&self) -> &FileTrajectoryMetadata {
        &self.metadata
    }
}

impl TrajectoryReader for SequentialFileTrajectoryReader {
    fn topology(&self) -> &Topology {
        match &self.inner {
            SequentialReaderInner::Xyz(reader) => reader.topology(),
            SequentialReaderInner::Dcd(reader) => reader.topology(),
            SequentialReaderInner::Trr(reader) => reader.topology(),
            SequentialReaderInner::Xtc(reader) => reader.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            SequentialReaderInner::Xyz(reader) => reader.shared_topology(),
            SequentialReaderInner::Dcd(reader) => reader.shared_topology(),
            SequentialReaderInner::Trr(reader) => reader.shared_topology(),
            SequentialReaderInner::Xtc(reader) => reader.shared_topology(),
        }
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        match &mut self.inner {
            SequentialReaderInner::Xyz(reader) => reader.read_next(destination),
            SequentialReaderInner::Dcd(reader) => reader.read_next(destination),
            SequentialReaderInner::Trr(reader) => {
                let read = reader.read_next(destination)?;
                self.metadata
                    .update_trr_precision(reader.first_header(), reader.precision_mixed());
                Ok(read)
            }
            SequentialReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

enum IndexedReaderInner {
    Xyz(xyz::IndexedXyzReader<BufReader<File>>),
    Dcd(dcd::IndexedDcdReader<BufReader<File>>),
    Trr(Box<trr::IndexedTrrReader<BufReader<File>>>),
    Xtc(xtc::IndexedXtcReader<BufReader<File>>),
}

/// Format-agnostic path-backed indexed reader retaining one file handle.
pub struct IndexedFileTrajectoryReader {
    inner: IndexedReaderInner,
    metadata: FileTrajectoryMetadata,
    open_report: TrajectoryOpenReport,
}

impl fmt::Debug for IndexedFileTrajectoryReader {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("IndexedFileTrajectoryReader")
            .field("metadata", &self.metadata)
            .field("open_report", &self.open_report)
            .finish_non_exhaustive()
    }
}

impl IndexedFileTrajectoryReader {
    /// Format detection and non-fatal facts recorded when this reader was opened.
    pub fn open_report(&self) -> &TrajectoryOpenReport {
        &self.open_report
    }

    /// Metadata verified for the whole file while building the index.
    pub fn metadata(&self) -> &FileTrajectoryMetadata {
        &self.metadata
    }
}

impl TrajectoryReader for IndexedFileTrajectoryReader {
    fn topology(&self) -> &Topology {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.topology(),
            IndexedReaderInner::Dcd(reader) => reader.topology(),
            IndexedReaderInner::Trr(reader) => reader.topology(),
            IndexedReaderInner::Xtc(reader) => reader.topology(),
        }
    }

    fn shared_topology(&self) -> Arc<Topology> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.shared_topology(),
            IndexedReaderInner::Dcd(reader) => reader.shared_topology(),
            IndexedReaderInner::Trr(reader) => reader.shared_topology(),
            IndexedReaderInner::Xtc(reader) => reader.shared_topology(),
        }
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_next(destination),
            IndexedReaderInner::Dcd(reader) => reader.read_next(destination),
            IndexedReaderInner::Trr(reader) => reader.read_next(destination),
            IndexedReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

impl SeekableTrajectoryReader for IndexedFileTrajectoryReader {
    fn frame_count(&self) -> Option<usize> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.frame_count(),
            IndexedReaderInner::Dcd(reader) => reader.frame_count(),
            IndexedReaderInner::Trr(reader) => reader.frame_count(),
            IndexedReaderInner::Xtc(reader) => reader.frame_count(),
        }
    }

    fn read_frame(
        &mut self,
        index: usize,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Dcd(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Trr(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Xtc(reader) => reader.read_frame(index, destination),
        }
    }
}

/// Reads an entire trajectory file into memory in file frame order.
///
/// Accepts an owned topology or a shared `Arc<Topology>` and detects the format
/// automatically. File coordinates must follow the topology's dense atom order;
/// counts and available format metadata are checked, but matching counts alone
/// cannot establish atom identity. All decoded frame fields are retained, and
/// a supplied shared topology retains its exact allocation.
///
/// This consumes the sequential reader through clean EOF. Any open or decoding
/// error is returned without publishing a partial trajectory. Memory usage grows
/// with the number of frames and atoms; use [`open_trajectory`] to process large
/// files one frame at a time or inspect file metadata and opening diagnostics.
///
/// `options` sets the format policies and limits; frame limits cause an error
/// rather than truncating the result. These codec limits do not impose a total
/// resident-memory budget for the loaded trajectory. No random-access index is
/// built.
///
/// # Errors
///
/// Returns the errors of [`open_trajectory`] and the first decoding,
/// validation or [`TrajectoryIoErrorKind::Io`] error of any frame.
///
/// [`TrajectoryIoErrorKind::Io`]: crate::TrajectoryIoErrorKind::Io
///
/// ```no_run
/// use kekule::mmcif;
/// use kekule_traj::io::read_trajectory;
///
/// let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?;
/// let topology = document.interpret(Default::default())?.into_topology();
/// let trajectory = read_trajectory("trajectory.xtc", topology, Default::default())?;
/// println!("{} frames, {} atoms", trajectory.len(), trajectory.topology().atom_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn read_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<Trajectory, TrajectoryIoError> {
    let mut reader = open_trajectory(path, topology, options)?;
    let mut buffer = reader.frame_buffer();
    let mut writer = MemoryTrajectoryWriter::new(reader.shared_topology());
    while reader.read_next(&mut buffer)? {
        writer.write_frame(buffer.frame_view())?;
    }
    Ok(writer.into_trajectory())
}

/// Opens one fast sequential path-backed trajectory reader without loading all frames.
///
/// Accepts an owned topology or a shared `Arc<Topology>`. File coordinates must
/// follow its dense atom order; counts and available format metadata are checked
/// automatically. Matching counts alone cannot establish atom identity.
/// `options` sets the format policies and limits. For a fully loaded in-memory
/// trajectory, use [`read_trajectory`] instead.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::Io`] when the file cannot be opened or
/// read, the detection errors of [`super::detect_trajectory_format`], and the
/// header errors of the selected codec's reader, such as
/// [`TrajectoryIoErrorKind::InconsistentAtomCount`] when the file's atom count
/// differs from the topology's.
///
/// [`TrajectoryIoErrorKind::Io`]: crate::TrajectoryIoErrorKind::Io
/// [`TrajectoryIoErrorKind::InconsistentAtomCount`]:
///     crate::TrajectoryIoErrorKind::InconsistentAtomCount
///
/// ```no_run
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::{io::open_trajectory, TrajectoryReader};
///
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Topology::from_molecule(molecule)?;
/// let mut reader = open_trajectory("ethane.xyz", topology, Default::default())?;
/// let mut frame = reader.frame_buffer();
/// while reader.read_next(&mut frame)? {
///     let model = frame.as_model_view();
///     assert_eq!(model.atom_count(), 2);
/// }
/// println!("{:?}", reader.open_report().selected_format());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn open_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<SequentialFileTrajectoryReader, TrajectoryIoError> {
    let topology = topology.into();
    let path = path.as_ref();
    let label = path.display().to_string();
    let file = File::open(path)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
    let mut reader = BufReader::new(file);
    let detection = detect::select_format(
        &mut reader,
        path,
        options.format_hint,
        &options.limits,
        &label,
    )?;
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let mut reader = xyz::XyzReader::new(
                reader,
                topology,
                options
                    .xyz
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            reader.validate_first_frame()?;
            let atom_count = reader.topology().atom_count();
            (
                SequentialReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, None),
                vec!["XYZ coordinates use the configured explicit/default length unit".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::DcdReader::new(
                reader,
                topology,
                options
                    .dcd
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), None);
            (
                SequentialReaderInner::Dcd(reader),
                metadata,
                vec![
                    "DCD steps are preserved from ISTART/NSAVC; time follows the explicit policy"
                        .into(),
                ],
            )
        }
        TrajectoryFormat::Trr => {
            let reader = trr::TrrReader::new(
                reader,
                topology,
                options
                    .trr
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let metadata = FileTrajectoryMetadata::trr(reader.first_header(), None, false);
            (
                SequentialReaderInner::Trr(Box::new(reader)),
                metadata,
                vec![
                    "TRR uses XDR big-endian scalars and explicit lambda preservation policy"
                        .into(),
                ],
            )
        }
        TrajectoryFormat::Xtc => {
            let cell_policy = options.xtc.cell_policy();
            let xtc_options = options.xtc;
            let reader = xtc::XtcReader::new(
                reader,
                topology,
                xtc_options
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let metadata = FileTrajectoryMetadata::xtc(reader.first_info(), cell_policy, None);
            (
                SequentialReaderInner::Xtc(reader),
                metadata,
                vec!["XTC decoding uses a bounded checked reader".into()],
            )
        }
    };
    let report = TrajectoryOpenReport {
        selected_format: detection.format,
        detection_evidence: detection.evidence,
        notes,
    };
    Ok(SequentialFileTrajectoryReader {
        inner,
        metadata,
        open_report: report,
    })
}

/// Opens one fully verified indexed path-backed trajectory reader.
///
/// Uses the topology order contract of [`open_trajectory`]: file coordinate
/// index `i` is topology dense atom index `i`, and readers check counts and
/// available metadata, including XYZ element order. Opening scans and validates
/// the entire file and stores bounded frame offsets for random access.
/// `options` sets the format policies and limits.
///
/// # Errors
///
/// Returns the errors of [`open_trajectory`], the decoding error of the first
/// invalid frame, and [`TrajectoryIoErrorKind::LimitExceeded`] when the index
/// outgrows the limits.
///
/// [`TrajectoryIoErrorKind::LimitExceeded`]: crate::TrajectoryIoErrorKind::LimitExceeded
pub fn open_indexed_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    options: TrajectoryOpenOptions,
) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
    let topology = topology.into();
    let path = path.as_ref();
    let label = path.display().to_string();
    let file = File::open(path)
        .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
    let mut reader = BufReader::new(file);
    let detection = detect::select_format(
        &mut reader,
        path,
        options.format_hint,
        &options.limits,
        &label,
    )?;
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let reader = xyz::IndexedXyzReader::new(
                reader,
                topology,
                options
                    .xyz
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let count = reader.frame_count().unwrap_or(0);
            let atom_count = reader.topology().atom_count();
            (
                IndexedReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, Some(count)),
                vec!["XYZ index verified every complete frame".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::IndexedDcdReader::new(
                reader,
                topology,
                options
                    .dcd
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), Some(count));
            (
                IndexedReaderInner::Dcd(reader),
                metadata,
                vec!["DCD index verified record markers and the declared frame count".into()],
            )
        }
        TrajectoryFormat::Trr => {
            let reader = trr::IndexedTrrReader::new(
                reader,
                topology,
                options
                    .trr
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::trr(
                reader.first_header(),
                Some(count),
                reader.precision_mixed(),
            );
            (
                IndexedReaderInner::Trr(Box::new(reader)),
                metadata,
                vec!["TRR index verified every XDR frame and payload block".into()],
            )
        }
        TrajectoryFormat::Xtc => {
            let cell_policy = options.xtc.cell_policy();
            let xtc_options = options.xtc;
            let reader = xtc::IndexedXtcReader::new(
                reader,
                topology,
                xtc_options
                    .with_limits(options.limits)
                    .with_source_label(label),
            )?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata =
                FileTrajectoryMetadata::xtc(reader.first_info(), cell_policy, Some(count));
            (
                IndexedReaderInner::Xtc(reader),
                metadata,
                vec!["XTC index fully decoded and verified every compressed frame".into()],
            )
        }
    };
    let report = TrajectoryOpenReport {
        selected_format: detection.format,
        detection_evidence: detection.evidence,
        notes,
    };
    Ok(IndexedFileTrajectoryReader {
        inner,
        metadata,
        open_report: report,
    })
}
