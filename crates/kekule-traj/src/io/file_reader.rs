use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

use kekule::topology::Topology;

use crate::{
    FrameBuffer, SeekableTrajectoryReader, TrajectoryError, TrajectoryFormat,
    TrajectoryIoOperation, TrajectoryReader,
};

use super::{
    dcd, detect, io_context, trr, xtc, xyz, FileTrajectoryMetadata, TrajectoryOpenOptions,
    TrajectoryOpenReport, TrajectoryTopologyBinding,
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
}

impl SequentialFileTrajectoryReader {
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

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
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
}

impl IndexedFileTrajectoryReader {
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

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_next(destination),
            IndexedReaderInner::Dcd(reader) => reader.read_next(destination),
            IndexedReaderInner::Trr(reader) => reader.read_next(destination),
            IndexedReaderInner::Xtc(reader) => reader.read_next(destination),
        }
    }
}

impl SeekableTrajectoryReader for IndexedFileTrajectoryReader {
    fn frame_count(&self) -> Option<u64> {
        match &self.inner {
            IndexedReaderInner::Xyz(reader) => reader.frame_count(),
            IndexedReaderInner::Dcd(reader) => reader.frame_count(),
            IndexedReaderInner::Trr(reader) => reader.frame_count(),
            IndexedReaderInner::Xtc(reader) => reader.frame_count(),
        }
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        match &mut self.inner {
            IndexedReaderInner::Xyz(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Dcd(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Trr(reader) => reader.read_frame(index, destination),
            IndexedReaderInner::Xtc(reader) => reader.read_frame(index, destination),
        }
    }
}

/// Opens one fast sequential path-backed trajectory reader.
pub fn open_trajectory(
    path: impl AsRef<Path>,
    binding: TrajectoryTopologyBinding,
    options: TrajectoryOpenOptions,
) -> Result<(SequentialFileTrajectoryReader, TrajectoryOpenReport), TrajectoryError> {
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
    let order_kind = binding.atom_order().kind();
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let mut reader =
                xyz::XyzReader::new(reader, binding, options.xyz, options.limits, label)?;
            reader.validate_first_frame()?;
            let atom_count = reader.topology().atom_count();
            (
                SequentialReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, None),
                vec!["XYZ coordinates use the configured explicit/default length unit".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::DcdReader::new(reader, binding, options.dcd, options.limits, label)?;
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
            let reader = trr::TrrReader::new(reader, binding, options.trr, options.limits, label)?;
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
            let xtc_options = options.xtc;
            let reader = xtc::XtcReader::new(reader, binding, xtc_options, options.limits, label)?;
            let metadata =
                FileTrajectoryMetadata::xtc(reader.first_info(), xtc_options.cell_policy(), None);
            (
                SequentialReaderInner::Xtc(reader),
                metadata,
                vec![
                    "XTC decoding uses a bounded checked reader; molly is confined to writing"
                        .into(),
                ],
            )
        }
    };
    let report = TrajectoryOpenReport {
        selected_format: detection.format,
        detection_evidence: detection.evidence,
        atom_order_evidence: order_kind,
        notes,
    };
    Ok((SequentialFileTrajectoryReader { inner, metadata }, report))
}

/// Opens one fully verified indexed path-backed trajectory reader.
pub fn open_indexed_trajectory(
    path: impl AsRef<Path>,
    binding: TrajectoryTopologyBinding,
    options: TrajectoryOpenOptions,
) -> Result<(IndexedFileTrajectoryReader, TrajectoryOpenReport), TrajectoryError> {
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
    let order_kind = binding.atom_order().kind();
    let (inner, metadata, notes) = match detection.format {
        TrajectoryFormat::Xyz => {
            let reader = xyz::XyzReader::new(reader, binding, options.xyz, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let atom_count = reader.topology().atom_count();
            (
                IndexedReaderInner::Xyz(reader),
                FileTrajectoryMetadata::xyz(atom_count, Some(count)),
                vec!["XYZ index verified every complete frame".into()],
            )
        }
        TrajectoryFormat::Dcd => {
            let reader = dcd::DcdReader::new(reader, binding, options.dcd, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::dcd(reader.header(), Some(count));
            (
                IndexedReaderInner::Dcd(reader),
                metadata,
                vec!["DCD index verified record markers and the declared frame count".into()],
            )
        }
        TrajectoryFormat::Trr => {
            let reader = trr::TrrReader::new(reader, binding, options.trr, options.limits, label)?;
            let reader = reader.to_indexed()?;
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
            let xtc_options = options.xtc;
            let reader = xtc::XtcReader::new(reader, binding, xtc_options, options.limits, label)?;
            let reader = reader.to_indexed()?;
            let count = reader.frame_count().unwrap_or(0);
            let metadata = FileTrajectoryMetadata::xtc(
                reader.first_info(),
                xtc_options.cell_policy(),
                Some(count),
            );
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
        atom_order_evidence: order_kind,
        notes,
    };
    Ok((IndexedFileTrajectoryReader { inner, metadata }, report))
}
