use super::*;
use crate::traj::trajectory::reason::{DcdField, DcdRecord, Part, Reason};

/// Sequential DCD reader over one seekable stream.
pub struct DcdReader<R> {
    pub(super) reader: R,
    pub(super) topology: Arc<Topology>,
    pub(super) options: DcdReadOptions,
    pub(super) source: ReadSource,
    pub(super) header: DcdHeader,
    pub(super) free_indices: Vec<usize>,
    pub(super) fixed_indices: Vec<usize>,
    pub(super) fixed_reference: Vec<Point3>,
    pub(super) positions: Vec<Point3>,
    pub(super) record: Vec<u8>,
    pub(super) frame_cursor: usize,
}

impl<R> fmt::Debug for DcdReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DcdReader")
            .field("atom_count", &self.topology.atom_count())
            .field("options", &self.options)
            .field("source_label", &self.source.label)
            .field("frame_cursor", &self.frame_cursor)
            .finish_non_exhaustive()
    }
}

#[derive(Debug, Clone, Copy)]
pub(crate) struct DcdDecodedFrame {
    pub(super) cell: Option<PeriodicCell>,
    pub(super) step: u64,
    pub(super) time: Option<f64>,
}

impl<R: Read + Seek> DcdReader<R> {
    /// Opens a stream whose coordinates follow the supplied topology's dense order.
    ///
    /// Accepts an owned topology or a shared `Arc<Topology>`. Counts and available
    /// format metadata are checked automatically; matching counts alone cannot
    /// establish atom identity. Limits and diagnostic labels are set in `options`.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidHeader`],
    /// [`TrajectoryIoErrorKind::UnsupportedVariant`],
    /// [`TrajectoryIoErrorKind::TruncatedRecord`] or
    /// [`TrajectoryIoErrorKind::RecordMarkerMismatch`] for a malformed or
    /// unsupported header, [`TrajectoryIoErrorKind::InconsistentAtomCount`]
    /// when its atom count differs from the topology's,
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when it exceeds the limits, and
    /// [`TrajectoryIoErrorKind::Io`] when reading fails.
    pub fn new(
        reader: R,
        topology: impl Into<Arc<Topology>>,
        options: impl Into<TrajectoryOpenOptions>,
    ) -> Result<Self, TrajectoryIoError> {
        let options = options.into();
        let source = options.read_source(FORMAT);
        Self::open(reader, topology.into(), options.dcd, source)
    }

    /// Opens the stream like [`Self::new`] with the format options and
    /// source taken apart.
    pub(crate) fn open(
        mut reader: R,
        topology: Arc<Topology>,
        options: DcdReadOptions,
        source: ReadSource,
    ) -> Result<Self, TrajectoryIoError> {
        let source_label = &source.label;
        let limits = &source.limits;
        let atom_count = topology.atom_count();
        validate_atom_count(FORMAT, atom_count, i32::MAX as usize, limits, source_label)?;
        let start = reader
            .stream_position()
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, source_label, error))?;
        let mut marker = [0_u8; 4];
        read_exact(
            FORMAT,
            &mut reader,
            &mut marker,
            TrajectoryIoOperation::ReadHeader,
            source_label,
            None,
            Reason::Missing(Part::RecordMarker),
        )?;
        let endian = if u32::from_le_bytes(marker) == HEADER_BYTES as u32 {
            DcdEndian::Little
        } else if u32::from_be_bytes(marker) == HEADER_BYTES as u32 {
            DcdEndian::Big
        } else {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                Reason::FirstMarkerNot84,
            ));
        };
        reader
            .seek(SeekFrom::Start(start))
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, source_label, error))?;

        let mut record = Vec::new();
        read_record(
            &mut reader,
            endian,
            &mut record,
            Some(0),
            limits,
            source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
            None,
        )?;
        if record.len() != HEADER_BYTES || &record[..4] != b"CORD" {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidHeader,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                Reason::NotCordHeader,
            ));
        }
        let controls = parse_controls(&record[4..], endian, source_label)?;
        if controls[11] != 0 || controls[12] != 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                Reason::FourDimensional,
            ));
        }
        let declared_frames = nonnegative_u64(controls[0], DcdField::Nset, source_label)? as usize;
        let start_step = nonnegative_u64(controls[1], DcdField::Istart, source_label)?;
        let step_interval = positive_u64(controls[2], DcdField::Nsavc, source_label)?;
        let fixed_count = nonnegative_usize(controls[8], DcdField::Namnf, source_label)?;
        let delta = endian.f32(record[40..44].try_into().expect("header slice"));
        if !delta.is_finite() || delta < 0.0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                Reason::InvalidDelta(delta),
            ));
        }
        let has_cell = controls[10] != 0;
        let charmm_version = controls[19];
        if charmm_version <= 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::UnsupportedVariant,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                Reason::UnsupportedDcdLayout,
            ));
        }

        read_record(
            &mut reader,
            endian,
            &mut record,
            Some(0),
            limits,
            source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
            None,
        )?;
        validate_title_record(&record, endian, source_label)?;
        read_record(
            &mut reader,
            endian,
            &mut record,
            Some(0),
            limits,
            source_label,
            TrajectoryIoOperation::ReadHeader,
            None,
            false,
            None,
        )?;
        if record.len() != 4 {
            return Err(header_error(
                FORMAT,
                source_label,
                Reason::RecordSize {
                    record: DcdRecord::AtomCount,
                    bytes: record.len(),
                    expected: 4,
                },
            ));
        }
        let file_atoms = nonnegative_usize(
            endian.i32(record[..4].try_into().expect("atom-count slice")),
            DcdField::Natom,
            source_label,
        )?;
        validate_atom_count(FORMAT, file_atoms, i32::MAX as usize, limits, source_label)?;
        if file_atoms != atom_count {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentAtomCount,
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
            )
            .with_source_label(source_label)
            .with_counts(atom_count as u64, file_atoms as u64));
        }
        if fixed_count > atom_count {
            return Err(header_error(
                FORMAT,
                source_label,
                Reason::FixedAtomsExceed {
                    fixed: fixed_count,
                    atoms: atom_count,
                },
            ));
        }
        let free_count = atom_count - fixed_count;
        let mut free_indices = Vec::new();
        let mut fixed_indices = Vec::new();
        let reserve_error = |reason: Reason| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::Open,
                source_label,
                None,
                reason,
            )
        };
        if fixed_count > 0 {
            read_record(
                &mut reader,
                endian,
                &mut record,
                Some(0),
                limits,
                source_label,
                TrajectoryIoOperation::ReadHeader,
                None,
                false,
                None,
            )?;
            let expected = checked_coordinate_bytes(free_count, source_label)?;
            if record.len() != expected {
                return Err(header_error(
                    FORMAT,
                    source_label,
                    Reason::RecordSize {
                        record: DcdRecord::FreeAtomIndex,
                        bytes: record.len(),
                        expected,
                    },
                ));
            }
            reserve_scratch(
                &mut free_indices,
                free_count,
                scratch_bytes(&record),
                limits.max_scratch_bytes(),
                reserve_error,
            )?;
            let mut seen = Vec::new();
            reserve_scratch(
                &mut seen,
                atom_count,
                scratch_bytes(&record)
                    .and_then(|bytes| bytes.checked_add(scratch_bytes(&free_indices)?)),
                limits.max_scratch_bytes(),
                reserve_error,
            )?;
            seen.resize(atom_count, false);
            for chunk in record.as_chunks::<4>().0 {
                let one_based = endian.i32(*chunk);
                if one_based <= 0 {
                    return Err(header_error(
                        FORMAT,
                        source_label,
                        Reason::FreeAtomIndexNotPositive(one_based),
                    ));
                }
                let index = usize::try_from(one_based - 1).map_err(|_| {
                    header_error(
                        FORMAT,
                        source_label,
                        Reason::AddressRange(Part::FreeAtomIndex),
                    )
                })?;
                if index >= atom_count || seen[index] {
                    return Err(header_error(
                        FORMAT,
                        source_label,
                        Reason::FreeAtomIndexInvalid {
                            index: one_based,
                            atoms: atom_count,
                        },
                    ));
                }
                seen[index] = true;
                free_indices.push(index);
            }
            reserve_scratch(
                &mut fixed_indices,
                fixed_count,
                scratch_bytes(&record)
                    .and_then(|bytes| bytes.checked_add(scratch_bytes(&free_indices)?))
                    .and_then(|bytes| bytes.checked_add(scratch_bytes(&seen)?)),
                limits.max_scratch_bytes(),
                reserve_error,
            )?;
            fixed_indices.extend(
                seen.iter()
                    .enumerate()
                    .filter_map(|(index, free)| (!free).then_some(index)),
            );
        }
        let index_bytes = scratch_bytes(&free_indices)
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&fixed_indices)?));
        let mut fixed_reference = Vec::new();
        reserve_scratch(
            &mut fixed_reference,
            fixed_count,
            index_bytes.and_then(|bytes| bytes.checked_add(scratch_bytes(&record)?)),
            limits.max_scratch_bytes(),
            reserve_error,
        )?;
        let mut positions = Vec::new();
        reserve_scratch(
            &mut positions,
            atom_count,
            index_bytes
                .and_then(|bytes| bytes.checked_add(scratch_bytes(&fixed_reference)?))
                .and_then(|bytes| bytes.checked_add(scratch_bytes(&record)?)),
            limits.max_scratch_bytes(),
            reserve_error,
        )?;
        let data_start = reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadHeader,
                Some(FORMAT),
                source_label,
                error,
            )
        })?;
        let header = DcdHeader {
            endian,
            declared_frames,
            start_step,
            step_interval,
            delta,
            fixed_count,
            atom_count,
            has_cell,
            charmm_version,
            header_start: start,
            data_start,
        };
        Ok(Self {
            reader,
            topology,
            options,
            source,
            header,
            free_indices,
            fixed_indices,
            fixed_reference,
            positions,
            record,
            frame_cursor: 0,
        })
    }

    /// The topology every frame belongs to.
    pub fn topology(&self) -> &Topology {
        self.topology.as_ref()
    }

    pub(crate) fn header(&self) -> &DcdHeader {
        &self.header
    }

    pub(super) fn parse_next(
        &mut self,
        capture_fixed_reference: bool,
    ) -> Result<Option<DcdDecodedFrame>, TrajectoryIoError> {
        // The position array holds one point per atom, whether it is kept or
        // allocated by this read.
        let other_scratch_bytes = filled_bytes(&self.positions, self.header.atom_count)
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.free_indices)?))
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.fixed_indices)?))
            .and_then(|bytes| bytes.checked_add(scratch_bytes(&self.fixed_reference)?));
        let frame_start = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
                &self.source.label,
                error,
            )
        })?;
        // Frames are counted in the stream, whatever the header's `NSET`
        // says: the stream ends cleanly before a frame's first record.
        if self.frame_cursor >= self.source.limits.max_frames() {
            if probe_seekable_eof(
                &mut self.reader,
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Dcd,
                &self.source.label,
            )? {
                return Ok(None);
            }
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                Some(self.frame_cursor),
                None,
            )
            .exceeding(
                Part::Trajectory,
                LimitExceeded::first_beyond(Resource::Frames, self.source.limits.max_frames()),
            ));
        }
        let mut cell = None;
        if self.header.has_cell {
            if !read_record(
                &mut self.reader,
                self.header.endian,
                &mut self.record,
                other_scratch_bytes,
                &self.source.limits,
                &self.source.label,
                TrajectoryIoOperation::ReadFrame,
                Some(self.frame_cursor),
                true,
                Some((CELL_BYTES, Reason::UnitCellValues)),
            )? {
                return Ok(None);
            }
            cell = Some(decode_cell(
                &self.record,
                self.header.endian,
                &self.source.label,
                self.frame_cursor,
            )?);
        }

        let coordinate_indices = if self.frame_cursor == 0 || self.header.fixed_count == 0 {
            // The first record writes every atom once. The indexed reader
            // hands the buffer over, so the next read allocates its own on
            // the thread that reads.
            self.positions.clear();
            self.positions.reserve_exact(self.header.atom_count);
            None
        } else {
            // A prior random read may have overwritten the shared position
            // scratch, or handed it over. Every partial frame restores its
            // immutable fixed atoms.
            self.positions
                .resize(self.header.atom_count, Point3::new(0.0, 0.0, 0.0));
            for (&index, &point) in self.fixed_indices.iter().zip(&self.fixed_reference) {
                self.positions[index] = point;
            }
            Some(self.free_indices.as_slice())
        };
        let coordinate_count = coordinate_indices.map_or(self.header.atom_count, <[usize]>::len);
        let coordinate_bytes = checked_coordinate_bytes(coordinate_count, &self.source.label)?;
        if coordinate_bytes as u64 > self.source.limits.max_frame_bytes() {
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                Some(self.frame_cursor),
                None,
            )
            .exceeding(
                Part::CoordinateRecord,
                exceeded_bytes(
                    self.source.limits.max_frame_bytes(),
                    coordinate_bytes as u64,
                ),
            ));
        }
        for axis in 0..3 {
            let clean_frame_eof = !self.header.has_cell && axis == 0;
            if !read_record(
                &mut self.reader,
                self.header.endian,
                &mut self.record,
                other_scratch_bytes,
                &self.source.limits,
                &self.source.label,
                TrajectoryIoOperation::ReadFrame,
                Some(self.frame_cursor),
                clean_frame_eof,
                Some((coordinate_bytes, Reason::CoordinateRecordSize)),
            )? {
                return Ok(None);
            }
            for (value_index, chunk) in self.record.as_chunks::<4>().0.iter().enumerate() {
                let value = f64::from(self.header.endian.f32(*chunk));
                if !value.is_finite() {
                    return Err(frame_error(
                        FORMAT,
                        TrajectoryIoErrorKind::InvalidFrame,
                        TrajectoryIoOperation::ReadFrame,
                        &self.source.label,
                        self.frame_cursor,
                        Reason::NotFinite(Part::Coordinate),
                    ));
                }
                let atom_index =
                    coordinate_indices.map_or(value_index, |indices| indices[value_index]);
                match axis {
                    0 if coordinate_indices.is_none() => {
                        self.positions.push(Point3::new(value, 0.0, 0.0));
                    }
                    0 => self.positions[atom_index].x = value,
                    1 => self.positions[atom_index].y = value,
                    _ => self.positions[atom_index].z = value,
                }
            }
        }
        let frame_end = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
                &self.source.label,
                error,
            )
        })?;
        let frame_bytes = frame_end.checked_sub(frame_start).ok_or_else(|| {
            frame_error(
                FORMAT,
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                self.frame_cursor,
                Reason::FrameEndsBeforeStart,
            )
        })?;
        if frame_bytes > self.source.limits.max_frame_bytes() {
            return Err(limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                Some(self.frame_cursor),
                None,
            )
            .exceeding(
                Part::Frame,
                exceeded_bytes(self.source.limits.max_frame_bytes(), frame_bytes),
            ));
        }
        if capture_fixed_reference && self.frame_cursor == 0 && self.header.fixed_count > 0 {
            self.fixed_reference.clear();
            self.fixed_reference.extend(
                self.fixed_indices
                    .iter()
                    .map(|&index| self.positions[index]),
            );
        }
        let (step, time) = self.step_and_time(self.frame_cursor)?;
        self.frame_cursor = self.frame_cursor.checked_add(1).ok_or_else(|| {
            limit_error(
                FORMAT,
                TrajectoryIoOperation::ReadFrame,
                &self.source.label,
                Some(self.frame_cursor),
                Reason::AddressRange(Part::FrameNumber),
            )
        })?;
        Ok(Some(DcdDecodedFrame { cell, step, time }))
    }

    /// The step of frame `frame` on the `ISTART`/`NSAVC` sequence and, under
    /// [`DcdTimePolicy::HeaderDelta`], its time in the policy's unit.
    pub(super) fn step_and_time(
        &self,
        frame: usize,
    ) -> Result<(u64, Option<f64>), TrajectoryIoError> {
        let step = self
            .header
            .start_step
            .checked_add(
                (frame as u64)
                    .checked_mul(self.header.step_interval)
                    .ok_or_else(|| {
                        frame_error(
                            FORMAT,
                            TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                            TrajectoryIoOperation::ReadFrame,
                            &self.source.label,
                            frame,
                            Reason::StepRange64,
                        )
                    })?,
            )
            .ok_or_else(|| {
                frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::NegativeOrUnrepresentableStep,
                    TrajectoryIoOperation::ReadFrame,
                    &self.source.label,
                    frame,
                    Reason::StepRange64,
                )
            })?;
        let time = if let DcdTimePolicy::HeaderDelta { unit } = self.options.time_policy {
            let time = (step as f64) * f64::from(self.header.delta);
            if !Quantity::new(time, unit)
                .value_in(CANONICAL_TIME_UNIT)
                .is_ok_and(f64::is_finite)
            {
                return Err(frame_error(
                    FORMAT,
                    TrajectoryIoErrorKind::InvalidFrame,
                    TrajectoryIoOperation::ReadFrame,
                    &self.source.label,
                    frame,
                    Reason::NotFinite(Part::FrameTime),
                ));
            }
            Some(time)
        } else {
            None
        };
        Ok((step, time))
    }

    /// The step and canonical time of frame `frame`, as an index records them.
    pub(super) fn stamp(&self, frame: usize) -> Result<FrameStamp, TrajectoryIoError> {
        let (step, time) = self.step_and_time(frame)?;
        let unit = header_delta_unit(self.options.time_policy);
        Ok(FrameStamp {
            time: time
                .and_then(|time| Quantity::new(time, unit).value_in(CANONICAL_TIME_UNIT).ok()),
            step: Some(step),
        })
    }

    pub(super) fn publish(
        &self,
        positions: &[Point3],
        decoded: DcdDecodedFrame,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryIoError> {
        let mut data = FrameData::new(Quantity::new(positions, ANGSTROM)).with_step(decoded.step);
        if let Some(cell) = decoded.cell {
            data = data.with_cell(cell);
        }
        if let (Some(time), DcdTimePolicy::HeaderDelta { unit }) =
            (decoded.time, self.options.time_policy)
        {
            data = data.with_time(Quantity::new(time, unit));
        }
        destination.replace_from_data(data).map_err(Into::into)
    }

    pub(super) fn position(
        &mut self,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader
            .stream_position()
            .map_err(|error| io_context(operation, Some(FORMAT), &self.source.label, error))
    }

    pub(super) fn seek_to(
        &mut self,
        offset: SeekFrom,
        operation: TrajectoryIoOperation,
    ) -> Result<u64, TrajectoryIoError> {
        self.reader
            .seek(offset)
            .map_err(|error| io_context(operation, Some(FORMAT), &self.source.label, error))
    }

    /// Reads the header's `NSET` again, which a writer appending frames may
    /// have updated, and returns it; the stream position is kept.
    pub(super) fn reread_declared_frames(&mut self) -> Result<usize, TrajectoryIoError> {
        let position = self.position(TrajectoryIoOperation::ReadHeader)?;
        let count_offset = self.header.header_start.checked_add(8).ok_or_else(|| {
            header_error(
                FORMAT,
                &self.source.label,
                Reason::Range64(Part::HeaderOffset),
            )
        })?;
        self.seek_to(
            SeekFrom::Start(count_offset),
            TrajectoryIoOperation::ReadHeader,
        )?;
        let mut count = [0_u8; 4];
        let read = read_exact(
            FORMAT,
            &mut self.reader,
            &mut count,
            TrajectoryIoOperation::ReadHeader,
            &self.source.label,
            None,
            Reason::Truncated(Part::Nset),
        );
        self.seek_to(SeekFrom::Start(position), TrajectoryIoOperation::ReadHeader)?;
        read?;
        let declared = nonnegative_u64(
            self.header.endian.i32(count),
            DcdField::Nset,
            &self.source.label,
        )?;
        Ok(declared as usize)
    }
}

impl<R: Read + Seek> TrajectoryReader for DcdReader<R> {
    fn topology(&self) -> &Topology {
        self.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryIoError> {
        if !self.topology().shares_layout(destination.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        let offset = self.reader.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(FORMAT),
                &self.source.label,
                error,
            )
        })?;
        let decoded = self
            .parse_next(true)
            .map_err(|error| frame_offset_context(error, self.frame_cursor, offset))?;
        let Some(decoded) = decoded else {
            return Ok(false);
        };
        self.publish(&self.positions, decoded, destination)?;
        Ok(true)
    }
}

pub(super) fn parse_controls(
    bytes: &[u8],
    endian: DcdEndian,
    source_label: &str,
) -> Result<[i32; 20], TrajectoryIoError> {
    if bytes.len() != 80 {
        return Err(header_error(
            FORMAT,
            source_label,
            Reason::RecordSize {
                record: DcdRecord::Control,
                bytes: bytes.len(),
                expected: 80,
            },
        ));
    }
    let mut controls = [0_i32; 20];
    for (target, chunk) in controls.iter_mut().zip(bytes.as_chunks::<4>().0.iter()) {
        *target = endian.i32(*chunk);
    }
    Ok(controls)
}

pub(super) fn validate_title_record(
    record: &[u8],
    endian: DcdEndian,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    if record.len() < 4 {
        return Err(header_error(
            FORMAT,
            source_label,
            Reason::TitleRecordShort(record.len()),
        ));
    }
    let count = nonnegative_usize(
        endian.i32(record[..4].try_into().expect("title count")),
        DcdField::Ntitle,
        source_label,
    )?;
    let expected = count
        .checked_mul(80)
        .and_then(|bytes| bytes.checked_add(4))
        .ok_or_else(|| {
            header_error(FORMAT, source_label, Reason::AddressRange(Part::TitleCount))
        })?;
    if record.len() != expected {
        return Err(header_error(
            FORMAT,
            source_label,
            Reason::TitleRecordSize {
                bytes: record.len(),
                titles: count,
                expected,
            },
        ));
    }
    Ok(())
}

pub(super) fn decode_cell(
    record: &[u8],
    endian: DcdEndian,
    source_label: &str,
    frame: usize,
) -> Result<PeriodicCell, TrajectoryIoError> {
    let mut values = [0_f64; 6];
    for (target, chunk) in values.iter_mut().zip(record.as_chunks::<8>().0.iter()) {
        *target = endian.f64(*chunk);
    }
    if values.iter().any(|value| !value.is_finite()) {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            Reason::CellNotFinite,
        ));
    }
    let [a, gamma_raw, b, beta_raw, alpha_raw, c] = values;
    if a <= 0.0 || b <= 0.0 || c <= 0.0 {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            Reason::CellLengths { a, b, c },
        ));
    }
    let cosine_variant = [alpha_raw, beta_raw, gamma_raw]
        .into_iter()
        .all(|value| (-1.0..=1.0).contains(&value));
    let angle = |value: f64| {
        if cosine_variant {
            value.acos()
        } else {
            value.to_radians()
        }
    };
    let alpha = angle(alpha_raw);
    let beta = angle(beta_raw);
    let gamma = angle(gamma_raw);
    let sin_gamma = gamma.sin();
    if !alpha.is_finite()
        || !beta.is_finite()
        || !gamma.is_finite()
        || sin_gamma.abs() <= f64::EPSILON
    {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            Reason::DegenerateCellAngles,
        ));
    }
    let ax = Vector3::new(a, 0.0, 0.0);
    let bx = Vector3::new(b * gamma.cos(), b * sin_gamma, 0.0);
    let cx = c * beta.cos();
    let cy = c * (alpha.cos() - beta.cos() * gamma.cos()) / sin_gamma;
    let cz_squared = c * c - cx * cx - cy * cy;
    if !cz_squared.is_finite() || cz_squared <= 0.0 {
        return Err(frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            Reason::DegenerateCellVectors,
        ));
    }
    PeriodicCell::new(
        Quantity::new([ax, bx, Vector3::new(cx, cy, cz_squared.sqrt())], ANGSTROM),
        [true; 3],
    )
    .map_err(|error| {
        frame_error(
            FORMAT,
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::ReadFrame,
            source_label,
            frame,
            Reason::InvalidDcdCell,
        )
        .with_source(error)
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn read_record<R: Read>(
    reader: &mut R,
    endian: DcdEndian,
    scratch: &mut Vec<u8>,
    other_scratch_bytes: Option<usize>,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    operation: TrajectoryIoOperation,
    frame: Option<usize>,
    clean_eof: bool,
    expected: Option<(usize, Reason)>,
) -> Result<bool, TrajectoryIoError> {
    let mut marker = [0_u8; 4];
    match reader.read(&mut marker[..1]) {
        Ok(0) if clean_eof => return Ok(false),
        Ok(0) => {
            return Err(stream_error(
                FORMAT,
                TrajectoryIoErrorKind::TruncatedRecord,
                operation,
                source_label,
                frame,
                Reason::Missing(Part::RecordMarker),
            ))
        }
        Ok(_) => {}
        Err(error) => return Err(io_context(operation, Some(FORMAT), source_label, error)),
    }
    read_exact(
        FORMAT,
        reader,
        &mut marker[1..],
        operation,
        source_label,
        frame,
        Reason::Incomplete(Part::RecordMarker),
    )?;
    let size = usize::try_from(endian.u32(marker)).map_err(|_| {
        limit_error(
            FORMAT,
            operation,
            source_label,
            frame,
            Reason::AddressRange(Part::RecordLength),
        )
    })?;
    // A record whose size the layout fixes is checked before its payload is
    // read, so that a corrupt marker is never taken for a record still being
    // written.
    if let Some((expected, reason)) = expected.filter(|(expected, _)| *expected != size) {
        let mut context =
            TrajectoryIoError::codec(TrajectoryIoErrorKind::InvalidFrame, operation, Some(FORMAT))
                .with_source_label(source_label)
                .with_counts(expected as u64, size as u64)
                .with_reason(reason);
        if let Some(frame) = frame {
            context = context.with_frame(frame);
        }
        return Err(context);
    }
    if size as u64 > limits.max_record_bytes() {
        return Err(
            limit_error(FORMAT, operation, source_label, frame, None).exceeding(
                Part::Record,
                exceeded_bytes(limits.max_record_bytes(), size as u64),
            ),
        );
    }
    if size > limits.max_scratch_bytes() {
        return Err(
            limit_error(FORMAT, operation, source_label, frame, None).exceeding(
                Part::Record,
                LimitExceeded::new(Resource::BufferBytes, limits.max_scratch_bytes(), size),
            ),
        );
    }
    reserve_scratch(
        scratch,
        size,
        other_scratch_bytes,
        limits.max_scratch_bytes(),
        |reason| limit_error(FORMAT, operation, source_label, frame, reason),
    )?;
    scratch.resize(size, 0);
    read_exact(
        FORMAT,
        reader,
        scratch,
        operation,
        source_label,
        frame,
        Reason::Truncated(Part::RecordPayload),
    )?;
    let mut trailing = [0_u8; 4];
    read_exact(
        FORMAT,
        reader,
        &mut trailing,
        operation,
        source_label,
        frame,
        Reason::Missing(Part::TrailingRecordMarker),
    )?;
    if endian.u32(trailing) != size as u32 {
        let mut context = TrajectoryIoError::codec(
            TrajectoryIoErrorKind::RecordMarkerMismatch,
            operation,
            Some(FORMAT),
        )
        .with_source_label(source_label)
        .with_counts(size as u64, u64::from(endian.u32(trailing)))
        .with_reason(None);
        if let Some(frame) = frame {
            context = context.with_frame(frame);
        }
        return Err(context);
    }
    Ok(true)
}
