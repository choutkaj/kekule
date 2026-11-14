use super::*;

/// The topology content an XYZ file of `topology` leaves out; see
/// [`XyzWriter::report`].
pub(super) fn topology_omissions(topology: &Topology) -> WriteReport {
    let (mut hydrogens, mut charges, mut isotopes, mut radicals, mut maps, mut sites) =
        (0, 0, 0, 0, 0, 0);
    for &atom_id in topology.atom_ids() {
        let view = topology
            .atom(atom_id)
            .expect("dense topology atom IDs name topology atoms");
        let atom = view.atom();
        hydrogens += usize::from(view.implicit_hydrogens() != Some(0));
        charges += usize::from(atom.formal_charge() != 0);
        isotopes += usize::from(atom.isotope().is_some());
        radicals += usize::from(atom.radical().is_some());
        maps += usize::from(atom.atom_map().is_some());
        sites += usize::from(view.atom_site().is_some_and(|site| {
            let metadata = site.metadata();
            metadata.label_atom_id().is_some() || metadata.auth_atom_id().is_some()
        }));
    }
    let stereo = topology
        .instances()
        .filter(|instance| {
            let molecule = instance.definition().molecule();
            molecule.stereo_elements().next().is_some() || molecule.stereo_groups().next().is_some()
        })
        .count();
    [
        WriteOmission::Bonds {
            bonds: topology.bond_count(),
        },
        WriteOmission::ImplicitHydrogens { atoms: hydrogens },
        WriteOmission::FormalCharges { atoms: charges },
        WriteOmission::Isotopes { atoms: isotopes },
        WriteOmission::Radicals { atoms: radicals },
        WriteOmission::AtomMaps { atoms: maps },
        WriteOmission::AtomSites { atoms: sites },
        WriteOmission::Stereo { molecules: stereo },
    ]
    .into_iter()
    .collect()
}

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
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when the length
    /// unit is not a length or the comment contains a line break,
    /// [`TrajectoryIoErrorKind::InvalidPrecision`] for more than 17 decimal
    /// places, and [`TrajectoryIoErrorKind::LimitExceeded`] for a comment over
    /// 1 MiB or a topology without atoms.
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
        options
            .length_unit
            .conversion_factor_to(CANONICAL_LENGTH_UNIT)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::Open,
                    Some(FORMAT),
                    &source_label,
                    "The XYZ length unit is not a length.",
                )
                .with_source(error)
            })?;
        if options.decimal_places > 17 {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidPrecision,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                format!(
                    "The decimal precision is {} places, more than 17.",
                    options.decimal_places
                ),
            ));
        }
        if options.comment.contains(['\r', '\n']) {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                "The comment contains a line break.",
            ));
        }
        if options.comment.len() > MAX_WRITER_COMMENT_BYTES {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                format!(
                    "The comment has {} bytes, more than {MAX_WRITER_COMMENT_BYTES}.",
                    options.comment.len()
                ),
            ));
        }
        if topology.atom_count() == 0 {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                &source_label,
                "The topology has no atoms.",
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

    /// What the stream written so far leaves out: XYZ carries each atom's
    /// element and position only.
    ///
    /// From construction, the report holds the topology content every frame
    /// leaves out: bonds; atoms whose implicit-hydrogen count is nonzero, or
    /// inferred and not yet perceived, since XYZ readers fix every count at
    /// zero; atoms with a nonzero formal charge, an isotope, radical
    /// electrons, or an atom-map number; atoms whose hierarchy atom site has
    /// an atom name, which is lost with the site's residue and chain; and
    /// molecule instances with represented stereo. Sites without a name,
    /// such as the synthetic `UNL` residues `kekule::xyz` gives an
    /// interpreted file, are not counted.
    ///
    /// Each frame or structure set written then adds what it left out, counted
    /// over every frame written: state omitted under
    /// [`XyzWriteOptions::cell`], [`XyzWriteOptions::geometry_metadata`], and
    /// [`XyzWriteOptions::properties`], and geometries' occupancies and B
    /// factors. The writer writes the same text whatever it reports.
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
        let factor = self.output_factor()?;
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
    /// Writes one frame per geometry of `set`, in geometry order, with the
    /// geometry's positions.
    ///
    /// The whole set is checked before its first frame is written. The
    /// geometries' places, velocities, and forces follow
    /// [`XyzWriteOptions::geometry_metadata`], their atom, bond, and own
    /// properties [`XyzWriteOptions::properties`], and their cells
    /// [`XyzWriteOptions::cell`]; occupancies and B factors are always left
    /// out. [`Self::report`] gains everything left out.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::TopologyMismatch`] when the set does
    /// not share the writer's topology layout;
    /// [`TrajectoryIoErrorKind::InvalidFrame`] for an empty set, for a
    /// coordinate that is not finite after conversion, or after an earlier
    /// write left a partial frame; [`TrajectoryIoErrorKind::UnsupportedField`]
    /// naming the [`GeometryField`] for set or state properties and for a
    /// cell, geometry property, or geometry state the options do not omit, with
    /// [`TrajectoryIoError::frame`] the frame of the first geometry concerned;
    /// and [`TrajectoryIoErrorKind::Io`] when the stream fails.
    pub fn write_structure_set(&mut self, set: &StructureSet) -> Result<(), TrajectoryIoError> {
        require_no_partial_frame(
            self.partial,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            &self.source_label,
        )?;
        if !self.topology().shares_layout(set.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        if set.is_empty() {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                "The structure set has no geometries.",
            ));
        }
        // The frame a geometry is written as, for errors.
        let frame = |geometry: usize| self.frame_count.saturating_add(geometry);
        if !set.properties().is_empty() {
            return Err(writer_field_error(
                TrajectoryFormat::Xyz,
                &self.source_label,
                frame(0),
                GeometryField::SetProperties,
            ));
        }
        let factor = self.output_factor()?;
        let mut omitted = Vec::new();
        for (index, geometry) in set.iter().enumerate() {
            let view = geometry.as_structure_view();
            // The set-level check below decides about the rest.
            let state = self
                .options
                .state_policies(GEOMETRY_CARRIES)
                .check(view)
                .map_err(|field| {
                    writer_field_error(
                        TrajectoryFormat::Xyz,
                        &self.source_label,
                        frame(index),
                        field,
                    )
                })?;
            self.check_positions(view.positions(), factor)?;
            omitted.extend(state);
        }
        let mut state = Omissions::default();
        record_multi_model_structure_set(set, self.options.geometry_metadata, &mut state).map_err(
            |(geometry, field)| {
                writer_field_error(
                    TrajectoryFormat::Xyz,
                    &self.source_label,
                    frame(geometry - 1),
                    field,
                )
            },
        )?;
        for geometry in set {
            self.write_positions(geometry.as_structure_view().positions(), factor)?;
        }
        let count = |present: fn(StructureView<'_>) -> bool| {
            set.iter()
                .filter(|geometry| present(geometry.as_structure_view()))
                .count()
        };
        omitted.extend([
            WriteOmission::Occupancies {
                structures: count(|view| {
                    view.occupancies()
                        .is_some_and(|values| values.iter().any(Option::is_some))
                }),
            },
            WriteOmission::BFactors {
                structures: count(|view| {
                    view.b_factors()
                        .is_some_and(|values| values.value().iter().any(Option::is_some))
                }),
            },
        ]);
        self.record(
            state
                .into_report()
                .omissions()
                .iter()
                .copied()
                .chain(omitted),
        );
        Ok(())
    }

    /// Adds `omissions` to [`Self::report`].
    pub(super) fn record(&mut self, omissions: impl IntoIterator<Item = WriteOmission>) {
        record(&mut self.report, omissions);
    }

    /// The factor from canonical lengths to the output unit.
    pub(super) fn output_factor(&self) -> Result<f64, TrajectoryIoError> {
        CANONICAL_LENGTH_UNIT
            .conversion_factor_to(self.options.length_unit)
            .map_err(|error| {
                codec_context(
                    TrajectoryIoErrorKind::InconsistentMetadata,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    "The XYZ output unit is not a length.",
                )
                .with_source(error)
            })
    }

    /// Fails unless every coordinate stays finite after conversion.
    pub(super) fn check_positions(
        &self,
        positions: &Positions,
        factor: f64,
    ) -> Result<(), TrajectoryIoError> {
        for point in positions.values().value().iter() {
            if !Point3::new(point.x * factor, point.y * factor, point.z * factor).is_finite() {
                return Err(codec_context(
                    TrajectoryIoErrorKind::InvalidFrame,
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    "A coordinate is not finite in the output unit.",
                ));
            }
        }
        Ok(())
    }

    /// Writes one frame of checked positions.
    pub(super) fn write_positions(
        &mut self,
        positions: &Positions,
        factor: f64,
    ) -> Result<(), TrajectoryIoError> {
        let positions = positions.values();
        self.partial = true;
        writeln!(self.writer, "{}", self.topology.atom_count()).map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })?;
        writeln!(self.writer, "{}", self.options.comment).map_err(|error| {
            io_context(
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                error,
            )
        })?;
        for (atom_index, point) in positions.value().iter().enumerate() {
            let atom_id = self.topology.atom_ids()[atom_index];
            let element = self
                .topology
                .atom(atom_id)
                .expect("dense topology atom IDs name topology atoms");
            writeln!(
                self.writer,
                "{} {x:.precision$} {y:.precision$} {z:.precision$}",
                element.element().symbol(),
                x = point.x * factor,
                y = point.y * factor,
                z = point.z * factor,
                precision = self.options.decimal_places,
            )
            .map_err(|error| {
                io_context(
                    TrajectoryIoOperation::WriteFrame,
                    Some(FORMAT),
                    &self.source_label,
                    error,
                )
            })?;
        }
        self.partial = false;
        self.frame_count = self.frame_count.checked_add(1).ok_or_else(|| {
            codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::WriteFrame,
                Some(FORMAT),
                &self.source_label,
                "The frame count exceeds this platform's address range.",
            )
        })?;
        Ok(())
    }
}
