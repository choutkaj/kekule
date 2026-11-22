//! XYZ writing: one frame per structure or per geometry of a structure set.

use std::io::Write;

use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, GeometryField, LimitExceeded, Resource,
    UnrepresentablePolicy, WriteOmission, WriteReport,
};
use crate::io::write_policy::{
    check_geometry_state, geometries_with_cells, record_multi_model_structure_set, set_fields,
    unrepresentable_state, Carries, Omissions, StatePolicies,
};
use crate::math::Point3;
use crate::structure::{
    AsStructureView, GeometryView, Positions, Structure, StructureSet, StructureView,
};
use crate::topology::Topology;
use crate::units::{Unit, ANGSTROM, CANONICAL_LENGTH_UNIT};

const MAX_WRITER_COMMENT_BYTES: usize = 1_048_576;

/// What a structure written by [`write_xyz_to`] carries beyond positions: its
/// occupancies and B factors, which are left out and reported under every
/// policy.
const STRUCTURE_CARRIES: Carries = Carries {
    sites: true,
    ..Carries::POSITIONS
};

/// What [`write_xyz_to`] checks of each geometry of a set: its cell and
/// properties. The set checks the geometries' other state.
const GEOMETRY_CARRIES: Carries = Carries {
    place: true,
    time: true,
    step: true,
    velocities: true,
    forces: true,
    sites: true,
    ..Carries::POSITIONS
};

/// The topology content an XYZ file of `topology` leaves out; see
/// [`XyzWriter::report`](crate::traj::io::xyz::XyzWriter::report).
pub(crate) fn topology_omissions(topology: &Topology) -> WriteReport {
    let (mut hydrogens, mut charges, mut isotopes, mut radicals, mut maps, mut sites) =
        (0, 0, 0, 0, 0, 0);
    let mut first_charged = None;
    for &atom_id in topology.atom_ids() {
        let view = topology
            .atom(atom_id)
            .expect("dense topology atom IDs name topology atoms");
        let atom = view.atom();
        hydrogens += usize::from(view.implicit_hydrogens() != Some(0));
        if atom.formal_charge() != 0 {
            charges += 1;
            first_charged.get_or_insert(atom_id);
        }
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
        WriteOmission::FormalCharges {
            atoms: charges,
            first: first_charged,
        },
        WriteOmission::Isotopes { atoms: isotopes },
        WriteOmission::Radicals { atoms: radicals },
        WriteOmission::AtomMaps { atoms: maps },
        WriteOmission::AtomSites { atoms: sites },
        WriteOmission::Stereo { molecules: stereo },
    ]
    .into_iter()
    .collect()
}

/// Options of [`write()`](crate::xyz::write) and of the trajectory codec
/// [`XyzWriter`](crate::traj::io::xyz::XyzWriter).
///
/// The default fails on every cell, frame state, and structure set state XYZ
/// cannot carry.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzWriteOptions {
    length_unit: Unit,
    decimal_places: usize,
    comment: String,
    cell: UnrepresentablePolicy,
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
}

impl Default for XyzWriteOptions {
    fn default() -> Self {
        Self {
            length_unit: ANGSTROM,
            decimal_places: 8,
            comment: "written by kekule".into(),
            cell: UnrepresentablePolicy::Fail,
            geometry_metadata: UnrepresentablePolicy::Fail,
            properties: UnrepresentablePolicy::Fail,
        }
    }
}

/// A write-stage XYZ error without a place.
fn fail(kind: DiagnosticKind, message: impl Into<String>) -> FormatError {
    FormatError::write(Format::Xyz, kind, message)
}

impl XyzWriteOptions {
    /// Sets [`Self::length_unit`].
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] of kind [`DiagnosticKind::InvalidValue`] when
    /// `length_unit` is not a unit of length.
    pub fn with_length_unit(mut self, length_unit: Unit) -> Result<Self, FormatError> {
        CANONICAL_LENGTH_UNIT
            .conversion_factor_to(length_unit)
            .map_err(|error| {
                fail(
                    DiagnosticKind::InvalidValue,
                    format!("The unit {length_unit} is not a length unit."),
                )
                .with_source(error)
            })?;
        self.length_unit = length_unit;
        Ok(self)
    }

    /// Sets [`Self::decimal_places`].
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] of kind [`DiagnosticKind::InvalidValue`] for
    /// more than 17 places.
    pub fn with_decimal_places(mut self, decimal_places: usize) -> Result<Self, FormatError> {
        if decimal_places > 17 {
            return Err(fail(
                DiagnosticKind::InvalidValue,
                format!("The decimal precision is {decimal_places} places, more than 17."),
            ));
        }
        self.decimal_places = decimal_places;
        Ok(self)
    }

    /// Sets [`Self::comment`].
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] of kind [`DiagnosticKind::InvalidValue`] when
    /// `comment` contains a line break, and of kind [`DiagnosticKind::Limit`]
    /// when it exceeds 1 MiB.
    pub fn with_comment(mut self, comment: impl Into<String>) -> Result<Self, FormatError> {
        let comment = comment.into();
        if comment.contains(['\r', '\n']) {
            return Err(fail(
                DiagnosticKind::InvalidValue,
                "The comment contains a line break.",
            ));
        }
        if comment.len() > MAX_WRITER_COMMENT_BYTES {
            let limit = LimitExceeded::new(
                Resource::ValueBytes,
                MAX_WRITER_COMMENT_BYTES,
                comment.len(),
            );
            return Err(fail(
                DiagnosticKind::Limit,
                format!("The comment exceeds its {}.", limit.describe()),
            )
            .exceeding(limit));
        }
        self.comment = comment;
        Ok(self)
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub const fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// Sets [`Self::geometry_metadata`].
    #[must_use]
    pub const fn with_geometry_metadata(mut self, policy: UnrepresentablePolicy) -> Self {
        self.geometry_metadata = policy;
        self
    }

    /// Sets [`Self::properties`].
    #[must_use]
    pub const fn with_properties(mut self, policy: UnrepresentablePolicy) -> Self {
        self.properties = policy;
        self
    }

    /// The length unit coordinates are written in; angstrom by default.
    pub const fn length_unit(&self) -> Unit {
        self.length_unit
    }

    /// Fixed decimal places per coordinate, at most 17; 8 by default.
    pub const fn decimal_places(&self) -> usize {
        self.decimal_places
    }

    /// The comment line written in every frame, without line breaks and at
    /// most 1 MiB; `"written by kekule"` by default.
    pub fn comment(&self) -> &str {
        &self.comment
    }

    /// The periodic cell of a structure, geometry, or frame, which XYZ cannot
    /// carry; by default one with a cell fails to write.
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// The state of a geometry or frame beyond positions, cells, and
    /// properties, which XYZ cannot carry: its time, step, or place on an
    /// order, velocities, and forces, and a frame's occupancies and B
    /// factors. By default one with any of them fails to write, and so does
    /// a set with a time order without times or steps. [`write()`](crate::xyz::write)
    /// leaves the occupancies and B factors of structures and geometries out
    /// under every policy.
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// The atom, bond, and own properties of a structure, geometry, or frame,
    /// which XYZ cannot carry; by default one with any of them fails to
    /// write. So do a structure set's own properties and columns, which an
    /// omission counts once per geometry.
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// The policies of the shared check for a format that carries `carries`.
    pub(crate) const fn state_policies(&self, carries: Carries) -> StatePolicies {
        StatePolicies {
            carries,
            geometry_metadata: self.geometry_metadata,
            properties: self.properties,
            cell: self.cell,
        }
    }

    /// The factor from canonical lengths to [`Self::length_unit`].
    pub(crate) fn output_factor(&self) -> f64 {
        CANONICAL_LENGTH_UNIT
            .conversion_factor_to(self.length_unit)
            .expect("the XYZ output unit is checked to be a length when set")
    }
}

/// What [`write()`](crate::xyz::write) writes: one structure as one frame, or
/// a structure set as one frame per geometry, in geometry order.
///
/// A file of several frames reads back as a structure set without an order,
/// so a set's order is lost: a geometry's place in it follows
/// [`XyzWriteOptions::geometry_metadata`], with its velocities and forces.
#[derive(Debug, Clone, Copy)]
pub struct XyzSource<'a>(Content<'a>);

#[derive(Debug, Clone, Copy)]
enum Content<'a> {
    Structure(StructureView<'a>),
    StructureSet(&'a StructureSet),
}

impl<'a> XyzSource<'a> {
    /// One structure, or any borrowed structure view such as a geometry of a
    /// set, as one frame.
    pub fn structure(structure: &'a (impl AsStructureView + ?Sized)) -> Self {
        Self(Content::Structure(structure.as_structure_view()))
    }

    /// A structure set, one frame per geometry.
    pub fn structure_set(set: &'a StructureSet) -> Self {
        Self(Content::StructureSet(set))
    }
}

impl<'a> From<&'a Structure> for XyzSource<'a> {
    fn from(structure: &'a Structure) -> Self {
        Self::structure(structure)
    }
}

impl<'a> From<&'a StructureSet> for XyzSource<'a> {
    fn from(set: &'a StructureSet) -> Self {
        Self::structure_set(set)
    }
}

impl<'a, R: ?Sized> From<GeometryView<'a, R>> for XyzSource<'a> {
    fn from(geometry: GeometryView<'a, R>) -> Self {
        Self(Content::Structure(geometry.as_structure_view()))
    }
}

/// Whether every coordinate of `positions` stays finite when multiplied by
/// `factor`.
pub(crate) fn finite_in_unit(positions: &Positions, factor: f64) -> bool {
    positions
        .values()
        .value()
        .iter()
        .all(|point| Point3::new(point.x * factor, point.y * factor, point.z * factor).is_finite())
}

/// Writes one frame of `positions`, multiplied by `factor`, with the atom
/// count, the comment of `options`, and one element and three coordinates per
/// atom of `topology`.
pub(crate) fn write_frame(
    writer: &mut impl Write,
    topology: &Topology,
    positions: &Positions,
    factor: f64,
    options: &XyzWriteOptions,
) -> std::io::Result<()> {
    writeln!(writer, "{}", topology.atom_count())?;
    writeln!(writer, "{}", options.comment)?;
    for (point, &atom_id) in positions.values().value().iter().zip(topology.atom_ids()) {
        let atom = topology
            .atom(atom_id)
            .expect("dense topology atom IDs name topology atoms");
        writeln!(
            writer,
            "{} {x:.precision$} {y:.precision$} {z:.precision$}",
            atom.element().symbol(),
            x = point.x * factor,
            y = point.y * factor,
            z = point.z * factor,
            precision = options.decimal_places,
        )?;
    }
    Ok(())
}

/// The error for the `field` of geometry `geometry` (1-based) of a set,
/// which XYZ cannot carry.
fn unrepresentable(field: GeometryField, geometry: usize) -> FormatError {
    let error = fail(
        DiagnosticKind::Unsupported,
        unrepresentable_state(Format::Xyz, field, None),
    )
    .with_detail(FormatErrorDetail::UnrepresentableState(field));
    match field {
        GeometryField::Order | GeometryField::SetProperties | GeometryField::SetSampling => error,
        _ => error.in_record(geometry),
    }
}

/// Writes `source` to a string. See [`write_xyz_to`].
pub(crate) fn write_xyz<'a>(
    source: impl Into<XyzSource<'a>>,
    options: XyzWriteOptions,
) -> Result<(String, WriteReport), FormatError> {
    let mut output = Vec::new();
    let report = write_xyz_to(&mut output, source, options)?;
    Ok((
        String::from_utf8(output).expect("the XYZ writer emits UTF-8"),
        report,
    ))
}

/// Checks the whole source, then streams its frames into `writer`.
pub(crate) fn write_xyz_to<'a>(
    writer: &mut impl Write,
    source: impl Into<XyzSource<'a>>,
    options: XyzWriteOptions,
) -> Result<WriteReport, FormatError> {
    let source = source.into();
    let topology = match source.0 {
        Content::Structure(view) => view.topology(),
        Content::StructureSet(set) => set.topology(),
    };
    if topology.atom_count() == 0 {
        return Err(fail(
            DiagnosticKind::Unsupported,
            "The topology has no atoms.",
        ));
    }
    let factor = options.output_factor();
    let (views, mut omissions) = match source.0 {
        Content::Structure(view) => {
            let omitted =
                check_geometry_state(Format::Xyz, view, options.state_policies(STRUCTURE_CARRIES))?;
            check_finite(view, factor, None)?;
            (vec![view], omitted)
        }
        Content::StructureSet(set) => {
            if set.is_empty() {
                return Err(fail(
                    DiagnosticKind::Inconsistent,
                    "The structure set has no geometries.",
                ));
            }
            // The set's properties and columns fail once for the whole set.
            // Omitted, the check of each geometry below records them.
            if let Some(&field) = set_fields(set).first() {
                if options.properties == UnrepresentablePolicy::Fail {
                    // The columns fail at the first geometry with cells.
                    let geometry = geometries_with_cells(set).next().unwrap_or(0);
                    return Err(unrepresentable(field, geometry + 1));
                }
            }
            let mut omitted = Vec::new();
            let views = set
                .iter()
                .map(|geometry| geometry.as_structure_view())
                .collect::<Vec<_>>();
            // The set-level check below decides about the rest.
            for (record, &view) in (1..).zip(&views) {
                omitted.extend(
                    check_geometry_state(
                        Format::Xyz,
                        view,
                        options.state_policies(GEOMETRY_CARRIES),
                    )
                    .map_err(|error| error.in_record(record))?,
                );
                check_finite(view, factor, Some(record))?;
            }
            let mut state = Omissions::default();
            record_multi_model_structure_set(set, options.geometry_metadata, &mut state)
                .map_err(|(geometry, field)| unrepresentable(field, geometry))?;
            let mut state = state.into_report().omissions().to_vec();
            state.extend(omitted);
            (views, state)
        }
    };
    let count = |present: fn(StructureView<'_>) -> bool| {
        views.iter().filter(|view| present(**view)).count()
    };
    omissions.extend([
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
    for view in &views {
        write_frame(writer, topology, view.positions(), factor, &options)
            .map_err(|error| FormatError::sink(Format::Xyz, error))?;
    }
    Ok(topology_omissions(topology)
        .omissions()
        .iter()
        .copied()
        .chain(omissions)
        .collect())
}

/// Fails unless every coordinate of `view` stays finite in the output unit.
fn check_finite(
    view: StructureView<'_>,
    factor: f64,
    record: Option<usize>,
) -> Result<(), FormatError> {
    if finite_in_unit(view.positions(), factor) {
        return Ok(());
    }
    let error = fail(
        DiagnosticKind::Unsupported,
        "A coordinate is not finite in the output unit.",
    );
    Err(match record {
        Some(record) => error.in_record(record),
        None => error,
    })
}
