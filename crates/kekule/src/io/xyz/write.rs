//! XYZ writer options and the topology content XYZ leaves out.

use crate::diagnostics::{UnrepresentablePolicy, WriteOmission, WriteReport};
use crate::io::write_policy::{Carries, StatePolicies};
use crate::topology::Topology;
use crate::traj::io::codec_context;
use crate::traj::io::xyz::validate_length_unit;
use crate::traj::trajectory::reason::Reason;
use crate::traj::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
};
use crate::units::{Unit, ANGSTROM};

/// The format of the writer whose options these are.
const FORMAT: TrajectoryFormat = TrajectoryFormat::Xyz;

const MAX_WRITER_COMMENT_BYTES: usize = 1_048_576;

/// The label of errors in setting writer options.
const WRITER_OPTIONS: &str = "XYZ writer options";

/// The topology content an XYZ file of `topology` leaves out; see
/// [`XyzWriter::report`](crate::traj::io::xyz::XyzWriter::report).
pub(crate) fn topology_omissions(topology: &Topology) -> WriteReport {
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

/// Deterministic strict XYZ writer options.
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

impl XyzWriteOptions {
    /// Sets [`Self::length_unit`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when
    /// `length_unit` is not a length.
    pub fn with_length_unit(mut self, length_unit: Unit) -> Result<Self, TrajectoryIoError> {
        validate_length_unit(length_unit, WRITER_OPTIONS)?;
        self.length_unit = length_unit;
        Ok(self)
    }

    /// Sets [`Self::decimal_places`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidPrecision`] for more than 17
    /// places.
    pub fn with_decimal_places(mut self, decimal_places: usize) -> Result<Self, TrajectoryIoError> {
        if decimal_places > 17 {
            return Err(codec_context(
                TrajectoryIoErrorKind::InvalidPrecision,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                WRITER_OPTIONS,
                Reason::DecimalPlaces(decimal_places),
            ));
        }
        self.decimal_places = decimal_places;
        Ok(self)
    }

    /// Sets [`Self::comment`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when `comment`
    /// contains a line break, and [`TrajectoryIoErrorKind::LimitExceeded`]
    /// when it exceeds 1 MiB.
    pub fn with_comment(mut self, comment: impl Into<String>) -> Result<Self, TrajectoryIoError> {
        let comment = comment.into();
        if comment.contains(['\r', '\n']) {
            return Err(codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                WRITER_OPTIONS,
                Reason::CommentLineBreak,
            ));
        }
        if comment.len() > MAX_WRITER_COMMENT_BYTES {
            return Err(codec_context(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                WRITER_OPTIONS,
                Reason::CommentTooLong(comment.len()),
            ));
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

    /// The periodic cell of a frame or geometry, which XYZ cannot carry; by
    /// default a frame or geometry with a cell fails to write.
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// A frame's or geometry's state beyond positions, cells, and
    /// properties, which XYZ cannot carry: its time, step, or place on an
    /// order, velocities, forces, occupancies, and B-factors. By default a
    /// frame or geometry with any of them fails to write, and so does a set
    /// with an order or a time order without times or steps among them
    /// ([`XyzWriter::write_structure_set`](crate::traj::io::xyz::XyzWriter::write_structure_set), which leaves a geometry's
    /// occupancies and B-factors out under every policy).
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// A frame's or geometry's atom, bond, and own properties, which XYZ
    /// cannot carry; by default a frame or geometry with any of them fails
    /// to write.
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
}
