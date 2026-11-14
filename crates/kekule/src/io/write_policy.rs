//! Failures of the molecule-encoding writers (SMILES and Molfile CTABs).
use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, GeometryField, LimitExceeded,
};
use crate::structure::{Order, Place, SetContext, StructureKind, StructureSet, StructureView};

/// Constructors of write-stage [`FormatError`]s for one format, following the
/// writer failure taxonomy: cannot encode, missing perception, limit,
/// inconsistent source, and sink.
#[derive(Debug, Clone, Copy)]
pub(crate) struct WriteErrors(Format);

pub(crate) const SMILES_WRITE: WriteErrors = WriteErrors(Format::Smiles);
pub(crate) const MOLFILE_WRITE: WriteErrors = WriteErrors(Format::Molfile);

impl WriteErrors {
    /// The format, or the requested format version, cannot encode valid
    /// source content without losing or changing it.
    pub(crate) fn unsupported(self, message: impl Into<String>) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::Unsupported, message)
    }

    /// The source lacks perception the writer requires.
    pub(crate) fn missing_perception(self, message: impl Into<String>) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::MissingPerception, message)
    }

    /// A bounded writer calculation could not complete within `limit`;
    /// `subject` starts the sentence, such as `Canonical labeling`.
    pub(crate) fn resource_limit(self, subject: &str, limit: LimitExceeded) -> FormatError {
        FormatError::write(
            self.0,
            DiagnosticKind::Limit,
            format!("{subject} exceeded its {}.", limit.describe()),
        )
        .exceeding(limit)
    }

    /// The source is internally inconsistent or incomplete, for example a
    /// structure without positions for every atom.
    pub(crate) fn invalid_structure(self, message: impl Into<String>) -> FormatError {
        FormatError::write(self.0, DiagnosticKind::Inconsistent, message)
    }

    /// The output sink failed.
    pub(crate) fn io(self, error: std::io::Error) -> FormatError {
        FormatError::sink(self.0, error)
    }
}

/// What a structure writer does with state its format cannot represent, such
/// as a periodic cell in a Molfile. Representable state is always written.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum UnrepresentablePolicy {
    /// Fail rather than write a record that silently lacks the state.
    #[default]
    Fail,
    /// Write the record without the state.
    Omit,
}

/// What a structure writer left out of its output: content the caller chose
/// [`UnrepresentablePolicy::Omit`] for, and content a writer gives no such
/// choice for because its format can never carry it, such as the bonds of
/// the XYZ writer of [`crate::traj::io::xyz`].
///
/// A writer reports every omission it makes, once per kind, with the number
/// of entities it affects summed over everything one call writes. An empty
/// report means the output holds everything the source does that the format
/// can carry.
///
/// Collecting omissions into a report merges repeated kinds by summing their
/// counts and drops zero counts, as writers do.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct WriteReport {
    omissions: Vec<WriteOmission>,
}

impl WriteReport {
    /// The omissions, in the order the writer first made each kind.
    pub fn omissions(&self) -> &[WriteOmission] {
        &self.omissions
    }

    /// Whether the writer omitted nothing.
    pub fn is_empty(&self) -> bool {
        self.omissions.is_empty()
    }
}

impl FromIterator<WriteOmission> for WriteReport {
    fn from_iter<I: IntoIterator<Item = WriteOmission>>(omissions: I) -> Self {
        let mut collected = Omissions::default();
        for omission in omissions {
            collected.record(omission);
        }
        collected.into_report()
    }
}

/// One kind of content a writer omitted, with how much of it.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum WriteOmission {
    /// Implicit-hydrogen counts of this many atoms that readers of the format
    /// do not recover: PDB and mmCIF readers infer counts from valence in
    /// place of fixed ones, and XYZ readers fix every count at zero.
    ImplicitHydrogens {
        /// The atoms whose counts were left out.
        atoms: usize,
    },
    /// Nonzero formal charges of this many atoms.
    FormalCharges {
        /// The atoms written without their charges.
        atoms: usize,
    },
    /// Isotopes the format cannot name, of this many atoms.
    Isotopes {
        /// The atoms written as their element alone.
        atoms: usize,
    },
    /// Radical electrons of this many atoms.
    Radicals {
        /// The atoms written without their radicals.
        atoms: usize,
    },
    /// Atom-map numbers of this many atoms.
    AtomMaps {
        /// The atoms written without their atom maps.
        atoms: usize,
    },
    /// Represented stereo of this many molecule instances.
    Stereo {
        /// The molecule instances written without their stereo.
        molecules: usize,
    },
    /// Named hierarchy atom sites of this many atoms: their atom names and
    /// their residues and chains.
    AtomSites {
        /// The atoms written without their sites.
        atoms: usize,
    },
    /// Bonds the format's connectivity records cannot state, left out.
    Bonds {
        /// The bonds left out.
        bonds: usize,
    },
    /// Metal coordination records of interpretation reports, of this many
    /// records, that the output does not state: records that no longer
    /// match the written structure, and, in PDB, records that not every
    /// geometry or item states, that name a chain written as a segment
    /// identifier, whose partners are not one metal and one non-metal, or
    /// whose partner identifiers also name another atom. A record that
    /// several reports state counts once.
    Coordination {
        /// The records left out.
        records: usize,
    },
    /// Dative bonds between residues, of this many bonds, written as metal
    /// coordination (mmCIF `metalc`, PDB `CONECT`), which macromolecular
    /// readers record in their report rather than as bonds.
    CoordinationBonds {
        /// The bonds written so.
        bonds: usize,
    },
    /// Residues whose bonds differ from the component dictionary that
    /// readers apply to them, so readers will see the dictionary's bonds.
    DictionaryBonds {
        /// The residues affected.
        residues: usize,
    },
    /// A periodic cell the format cannot carry, of this many written
    /// structures or blocks.
    Cell {
        /// The structures or blocks written without their cell.
        structures: usize,
    },
    /// Crystallographic occupancies of this many written structures.
    Occupancies {
        /// The structures written without their occupancies.
        structures: usize,
    },
    /// Crystallographic B factors of this many written structures.
    BFactors {
        /// The structures written without their B factors.
        structures: usize,
    },
    /// Chain identifiers the format's chain field cannot hold, of this many
    /// chains, written in another field that readers taking chains from the
    /// chain field alone do not see: PDB segment identifiers
    /// ([`crate::pdb::PdbChainIdPolicy::SegmentId`]).
    ChainIds {
        /// The chains written so.
        chains: usize,
    },
    /// Geometry or trajectory frame state the format cannot carry.
    GeometryState {
        /// The state left out.
        field: GeometryField,
        /// The geometries or frames whose state was left out.
        geometries: usize,
    },
}

/// What a structure writer does with the state of a structure that its
/// format cannot carry: a geometry's, a structure set's, or a trajectory frame's
/// state, a structure's occupancies and B factors, and properties.
#[derive(Debug, Clone, Copy)]
pub(crate) struct StatePolicies {
    /// The policy for geometry state, and for occupancies and B factors.
    pub(crate) policy: UnrepresentablePolicy,
    /// The policy for per-atom, per-bond, and own structure properties and
    /// for the properties of a geometry's set.
    pub(crate) properties: UnrepresentablePolicy,
    /// Whether the format writes a geometry's step or coordinate place itself.
    pub(crate) writes_place: bool,
    /// Whether the format writes occupancies and B factors.
    pub(crate) writes_sites: bool,
}

impl StatePolicies {
    /// The policies of a format that writes coordinates and nothing else of
    /// a structure's state.
    pub(crate) const fn positions_only(
        policy: UnrepresentablePolicy,
        properties: UnrepresentablePolicy,
    ) -> Self {
        Self {
            policy,
            properties,
            writes_place: false,
            writes_sites: false,
        }
    }

    /// The policies of a format that also writes occupancies and B factors.
    pub(crate) const fn with_sites(
        policy: UnrepresentablePolicy,
        properties: UnrepresentablePolicy,
    ) -> Self {
        Self {
            writes_sites: true,
            ..Self::positions_only(policy, properties)
        }
    }

    /// The policy that decides about `field`.
    const fn of(self, field: GeometryField) -> UnrepresentablePolicy {
        match field {
            GeometryField::AtomProperties
            | GeometryField::BondProperties
            | GeometryField::Properties
            | GeometryField::SetProperties => self.properties,
            _ => self.policy,
        }
    }
}

/// The sentence for the `field` of a geometry that `format` cannot carry: of
/// the 1-based geometry `geometry` of a set, or of a lone geometry for `None`. The
/// order belongs to the whole set.
pub(crate) fn unrepresentable_state(
    format: impl std::fmt::Display,
    field: GeometryField,
    geometry: Option<usize>,
) -> String {
    match (field, geometry) {
        (GeometryField::Order, _) => format!("The geometry order has no {format} representation."),
        (GeometryField::SetProperties, _) => {
            format!("The structure set properties have no {format} representation.")
        }
        (field, Some(geometry)) => {
            let verb = verb(field);
            format!("The {field} of geometry {geometry} {verb} no {format} representation.")
        }
        (field, None) => {
            let verb = verb(field);
            format!("The {field} of the geometry {verb} no {format} representation.")
        }
    }
}

/// The sentence for the `field` of a lone structure that `format` cannot
/// carry: of the geometry for a geometry or a frame, of the structure otherwise.
fn unrepresentable_geometry_state(
    format: Format,
    field: GeometryField,
    kind: StructureKind<'_>,
) -> String {
    match kind {
        StructureKind::InSet(_) => unrepresentable_state(format, field, None),
        StructureKind::Standalone => {
            let verb = verb(field);
            format!("The {field} of the structure {verb} no {format} representation.")
        }
    }
}

/// The verb for `field` as the subject: "have" for a plural field, "has"
/// otherwise.
const fn verb(field: GeometryField) -> &'static str {
    match field {
        GeometryField::Velocities
        | GeometryField::Forces
        | GeometryField::AtomProperties
        | GeometryField::BondProperties
        | GeometryField::Properties
        | GeometryField::SetProperties
        | GeometryField::Occupancies
        | GeometryField::BFactors => "have",
        _ => "has",
    }
}

/// Collects a writer's omissions, merging repeated kinds by summing counts.
#[derive(Debug, Default)]
pub(crate) struct Omissions(Vec<WriteOmission>);

impl Omissions {
    /// Records `omission`, adding its count to an earlier one of its kind.
    /// Zero counts record nothing.
    pub(crate) fn record(&mut self, omission: WriteOmission) {
        use WriteOmission as O;
        let count = |omission: &WriteOmission| match *omission {
            O::ImplicitHydrogens { atoms }
            | O::FormalCharges { atoms }
            | O::Isotopes { atoms }
            | O::Radicals { atoms }
            | O::AtomMaps { atoms }
            | O::AtomSites { atoms } => atoms,
            O::Stereo { molecules } => molecules,
            O::Bonds { bonds } => bonds,
            O::Coordination { records } => records,
            O::CoordinationBonds { bonds } => bonds,
            O::DictionaryBonds { residues } => residues,
            O::ChainIds { chains } => chains,
            O::Cell { structures }
            | O::Occupancies { structures }
            | O::BFactors { structures }
            | O::GeometryState {
                geometries: structures,
                ..
            } => structures,
        };
        let added = count(&omission);
        if added == 0 {
            return;
        }
        let same_kind = |existing: &WriteOmission| match (existing, &omission) {
            (O::GeometryState { field: a, .. }, O::GeometryState { field: b, .. }) => a == b,
            _ => std::mem::discriminant(existing) == std::mem::discriminant(&omission),
        };
        let Some(existing) = self.0.iter_mut().find(|existing| same_kind(existing)) else {
            self.0.push(omission);
            return;
        };
        match existing {
            O::ImplicitHydrogens { atoms }
            | O::FormalCharges { atoms }
            | O::Isotopes { atoms }
            | O::Radicals { atoms }
            | O::AtomMaps { atoms }
            | O::AtomSites { atoms } => *atoms += added,
            O::Stereo { molecules } => *molecules += added,
            O::Bonds { bonds } => *bonds += added,
            O::Coordination { records } => *records += added,
            O::CoordinationBonds { bonds } => *bonds += added,
            O::DictionaryBonds { residues } => *residues += added,
            O::ChainIds { chains } => *chains += added,
            O::Cell { structures }
            | O::Occupancies { structures }
            | O::BFactors { structures }
            | O::GeometryState {
                geometries: structures,
                ..
            } => *structures += added,
        }
    }

    /// Records the state of the lone structure `view` that the policies
    /// omit, in the order of [`lost_geometry_state`]. Call it after
    /// [`check_geometry_state`] has passed.
    pub(crate) fn record_geometry_state(
        &mut self,
        view: StructureView<'_>,
        policies: StatePolicies,
    ) {
        for field in lost_geometry_state(view, policies) {
            if policies.of(field) == UnrepresentablePolicy::Omit {
                self.record(structure_omission(field, 1));
            }
        }
    }

    pub(crate) fn into_report(self) -> WriteReport {
        WriteReport { omissions: self.0 }
    }
}

/// The omission of `field` from `structures` structures: occupancies and B
/// factors have omissions of their own.
fn structure_omission(field: GeometryField, structures: usize) -> WriteOmission {
    match field {
        GeometryField::Occupancies => WriteOmission::Occupancies { structures },
        GeometryField::BFactors => WriteOmission::BFactors { structures },
        field => WriteOmission::GeometryState {
            field,
            geometries: structures,
        },
    }
}

/// The state of the lone structure `view` that a format with `policies`
/// cannot carry, whatever the policies decide about it, in checking order:
/// its step or coordinate place, unless the format writes places; its time,
/// MD step, velocities, and forces; the structure's occupancies and B
/// factors, unless the format writes them; and its per-atom, per-bond, and
/// own properties and those of its set.
fn lost_geometry_state(view: StructureView<'_>, policies: StatePolicies) -> Vec<GeometryField> {
    let mut fields = Vec::new();
    let context = match view.kind() {
        StructureKind::InSet(context) => Some(context),
        StructureKind::Standalone => None,
    };
    if let Some(context) = context {
        if !policies.writes_place {
            fields.extend(place_field(context));
        }
        fields.extend(frame_fields(context));
    }
    if !policies.writes_sites {
        fields.extend(site_fields(view));
    }
    fields.extend(property_fields(view));
    if context
        .and_then(SetContext::set_properties)
        .is_some_and(|properties| !properties.is_empty())
    {
        fields.push(GeometryField::SetProperties);
    }
    fields
}

/// The occupancies and B factors `view` records.
fn site_fields(view: StructureView<'_>) -> impl Iterator<Item = GeometryField> {
    [
        (
            view.occupancies()
                .is_some_and(|values| values.iter().any(Option::is_some)),
            GeometryField::Occupancies,
        ),
        (
            view.b_factors()
                .is_some_and(|values| values.value().iter().any(Option::is_some)),
            GeometryField::BFactors,
        ),
    ]
    .into_iter()
    .filter_map(|(present, field)| present.then_some(field))
}

/// The per-atom, per-bond, and own properties `view` has, in that order.
fn property_fields(view: StructureView<'_>) -> impl Iterator<Item = GeometryField> {
    let properties = view.properties();
    [
        (properties.atoms().has_data(), GeometryField::AtomProperties),
        (properties.bonds().has_data(), GeometryField::BondProperties),
        (!properties.owner().is_empty(), GeometryField::Properties),
    ]
    .into_iter()
    .filter_map(|(present, field)| present.then_some(field))
}

/// The field a geometry's step or coordinate place is written as, if it has
/// one.
fn place_field(context: SetContext<'_>) -> Option<GeometryField> {
    match context.place() {
        Some(Place::Step(_)) => Some(GeometryField::OrderStep),
        Some(Place::Coordinate(_)) => Some(GeometryField::OrderCoordinate),
        _ => None,
    }
}

/// A geometry's time, MD step, velocities, and forces that are present, in
/// that order.
fn frame_fields(context: SetContext<'_>) -> impl Iterator<Item = GeometryField> {
    let (time, step) = match context.place() {
        Some(Place::Time { time, step }) => (time.is_some(), step.is_some()),
        _ => (false, false),
    };
    [
        (time, GeometryField::Time),
        (step, GeometryField::Step),
        (context.velocities().is_some(), GeometryField::Velocities),
        (context.forces().is_some(), GeometryField::Forces),
    ]
    .into_iter()
    .filter_map(|(present, field)| present.then_some(field))
}

/// Fails for the first state of the lone structure `view` that `format`
/// cannot carry and the policies do not omit, in the order of
/// [`lost_geometry_state`].
pub(crate) fn check_geometry_state(
    format: Format,
    view: StructureView<'_>,
    policies: StatePolicies,
) -> Result<(), FormatError> {
    let lost = lost_geometry_state(view, policies);
    let field = lost
        .iter()
        .copied()
        .find(|&field| policies.of(field) == UnrepresentablePolicy::Fail);
    match field {
        None => Ok(()),
        Some(field) => Err(FormatError::write(
            format,
            DiagnosticKind::Unsupported,
            unrepresentable_geometry_state(format, field, view.kind()),
        )
        .with_detail(FormatErrorDetail::UnrepresentableState(field))),
    }
}

/// Records or rejects the properties of `set` and of its geometries that a
/// multi-model format cannot carry: the set's own properties, then each
/// geometry's per-atom, per-bond, and own properties. With
/// [`UnrepresentablePolicy::Fail`], returns the 1-based number of the first
/// geometry with such properties and the field, geometry 1 for the set's.
/// With [`UnrepresentablePolicy::Omit`], records the omissions, the set's
/// counting every geometry.
pub(crate) fn record_structure_set_properties(
    set: &StructureSet,
    policy: UnrepresentablePolicy,
    omissions: &mut Omissions,
) -> Result<(), (usize, GeometryField)> {
    let set_properties = !set.properties().is_empty() && !set.is_empty();
    if policy == UnrepresentablePolicy::Fail {
        if set_properties {
            return Err((1, GeometryField::SetProperties));
        }
        for (index, geometry) in set.iter().enumerate() {
            if let Some(field) = property_fields(geometry.as_structure_view()).next() {
                return Err((index + 1, field));
            }
        }
        return Ok(());
    }
    if set_properties {
        omissions.record(WriteOmission::GeometryState {
            field: GeometryField::SetProperties,
            geometries: set.len(),
        });
    }
    for geometry in set {
        for field in property_fields(geometry.as_structure_view()) {
            omissions.record(structure_omission(field, 1));
        }
    }
    Ok(())
}

/// Records or rejects the state of `set` that a multi-model format cannot
/// carry. Its readers give geometries without an order, so the geometries'
/// places, velocities, and forces, and a time order without times or steps,
/// are lost.
///
/// With [`UnrepresentablePolicy::Fail`], returns the 1-based number of the
/// first geometry with such state and the field, checked for each geometry in
/// the order above. The time order is the first geometry's and is checked
/// last. With [`UnrepresentablePolicy::Omit`], records the omissions.
pub(crate) fn record_multi_model_structure_set(
    set: &StructureSet,
    policy: UnrepresentablePolicy,
    omissions: &mut Omissions,
) -> Result<(), (usize, GeometryField)> {
    let contexts = set.iter().map(|view| SetContext::of_view(view.as_dyn()));
    // The state of one geometry that is lost, in checking order.
    let lost = |context: SetContext<'_>| {
        let mut fields = Vec::new();
        fields.extend(place_field(context));
        fields.extend(frame_fields(context));
        fields
    };
    let bare_order = matches!(
        set.order(),
        Some(Order::Time(axis)) if axis.times().is_none() && axis.steps().is_none()
    );
    if policy == UnrepresentablePolicy::Fail {
        for (index, context) in contexts.enumerate() {
            if let Some(&field) = lost(context).first() {
                return Err((index + 1, field));
            }
        }
        if bare_order && !set.is_empty() {
            return Err((1, GeometryField::Order));
        }
        return Ok(());
    }
    let count = |field: GeometryField| {
        contexts
            .clone()
            .filter(|context| lost(*context).contains(&field))
            .count()
    };
    for (field, geometries) in [
        (GeometryField::Order, if bare_order { set.len() } else { 0 }),
        (GeometryField::Time, count(GeometryField::Time)),
        (GeometryField::Step, count(GeometryField::Step)),
        (GeometryField::OrderStep, count(GeometryField::OrderStep)),
        (
            GeometryField::OrderCoordinate,
            count(GeometryField::OrderCoordinate),
        ),
        (GeometryField::Velocities, count(GeometryField::Velocities)),
        (GeometryField::Forces, count(GeometryField::Forces)),
    ] {
        omissions.record(WriteOmission::GeometryState { field, geometries });
    }
    Ok(())
}
