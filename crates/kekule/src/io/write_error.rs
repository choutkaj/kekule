//! Failures of the molecule-encoding writers (SMILES and Molfile CTABs).
use crate::diagnostics::{
    DiagnosticKind, Format, FormatError, FormatErrorDetail, LimitExceeded, StructureField,
};
use crate::structure::{
    Distribution, Ensemble, Place, SampleContext, SampleSet, Series, SeriesAxis, StructureKind,
};

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
    /// model without positions for every atom.
    pub(crate) fn invalid_model(self, message: impl Into<String>) -> FormatError {
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
    /// sample or item states, that name a chain written as a segment
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
    /// Sample set sample, ensemble sample, trajectory frame, or series item
    /// state the format cannot carry.
    StructureState {
        /// The state left out.
        field: StructureField,
        /// The samples, frames, or items whose state was left out.
        structures: usize,
    },
}

/// What a structure writer does with the state of an ensemble sample,
/// trajectory frame, or series item that its format cannot carry.
#[derive(Debug, Clone, Copy)]
pub(crate) struct StatePolicies {
    /// `None` for a format that writes sample state itself.
    pub(crate) sample: Option<UnrepresentablePolicy>,
    pub(crate) frame: UnrepresentablePolicy,
    /// `None` for a format that writes an item's step or coordinate itself.
    pub(crate) item: Option<UnrepresentablePolicy>,
    /// The policy for all state of a sample of a [`SampleSet`], in place of
    /// the three above.
    pub(crate) set: UnrepresentablePolicy,
}

impl StatePolicies {
    /// The policies for a sample with `context`: the set policy for every
    /// state the format cannot carry when the sample belongs to a set.
    fn of(self, context: SampleContext<'_>) -> Self {
        if !context.is_of_set() {
            return self;
        }
        Self {
            sample: self.sample.map(|_| self.set),
            frame: self.set,
            item: self.item.map(|_| self.set),
            set: self.set,
        }
    }
}

/// The item state a series' items have that a format without series
/// cannot carry: their steps or their coordinates.
const fn item_field(axis: SeriesAxis) -> StructureField {
    match axis {
        SeriesAxis::Step => StructureField::ItemStep,
        SeriesAxis::Coordinate(_) => StructureField::ItemCoordinate,
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
            | O::StructureState { structures, .. } => structures,
        };
        let added = count(&omission);
        if added == 0 {
            return;
        }
        let same_kind = |existing: &WriteOmission| match (existing, &omission) {
            (O::StructureState { field: a, .. }, O::StructureState { field: b, .. }) => a == b,
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
            | O::StructureState { structures, .. } => *structures += added,
        }
    }

    /// Records the sample state `kind` has that the policies omit: its
    /// membership, its step or coordinate, and its time, step, velocities,
    /// and forces. Call it after [`check_structure_state`] has passed.
    pub(crate) fn record_structure_state(
        &mut self,
        kind: StructureKind<'_>,
        policies: StatePolicies,
    ) {
        let StructureKind::Sample(context) = kind else {
            return;
        };
        let policies = policies.of(context);
        let omit = Some(UnrepresentablePolicy::Omit);
        let mut fields = Vec::new();
        if outside_distribution(context, policies) && policies.set == UnrepresentablePolicy::Omit {
            fields.push(StructureField::Membership);
        }
        if let Some(membership) = context.membership().filter(|_| policies.sample == omit) {
            // A lone record cannot carry a weight relative to other
            // samples, so its weight is always left out.
            fields.push(StructureField::Weight);
            if membership.state_energy().is_some() {
                fields.push(StructureField::Energy);
            }
            if membership.temperature().is_some() {
                fields.push(StructureField::Temperature);
            }
        }
        if policies.item == omit {
            fields.extend(place_field(context));
        }
        if policies.frame == UnrepresentablePolicy::Omit {
            fields.extend(frame_fields(context));
        }
        for field in fields {
            self.record(WriteOmission::StructureState {
                field,
                structures: 1,
            });
        }
    }

    pub(crate) fn into_report(self) -> WriteReport {
        WriteReport { omissions: self.0 }
    }
}

/// Whether a format that writes distribution membership itself writes a
/// sample of a set whose distribution leaves it out, which its records
/// cannot tell from a member.
fn outside_distribution(context: SampleContext<'_>, policies: StatePolicies) -> bool {
    policies.sample.is_none() && context.is_outside_distribution()
}

/// The field a sample's step or coordinate place is written as, if it has
/// one.
fn place_field(context: SampleContext<'_>) -> Option<StructureField> {
    match context.place() {
        Some(Place::Step(_)) => Some(StructureField::ItemStep),
        Some(Place::Coordinate(_)) => Some(StructureField::ItemCoordinate),
        _ => None,
    }
}

/// A sample's time, MD step, velocities, and forces that are present, in
/// that order.
fn frame_fields(context: SampleContext<'_>) -> impl Iterator<Item = StructureField> {
    let (time, step) = match context.place() {
        Some(Place::Time { time, step }) => (time.is_some(), step.is_some()),
        _ => (false, false),
    };
    [
        (time, StructureField::Time),
        (step, StructureField::Step),
        (context.velocities().is_some(), StructureField::Velocities),
        (context.forces().is_some(), StructureField::Forces),
    ]
    .into_iter()
    .filter_map(|(present, field)| present.then_some(field))
}

/// Fails for sample state that `format` cannot carry, unless the caller
/// chose to omit it: a membership, or for a format that writes memberships,
/// a sample outside its set's distribution; a step or coordinate place; and
/// a time, MD step, velocities, or forces, checked in that order.
pub(crate) fn check_structure_state(
    format: Format,
    kind: StructureKind<'_>,
    policies: StatePolicies,
) -> Result<(), FormatError> {
    let StructureKind::Sample(context) = kind else {
        return Ok(());
    };
    let policies = policies.of(context);
    let fail = Some(UnrepresentablePolicy::Fail);
    let membership = context
        .membership()
        .filter(|_| policies.sample == fail)
        .map(|membership| {
            // A lone record cannot carry a weight relative to other samples.
            if membership.state_energy().is_some() {
                StructureField::Energy
            } else {
                StructureField::Weight
            }
        })
        .or_else(|| {
            (outside_distribution(context, policies) && policies.set == UnrepresentablePolicy::Fail)
                .then_some(StructureField::Membership)
        });
    let field = membership
        .or_else(|| place_field(context).filter(|_| policies.item == fail))
        .or_else(|| {
            frame_fields(context)
                .next()
                .filter(|_| policies.frame == UnrepresentablePolicy::Fail)
        });
    match field {
        None => Ok(()),
        Some(field) => Err(FormatError::write(
            format,
            DiagnosticKind::Unsupported,
            format!("{format} cannot represent the {field}."),
        )
        .with_detail(FormatErrorDetail::UnrepresentableState(field))),
    }
}

/// Records or rejects the state of `ensemble` that a multi-model format
/// cannot carry. Its readers give one state of equally weighted samples
/// without an energy or temperature, so only that survives a round trip.
///
/// With [`UnrepresentablePolicy::Fail`], returns the 1-based number of the
/// first sample with such state and the field; the temperature is the first
/// sample's. With [`UnrepresentablePolicy::Omit`], records the omissions.
pub(crate) fn record_multi_model_ensemble(
    ensemble: &Ensemble,
    policy: UnrepresentablePolicy,
    omissions: &mut Omissions,
) -> Result<(), (usize, StructureField)> {
    let weight = ensemble.get(0).map(|sample| sample.weight());
    if policy == UnrepresentablePolicy::Fail {
        for sample in ensemble.iter() {
            let field = if sample.state().energy().is_some() {
                StructureField::Energy
            } else if Some(sample.weight()) != weight {
                StructureField::Weight
            } else if sample.state().index() != 0 {
                StructureField::State
            } else {
                continue;
            };
            return Err((sample.index() + 1, field));
        }
        if ensemble.temperature().is_some() {
            return Err((1, StructureField::Temperature));
        }
        return Ok(());
    }
    let all = ensemble.len();
    let unequal = ensemble
        .iter()
        .any(|sample| Some(sample.weight()) != weight);
    for (field, structures) in [
        (StructureField::Weight, if unequal { all } else { 0 }),
        (
            StructureField::Energy,
            ensemble
                .iter()
                .filter(|sample| sample.state().energy().is_some())
                .count(),
        ),
        (
            StructureField::State,
            if ensemble.state_count() > 1 { all } else { 0 },
        ),
        (
            StructureField::Temperature,
            if ensemble.temperature().is_some() {
                all
            } else {
                0
            },
        ),
    ] {
        omissions.record(WriteOmission::StructureState { field, structures });
    }
    Ok(())
}

/// Records or rejects the steps or coordinates of the items of `series`,
/// which a multi-model format without series cannot carry.
///
/// With [`UnrepresentablePolicy::Fail`], returns the 1-based number of the
/// first item, every item having one, and the field. With
/// [`UnrepresentablePolicy::Omit`], records the omission for every item.
pub(crate) fn record_multi_model_series(
    series: &Series,
    policy: UnrepresentablePolicy,
    omissions: &mut Omissions,
) -> Result<(), (usize, StructureField)> {
    let field = item_field(series.axis());
    if policy == UnrepresentablePolicy::Fail && !series.is_empty() {
        return Err((1, field));
    }
    omissions.record(WriteOmission::StructureState {
        field,
        structures: series.len(),
    });
    Ok(())
}

/// Records or rejects the state of `set` that a multi-model format cannot
/// carry. Its readers give one state of equally weighted samples without an
/// energy, a temperature, or an order, so only that survives a round trip:
/// a sample outside the distribution, unequal weights, more than one state,
/// state energies, the temperature, the samples' places, velocities, and
/// forces are lost.
///
/// With [`UnrepresentablePolicy::Fail`], returns the 1-based number of the
/// first sample with such state and the field, checked for each sample in
/// the order above; the temperature is the first sample's and is checked
/// last. With [`UnrepresentablePolicy::Omit`], records the omissions.
pub(crate) fn record_multi_model_sample_set(
    set: &SampleSet,
    policy: UnrepresentablePolicy,
    omissions: &mut Omissions,
) -> Result<(), (usize, StructureField)> {
    let distribution = set.distribution();
    let contexts = (0..set.len()).map(|index| SampleContext::of_set(set, index));
    let weight = contexts
        .clone()
        .find_map(|context| context.membership())
        .map(|membership| membership.weight());
    // The state of one sample that is lost, in checking order.
    let lost = |context: SampleContext<'_>| {
        let membership = context.membership();
        let mut fields = Vec::new();
        match membership {
            None if distribution.is_some() => fields.push(StructureField::Membership),
            None => {}
            Some(membership) => {
                if membership.state_energy().is_some() {
                    fields.push(StructureField::Energy);
                }
                if Some(membership.weight()) != weight {
                    fields.push(StructureField::Weight);
                }
                if membership.state() != 0 {
                    fields.push(StructureField::State);
                }
            }
        }
        fields.extend(place_field(context));
        let (time, step) = match context.place() {
            Some(Place::Time { time, step }) => (time.is_some(), step.is_some()),
            _ => (false, false),
        };
        fields.extend(
            [
                (time, StructureField::Time),
                (step, StructureField::Step),
                (context.velocities().is_some(), StructureField::Velocities),
                (context.forces().is_some(), StructureField::Forces),
            ]
            .into_iter()
            .filter_map(|(present, field)| present.then_some(field)),
        );
        fields
    };
    let temperature = distribution.and_then(|distribution| distribution.temperature());
    if policy == UnrepresentablePolicy::Fail {
        for (index, context) in contexts.enumerate() {
            if let Some(&field) = lost(context).first() {
                return Err((index + 1, field));
            }
        }
        if temperature.is_some() {
            return Err((1, StructureField::Temperature));
        }
        return Ok(());
    }
    let members = distribution.map_or(0, Distribution::member_count);
    let unequal = contexts
        .clone()
        .filter_map(|context| context.membership())
        .any(|membership| Some(membership.weight()) != weight);
    let count = |field: StructureField| {
        contexts
            .clone()
            .filter(|context| lost(*context).contains(&field))
            .count()
    };
    let states = distribution.map_or(0, Distribution::state_count);
    for (field, structures) in [
        (StructureField::Weight, if unequal { members } else { 0 }),
        (StructureField::Energy, count(StructureField::Energy)),
        (StructureField::State, if states > 1 { members } else { 0 }),
        (
            StructureField::Temperature,
            if temperature.is_some() { members } else { 0 },
        ),
        (
            StructureField::Membership,
            count(StructureField::Membership),
        ),
        (StructureField::Time, count(StructureField::Time)),
        (StructureField::Step, count(StructureField::Step)),
        (StructureField::ItemStep, count(StructureField::ItemStep)),
        (
            StructureField::ItemCoordinate,
            count(StructureField::ItemCoordinate),
        ),
        (
            StructureField::Velocities,
            count(StructureField::Velocities),
        ),
        (StructureField::Forces, count(StructureField::Forces)),
    ] {
        omissions.record(WriteOmission::StructureState { field, structures });
    }
    Ok(())
}

/// The first ensemble sample state that a file of equally weighted structures
/// in one state cannot carry; see [`check_multi_model_ensemble`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct UnrepresentableSampleState {
    sample: usize,
    field: StructureField,
}

impl UnrepresentableSampleState {
    /// The 0-based index of the sample. The temperature is reported for the
    /// first sample.
    pub const fn sample(&self) -> usize {
        self.sample
    }

    /// The state the sample has: [`StructureField::Energy`],
    /// [`StructureField::Weight`], [`StructureField::State`] or
    /// [`StructureField::Temperature`].
    pub const fn field(&self) -> StructureField {
        self.field
    }
}

/// Checks the sample state of `ensemble` that a file of equally weighted
/// structures in one state cannot carry, by the rules of the PDB and mmCIF
/// ensemble writers, for writers of other formats that share them, such as
/// the XYZ writer of [`crate::traj::io::xyz`].
///
/// Such a file reads back as one state of equally weighted samples without
/// an energy or temperature. Unequal weights, more than one state, any state
/// energy, and a temperature are therefore lost: with
/// [`UnrepresentablePolicy::Omit`] the returned report counts the samples
/// each one affects (every sample for weights, states and the temperature),
/// and with [`UnrepresentablePolicy::Fail`] the report is empty when there is
/// nothing to lose.
///
/// # Errors
///
/// With [`UnrepresentablePolicy::Fail`], returns the first sample, in sample
/// order, whose state energy, weight (unequal to the first sample's) or
/// state (other than the first) would be lost, checked in that order for
/// each sample, and otherwise the temperature.
pub fn check_multi_model_ensemble(
    ensemble: &Ensemble,
    policy: UnrepresentablePolicy,
) -> Result<WriteReport, UnrepresentableSampleState> {
    let mut omissions = Omissions::default();
    record_multi_model_ensemble(ensemble, policy, &mut omissions).map_err(|(sample, field)| {
        UnrepresentableSampleState {
            sample: sample - 1,
            field,
        }
    })?;
    Ok(omissions.into_report())
}
