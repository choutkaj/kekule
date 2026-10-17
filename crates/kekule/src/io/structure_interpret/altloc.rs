//! Alternate locations are source alternatives, not additional coordinate models.
//!
//! Selection keeps unlabelled shared atoms and chooses whole residues (or explicit
//! groups). Occupancy ranking is a deterministic representative-selection heuristic,
//! not a joint probability or a guarantee of compatibility between unrelated groups.

mod select;
pub(super) mod source;

use std::collections::BTreeMap;

use super::atom_site::AtomRow;

pub(super) use select::{inventory, select_alt_locations};

/// Policy for reducing source alternate locations to one coordinate set.
///
/// The default ranks complete residue alternatives by their mean labelled-atom
/// occupancy. Missing occupancies contribute zero and are reported; shared atoms
/// do not affect ranking. Explicit groups use the mean of their residue scores.
/// Ties choose the lexicographically smallest label and are reported. Incomplete
/// alternatives are reported and excluded, never filled from another label.
/// Different chemical components at one residue position are mutually exclusive.
///
/// Supplied `_atom_sites_alt_gen` configurations constrain selection across the
/// structure. Only declared combinations are considered, without enumeration or
/// occupancy-derived ensemble weights. In their absence, residues are independent
/// unless grouped explicitly. Source coordinate models are always separate scopes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AltLocPolicy {
    HighestOccupancy,
    /// Require this label at every residue with alternatives; retain shared atoms.
    SelectLabel(String),
    /// Prefer this label, falling back to occupancy ranking where unavailable.
    PreferLabel(String),
    /// Reject any labelled alternate site, even when only one label is present.
    ErrorOnAlternateLocations,
    Configured(AltLocSelection),
}

/// Default choice for residues not covered by an exact override.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum AltLocPreference {
    #[default]
    HighestOccupancy,
    SelectLabel(String),
    PreferLabel(String),
    ErrorOnAlternateLocations,
}

/// Exact choices and optional correlation constraints for one source block.
///
/// Obtain source identities from [`super::super::MmcifBlock::alternate_locations`].
/// Overrides take precedence over the default preference. Every override and group member must
/// identify a residue with alternatives. Unknown identities and conflicting group
/// choices are errors, including those in unselected source coordinate models.
///
/// ```
/// use kekule::diagnostics::FormatError;
/// use kekule::macromolecular::{AltLocPolicy, AltLocPreference, AltLocSelection,
///     StructureInterpretOptions, StructureInterpretation};
/// use kekule::mmcif::MmcifBlock;
/// # fn choose(block: &MmcifBlock) -> Result<StructureInterpretation, FormatError> {
/// let mut selection =
///     AltLocSelection::default().with_default_preference(AltLocPreference::PreferLabel("A".into()));
/// // Override one source residue, using its exact block-scoped identity.
/// if let Some(residue) = block.alternate_locations()?.into_iter()
///     .find(|residue| residue.labels().iter().any(|label| label == "B"))
/// {
///     selection = selection.with_override(residue.residue().clone(), "B");
/// }
/// block.interpret(
///     StructureInterpretOptions::default()
///         .with_altloc_policy(AltLocPolicy::Configured(selection)),
/// )
/// # }
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AltLocSelection {
    pub(crate) default: AltLocPreference,
    pub(crate) overrides: BTreeMap<SourceResidueId, String>,
    pub(crate) groups: Vec<Vec<SourceResidueId>>,
    pub(crate) source_conformation: Option<String>,
}

impl AltLocSelection {
    /// The preference for residues without an override.
    pub fn default_preference(&self) -> &AltLocPreference {
        &self.default
    }

    /// Sets [`Self::default_preference`].
    #[must_use]
    pub fn with_default_preference(mut self, default: AltLocPreference) -> Self {
        self.default = default;
        self
    }

    /// Exact label choices that take precedence over the default preference.
    pub fn overrides(&self) -> &BTreeMap<SourceResidueId, String> {
        &self.overrides
    }

    /// Chooses `label` for one source residue, replacing an earlier override.
    #[must_use]
    pub fn with_override(mut self, residue: SourceResidueId, label: impl Into<String>) -> Self {
        self.overrides.insert(residue, label.into());
        self
    }

    /// Replaces every override.
    #[must_use]
    pub fn with_overrides(mut self, overrides: BTreeMap<SourceResidueId, String>) -> Self {
        self.overrides = overrides;
        self
    }

    /// Residues constrained to use the same label. Groups must contain at least
    /// two distinct residues in one source model; overlapping groups are united.
    pub fn groups(&self) -> &[Vec<SourceResidueId>] {
        &self.groups
    }

    /// Sets [`Self::groups`].
    #[must_use]
    pub fn with_groups(mut self, groups: Vec<Vec<SourceResidueId>>) -> Self {
        self.groups = groups;
        self
    }

    /// Select one `_atom_sites_alt_gen.ens_id` combination explicitly.
    /// This identifier does not create or select an `Ensemble` member.
    pub fn source_conformation(&self) -> Option<&str> {
        self.source_conformation.as_deref()
    }

    /// Sets [`Self::source_conformation`].
    #[must_use]
    pub fn with_source_conformation(mut self, source_conformation: Option<String>) -> Self {
        self.source_conformation = source_conformation;
        self
    }
}

/// A source-scoped residue identity, independent of selected chemistry and rows.
///
/// `asym_id` uses the label namespace, falling back to author only when absent;
/// a PDB entry has only the author namespace. Position similarly prefers label
/// sequence, then author sequence, then source occurrence. Insertion codes and
/// source coordinate models remain significant.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SourceResidueId {
    pub(crate) model_id: String,
    pub(crate) asym_id: String,
    pub(crate) position: SourceResiduePosition,
    pub(crate) insertion_code: Option<String>,
}

impl SourceResidueId {
    /// A residue identity from its source-scoped parts, for example to name
    /// an [`AltLocSelection`] override without reading the block first.
    pub fn from_parts(
        model_id: impl Into<String>,
        asym_id: impl Into<String>,
        position: SourceResiduePosition,
        insertion_code: Option<String>,
    ) -> Self {
        Self {
            model_id: model_id.into(),
            asym_id: asym_id.into(),
            position,
            insertion_code,
        }
    }

    pub fn model_id(&self) -> &str {
        &self.model_id
    }

    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }

    pub fn position(&self) -> &SourceResiduePosition {
        &self.position
    }

    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SourceResiduePosition {
    Label(i32),
    Author(String),
    /// An author sequence number that an earlier, separate residue of the same
    /// chain already uses, as in PDB files whose numbering wraps around;
    /// `repeat` counts the earlier residues with this number.
    RepeatedAuthor {
        sequence: String,
        repeat: usize,
    },
    /// Zero-based source occurrence when neither sequence namespace is available.
    Occurrence {
        component_id: String,
        index: usize,
    },
}

impl SourceResidueId {
    pub(super) fn of(row: &AtomRow) -> Self {
        Self::new(
            &row.model_id,
            &row.asym_id,
            row.label_seq_id,
            row.auth_seq_id.as_deref(),
            row.insertion_code.as_deref(),
            row.occurrence,
            &row.comp_id,
        )
    }

    /// The identity of a residue from its source fields. An occurrence counts
    /// repeats of an author number, or residues when no number is available.
    pub(in crate::io) fn new(
        model_id: &str,
        asym_id: &str,
        label_sequence_id: Option<i32>,
        author_sequence_id: Option<&str>,
        insertion_code: Option<&str>,
        occurrence: Option<usize>,
        component_id: &str,
    ) -> Self {
        let position = match (label_sequence_id, author_sequence_id, occurrence) {
            (Some(sequence), _, _) => SourceResiduePosition::Label(sequence),
            (None, Some(sequence), None) => SourceResiduePosition::Author(sequence.to_owned()),
            (None, Some(sequence), Some(repeat)) => SourceResiduePosition::RepeatedAuthor {
                sequence: sequence.to_owned(),
                repeat,
            },
            (None, None, occurrence) => SourceResiduePosition::Occurrence {
                component_id: component_id.to_owned(),
                index: occurrence.unwrap_or_default(),
            },
        };
        Self {
            model_id: model_id.to_owned(),
            asym_id: asym_id.to_owned(),
            position,
            insertion_code: insertion_code.map(str::to_owned),
        }
    }
}

/// Source alternatives available at a residue before any coordinate selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AltLocResidue {
    pub(crate) residue: SourceResidueId,
    pub(crate) labels: Vec<String>,
}

impl AltLocResidue {
    pub fn residue(&self) -> &SourceResidueId {
        &self.residue
    }

    /// Labels in lexical order. Shared unlabelled atoms are not a separate choice.
    pub fn labels(&self) -> &[String] {
        &self.labels
    }
}

/// Why an alternate configuration was selected for a residue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AltLocSelectionReason {
    HighestOccupancy,
    ExactSelection,
    PreferredLabel,
    PreferredLabelUnavailable,
    /// The source declared a combination containing the selected labels.
    SourceConformation {
        id: String,
    },
}

/// Source correspondence for one residue's alternate-location decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AltLocDecision {
    pub(crate) residue: SourceResidueId,
    pub(crate) available_labels: Vec<String>,
    pub(crate) selected_labels: Vec<String>,
    pub(crate) reason: AltLocSelectionReason,
    pub(crate) selected_source_lines: Vec<usize>,
    pub(crate) omitted_source_lines: Vec<usize>,
}

impl AltLocDecision {
    pub fn residue(&self) -> &SourceResidueId {
        &self.residue
    }

    pub fn available_labels(&self) -> &[String] {
        &self.available_labels
    }

    /// Usually one label; source-declared combinations may include several
    /// non-overlapping labelled subsets of one residue.
    pub fn selected_labels(&self) -> &[String] {
        &self.selected_labels
    }

    pub fn reason(&self) -> &AltLocSelectionReason {
        &self.reason
    }

    /// One-based source lines, including retained unlabelled shared atoms.
    pub fn selected_source_lines(&self) -> &[usize] {
        &self.selected_source_lines
    }

    pub fn omitted_source_lines(&self) -> &[usize] {
        &self.omitted_source_lines
    }
}
