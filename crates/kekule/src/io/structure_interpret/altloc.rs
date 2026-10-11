//! Alternate locations are source alternatives, not additional coordinate models.
//!
//! Selection keeps unlabelled shared atoms and chooses whole residues (or explicit
//! groups). Occupancy ranking is a deterministic representative-selection heuristic,
//! not a joint probability or a guarantee of compatibility between unrelated groups.

mod select;
mod source;

use std::collections::BTreeMap;

use super::atom_site::AtomRow;

pub(crate) use select::inventory;
pub(super) use select::select_alt_locations;

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
/// Overrides take precedence over `default`. Every override and group member must
/// identify a residue with alternatives. Unknown identities and conflicting group
/// choices are errors, including those in unselected source coordinate models.
///
/// ```
/// use kekule::mmcif::{MmcifBlock, AltLocPolicy, AltLocPreference,
///     AltLocSelection, StructureInterpretOptions, StructureInterpretation,
///     StructureInterpretError};
/// # fn choose(block: &MmcifBlock) -> Result<StructureInterpretation, StructureInterpretError> {
/// let mut selection = AltLocSelection {
///     default: AltLocPreference::PreferLabel("A".into()),
///     ..Default::default()
/// };
/// // Override one source residue, using its exact block-scoped identity.
/// if let Some(residue) = block.alternate_locations()?.into_iter()
///     .find(|residue| residue.labels.iter().any(|label| label == "B"))
/// {
///     selection.overrides.insert(residue.residue, "B".into());
/// }
/// block.interpret_with_options(StructureInterpretOptions {
///     altloc_policy: AltLocPolicy::Configured(selection),
///     ..Default::default()
/// })
/// # }
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AltLocSelection {
    pub default: AltLocPreference,
    pub overrides: BTreeMap<SourceResidueId, String>,
    /// Residues constrained to use the same label. Groups must contain at least
    /// two distinct residues in one source model; overlapping groups are united.
    pub groups: Vec<Vec<SourceResidueId>>,
    /// Select one `_atom_sites_alt_gen.ens_id` combination explicitly.
    /// This identifier does not create or select an `Ensemble` member.
    pub source_conformation: Option<String>,
}

/// A block-scoped residue identity, independent of selected chemistry and rows.
///
/// `asym_id` uses the label namespace, falling back to author only when absent.
/// Position similarly prefers label sequence, then author sequence, then source
/// occurrence. Insertion codes and source coordinate models remain significant.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SourceResidueId {
    pub model_id: String,
    pub asym_id: String,
    pub position: SourceResiduePosition,
    pub insertion_code: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SourceResiduePosition {
    Label(i32),
    Author(String),
    /// Zero-based source occurrence when neither sequence namespace is available.
    Occurrence {
        component_id: String,
        index: usize,
    },
}

impl SourceResidueId {
    pub(super) fn of(row: &AtomRow) -> Self {
        Self {
            model_id: row.model_id.clone(),
            asym_id: row.asym_id.clone(),
            position: if let Some(sequence) = row.label_seq_id {
                SourceResiduePosition::Label(sequence)
            } else if let Some(sequence) = &row.auth_seq_id {
                SourceResiduePosition::Author(sequence.clone())
            } else {
                SourceResiduePosition::Occurrence {
                    component_id: row.comp_id.clone(),
                    index: row
                        .occurrence
                        .expect("unsequenced residue has an occurrence"),
                }
            },
            insertion_code: row.insertion_code.clone(),
        }
    }
}

/// Source alternatives available at a residue before any coordinate selection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AltLocResidue {
    pub residue: SourceResidueId,
    /// Labels in lexical order. Shared unlabelled atoms are not a separate choice.
    pub labels: Vec<String>,
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
    pub residue: SourceResidueId,
    pub available_labels: Vec<String>,
    /// Usually one label; source-declared combinations may include several
    /// non-overlapping labelled subsets of one residue.
    pub selected_labels: Vec<String>,
    pub reason: AltLocSelectionReason,
    /// One-based source lines, including retained unlabelled shared atoms.
    pub selected_source_lines: Vec<usize>,
    pub omitted_source_lines: Vec<usize>,
}
