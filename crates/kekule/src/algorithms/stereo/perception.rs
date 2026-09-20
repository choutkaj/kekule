use super::*;
use std::borrow::Cow;
use std::result::Result;

/// Bounds for exact graph-symmetry searches used by stereo perception.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StereoPerceptionOptions {
    /// Maximum analysis vertices, including non-graph hydrogens.
    pub max_atoms: usize,
    /// Maximum trial atom mappings across all sites and refinement rounds.
    pub max_search_states: usize,
}

impl Default for StereoPerceptionOptions {
    fn default() -> Self {
        Self {
            max_atoms: 4096,
            max_search_states: 1_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum StereoPerceptionError {
    InvalidOptions(&'static str),
    Perception(crate::chemistry::PerceptionError),
    InvalidStereo(StereoValidationError),
    ResourceLimit {
        resource: &'static str,
        observed: usize,
        limit: usize,
    },
}

impl fmt::Display for StereoPerceptionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => {
                write!(f, "invalid stereo perception options: {message}")
            }
            Self::Perception(error) => error.fmt(f),
            Self::InvalidStereo(error) => error.fmt(f),
            Self::ResourceLimit {
                resource,
                observed,
                limit,
            } => write!(
                f,
                "stereo {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
        }
    }
}

impl std::error::Error for StereoPerceptionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Perception(error) => Some(error),
            Self::InvalidStereo(error) => Some(error),
            _ => None,
        }
    }
}

pub(super) fn prepared(mol: &Molecule) -> Result<Cow<'_, Molecule>, StereoPerceptionError> {
    if mol.perception().has_valence()
        && mol.perception().has_aromaticity()
        && mol.ring_set().is_some()
    {
        Ok(Cow::Borrowed(mol))
    } else {
        let mut perceived = mol.clone();
        perceived
            .perceive()
            .map_err(StereoPerceptionError::Perception)?;
        Ok(Cow::Owned(perceived))
    }
}

/// Finds stereo candidates with explicit search bounds.
///
/// Complete installed valence/ring/aromaticity perception is reused; otherwise the
/// default profile is computed on a temporary copy. Graph hydrogens and known
/// non-graph hydrogens are equivalent in the analysis. Atom maps are ignored.
///
/// Neutral, closed-shell three-coordinate nitrogen requires three single-bond
/// ligand directions, no conjugated attachment, and either a three-membered
/// ring or the RDKit-like bridgehead criterion on the selected ring set. This
/// models inversion eligibility, not kinetic stability. A hydrogen ligand is
/// treated identically whether declared, inferred, or represented as an atom.
/// Three-coordinate charged/radical nitrogen remains outside the lone-pair model.
///
/// A site is excluded only after proving an orientation-reversing automorphism
/// that preserves other stereo constraints. Unspecified sites are held fixed,
/// including their local orientation, so dependent ring stereo is retained.
/// Removed sites cease to constrain subsequent rounds. Specified absolute
/// centers can be exchanged only with matching configuration; this distinguishes
/// homomorphic and enantiomorphic ligands. Enhanced-group and axial assertions
/// are conservatively held fixed pending their complete dependency analysis.
/// Exhaustion is an error, never evidence of equivalence or stereogenicity.
pub fn detect_stereo_candidates_with_options(
    mol: &Molecule,
    options: StereoPerceptionOptions,
) -> Result<Vec<StereoCandidate>, StereoPerceptionError> {
    if options.max_atoms == 0 || options.max_search_states == 0 {
        return Err(StereoPerceptionError::InvalidOptions(
            "limits must be positive",
        ));
    }
    if mol.atom_count() > options.max_atoms {
        return Err(StereoPerceptionError::ResourceLimit {
            resource: "analysis atoms",
            observed: mol.atom_count(),
            limit: options.max_atoms,
        });
    }
    validate_stereo(mol).map_err(StereoPerceptionError::InvalidStereo)?;
    let perceived = prepared(mol)?;
    let mol = perceived.as_ref();
    let mut candidates = tetrahedral_candidates(mol);
    candidates.extend(double_bond_candidates(mol));
    super::symmetry::filter(mol, candidates, options)
}

/// Stereo assertions removed by explicit chemical cleanup.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StereoCleanupReport {
    pub removed_elements: Vec<StereoElementId>,
    /// Assertions whose local geometry is outside this cleanup model. These
    /// are preserved, not classified as chemically invalid from their absence
    /// in the candidate set.
    pub unclassified_elements: Vec<StereoElementId>,
}

/// Removes chemically ineligible or symmetry-degenerate tetrahedral/double-bond
/// assertions in an explicit, transactional editor operation.
///
/// Surviving assertions retain their carrier ordering and configuration. Group
/// membership is pruned by the normal checked editor operations. Unsupported
/// local geometries and axial assertions are reported as unclassified and
/// preserved. Default perception and parsing never invoke this operation.
/// Failure leaves the editor, its properties and installed perception untouched.
pub fn cleanup_stereo(
    editor: &mut MoleculeEditor,
    options: StereoPerceptionOptions,
) -> Result<StereoCleanupReport, StereoPerceptionError> {
    let candidates = detect_stereo_candidates_with_options(editor.working(), options)?;
    let perceived = prepared(editor.working())?;
    let mut unclassified_elements = Vec::new();
    let removed_elements: Vec<_> = editor
        .stereo_elements()
        .filter_map(|(id, element)| {
            let eligible = match &element.kind {
                StereoElementKind::Tetrahedral(stereo)
                    if unclassified_tetrahedral_geometry(&perceived, stereo.center) =>
                {
                    unclassified_elements.push(id);
                    true
                }
                StereoElementKind::Tetrahedral(stereo) => candidates.iter().any(|candidate| {
                    matches!(candidate,
                StereoCandidate::Tetrahedral { center, .. } if *center == stereo.center)
                }),
                StereoElementKind::DoubleBond(stereo) => candidates.iter().any(|candidate| {
                    matches!(candidate,
                StereoCandidate::DoubleBond { bond, .. } if *bond == stereo.bond)
                }),
                StereoElementKind::Axis(_) => {
                    unclassified_elements.push(id);
                    true
                }
            };
            (!eligible).then_some(id)
        })
        .collect();
    drop(perceived);
    if !removed_elements.is_empty() {
        let mut staged = editor.clone();
        for &id in &removed_elements {
            staged
                .remove_stereo_element(id)
                .expect("collected live stereo element");
        }
        *editor = staged;
    }
    Ok(StereoCleanupReport {
        removed_elements,
        unclassified_elements,
    })
}
