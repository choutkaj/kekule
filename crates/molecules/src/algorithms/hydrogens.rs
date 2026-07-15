use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::core::{
    Atom, AtomId, BondId, BondOrder, Element, Molecule, MoleculeError, StereoCarrier,
    StereoElementId, StereoElementKind,
};

use super::{perceive_valence_with_options, ValenceModel, ValenceOptions};

const DEFAULT_MAX_ADDED_HYDROGENS: usize = 1_000_000;

/// Controls which encoded hydrogens are materialized and bounds graph growth.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AddHydrogensOptions {
    /// Materialize only bracket-style hydrogen counts and leave perceived
    /// implicit hydrogens implicit.
    pub explicit_only: bool,
    /// Bound topology growth before any mutation is committed.
    pub max_added_hydrogens: usize,
}

impl Default for AddHydrogensOptions {
    fn default() -> Self {
        Self {
            explicit_only: false,
            max_added_hydrogens: DEFAULT_MAX_ADDED_HYDROGENS,
        }
    }
}

/// Identifies the count representation consumed for an added hydrogen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AddedHydrogenOrigin {
    ExplicitCount,
    Implicit,
}

/// Stable-ID mapping for one materialized hydrogen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AddedHydrogen {
    pub hydrogen: AtomId,
    pub parent: AtomId,
    pub origin: AddedHydrogenOrigin,
}

/// Complete mapping produced by a successful addition.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AddHydrogensReport {
    pub added: Vec<AddedHydrogen>,
}

/// Stable-ID mapping for one collapsed hydrogen.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RemovedHydrogen {
    pub hydrogen: AtomId,
    pub parent: AtomId,
}

/// Explains why a graph hydrogen could not be collapsed losslessly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RetainedHydrogenReason {
    Isotopic,
    Mapped,
    Charged,
    Radical,
    EncodedHydrogenCount,
    AtomProperties,
    NoHeavyAtomNeighbor,
    MultipleNeighbors,
    HydrogenNeighbor,
    NonSingleBond,
    BondProperties,
    StereoBondMark,
    UnsupportedStereoRole,
}

/// A hydrogen deliberately retained by conservative removal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RetainedHydrogen {
    pub hydrogen: AtomId,
    pub reason: RetainedHydrogenReason,
}

/// Encoded count chosen for a parent after graph-hydrogen collapse.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HydrogenCountAdjustment {
    pub parent: AtomId,
    pub removed_graph_hydrogens: usize,
    pub explicit_hydrogens: u8,
    pub implicit_hydrogens: u8,
}

/// Complete topology mapping and count reconstruction from removal.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RemoveHydrogensReport {
    pub removed: Vec<RemovedHydrogen>,
    pub retained: Vec<RetainedHydrogen>,
    pub adjustments: Vec<HydrogenCountAdjustment>,
}

/// Failure from planning or committing a transactional hydrogen transform.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum HydrogenNormalizationError {
    MissingValencePerception,
    ResourceLimit {
        requested_hydrogens: usize,
        limit: usize,
    },
    HydrogenCountOverflow {
        atom: AtomId,
        count: usize,
    },
    InconsistentStereoHydrogen {
        element: StereoElementId,
        atom: AtomId,
    },
    UnsupportedImplicitAxisHydrogen {
        element: StereoElementId,
    },
    HydrogenCountNotPreserved {
        atom: AtomId,
        expected: usize,
        actual: usize,
    },
    Molecule(MoleculeError),
}

impl fmt::Display for HydrogenNormalizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingValencePerception => write!(
                f,
                "hydrogen normalization requires current valence perception"
            ),
            Self::ResourceLimit {
                requested_hydrogens,
                limit,
            } => write!(
                f,
                "adding {requested_hydrogens} hydrogens exceeds the configured limit of {limit}"
            ),
            Self::HydrogenCountOverflow { atom, count } => write!(
                f,
                "atom {atom} requires {count} encoded hydrogens, which exceeds the supported count"
            ),
            Self::InconsistentStereoHydrogen { element, atom } => write!(
                f,
                "stereo element {element} has an inconsistent implicit hydrogen at atom {atom}"
            ),
            Self::UnsupportedImplicitAxisHydrogen { element } => write!(
                f,
                "axis stereo element {element} cannot materialize an implicit hydrogen carrier"
            ),
            Self::HydrogenCountNotPreserved {
                atom,
                expected,
                actual,
            } => write!(
                f,
                "hydrogen normalization changed the total hydrogen count at atom {atom}: expected {expected}, got {actual}"
            ),
            Self::Molecule(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for HydrogenNormalizationError {}

impl From<MoleculeError> for HydrogenNormalizationError {
    fn from(error: MoleculeError) -> Self {
        Self::Molecule(error)
    }
}

pub(crate) fn add_hydrogens_to_molecule(
    molecule: &mut Molecule,
    options: AddHydrogensOptions,
) -> Result<AddHydrogensReport, HydrogenNormalizationError> {
    if !options.explicit_only && !molecule.perception().has_valence() {
        return Err(HydrogenNormalizationError::MissingValencePerception);
    }

    let plan = molecule
        .atoms()
        .map(|(parent, atom)| {
            let explicit = usize::from(atom.explicit_hydrogens);
            let implicit = if options.explicit_only {
                0
            } else {
                usize::from(molecule.implicit_hydrogens(parent)?.unwrap_or(0))
            };
            Ok((parent, explicit, implicit))
        })
        .collect::<Result<Vec<_>, MoleculeError>>()?;
    let requested_hydrogens = plan
        .iter()
        .map(|(_, explicit, implicit)| explicit + implicit)
        .sum::<usize>();
    if requested_hydrogens > options.max_added_hydrogens {
        return Err(HydrogenNormalizationError::ResourceLimit {
            requested_hydrogens,
            limit: options.max_added_hydrogens,
        });
    }

    validate_materialized_stereo_hydrogens(molecule, &plan, options.explicit_only)?;

    let mut staged = molecule.clone();
    let hydrogen = Element::from_atomic_number(1).expect("hydrogen is a periodic-table element");
    let mut report = AddHydrogensReport::default();
    let mut added_by_parent = BTreeMap::<AtomId, Vec<AtomId>>::new();

    for (parent, explicit, implicit) in plan {
        for origin in std::iter::repeat_n(AddedHydrogenOrigin::ExplicitCount, explicit)
            .chain(std::iter::repeat_n(AddedHydrogenOrigin::Implicit, implicit))
        {
            let mut atom = Atom::new(hydrogen);
            atom.no_implicit_hydrogens = true;
            let hydrogen_id = staged.add_atom(atom);
            staged.add_bond(parent, hydrogen_id, BondOrder::Single)?;
            added_by_parent.entry(parent).or_default().push(hydrogen_id);
            report.added.push(AddedHydrogen {
                hydrogen: hydrogen_id,
                parent,
                origin,
            });
        }
        if explicit > 0 {
            staged.atom_mut(parent)?.explicit_hydrogens = 0;
        }
    }

    materialize_stereo_hydrogen_carriers(&mut staged, &added_by_parent)?;
    *molecule = staged;
    Ok(report)
}

pub(crate) fn remove_hydrogens_from_molecule(
    molecule: &mut Molecule,
) -> Result<RemoveHydrogensReport, HydrogenNormalizationError> {
    if !molecule.perception().has_valence() {
        return Err(HydrogenNormalizationError::MissingValencePerception);
    }

    let mut report = RemoveHydrogensReport::default();
    let mut candidates = Vec::<(AtomId, AtomId, BondId)>::new();
    for (hydrogen, atom) in molecule
        .atoms()
        .filter(|(_, atom)| atom.element.atomic_number() == 1)
    {
        match removable_hydrogen(molecule, hydrogen, atom)? {
            Ok((parent, bond)) => candidates.push((hydrogen, parent, bond)),
            Err(reason) => report.retained.push(RetainedHydrogen { hydrogen, reason }),
        }
    }

    let candidate_ids = candidates
        .iter()
        .map(|(hydrogen, _, _)| *hydrogen)
        .collect::<BTreeSet<_>>();
    let mut removable = BTreeSet::new();
    for (hydrogen, parent, bond) in &candidates {
        let stereo_is_collapsible =
            stereo_hydrogen_is_collapsible(molecule, *hydrogen, *parent, &candidate_ids);
        let source_mark_is_safe = molecule.stereo_bond_mark(*bond).is_none()
            || stereo_hydrogen_has_collapsible_role(molecule, *hydrogen, *parent);
        if stereo_is_collapsible && source_mark_is_safe {
            removable.insert(*hydrogen);
        } else {
            report.retained.push(RetainedHydrogen {
                hydrogen: *hydrogen,
                reason: if source_mark_is_safe {
                    RetainedHydrogenReason::UnsupportedStereoRole
                } else {
                    RetainedHydrogenReason::StereoBondMark
                },
            });
        }
    }

    let mut by_parent = BTreeMap::<AtomId, Vec<AtomId>>::new();
    for (hydrogen, parent, _) in candidates {
        if removable.contains(&hydrogen) {
            by_parent.entry(parent).or_default().push(hydrogen);
        }
    }

    let mut expected_totals = BTreeMap::<AtomId, usize>::new();
    for (parent, hydrogens) in &by_parent {
        let atom = molecule.atom(*parent)?;
        let implicit = usize::from(molecule.implicit_hydrogens(*parent)?.unwrap_or(0));
        expected_totals.insert(
            *parent,
            usize::from(atom.explicit_hydrogens) + implicit + hydrogens.len(),
        );
    }

    let stereo_count_parents = stereo_hydrogen_parents(molecule, &removable);
    let mut staged = molecule.clone();
    collapse_stereo_hydrogen_carriers(&mut staged, &removable)?;
    for (parent, hydrogens) in &by_parent {
        for hydrogen in hydrogens {
            staged.delete_atom(*hydrogen)?;
            report.removed.push(RemovedHydrogen {
                hydrogen: *hydrogen,
                parent: *parent,
            });
        }
    }

    adjust_collapsed_hydrogen_counts(
        &mut staged,
        &expected_totals,
        &by_parent,
        &stereo_count_parents,
    )?;
    report.adjustments = verify_collapsed_hydrogen_counts(&staged, &expected_totals, &by_parent)?;
    report.removed.sort_by_key(|entry| entry.hydrogen);
    report.retained.sort_by_key(|entry| entry.hydrogen);
    report.adjustments.sort_by_key(|entry| entry.parent);
    *molecule = staged;
    Ok(report)
}

fn validate_materialized_stereo_hydrogens(
    molecule: &Molecule,
    plan: &[(AtomId, usize, usize)],
    explicit_only: bool,
) -> Result<(), HydrogenNormalizationError> {
    let totals = plan
        .iter()
        .map(|(atom, explicit, implicit)| (*atom, explicit + implicit))
        .collect::<BTreeMap<_, _>>();
    for (element_id, element) in molecule.stereo_elements() {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                let implicit = stereo
                    .carriers
                    .iter()
                    .filter(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen))
                    .count();
                validate_implicit_stereo_count(
                    element_id,
                    stereo.center,
                    implicit,
                    totals.get(&stereo.center).copied().unwrap_or(0),
                    explicit_only,
                )?;
            }
            StereoElementKind::DoubleBond(stereo) => {
                for (parent, carrier) in [
                    (stereo.left, stereo.left_carrier),
                    (stereo.right, stereo.right_carrier),
                ] {
                    let implicit = usize::from(matches!(carrier, StereoCarrier::ImplicitHydrogen));
                    validate_implicit_stereo_count(
                        element_id,
                        parent,
                        implicit,
                        totals.get(&parent).copied().unwrap_or(0),
                        explicit_only,
                    )?;
                }
            }
            StereoElementKind::Axis(stereo) => {
                if stereo
                    .carriers
                    .iter()
                    .any(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen))
                {
                    return Err(
                        HydrogenNormalizationError::UnsupportedImplicitAxisHydrogen {
                            element: element_id,
                        },
                    );
                }
            }
        }
    }
    Ok(())
}

fn validate_implicit_stereo_count(
    element: StereoElementId,
    atom: AtomId,
    implicit_carriers: usize,
    added_hydrogens: usize,
    explicit_only: bool,
) -> Result<(), HydrogenNormalizationError> {
    if implicit_carriers > 1
        || (implicit_carriers == 1
            && ((!explicit_only && added_hydrogens != 1) || (explicit_only && added_hydrogens > 1)))
    {
        return Err(HydrogenNormalizationError::InconsistentStereoHydrogen { element, atom });
    }
    Ok(())
}

fn materialize_stereo_hydrogen_carriers(
    molecule: &mut Molecule,
    added_by_parent: &BTreeMap<AtomId, Vec<AtomId>>,
) -> Result<(), HydrogenNormalizationError> {
    let replacements = molecule
        .stereo_elements()
        .filter_map(|(id, element)| {
            let mut replacement = element.clone();
            let mut changed = false;
            match &mut replacement.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    if let Some(hydrogen) = added_by_parent
                        .get(&stereo.center)
                        .and_then(|hydrogens| hydrogens.first())
                    {
                        for carrier in &mut stereo.carriers {
                            if matches!(carrier, StereoCarrier::ImplicitHydrogen) {
                                *carrier = StereoCarrier::Atom(*hydrogen);
                                changed = true;
                            }
                        }
                    }
                }
                StereoElementKind::DoubleBond(stereo) => {
                    if matches!(stereo.left_carrier, StereoCarrier::ImplicitHydrogen) {
                        if let Some(hydrogen) = added_by_parent
                            .get(&stereo.left)
                            .and_then(|hydrogens| hydrogens.first())
                        {
                            stereo.left_carrier = StereoCarrier::Atom(*hydrogen);
                            changed = true;
                        }
                    }
                    if matches!(stereo.right_carrier, StereoCarrier::ImplicitHydrogen) {
                        if let Some(hydrogen) = added_by_parent
                            .get(&stereo.right)
                            .and_then(|hydrogens| hydrogens.first())
                        {
                            stereo.right_carrier = StereoCarrier::Atom(*hydrogen);
                            changed = true;
                        }
                    }
                }
                StereoElementKind::Axis(_) => {}
            }
            changed.then_some((id, replacement))
        })
        .collect::<Vec<_>>();
    for (id, replacement) in replacements {
        molecule.replace_stereo_element(id, replacement)?;
    }
    Ok(())
}

fn removable_hydrogen(
    molecule: &Molecule,
    hydrogen: AtomId,
    atom: &Atom,
) -> Result<Result<(AtomId, BondId), RetainedHydrogenReason>, MoleculeError> {
    if atom.isotope.is_some() {
        return Ok(Err(RetainedHydrogenReason::Isotopic));
    }
    if atom.atom_map.is_some() {
        return Ok(Err(RetainedHydrogenReason::Mapped));
    }
    if atom.formal_charge != 0 {
        return Ok(Err(RetainedHydrogenReason::Charged));
    }
    if atom.radical.is_some() {
        return Ok(Err(RetainedHydrogenReason::Radical));
    }
    if atom.explicit_hydrogens != 0 {
        return Ok(Err(RetainedHydrogenReason::EncodedHydrogenCount));
    }
    if !atom.props.is_empty() {
        return Ok(Err(RetainedHydrogenReason::AtomProperties));
    }
    let incident = molecule.incident_bonds(hydrogen)?.collect::<Vec<_>>();
    if incident.is_empty() {
        return Ok(Err(RetainedHydrogenReason::NoHeavyAtomNeighbor));
    }
    if incident.len() != 1 {
        return Ok(Err(RetainedHydrogenReason::MultipleNeighbors));
    }
    let (bond_id, bond) = incident[0];
    let parent = bond.other_atom(hydrogen);
    if molecule.atom(parent)?.element.atomic_number() == 1 {
        return Ok(Err(RetainedHydrogenReason::HydrogenNeighbor));
    }
    if bond.order != BondOrder::Single {
        return Ok(Err(RetainedHydrogenReason::NonSingleBond));
    }
    if !bond.props.is_empty() {
        return Ok(Err(RetainedHydrogenReason::BondProperties));
    }
    Ok(Ok((parent, bond_id)))
}

fn stereo_hydrogen_has_collapsible_role(
    molecule: &Molecule,
    hydrogen: AtomId,
    parent: AtomId,
) -> bool {
    molecule
        .stereo_elements()
        .any(|(_, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                stereo.center == parent && stereo.carriers.contains(&StereoCarrier::Atom(hydrogen))
            }
            StereoElementKind::DoubleBond(stereo) => {
                (stereo.left == parent && stereo.left_carrier == StereoCarrier::Atom(hydrogen))
                    || (stereo.right == parent
                        && stereo.right_carrier == StereoCarrier::Atom(hydrogen))
            }
            StereoElementKind::Axis(_) => false,
        })
}

fn stereo_hydrogen_is_collapsible(
    molecule: &Molecule,
    hydrogen: AtomId,
    parent: AtomId,
    candidates: &BTreeSet<AtomId>,
) -> bool {
    molecule.stereo_elements().all(|(_, element)| match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            if stereo.center == hydrogen {
                return false;
            }
            let referenced = stereo
                .carriers
                .iter()
                .filter(|carrier| matches!(carrier, StereoCarrier::Atom(atom) if *atom == hydrogen))
                .count();
            referenced == 0
                || (referenced == 1
                    && stereo.center == parent
                    && !stereo
                        .carriers
                        .iter()
                        .any(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen))
                    && stereo
                        .carriers
                        .iter()
                        .filter(|carrier| {
                            matches!(carrier, StereoCarrier::Atom(atom) if candidates.contains(atom))
                        })
                        .count()
                        == 1)
        }
        StereoElementKind::DoubleBond(stereo) => {
            if stereo.left == hydrogen || stereo.right == hydrogen {
                return false;
            }
            let left = matches!(stereo.left_carrier, StereoCarrier::Atom(atom) if atom == hydrogen);
            let right = matches!(stereo.right_carrier, StereoCarrier::Atom(atom) if atom == hydrogen);
            match (left, right) {
                (false, false) => true,
                (true, false) => stereo.left == parent,
                (false, true) => stereo.right == parent,
                (true, true) => false,
            }
        }
        StereoElementKind::Axis(stereo) => !stereo
            .carriers
            .iter()
            .any(|carrier| matches!(carrier, StereoCarrier::Atom(atom) if *atom == hydrogen)),
    })
}

fn collapse_stereo_hydrogen_carriers(
    molecule: &mut Molecule,
    removable: &BTreeSet<AtomId>,
) -> Result<(), HydrogenNormalizationError> {
    let replacements = molecule
        .stereo_elements()
        .filter_map(|(id, element)| {
            let mut replacement = element.clone();
            let mut changed = false;
            match &mut replacement.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    for carrier in &mut stereo.carriers {
                        if matches!(carrier, StereoCarrier::Atom(atom) if removable.contains(atom)) {
                            *carrier = StereoCarrier::ImplicitHydrogen;
                            changed = true;
                        }
                    }
                }
                StereoElementKind::DoubleBond(stereo) => {
                    if matches!(stereo.left_carrier, StereoCarrier::Atom(atom) if removable.contains(&atom)) {
                        stereo.left_carrier = StereoCarrier::ImplicitHydrogen;
                        changed = true;
                    }
                    if matches!(stereo.right_carrier, StereoCarrier::Atom(atom) if removable.contains(&atom)) {
                        stereo.right_carrier = StereoCarrier::ImplicitHydrogen;
                        changed = true;
                    }
                }
                StereoElementKind::Axis(_) => {}
            }
            changed.then_some((id, replacement))
        })
        .collect::<Vec<_>>();
    for (id, replacement) in replacements {
        molecule.replace_stereo_element(id, replacement)?;
    }
    Ok(())
}

fn stereo_hydrogen_parents(molecule: &Molecule, removable: &BTreeSet<AtomId>) -> BTreeSet<AtomId> {
    let mut parents = BTreeSet::new();
    for (_, element) in molecule.stereo_elements() {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                if stereo.carriers.iter().any(
                    |carrier| matches!(carrier, StereoCarrier::Atom(atom) if removable.contains(atom)),
                ) {
                    parents.insert(stereo.center);
                }
            }
            StereoElementKind::DoubleBond(stereo) => {
                if matches!(stereo.left_carrier, StereoCarrier::Atom(atom) if removable.contains(&atom))
                {
                    parents.insert(stereo.left);
                }
                if matches!(stereo.right_carrier, StereoCarrier::Atom(atom) if removable.contains(&atom))
                {
                    parents.insert(stereo.right);
                }
            }
            StereoElementKind::Axis(_) => {}
        }
    }
    parents
}

fn adjust_collapsed_hydrogen_counts(
    molecule: &mut Molecule,
    expected_totals: &BTreeMap<AtomId, usize>,
    by_parent: &BTreeMap<AtomId, Vec<AtomId>>,
    stereo_count_parents: &BTreeSet<AtomId>,
) -> Result<(), HydrogenNormalizationError> {
    let mut probe = molecule.clone();
    let _ = perceive_valence_with_options(
        &mut probe,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
    for (parent, expected) in expected_totals {
        if stereo_count_parents.contains(parent) {
            let adjusted = u8::try_from(*expected).map_err(|_| {
                HydrogenNormalizationError::HydrogenCountOverflow {
                    atom: *parent,
                    count: *expected,
                }
            })?;
            molecule.atom_mut(*parent)?.explicit_hydrogens = adjusted;
            continue;
        }
        let explicit = usize::from(molecule.atom(*parent)?.explicit_hydrogens);
        let implicit = usize::from(probe.implicit_hydrogens(*parent)?.unwrap_or(0));
        let actual = explicit + implicit;
        if actual < *expected {
            let adjusted = explicit + (*expected - actual);
            let adjusted = u8::try_from(adjusted).map_err(|_| {
                HydrogenNormalizationError::HydrogenCountOverflow {
                    atom: *parent,
                    count: adjusted,
                }
            })?;
            molecule.atom_mut(*parent)?.explicit_hydrogens = adjusted;
        } else if actual > *expected {
            let adjusted = u8::try_from(*expected).map_err(|_| {
                HydrogenNormalizationError::HydrogenCountOverflow {
                    atom: *parent,
                    count: *expected,
                }
            })?;
            let mut atom = molecule.atom_mut(*parent)?;
            atom.explicit_hydrogens = adjusted;
            atom.no_implicit_hydrogens = true;
        }
        debug_assert!(!by_parent[parent].is_empty());
    }
    Ok(())
}

fn verify_collapsed_hydrogen_counts(
    molecule: &Molecule,
    expected_totals: &BTreeMap<AtomId, usize>,
    by_parent: &BTreeMap<AtomId, Vec<AtomId>>,
) -> Result<Vec<HydrogenCountAdjustment>, HydrogenNormalizationError> {
    let mut probe = molecule.clone();
    let _ = perceive_valence_with_options(
        &mut probe,
        ValenceModel::RdkitLike,
        ValenceOptions { strict: false },
    );
    let mut adjustments = Vec::new();
    for (parent, expected) in expected_totals {
        let explicit = molecule.atom(*parent)?.explicit_hydrogens;
        let implicit = probe.implicit_hydrogens(*parent)?.unwrap_or(0);
        let actual = usize::from(explicit) + usize::from(implicit);
        if actual != *expected {
            return Err(HydrogenNormalizationError::HydrogenCountNotPreserved {
                atom: *parent,
                expected: *expected,
                actual,
            });
        }
        adjustments.push(HydrogenCountAdjustment {
            parent: *parent,
            removed_graph_hydrogens: by_parent[parent].len(),
            explicit_hydrogens: explicit,
            implicit_hydrogens: implicit,
        });
    }
    Ok(adjustments)
}
