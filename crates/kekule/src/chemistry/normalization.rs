use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::algorithms::StereoValidationIssue;
use crate::core::{
    canonicalize_represented_chemistry, Atom, AtomId, BondId, BondOrder, Molecule, MoleculeError,
    StereoElementId,
};

use super::source_stereo::{
    normalize_source_stereo, SourceStereoBondMark, SourceStereoBondMarkKind,
};

const MAX_AROMATIC_LOCALIZATION_MATCHING_STATES: usize = 100_000;

/// Failure to publish Kekule's canonical represented chemistry.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NormalizationError {
    /// A meaning-preserving representation rewrite requires an unsupported charge.
    FormalChargeOutOfRange { atom: AtomId, charge: usize },
    /// Source-aromatic bonds cannot be localized under the fixed rules.
    InvalidAromaticRepresentation(AtomId),
    /// Imported aromatic localization exhausted its deterministic search budget.
    AromaticLocalizationLimit {
        atom: AtomId,
        examined_states: usize,
        limit: usize,
    },
    /// Source-declared stereo could not be represented canonically.
    SourceStereo(SourceStereoNormalizationError),
}

/// Successful sidecar output from representation normalization.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct NormalizationReport {
    /// Canonical stereo elements created from source bond marks.
    pub created_stereo_elements: Vec<StereoElementId>,
    /// Nonfatal source-representation diagnostics.
    pub warnings: Vec<NormalizationWarning>,
}

/// Nonfatal source-representation diagnostic emitted during normalization.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum NormalizationWarning {
    AmbiguousTetrahedralWedgeMarks { center: AtomId, mark_count: usize },
}

/// One concrete source-stereo normalization issue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SourceStereoNormalizationIssue {
    MissingSourceBondMark {
        bond: BondId,
    },
    InvalidSourceBondMarkEndpoint {
        bond: BondId,
        from: AtomId,
    },
    UnassembledTetrahedralBondMark {
        bond: BondId,
        kind: SourceStereoBondMarkKind,
    },
    AmbiguousDirectionalBondMarks {
        double_bond: BondId,
        endpoint: AtomId,
        mark_count: usize,
    },
    UnpairedDirectionalBondMark {
        bond: BondId,
    },
    UnsupportedSourceBondMark {
        bond: BondId,
        kind: SourceStereoBondMarkKind,
    },
    InvalidStereo(StereoValidationIssue),
    CouldNotCreateStereoElement(MoleculeError),
}

/// Collected failure to canonicalize source-declared stereochemistry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceStereoNormalizationError {
    pub issues: Vec<SourceStereoNormalizationIssue>,
}

impl fmt::Display for SourceStereoNormalizationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "source-stereo canonicalization reported {} issue(s): {:?}",
            self.issues.len(),
            self.issues
        )
    }
}

impl std::error::Error for SourceStereoNormalizationError {}

impl fmt::Display for NormalizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "normalizing atom {atom} requires formal charge +{charge}, which is outside the supported range"
            ),
            Self::InvalidAromaticRepresentation(atom) => {
                write!(f, "invalid imported aromatic representation at atom {atom}")
            }
            Self::AromaticLocalizationLimit {
                atom,
                examined_states,
                limit,
            } => write!(
                f,
                "imported aromatic localization limit exceeded at atom {atom}: examined {examined_states} matching states, limit {limit}"
            ),
            Self::SourceStereo(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for NormalizationError {}

impl NormalizationError {
    pub(crate) fn atom_location_hint(&self) -> Option<AtomId> {
        match self {
            Self::FormalChargeOutOfRange { atom, .. }
            | Self::InvalidAromaticRepresentation(atom)
            | Self::AromaticLocalizationLimit { atom, .. } => Some(*atom),
            Self::SourceStereo(error) => error.atom_location_hint(),
        }
    }

    pub(crate) fn bond_location_hint(&self) -> Option<BondId> {
        match self {
            Self::SourceStereo(error) => error.bond_location_hint(),
            _ => None,
        }
    }
}

impl SourceStereoNormalizationError {
    fn atom_location_hint(&self) -> Option<AtomId> {
        self.issues.iter().find_map(|issue| match issue {
            SourceStereoNormalizationIssue::InvalidSourceBondMarkEndpoint { from, .. } => {
                Some(*from)
            }
            SourceStereoNormalizationIssue::AmbiguousDirectionalBondMarks { endpoint, .. } => {
                Some(*endpoint)
            }
            SourceStereoNormalizationIssue::InvalidStereo(
                StereoValidationIssue::MissingStereoAtom { atom, .. }
                | StereoValidationIssue::InvalidTetrahedralCarrierCount { center: atom, .. }
                | StereoValidationIssue::DuplicateTetrahedralCarrier { center: atom, .. }
                | StereoValidationIssue::TetrahedralCarrierNotAdjacent { center: atom, .. }
                | StereoValidationIssue::DoubleBondCarrierIsFocusAtom { endpoint: atom, .. }
                | StereoValidationIssue::DoubleBondCarrierNotAdjacent { endpoint: atom, .. }
                | StereoValidationIssue::UnsupportedDoubleBondCarrier { endpoint: atom, .. },
            ) => Some(*atom),
            _ => None,
        })
    }

    fn bond_location_hint(&self) -> Option<BondId> {
        self.issues.iter().find_map(|issue| match issue {
            SourceStereoNormalizationIssue::MissingSourceBondMark { bond }
            | SourceStereoNormalizationIssue::InvalidSourceBondMarkEndpoint { bond, .. }
            | SourceStereoNormalizationIssue::UnassembledTetrahedralBondMark { bond, .. }
            | SourceStereoNormalizationIssue::UnpairedDirectionalBondMark { bond }
            | SourceStereoNormalizationIssue::UnsupportedSourceBondMark { bond, .. } => Some(*bond),
            SourceStereoNormalizationIssue::AmbiguousDirectionalBondMarks {
                double_bond, ..
            } => Some(*double_bond),
            SourceStereoNormalizationIssue::InvalidStereo(
                StereoValidationIssue::MissingStereoBond { bond, .. }
                | StereoValidationIssue::InvalidDoubleBondOrder { bond, .. }
                | StereoValidationIssue::DoubleBondFocusMismatch { bond, .. }
                | StereoValidationIssue::InvalidAxisCarrierCount { axis: bond, .. }
                | StereoValidationIssue::AxisCarrierIsFocusAtom { axis: bond, .. }
                | StereoValidationIssue::AxisCarrierNotAdjacent { axis: bond, .. }
                | StereoValidationIssue::UnsupportedAxisCarrier { axis: bond, .. },
            ) => Some(*bond),
            _ => None,
        })
    }
}

/// Canonicalize an interpreter-owned staging graph before publication.
///
/// Callers must discard the staging graph on failure. Successful
/// canonicalization clears all derived perception state so interpretation
/// publishes represented chemistry only.
pub(crate) fn canonicalize_molecule_for_publication(
    molecule: &mut Molecule,
    geometry: Option<&dyn crate::chemistry::AtomPositionSource>,
    source_stereo: &[SourceStereoBondMark],
) -> Result<NormalizationReport, NormalizationError> {
    canonicalize_represented_chemistry(molecule).map_err(|error| {
        NormalizationError::FormalChargeOutOfRange {
            atom: error.atom,
            charge: error.charge,
        }
    })?;
    // Source-stereo normalization must not observe arbitrary installed
    // perception. Representation rewrites above already invalidate it
    // conceptually, so clear it before decoding any source marks.
    molecule.clear_perception();
    let report = normalize_source_stereo(molecule, geometry, source_stereo)
        .map_err(NormalizationError::SourceStereo)?;
    // Adding represented stereo invalidates only stereo-derived state. The
    // publication contract clears the complete perception state.
    molecule.clear_perception();
    Ok(report)
}

/// Localize source-aromatic bonds before an interpreter publishes a molecule.
///
/// Interpreters stage each source-aromatic edge as a single bond and identify
/// it here by ID. This keeps source syntax outside [`BondOrder`] while reusing
/// the deterministic representation-localization rules.
pub(crate) fn localize_source_aromatic_bonds(
    molecule: &mut Molecule,
    aromatic_bonds: &BTreeSet<BondId>,
) -> Result<(), NormalizationError> {
    for component in source_aromatic_bond_components(molecule, aromatic_bonds) {
        if !try_localize_aromatic_component_with_limit(
            molecule,
            &component,
            aromatic_bonds,
            MAX_AROMATIC_LOCALIZATION_MATCHING_STATES,
        )? {
            return Err(NormalizationError::InvalidAromaticRepresentation(
                component[0],
            ));
        }
    }
    debug_assert!(aromatic_bonds.iter().all(|bond_id| molecule
        .bond(*bond_id)
        .is_ok_and(|bond| matches!(bond.order, BondOrder::Single | BondOrder::Double))));
    Ok(())
}

/// Localize one imported aromatic-bond component using only represented state.
///
/// Aromatic source bonds contribute one unit of baseline valence. The fixed
/// rules below reserve source-form implicit hydrogens before matching atoms
/// that require one additional bond-order unit. They intentionally do not
/// consult installed perception or any selectable chemical model.
fn try_localize_aromatic_component_with_limit(
    molecule: &mut Molecule,
    component: &[AtomId],
    aromatic_bonds: &BTreeSet<BondId>,
    max_matching_states: usize,
) -> Result<bool, NormalizationError> {
    let component_atoms = component.iter().copied().collect::<BTreeSet<_>>();
    let mut demand = BTreeSet::new();
    for atom_id in component {
        let Ok(atom) = molecule.atom(*atom_id) else {
            return Ok(false);
        };
        let baseline_bond_valence = molecule
            .incident_bonds(*atom_id)
            .ok()
            .into_iter()
            .flatten()
            .map(|(_, bond)| represented_bond_valence(bond.order))
            .sum::<usize>();
        let Some(target_valence) =
            aromatic_localization_target_valence(atom, baseline_bond_valence)
        else {
            return Ok(false);
        };
        let explicit_valence =
            baseline_bond_valence.saturating_add(usize::from(atom.hydrogens.explicit_count()));
        let implicit_hydrogens = source_aromatic_implicit_hydrogens(atom, explicit_valence);
        let occupied_valence = explicit_valence
            .saturating_add(implicit_hydrogens)
            .saturating_add(usize::from(
                atom.radical
                    .map_or(0, |radical| radical.unpaired_electron_count()),
            ));
        match target_valence.checked_sub(occupied_valence) {
            Some(0) => {}
            Some(1) => {
                demand.insert(*atom_id);
            }
            _ => return Ok(false),
        }
    }
    if demand.len() % 2 != 0 {
        return Ok(false);
    }

    let mut adjacency = BTreeMap::<AtomId, Vec<(AtomId, BondId)>>::new();
    for (bond_id, bond) in molecule.bonds().filter(|(bond_id, bond)| {
        aromatic_bonds.contains(bond_id)
            && component_atoms.contains(&bond.a())
            && component_atoms.contains(&bond.b())
    }) {
        if demand.contains(&bond.a()) && demand.contains(&bond.b()) {
            adjacency
                .entry(bond.a())
                .or_default()
                .push((bond.b(), bond_id));
            adjacency
                .entry(bond.b())
                .or_default()
                .push((bond.a(), bond_id));
        }
    }
    for neighbors in adjacency.values_mut() {
        neighbors.sort_unstable();
    }

    let mut stack = vec![(demand, Vec::<BondId>::new())];
    let mut examined_states = 0usize;
    let selected_double_bonds = loop {
        let Some((unmatched, selected)) = stack.pop() else {
            return Ok(false);
        };
        examined_states = examined_states.saturating_add(1);
        if examined_states > max_matching_states {
            return Err(NormalizationError::AromaticLocalizationLimit {
                atom: component[0],
                examined_states,
                limit: max_matching_states,
            });
        }
        if unmatched.is_empty() {
            break selected;
        }
        let Some(atom_id) = unmatched.iter().copied().min_by_key(|atom_id| {
            adjacency
                .get(atom_id)
                .map(|neighbors| {
                    neighbors
                        .iter()
                        .filter(|(neighbor, _)| unmatched.contains(neighbor))
                        .count()
                })
                .unwrap_or(0)
        }) else {
            return Ok(false);
        };
        let candidates = adjacency
            .get(&atom_id)
            .into_iter()
            .flatten()
            .filter(|(neighbor, _)| unmatched.contains(neighbor))
            .copied()
            .collect::<Vec<_>>();
        for (neighbor, bond_id) in candidates.into_iter().rev() {
            let mut next_unmatched = unmatched.clone();
            next_unmatched.remove(&atom_id);
            next_unmatched.remove(&neighbor);
            let mut next_selected = selected.clone();
            next_selected.push(bond_id);
            stack.push((next_unmatched, next_selected));
        }
    };

    let selected_double_bonds = selected_double_bonds.into_iter().collect::<BTreeSet<_>>();
    for (bond_id, bond) in (0..=u32::MAX)
        .zip(molecule.graph.bonds.iter_mut())
        .filter_map(|(raw, bond)| bond.as_mut().map(|bond| (BondId::new(raw), bond)))
    {
        if aromatic_bonds.contains(&bond_id)
            && component_atoms.contains(&bond.a())
            && component_atoms.contains(&bond.b())
        {
            bond.order = if selected_double_bonds.contains(&bond_id) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
        }
    }
    Ok(true)
}

fn source_aromatic_bond_components(
    molecule: &Molecule,
    aromatic_bonds: &BTreeSet<BondId>,
) -> Vec<Vec<AtomId>> {
    let mut adjacency = BTreeMap::<AtomId, Vec<AtomId>>::new();
    for (_, bond) in molecule
        .bonds()
        .filter(|(bond_id, _)| aromatic_bonds.contains(bond_id))
    {
        adjacency.entry(bond.a()).or_default().push(bond.b());
        adjacency.entry(bond.b()).or_default().push(bond.a());
    }
    for neighbors in adjacency.values_mut() {
        neighbors.sort_unstable();
    }

    let mut components = Vec::new();
    let mut visited = BTreeSet::new();
    for start in adjacency.keys().copied() {
        if !visited.insert(start) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        while let Some(atom_id) = stack.pop() {
            component.push(atom_id);
            if let Some(neighbors) = adjacency.get(&atom_id) {
                for neighbor in neighbors.iter().rev().copied() {
                    if visited.insert(neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}

fn represented_bond_valence(order: BondOrder) -> usize {
    match order {
        BondOrder::Zero | BondOrder::Dative => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
    }
}

fn source_aromatic_implicit_hydrogens(atom: &Atom, explicit_valence: usize) -> usize {
    if !atom.hydrogens.allows_implicit() {
        return 0;
    }
    let target = match atom.element.symbol() {
        "B" | "C" => 3,
        "N" | "O" | "S" | "Se" | "Te" => {
            if atom.hydrogens.explicit_count() > 0 || atom.formal_charge > 0 {
                3
            } else {
                2
            }
        }
        "P" => explicit_valence,
        _ => return 0,
    };
    target.saturating_sub(explicit_valence)
}

fn aromatic_localization_target_valence(
    atom: &Atom,
    baseline_bond_valence: usize,
) -> Option<usize> {
    let target = match (atom.element.symbol(), atom.formal_charge) {
        ("B", -1) => 4,
        ("B", 0) => 3,
        ("B", 1) => 2,
        ("C", -1 | 1) => 3,
        ("C", 0)
            if !atom.hydrogens.allows_implicit()
                && atom.hydrogens.explicit_count() == 0
                && baseline_bond_valence == 2 =>
        {
            3
        }
        ("C", 0) => 4,
        ("N" | "P", -1) => 2,
        ("N" | "P", 0) => 3,
        ("N" | "P", 1) => 4,
        ("O" | "S" | "Se" | "Te", -1) => 1,
        ("O" | "S" | "Se" | "Te", 0) => 2,
        ("O" | "S" | "Se" | "Te", 1) => 3,
        _ => return None,
    };
    Some(target)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aromatic_localization_limit_is_structured_and_transactional() {
        let mut molecule = crate::core::MoleculeEditor::new();
        let atoms = (0..6)
            .map(|_| {
                molecule
                    .add_atom(Atom::new(crate::core::Element::from_symbol("C").unwrap()))
                    .unwrap()
            })
            .collect::<Vec<_>>();
        let aromatic_bonds = (0..6)
            .map(|index| {
                molecule
                    .add_bond(atoms[index], atoms[(index + 1) % 6], BondOrder::Single)
                    .unwrap()
            })
            .collect::<BTreeSet<_>>();
        let component = source_aromatic_bond_components(molecule.working(), &aromatic_bonds)
            .into_iter()
            .next()
            .expect("imported aromatic component");
        let before = molecule.clone();

        let error = try_localize_aromatic_component_with_limit(
            molecule.working_mut(),
            &component,
            &aromatic_bonds,
            0,
        )
        .expect_err("zero matching budget should fail structurally");

        assert!(matches!(
            error,
            NormalizationError::AromaticLocalizationLimit {
                examined_states: 1,
                limit: 0,
                ..
            }
        ));
        assert_eq!(molecule, before);
    }
}
