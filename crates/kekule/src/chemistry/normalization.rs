use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::algorithms::StereoValidationIssue;
use crate::core::{
    Atom, AtomId, BondId, BondOrder, Molecule, MoleculeError, StereoBondMarkKind, StereoElementId,
};

use super::source_stereo::normalize_source_stereo;

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
    UnassembledTetrahedralBondMark {
        bond: BondId,
        kind: StereoBondMarkKind,
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
        kind: StereoBondMarkKind,
    },
    InvalidStereo(StereoValidationIssue),
    CouldNotCreateStereoElement(MoleculeError),
    CouldNotConsumeSourceBondMark(MoleculeError),
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
            "source-stereo normalization reported {} issue(s)",
            self.issues.len()
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

/// Deterministically normalize represented chemistry into Kekule's canonical form.
///
/// Normalization is transactional and idempotent. A successful call clears the
/// complete installed perception state, including when the primary
/// representation was already normalized.
pub fn normalize_molecule(
    molecule: &mut Molecule,
) -> Result<NormalizationReport, NormalizationError> {
    let mut staged = molecule.clone();
    let report = normalize_molecule_in_place(&mut staged)?;
    *molecule = staged;
    Ok(report)
}

pub(crate) fn normalize_molecule_in_place(
    molecule: &mut Molecule,
) -> Result<NormalizationReport, NormalizationError> {
    normalize_hypervalent_oxo_halides(molecule)?;
    // Source-stereo normalization must not observe arbitrary installed
    // perception. Representation rewrites above already invalidate it
    // conceptually, so clear it before decoding any source marks.
    molecule.invalidate_topology();
    let report = normalize_source_stereo(molecule).map_err(NormalizationError::SourceStereo)?;
    // Adding represented stereo invalidates only stereo-derived state. The
    // normalization publication contract clears the complete perception state.
    molecule.invalidate_topology();
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
            baseline_bond_valence.saturating_add(usize::from(atom.explicit_hydrogens));
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
        .zip(molecule.bonds.iter_mut())
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
    if atom.no_implicit_hydrogens {
        return 0;
    }
    let target = match atom.element.symbol() {
        "B" | "C" => 3,
        "N" | "O" | "S" | "Se" | "Te" => {
            if atom.explicit_hydrogens > 0 || atom.formal_charge > 0 {
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
            if atom.no_implicit_hydrogens
                && atom.explicit_hydrogens == 0
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

fn normalize_hypervalent_oxo_halides(molecule: &mut Molecule) -> Result<(), NormalizationError> {
    let halogens = molecule
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.formal_charge == 0
                && matches!(atom.element.symbol(), "Cl" | "Br" | "I")
                && has_terminal_single_bond_oxygen_neighbor(molecule, atom_id))
            .then_some(atom_id)
        })
        .collect::<Vec<_>>();

    for atom_id in halogens {
        let oxo_bonds = oxo_bonds_to_neutral_oxygen(molecule, atom_id);
        if oxo_bonds.is_empty() {
            continue;
        }
        let charge = oxo_bonds.len();
        let formal_charge =
            i8::try_from(charge).map_err(|_| NormalizationError::FormalChargeOutOfRange {
                atom: atom_id,
                charge,
            })?;

        if let Some(atom) = molecule.atoms[atom_id.index()].as_mut() {
            atom.formal_charge = formal_charge;
        }
        for (oxygen_id, bond_id) in oxo_bonds {
            if let Some(atom) = molecule.atoms[oxygen_id.index()].as_mut() {
                atom.formal_charge = -1;
            }
            if let Some(bond) = molecule.bonds[bond_id.index()].as_mut() {
                bond.order = BondOrder::Single;
            }
        }
    }
    Ok(())
}

fn has_terminal_single_bond_oxygen_neighbor(molecule: &Molecule, atom_id: AtomId) -> bool {
    molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(_, bond)| {
            let oxygen_id = bond.other_atom(atom_id);
            bond.order == BondOrder::Single
                && molecule
                    .atom(oxygen_id)
                    .is_ok_and(|neighbor| neighbor.element.symbol() == "O")
                && molecule.incident_bonds(oxygen_id).is_ok_and(|mut bonds| {
                    bonds.all(|(_, oxygen_bond)| {
                        let neighbor_id = oxygen_bond.other_atom(oxygen_id);
                        neighbor_id == atom_id
                            || molecule
                                .atom(neighbor_id)
                                .is_ok_and(|neighbor| neighbor.element.symbol() == "H")
                    })
                })
        })
}

fn oxo_bonds_to_neutral_oxygen(molecule: &Molecule, atom_id: AtomId) -> Vec<(AtomId, BondId)> {
    molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(bond_id, bond)| {
            if bond.order != BondOrder::Double {
                return None;
            }
            let oxygen_id = bond.other_atom(atom_id);
            let oxygen = molecule.atom(oxygen_id).ok()?;
            (oxygen.element.symbol() == "O" && oxygen.formal_charge == 0)
                .then_some((oxygen_id, bond_id))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aromatic_localization_limit_is_structured_and_transactional() {
        let mut molecule = Molecule::new();
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
        let component = source_aromatic_bond_components(&molecule, &aromatic_bonds)
            .into_iter()
            .next()
            .expect("imported aromatic component");
        let before = molecule.clone();

        let error = try_localize_aromatic_component_with_limit(
            &mut molecule,
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
