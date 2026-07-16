use std::fmt;

use crate::algorithms::*;
use crate::core::*;
use crate::small::model::SmallMolecule;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeOptions {
    pub perceive_valence: bool,
    pub perceive_rings: bool,
    pub perceive_aromaticity: bool,
    pub perceive_stereo: bool,
}

impl Default for SanitizeOptions {
    fn default() -> Self {
        Self {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: true,
            perceive_stereo: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SanitizeReport {
    pub valence: Option<ValenceReport>,
    pub ring_count: Option<usize>,
    pub stereo: Option<StereoPerceptionReport>,
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SanitizeError {
    Valence(ValenceReport),
    Rings(RingPerceptionError),
    Aromaticity(AromaticityError),
    Stereo(StereoPerceptionReport),
}

impl fmt::Display for SanitizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Valence(report) => write!(
                f,
                "valence perception reported {} issue(s)",
                report.issues.len()
            ),
            Self::Rings(error) => write!(f, "{error}"),
            Self::Aromaticity(error) => write!(f, "{error}"),
            Self::Stereo(report) => write!(
                f,
                "stereo perception reported {} issue(s)",
                report.issues.len()
            ),
        }
    }
}

impl std::error::Error for SanitizeError {}

pub fn sanitize_small_molecule(
    molecule: &mut SmallMolecule,
    options: SanitizeOptions,
) -> std::result::Result<SanitizeReport, SanitizeError> {
    sanitize_small_molecule_with_ring_options(molecule, options, RingPerceptionOptions::default())
}

pub fn sanitize_small_molecule_with_ring_options(
    molecule: &mut SmallMolecule,
    options: SanitizeOptions,
    ring_options: RingPerceptionOptions,
) -> std::result::Result<SanitizeReport, SanitizeError> {
    let mut staged = molecule.clone();
    normalize_sanitize_charges(staged.graph_mut_raw());
    prepare_sanitize_states(staged.graph_mut_raw(), options);
    let valence = if options.perceive_valence {
        let report = perceive_valence(staged.graph_mut_raw(), ValenceModel::RdkitLike);
        if !report.is_ok() {
            return Err(SanitizeError::Valence(report));
        }
        Some(report)
    } else {
        None
    };
    let ring_count = if options.perceive_rings {
        Some(
            perceive_ring_set_with_options(staged.graph_mut_raw(), ring_options)
                .map_err(SanitizeError::Rings)?
                .len(),
        )
    } else {
        None
    };
    if options.perceive_aromaticity {
        perceive_aromaticity_with_ring_options(
            staged.graph_mut_raw(),
            AromaticityModel::RdkitLike,
            ring_options,
        )
        .map_err(SanitizeError::Aromaticity)?;
        if options.perceive_valence {
            normalize_aromatic_nitrogen_hydrogens(staged.graph_mut_raw());
        }
        if !options.perceive_rings {
            staged.graph_mut_raw().discard_ring_results();
        }
    }
    let stereo = if options.perceive_stereo {
        let report = perceive_stereo_with_options(
            staged.graph_mut_raw(),
            StereoPerceptionOptions {
                assign_coordinates: false,
                ..StereoPerceptionOptions::default()
            },
        );
        if stereo_report_has_fatal_issues(&report) {
            return Err(SanitizeError::Stereo(report));
        }
        Some(report)
    } else {
        None
    };
    *molecule = staged;
    Ok(SanitizeReport {
        valence,
        ring_count,
        stereo,
    })
}

fn stereo_report_has_fatal_issues(report: &StereoPerceptionReport) -> bool {
    report.issues.iter().any(|issue| {
        !matches!(
            issue,
            StereoPerceptionIssue::AmbiguousTetrahedralWedgeMarks { .. }
        )
    })
}

fn prepare_sanitize_states(mol: &mut Molecule, options: SanitizeOptions) {
    if !options.perceive_valence {
        mol.clear_valence();
    }
    if !options.perceive_rings {
        mol.clear_rings();
    }
    if !options.perceive_aromaticity {
        mol.clear_aromaticity();
    }
    if !options.perceive_stereo {
        mol.clear_cip_descriptors();
    }
}

fn normalize_sanitize_charges(mol: &mut Molecule) {
    normalize_hypervalent_oxo_halides(mol);
}

fn normalize_aromatic_nitrogen_hydrogens(mol: &mut Molecule) {
    let nitrogens = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            let aromatic = mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true);
            let implicit = mol.implicit_hydrogens(atom_id).ok().flatten().unwrap_or(0);
            (atom.element.symbol() == "N"
                && aromatic
                && atom.formal_charge == 0
                && atom.explicit_hydrogens.saturating_add(implicit) == 1)
                .then_some(atom_id)
        })
        .collect::<Vec<_>>();
    for atom_id in nitrogens {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.explicit_hydrogens = 1;
            atom.no_implicit_hydrogens = false;
        }
        mol.set_implicit_hydrogens(atom_id, 0);
    }
}

fn normalize_hypervalent_oxo_halides(mol: &mut Molecule) {
    let halogens = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.formal_charge == 0
                && matches!(atom.element.symbol(), "Cl" | "Br" | "I")
                && has_terminal_single_bond_oxygen_neighbor(mol, atom_id))
            .then_some(atom_id)
        })
        .collect::<Vec<_>>();

    let mut changed = false;
    for atom_id in halogens {
        let oxo_bonds = oxo_bonds_to_neutral_oxygen(mol, atom_id);
        if oxo_bonds.is_empty() {
            continue;
        };
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.formal_charge = atom
                .formal_charge
                .saturating_add(i8::try_from(oxo_bonds.len()).unwrap_or(i8::MAX));
            changed = true;
        }
        for (oxygen_id, bond_id) in oxo_bonds {
            if let Some(atom) = mol.atoms[oxygen_id.index()].as_mut() {
                atom.formal_charge = -1;
                changed = true;
            }
            if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                bond.order = BondOrder::Single;
                changed = true;
            }
        }
    }
    if changed {
        mol.invalidate_topology();
    }
}

fn has_terminal_single_bond_oxygen_neighbor(mol: &Molecule, atom_id: AtomId) -> bool {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(_, bond)| {
            let oxygen_id = bond.other_atom(atom_id);
            bond.order == BondOrder::Single
                && mol
                    .atom(oxygen_id)
                    .is_ok_and(|neighbor| neighbor.element.symbol() == "O")
                && mol.incident_bonds(oxygen_id).is_ok_and(|mut bonds| {
                    bonds.all(|(_, oxygen_bond)| {
                        let neighbor_id = oxygen_bond.other_atom(oxygen_id);
                        neighbor_id == atom_id
                            || mol
                                .atom(neighbor_id)
                                .is_ok_and(|neighbor| neighbor.element.symbol() == "H")
                    })
                })
        })
}

fn oxo_bonds_to_neutral_oxygen(mol: &Molecule, atom_id: AtomId) -> Vec<(AtomId, BondId)> {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(bond_id, bond)| {
            if !matches!(bond.order, BondOrder::Double) {
                return None;
            }
            let oxygen_id = if bond.a == atom_id { bond.b } else { bond.a };
            let oxygen = mol.atom(oxygen_id).ok()?;
            (oxygen.element.symbol() == "O" && oxygen.formal_charge == 0)
                .then_some((oxygen_id, bond_id))
        })
        .collect()
}
