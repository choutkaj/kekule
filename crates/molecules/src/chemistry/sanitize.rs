use std::fmt;

use crate::algorithms::*;
use crate::core::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeOptions {
    pub perceive_valence: bool,
    pub perceive_rings: bool,
    pub perceive_aromaticity: bool,
}

impl Default for SanitizeOptions {
    fn default() -> Self {
        Self {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SanitizeReport {
    pub valence: Option<ValenceReport>,
    pub ring_count: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SanitizeError {
    Valence(ValenceReport),
    Rings(RingPerceptionError),
    Aromaticity(AromaticityError),
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
    normalize_sanitize_charges(&mut staged.mol);
    let skipped_ring_state =
        (!options.perceive_rings).then(|| invalidate(staged.mol.perception.rings));
    prepare_sanitize_states(&mut staged.mol, options);
    let valence = if options.perceive_valence {
        let report = perceive_valence(&mut staged.mol, ValenceModel::RdkitLike);
        if !report.is_ok() {
            return Err(SanitizeError::Valence(report));
        }
        Some(report)
    } else {
        None
    };
    let ring_count = if options.perceive_rings {
        Some(
            perceive_ring_set_with_options(&mut staged.mol, ring_options)
                .map_err(SanitizeError::Rings)?
                .len(),
        )
    } else {
        None
    };
    if options.perceive_aromaticity {
        perceive_aromaticity_with_ring_options(
            &mut staged.mol,
            AromaticityModel::RdkitLike,
            ring_options,
        )
        .map_err(SanitizeError::Aromaticity)?;
        if let Some(state) = skipped_ring_state {
            staged.mol.perception.rings = state;
            staged.mol.ring_membership = None;
            staged.mol.ring_set = None;
        }
    }
    *molecule = staged;
    Ok(SanitizeReport {
        valence,
        ring_count,
    })
}

fn prepare_sanitize_states(mol: &mut Molecule, options: SanitizeOptions) {
    if !options.perceive_valence {
        mol.perception.valence = invalidate(mol.perception.valence);
    }
    if !options.perceive_rings {
        mol.perception.rings = invalidate(mol.perception.rings);
        mol.ring_membership = None;
        mol.ring_set = None;
    }
    if !options.perceive_aromaticity {
        mol.perception.aromaticity = invalidate(mol.perception.aromaticity);
    }
}

fn normalize_sanitize_charges(mol: &mut Molecule) {
    normalize_hypervalent_oxo_halides(mol);
}

fn normalize_hypervalent_oxo_halides(mol: &mut Molecule) {
    let halogens = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.formal_charge == 0
                && matches!(atom.element.symbol(), "Cl" | "Br" | "I")
                && explicit_valence(mol, atom_id) > 1)
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
                bond.aromatic = false;
                changed = true;
            }
        }
    }
    if changed {
        mol.invalidate_topology();
    }
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
