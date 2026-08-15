use std::fmt;

use crate::algorithms::*;
use crate::core::*;
use crate::small::model::SmallMolecule;

use super::normalization::{normalize_molecule, NormalizationError};

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
    pub stereo: Option<StereoPerceptionReport>,
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SanitizeError {
    Normalization(NormalizationError),
    Valence(ValenceError),
    Rings(RingPerceptionError),
    Aromaticity(AromaticityError),
    Stereo(StereoPerceptionError),
}

impl fmt::Display for SanitizeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Normalization(error) => write!(f, "{error}"),
            Self::Valence(error) => write!(f, "{error}"),
            Self::Rings(error) => write!(f, "{error}"),
            Self::Aromaticity(error) => write!(f, "{error}"),
            Self::Stereo(error) => write!(f, "{error}"),
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
    normalize_molecule(staged.graph_mut_raw()).map_err(SanitizeError::Normalization)?;
    prepare_sanitize_states(staged.graph_mut_raw(), options);
    if options.perceive_valence {
        perceive_valence(staged.graph_mut_raw(), ValenceModel::RdkitLike)
            .map_err(SanitizeError::Valence)?;
    }
    if options.perceive_rings {
        perceive_ring_set_with_options(staged.graph_mut_raw(), ring_options)
            .map_err(SanitizeError::Rings)?;
    }
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
        Some(
            perceive_stereo_with_options(
                staged.graph_mut_raw(),
                StereoPerceptionOptions {
                    assign_coordinates: false,
                    ..StereoPerceptionOptions::default()
                },
            )
            .map_err(SanitizeError::Stereo)?,
        )
    } else {
        None
    };
    *molecule = staged;
    Ok(SanitizeReport { stereo })
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
