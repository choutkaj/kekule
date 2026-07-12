mod mmcif;
mod mmcif_document;
mod mmcif_interpret;
mod smiles;
mod v2000;
mod v3000;

use crate::algorithms::{allowed_valences, explicit_valence};
use crate::core::{BondOrder, Molecule, StereoBondMarkKind};

pub use mmcif::*;
pub use mmcif_document::*;
pub use mmcif_interpret::*;
pub use smiles::*;
pub use v2000::*;
pub use v3000::*;

fn preserve_molfile_tetrahedral_hydrogens(mol: &mut Molecule) {
    let mut centers = mol
        .stereo_bond_marks()
        .filter(|mark| {
            matches!(
                mark.kind,
                StereoBondMarkKind::WedgeUp
                    | StereoBondMarkKind::WedgeDown
                    | StereoBondMarkKind::WedgeEither
            )
        })
        .filter_map(|mark| mol.bond(mark.bond).ok().map(|bond| bond.a()))
        .collect::<Vec<_>>();
    centers.sort_unstable();
    centers.dedup();

    for center in centers {
        let Ok(incident) = mol.incident_bonds(center) else {
            continue;
        };
        let incident = incident.collect::<Vec<_>>();
        if incident.len() != 3
            || incident
                .iter()
                .any(|(_, bond)| bond.order != BondOrder::Single)
        {
            continue;
        }
        let Ok(atom) = mol.atom(center) else {
            continue;
        };
        if atom.explicit_hydrogens != 0 || atom.no_implicit_hydrogens {
            continue;
        }
        let explicit = explicit_valence(mol, center);
        let implied_hydrogens = allowed_valences(atom)
            .and_then(|allowed| allowed.iter().copied().find(|target| *target >= explicit))
            .map(|target| target.saturating_sub(explicit))
            .unwrap_or(0);
        if implied_hydrogens == 1 {
            mol.atom_mut(center)
                .expect("existing Molfile atom")
                .explicit_hydrogens = 1;
        }
    }
}
