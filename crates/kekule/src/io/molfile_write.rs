use std::collections::BTreeMap;

use crate::chemistry::{project_molfile_stereo_bond_marks, SourceStereoBondMarkKind};
use crate::core::{Atom, AtomId, BondId, BondOrder, Molecule};
use crate::geometry::Point3;
use crate::structure::ModelView;
use crate::topology::MoleculeInstanceId;
use crate::units::ANGSTROM;

use super::MolWriteError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum ComponentKey {
    Molecule,
    Instance(MoleculeInstanceId),
}

#[derive(Debug, Clone, Copy)]
pub(super) struct MolfileAtom<'a> {
    pub(super) molecule: &'a Molecule,
    pub(super) id: AtomId,
    pub(super) atom: &'a Atom,
    pub(super) position: Point3,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct MolfileBond {
    pub(super) order: BondOrder,
    pub(super) from: u64,
    pub(super) to: u64,
    pub(super) stereo: Option<SourceStereoBondMarkKind>,
}

#[derive(Debug)]
pub(super) struct MolfileRecord<'a> {
    pub(super) atoms: Vec<MolfileAtom<'a>>,
    pub(super) bonds: Vec<MolfileBond>,
}

impl<'a> MolfileRecord<'a> {
    pub(super) fn molecule(molecule: &'a Molecule) -> Result<Self, MolWriteError> {
        let atoms = molecule
            .atoms()
            .map(|(id, atom)| MolfileAtom {
                molecule,
                id,
                atom,
                position: Point3::default(),
            })
            .collect::<Vec<_>>();
        let indexes = molecule
            .atom_ids()
            .zip(1u64..)
            .map(|(atom, serial)| ((ComponentKey::Molecule, atom), serial))
            .collect::<BTreeMap<_, _>>();
        let projected = stereo_projections(molecule)?;
        let bonds = molecule
            .bonds()
            .map(|(id, bond)| {
                prepare_bond(
                    ComponentKey::Molecule,
                    id,
                    bond.order,
                    bond.endpoints(),
                    &projected,
                    &indexes,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { atoms, bonds })
    }

    pub(super) fn model(model: ModelView<'a>) -> Result<Self, MolWriteError> {
        let topology = model.topology();
        let mut atoms = Vec::with_capacity(topology.atom_count());
        let mut indexes = BTreeMap::new();
        for (serial, (qualified, atom)) in (1u64..).zip(topology.atoms()) {
            let occurrence = topology
                .molecule(qualified.molecule())
                .map_err(|error| MolWriteError::invalid_model(error.to_string()))?;
            let position = model
                .position(qualified)
                .map_err(|error| MolWriteError::invalid_model(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MolWriteError::invalid_model(error.to_string()))?;
            indexes.insert(
                (
                    ComponentKey::Instance(qualified.molecule()),
                    qualified.atom(),
                ),
                serial,
            );
            atoms.push(MolfileAtom {
                molecule: occurrence.molecule(),
                id: qualified.atom(),
                atom,
                position,
            });
        }

        let mut projections = BTreeMap::new();
        for occurrence in topology.molecules() {
            projections.insert(occurrence.id(), stereo_projections(occurrence.molecule())?);
        }
        let mut bonds = Vec::with_capacity(topology.bond_count());
        for (qualified, bond) in topology.bonds() {
            let projected = projections
                .get(&qualified.molecule())
                .expect("every topology instance has a stereo projection");
            bonds.push(prepare_bond(
                ComponentKey::Instance(qualified.molecule()),
                qualified.bond(),
                bond.order,
                bond.endpoints(),
                projected,
                &indexes,
            )?);
        }
        Ok(Self { atoms, bonds })
    }
}

fn prepare_bond(
    component: ComponentKey,
    id: BondId,
    order: BondOrder,
    endpoints: (AtomId, AtomId),
    projected: &BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>,
    indexes: &BTreeMap<(ComponentKey, AtomId), u64>,
) -> Result<MolfileBond, MolWriteError> {
    let projection = projected.get(&id).copied();
    let (from, to) = projection
        .map(|(from, _)| {
            let other = if endpoints.0 == from {
                endpoints.1
            } else {
                endpoints.0
            };
            (from, other)
        })
        .unwrap_or(endpoints);
    Ok(MolfileBond {
        order,
        from: *indexes
            .get(&(component, from))
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from atom table"))?,
        to: *indexes
            .get(&(component, to))
            .ok_or_else(|| MolWriteError::new("bond endpoint missing from atom table"))?,
        stereo: projection.map(|(_, kind)| kind),
    })
}

fn stereo_projections(
    molecule: &Molecule,
) -> Result<BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>, MolWriteError> {
    Ok(project_molfile_stereo_bond_marks(molecule)
        .map_err(MolWriteError::new)?
        .into_iter()
        .map(|(bond, projection)| (bond, (projection.from, projection.kind)))
        .collect())
}
