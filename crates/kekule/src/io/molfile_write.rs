use std::collections::BTreeMap;

use crate::chemistry::{
    project_molfile_stereo_bond_marks, AtomPositionSource, SourceStereoBondMarkKind,
};
use crate::core::{Atom, AtomId, BondId, BondOrder, Molecule, StereoElementKind, StereoGroupKind};
use crate::geometry::Point3;
use crate::structure::ModelView;
use crate::topology::{InstanceAtomId, MoleculeInstanceId};
use crate::units::ANGSTROM;

use super::structure_documents::molfile_stereo_group_members_at_atom;
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
    pub(super) stereo_groups: Vec<MolfileStereoGroup>,
}

#[derive(Debug)]
pub(super) struct MolfileStereoGroup {
    pub(super) kind: StereoGroupKind,
    pub(super) atoms: Vec<u64>,
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
        let projected = stereo_projections(molecule, None)?;
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
        let stereo_groups =
            project_stereo_groups(molecule, |atom| indexes[&(ComponentKey::Molecule, atom)])?;
        Ok(Self {
            atoms,
            bonds,
            stereo_groups,
        })
    }

    pub(super) fn model(model: ModelView<'a>) -> Result<Self, MolWriteError> {
        let topology = model.topology();
        let mut atoms = Vec::with_capacity(topology.atom_count());
        let mut indexes = BTreeMap::new();
        let mut positions = BTreeMap::new();
        for (serial, (qualified, atom)) in (1u64..).zip(topology.atoms()) {
            let occurrence = topology
                .molecule(qualified.molecule())
                .map_err(|error| MolWriteError::invalid_model(error.to_string()))?;
            let position = model
                .position(qualified)
                .map_err(|error| MolWriteError::invalid_model(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MolWriteError::invalid_model(error.to_string()))?;
            // Parse the emitted decimal representation once so validation and
            // rendering use exactly the same coordinates, including tie rounding.
            let [x, y, z] = [position.x, position.y, position.z].map(|value| {
                format!("{value:.4}")
                    .parse()
                    .expect("a formatted coordinate is a valid floating-point number")
            });
            let position = Point3::new(x, y, z);
            positions.insert(qualified, position);
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
            let geometry = ModelStereoPositions {
                positions: &positions,
                instance: occurrence.id(),
            };
            projections.insert(
                occurrence.id(),
                stereo_projections(occurrence.molecule(), Some(&geometry))?,
            );
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
        let mut stereo_groups = Vec::new();
        for occurrence in topology.molecules() {
            stereo_groups.extend(project_stereo_groups(occurrence.molecule(), |atom| {
                indexes[&(ComponentKey::Instance(occurrence.id()), atom)]
            })?);
        }
        Ok(Self {
            atoms,
            bonds,
            stereo_groups,
        })
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
    geometry: Option<&dyn AtomPositionSource>,
) -> Result<BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>, MolWriteError> {
    Ok(project_molfile_stereo_bond_marks(molecule, geometry)
        .map_err(MolWriteError::new)?
        .into_iter()
        .map(|(bond, projection)| (bond, (projection.from, projection.kind)))
        .collect())
}

fn project_stereo_groups(
    molecule: &Molecule,
    atom_serial: impl Fn(AtomId) -> u64,
) -> Result<Vec<MolfileStereoGroup>, MolWriteError> {
    molecule
        .stereo_groups()
        .map(|(_, group)| {
            if !matches!(
                group.kind,
                StereoGroupKind::Absolute | StereoGroupKind::And | StereoGroupKind::Or
            ) {
                return Err(MolWriteError::new(
                    "Molfile writer cannot encode this stereo group kind",
                ));
            }
            let mut projected = MolfileStereoGroup {
                kind: group.kind,
                atoms: Vec::new(),
            };
            for member in &group.members {
                let element = molecule
                    .stereo_element(*member)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                let candidates = match &element.kind {
                    StereoElementKind::Tetrahedral(stereo) => vec![stereo.center],
                    StereoElementKind::Axis(stereo) => {
                        let bond = molecule
                            .bond(stereo.axis)
                            .map_err(|error| MolWriteError::new(error.to_string()))?;
                        vec![bond.a(), bond.b()]
                    }
                    StereoElementKind::DoubleBond(_) => {
                        return Err(MolWriteError::new(
                            "Molfile writer cannot encode double-bond stereo group members",
                        ))
                    }
                };
                let atom = candidates
                    .into_iter()
                    .find(|atom| molfile_stereo_group_members_at_atom(molecule, *atom) == [*member])
                    .ok_or_else(|| {
                        MolWriteError::new(
                            "V3000 stereo group member has no unambiguous atom or axis endpoint",
                        )
                    })?;
                projected.atoms.push(atom_serial(atom));
            }
            projected.atoms.sort_unstable();
            Ok(projected)
        })
        .collect()
}

struct ModelStereoPositions<'a> {
    positions: &'a BTreeMap<InstanceAtomId, Point3>,
    instance: MoleculeInstanceId,
}

impl AtomPositionSource for ModelStereoPositions<'_> {
    fn position_value(&self, atom: AtomId) -> Option<Point3> {
        self.positions
            .get(&InstanceAtomId::new(self.instance, atom))
            .copied()
    }
}
