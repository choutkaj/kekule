//! Molfile record preparation: atom order, emitted coordinates, and stereo
//! projected against them.
use std::collections::BTreeMap;

use crate::chemistry::{
    project_molfile_stereo_bond_marks, AtomPositionSource, SourceStereoBondMarkKind,
};
use crate::core::{Atom, AtomId, BondId, BondOrder, Molecule, StereoElementKind, StereoGroupKind};
use crate::geometry::Point3;
use crate::structure::ModelView;
use crate::topology::{InstanceAtomId, MoleculeInstanceId};
use crate::units::ANGSTROM;

use crate::io::structure_documents::{
    materialize_molfile_stereo_hydrogens, molfile_coordinate_double_bonds,
    molfile_coordinate_tetrahedra, molfile_coordinates_are_flat,
    molfile_stereo_group_members_at_atom,
};
use crate::io::{MoleculeWriteError, MolfileVersion};

const PROGRAM: &str = "kekule";

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
    /// The CTfile chiral flag. Readers following the BIOVIA convention, and
    /// RDKit's `TranslateChiralFlagToStereoGroups`, read an unset flag as
    /// putting every ungrouped center in one racemic AND group, so the flag
    /// is set whenever the source asserts an absolute configuration.
    pub(super) chiral: bool,
}

#[derive(Debug)]
pub(super) struct MolfileStereoGroup {
    pub(super) kind: StereoGroupKind,
    pub(super) atoms: Vec<u64>,
}

impl<'a> MolfileRecord<'a> {
    /// Writes the three header lines. The program line follows the CTfile
    /// `IIPPPPPPPPMMDDYYHHmmdd` layout with blank initials and timestamp, so
    /// output stays deterministic and the dimensional code lands in columns
    /// 21-22. Readers such as RDKit take any nonzero z as 3D whatever this
    /// code says, so it is `2D` exactly when every emitted z is zero and `3D`
    /// otherwise.
    pub(super) fn push_header(&self, out: &mut String, title: &str) {
        let dimension = if self.atoms.iter().all(|atom| atom.position.z == 0.0) {
            "2D"
        } else {
            "3D"
        };
        out.push_str(&format!("{title}\n  {PROGRAM:<8}{:10}{dimension}\n\n", ""));
    }

    pub(super) fn molecule(molecule: &'a Molecule) -> Result<Self, MoleculeWriteError> {
        let atoms = molecule
            .atoms()
            .map(|(id, atom)| MolfileAtom {
                molecule,
                id,
                atom,
                position: Point3::default(),
            })
            .collect::<Vec<_>>();
        reject_zero_atom_maps(&atoms)?;
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
            .collect::<Vec<_>>();
        let stereo_groups =
            project_stereo_groups(molecule, |atom| indexes[&(ComponentKey::Molecule, atom)])?;
        Ok(Self {
            atoms,
            bonds,
            stereo_groups,
            chiral: asserts_absolute_stereo(molecule),
        })
    }

    pub(super) fn model(
        model: ModelView<'a>,
        version: MolfileVersion,
    ) -> Result<Self, MoleculeWriteError> {
        let topology = model.topology();
        let mut atoms = Vec::with_capacity(topology.atom_count());
        let mut indexes = BTreeMap::new();
        let mut positions = BTreeMap::new();
        for (serial, view) in (1u64..).zip(topology.atoms()) {
            let qualified = view.id();
            let atom = view.atom();
            let occurrence = view.molecule();
            let position = model
                .position(qualified)
                .map_err(|error| MoleculeWriteError::invalid_model(error.to_string()))?
                .value_in(ANGSTROM)
                .map_err(|error| MoleculeWriteError::invalid_model(error.to_string()))?;
            // V2000 has four-place fixed-width fields. V3000 emits round-trip
            // decimal values. Stereo projection must use each format's actual
            // emitted coordinates, including V2000 tie rounding.
            let position = match version {
                MolfileVersion::V2000 => {
                    let [x, y, z] =
                        [position.x, position.y, position.z].map(super::v2000::emitted_coordinate);
                    Point3::new(x, y, z)
                }
                MolfileVersion::V3000 => position,
            };
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

        reject_zero_atom_maps(&atoms)?;
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
        for view in topology.bonds() {
            let qualified = view.id();
            let bond = view.bond();
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
            ));
        }
        // Readers such as RDKit take any nonzero z as 3D, then derive stereo
        // from coordinates and ignore bond marks. A flat drawing whose marks
        // carry stereo is emitted on z = 0, which moves no mark because XY is
        // unchanged. Every other record keeps its coordinates.
        let points = atoms.iter().map(|atom| atom.position).collect::<Vec<_>>();
        if molfile_coordinates_are_flat(&points)
            && bonds.iter().any(|bond| {
                bond.stereo
                    .is_some_and(|kind| kind != SourceStereoBondMarkKind::MolfileDoubleBondGeometry)
            })
        {
            for atom in &mut atoms {
                atom.position.z = 0.0;
            }
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
            chiral: topology
                .molecules()
                .any(|occurrence| asserts_absolute_stereo(occurrence.molecule())),
        })
    }
}

/// Both Molfile versions write and read atom map 0 as no atom map.
fn reject_zero_atom_maps(atoms: &[MolfileAtom<'_>]) -> Result<(), MoleculeWriteError> {
    match atoms.iter().position(|atom| atom.atom.atom_map == Some(0)) {
        Some(index) => Err(MoleculeWriteError::unsupported(format!(
            "Molfile cannot encode atom map 0 on atom {}, which it reads as no atom map",
            index + 1
        ))),
        None => Ok(()),
    }
}

/// Whether a specified tetrahedral or axis configuration is absolute: either
/// ungrouped, which Kekule reads as absolute, or in an ABS group.
fn asserts_absolute_stereo(molecule: &Molecule) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        element.is_specified()
            && !matches!(element.kind, StereoElementKind::DoubleBond(_))
            && element.group.is_none_or(|group| {
                molecule
                    .stereo_group(group)
                    .expect("a grouped stereo element names a live group")
                    .kind
                    == StereoGroupKind::Absolute
            })
    })
}

fn prepare_bond(
    component: ComponentKey,
    id: BondId,
    order: BondOrder,
    endpoints: (AtomId, AtomId),
    projected: &BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>,
    indexes: &BTreeMap<(ComponentKey, AtomId), u64>,
) -> MolfileBond {
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
    let serial = |atom| indexes[&(component, atom)];
    MolfileBond {
        order,
        from: serial(from),
        to: serial(to),
        stereo: projection.map(|(_, kind)| kind),
    }
}

fn stereo_projections(
    molecule: &Molecule,
    geometry: Option<&dyn AtomPositionSource>,
) -> Result<BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>, MoleculeWriteError> {
    // A 3D Molfile implicitly asserts configurations even without wedges.
    // Add unknown annotations in detached writer state where the source graph
    // makes no assertion, so writing a conformation cannot invent chemistry.
    let mut staged;
    let molecule = if let Some(geometry) = geometry {
        let mut warnings = Vec::new();
        let inferred = molfile_coordinate_tetrahedra(molecule, geometry, &mut warnings)
            .map_err(|error| MoleculeWriteError::invalid_model(error.to_string()))?;
        if !warnings.is_empty() {
            return Err(MoleculeWriteError::unsupported(format!(
                "cannot verify 3D stereo preservation: {warnings:?}"
            )));
        }
        if inferred.is_empty() {
            molecule
        } else {
            staged = molecule.clone();
            for mut element in inferred {
                if let StereoElementKind::Tetrahedral(stereo) = &mut element.kind {
                    materialize_molfile_stereo_hydrogens(&mut staged, [stereo.center]);
                    stereo.orientation = None;
                }
                staged
                    .add_stereo_element(element)
                    .map_err(|error| MoleculeWriteError::invalid_model(error.to_string()))?;
            }
            &staged
        }
    } else {
        molecule
    };
    // The reader keeps only the double-bond configurations a conformation's
    // coordinates specify, so planning verifies against the same set.
    let coordinate_double_bonds = geometry
        .map(|geometry| molfile_coordinate_double_bonds(molecule, geometry))
        .transpose()
        .map_err(|error| MoleculeWriteError::invalid_model(error.to_string()))?
        .flatten();
    Ok(
        project_molfile_stereo_bond_marks(molecule, geometry, coordinate_double_bonds.as_ref())
            .map_err(MoleculeWriteError::unsupported)?
            .into_iter()
            .map(|(bond, projection)| (bond, (projection.from, projection.kind)))
            .collect(),
    )
}

fn project_stereo_groups(
    molecule: &Molecule,
    atom_serial: impl Fn(AtomId) -> u64,
) -> Result<Vec<MolfileStereoGroup>, MoleculeWriteError> {
    molecule
        .stereo_groups()
        .map(|(_, group)| {
            if !matches!(
                group.kind,
                StereoGroupKind::Absolute | StereoGroupKind::And | StereoGroupKind::Or
            ) {
                return Err(MoleculeWriteError::unsupported(
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
                    .expect("a stereo group member is a live stereo element");
                let candidates = match &element.kind {
                    StereoElementKind::Tetrahedral(stereo) => vec![stereo.center],
                    StereoElementKind::Axis(stereo) => {
                        let bond = molecule
                            .bond(stereo.axis)
                            .expect("an axis stereo element names a live bond");
                        vec![bond.a(), bond.b()]
                    }
                    StereoElementKind::DoubleBond(_) => {
                        return Err(MoleculeWriteError::unsupported(
                            "Molfile writer cannot encode double-bond stereo group members",
                        ))
                    }
                };
                let atom = candidates
                    .into_iter()
                    .find(|atom| molfile_stereo_group_members_at_atom(molecule, *atom) == [*member])
                    .ok_or_else(|| {
                        MoleculeWriteError::unsupported(
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
