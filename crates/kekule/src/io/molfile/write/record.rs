//! Molfile record preparation: atom order, emitted coordinates, and stereo
//! projected against them.
use std::collections::BTreeMap;

use crate::algorithms::{
    is_metal, onium_hydrogens, rdkit_inferred_hydrogen_count, rdkit_reader_hydrogens,
};
use crate::io::{project_molfile_stereo_bond_marks, AtomPositionSource, SourceStereoBondMarkKind};
use crate::math::Point3;
use crate::molecule::{Atom, AtomId, BondId, BondOrder, ImplicitHydrogens, Molecule};
use crate::stereo::{StereoElementKind, StereoGroupKind};
use crate::structure::StructureView;
use crate::topology::{InstanceAtomId, MoleculeInstanceId};
use crate::units::ANGSTROM;

use crate::diagnostics::{Capitalized, FormatError};
use crate::io::molfile::interpret::{
    materialize_molfile_stereo_hydrogens, molfile_coordinate_double_bonds,
    molfile_coordinate_tetrahedra, molfile_coordinates_are_flat,
    molfile_stereo_group_members_at_atom,
};
use crate::io::write_policy::MOLFILE_WRITE;
use crate::io::MolfileVersion;

const PROGRAM: &str = "kekule";

/// The hydrogens a writer states, as a valence, for an atom with inferred
/// hydrogens that a reader would read otherwise: those of an onium centre,
/// which this crate's reader charges instead, and a count that differs from
/// RDKit's rule, which Molfile readers share, as the none of a metal atom.
pub(super) fn stated_inferred_hydrogens(
    molecule: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
) -> Option<u8> {
    onium_hydrogens(molecule, atom_id, atom).or_else(|| {
        let hydrogens = rdkit_inferred_hydrogen_count(molecule, atom_id, atom);
        (hydrogens != rdkit_reader_hydrogens(molecule, atom_id, atom)).then_some(hydrogens)
    })
}

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

    pub(super) fn molecule(molecule: &'a Molecule) -> Result<Self, FormatError> {
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
                    &|atom| is_metal_atom(molecule, atom),
                )
            })
            .collect::<Vec<_>>();
        check_dative_acceptors(&atoms, &bonds)?;
        let stereo_groups =
            project_stereo_groups(molecule, |atom| indexes[&(ComponentKey::Molecule, atom)])?;
        Ok(Self {
            atoms,
            bonds,
            stereo_groups,
            chiral: asserts_absolute_stereo(molecule),
        })
    }

    pub(super) fn structure(
        structure: StructureView<'a>,
        version: MolfileVersion,
    ) -> Result<Self, FormatError> {
        let topology = structure.topology();
        let mut atoms = Vec::with_capacity(topology.atom_count());
        let mut indexes = BTreeMap::new();
        let mut positions = BTreeMap::new();
        for (serial, view) in (1u64..).zip(topology.atoms()) {
            let qualified = view.id();
            let atom = view.atom();
            let instance = view.instance();
            let position = structure
                .position(qualified)
                .ok_or_else(|| {
                    MOLFILE_WRITE
                        .invalid_structure(format!("{} has no position.", Capitalized(qualified)))
                })?
                .value_in(ANGSTROM)
                .map_err(|error| {
                    MOLFILE_WRITE
                        .invalid_structure(format!("Cannot read the position of {qualified}."))
                        .with_source(error)
                })?;
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
                    ComponentKey::Instance(qualified.instance()),
                    qualified.atom(),
                ),
                serial,
            );
            atoms.push(MolfileAtom {
                molecule: instance.molecule(),
                id: qualified.atom(),
                atom,
                position,
            });
        }

        let mut projections = BTreeMap::new();
        for instance in topology.instances() {
            let geometry = StructureStereoPositions {
                positions: &positions,
                instance: instance.id(),
            };
            projections.insert(
                instance.id(),
                stereo_projections(instance.molecule(), Some(&geometry))?,
            );
        }
        let mut bonds = Vec::with_capacity(topology.bond_count());
        for view in topology.bonds() {
            let qualified = view.id();
            let bond = view.bond();
            let projected = projections
                .get(&qualified.instance())
                .expect("every topology instance has a stereo projection");
            bonds.push(prepare_bond(
                ComponentKey::Instance(qualified.instance()),
                qualified.bond(),
                bond.order,
                bond.endpoints(),
                projected,
                &indexes,
                &|atom| is_metal_atom(view.instance().molecule(), atom),
            ));
        }
        check_dative_acceptors(&atoms, &bonds)?;
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
        for instance in topology.instances() {
            stereo_groups.extend(project_stereo_groups(instance.molecule(), |atom| {
                indexes[&(ComponentKey::Instance(instance.id()), atom)]
            })?);
        }
        Ok(Self {
            atoms,
            bonds,
            stereo_groups,
            chiral: topology
                .instances()
                .any(|instance| asserts_absolute_stereo(instance.molecule())),
        })
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

fn is_metal_atom(molecule: &Molecule, atom: AtomId) -> bool {
    molecule
        .atom(atom)
        .is_some_and(|atom| is_metal(atom.element))
}

/// Readers such as RDKit count a dative bond toward the valence of its
/// second atom, the acceptor, and none of them gives that atom back the
/// hydrogens this crate's model infers with the bond not counted. A metal
/// takes no hydrogens, so a metal acceptor reads back unchanged.
fn check_dative_acceptors(
    atoms: &[MolfileAtom<'_>],
    bonds: &[MolfileBond],
) -> Result<(), FormatError> {
    for bond in bonds.iter().filter(|bond| bond.order == BondOrder::Dative) {
        let acceptor = &atoms[usize::try_from(bond.to - 1).expect("serials index atoms")];
        if is_metal(acceptor.atom.element) {
            continue;
        }
        let hydrogens = match acceptor.atom.hydrogens {
            ImplicitHydrogens::Fixed(count) => count,
            ImplicitHydrogens::Inferred => {
                rdkit_inferred_hydrogen_count(acceptor.molecule, acceptor.id, acceptor.atom)
            }
        };
        if hydrogens > 0 {
            return Err(MOLFILE_WRITE.unsupported(format!(
                "Molfile readers count a dative bond toward its acceptor, {}, so they would read \
                 fewer hydrogens on it. Add its hydrogens as atoms, or bond it another way.",
                acceptor.id
            )));
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn prepare_bond(
    component: ComponentKey,
    id: BondId,
    order: BondOrder,
    endpoints: (AtomId, AtomId),
    projected: &BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>,
    indexes: &BTreeMap<(ComponentKey, AtomId), u64>,
    is_metal: &dyn Fn(AtomId) -> bool,
) -> MolfileBond {
    let projection = projected.get(&id).copied();
    // A dative bond is written toward its metal end, which readers take as
    // the acceptor.
    let endpoints = if order == BondOrder::Dative && is_metal(endpoints.0) && !is_metal(endpoints.1)
    {
        (endpoints.1, endpoints.0)
    } else {
        endpoints
    };
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
) -> Result<BTreeMap<BondId, (AtomId, SourceStereoBondMarkKind)>, FormatError> {
    // A 3D Molfile implicitly asserts configurations even without wedges.
    // Add unknown annotations in detached writer state where the source graph
    // makes no assertion, so writing a conformation cannot invent chemistry.
    let staged;
    let molecule = if let Some(geometry) = geometry {
        let inferred = molfile_coordinate_tetrahedra(molecule, geometry)
            .map_err(|error| {
                MOLFILE_WRITE
                    .invalid_structure("Cannot perceive stereo from the 3D coordinates.")
                    .with_source(error)
            })?
            .map_err(|valence| {
                MOLFILE_WRITE
                    .unsupported(
                        "Cannot check that the 3D coordinates keep the stereo of the molecule.",
                    )
                    .with_source(valence)
            })?;
        if inferred.is_empty() {
            molecule
        } else {
            let mut editor = molecule.to_editor();
            for mut element in inferred {
                if let StereoElementKind::Tetrahedral(stereo) = &mut element.kind {
                    materialize_molfile_stereo_hydrogens(&mut editor, [stereo.center]);
                    stereo.orientation = None;
                }
                editor.add_stereo_element(element).map_err(|error| {
                    MOLFILE_WRITE
                        .invalid_structure("Cannot add the stereo of the 3D coordinates.")
                        .with_source(error)
                })?;
            }
            staged = editor.snapshot_keeping_ids();
            &staged
        }
    } else {
        molecule
    };
    // The reader keeps only the double-bond configurations the coordinates,
    // flat or not, specify, so planning verifies against the same set.
    let coordinate_double_bonds = geometry
        .map(|geometry| molfile_coordinate_double_bonds(molecule, geometry))
        .transpose()
        .map_err(|error| {
            MOLFILE_WRITE
                .invalid_structure("Cannot perceive double-bond stereo from the 3D coordinates.")
                .with_source(error)
        })?
        .flatten();
    Ok(
        project_molfile_stereo_bond_marks(molecule, geometry, coordinate_double_bonds.as_ref())?
            .into_iter()
            .map(|(bond, projection)| (bond, (projection.from, projection.kind)))
            .collect(),
    )
}

fn project_stereo_groups(
    molecule: &Molecule,
    atom_serial: impl Fn(AtomId) -> u64,
) -> Result<Vec<MolfileStereoGroup>, FormatError> {
    molecule
        .stereo_groups()
        .map(|(_, group)| {
            if !matches!(
                group.kind,
                StereoGroupKind::Absolute | StereoGroupKind::And | StereoGroupKind::Or
            ) {
                return Err(MOLFILE_WRITE
                    .unsupported("A Molfile cannot state racemic or relative stereo groups."));
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
                        return Err(MOLFILE_WRITE
                            .unsupported("A Molfile stereo group cannot hold double-bond stereo."))
                    }
                };
                let atom = candidates
                    .into_iter()
                    .find(|atom| molfile_stereo_group_members_at_atom(molecule, *atom) == [*member])
                    .ok_or_else(|| {
                        MOLFILE_WRITE
                            .unsupported("A stereo group member has no atom of its own to list.")
                    })?;
                projected.atoms.push(atom_serial(atom));
            }
            projected.atoms.sort_unstable();
            Ok(projected)
        })
        .collect()
}

struct StructureStereoPositions<'a> {
    positions: &'a BTreeMap<InstanceAtomId, Point3>,
    instance: MoleculeInstanceId,
}

impl AtomPositionSource for StructureStereoPositions<'_> {
    fn position_value(&self, atom: AtomId) -> Option<Point3> {
        self.positions
            .get(&InstanceAtomId::new(self.instance, atom))
            .copied()
    }
}
