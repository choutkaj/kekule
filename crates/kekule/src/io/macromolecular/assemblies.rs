//! Biological assemblies a source states, and the structures they make: mmCIF
//! `_pdbx_struct_assembly`, `_pdbx_struct_assembly_gen` and
//! `_pdbx_struct_oper_list`, and PDB `REMARK 350` biomolecules.
//!
//! An assembly is copies of chains, each placed by an operation: one stated
//! operator, or the product of several (`(1-60)(61-88)`). Operations apply
//! to the coordinates as they were read. An assembly the source states
//! incompletely, whose operators are not rotations, or that would take the
//! source's operator uses past 100,000 is left out of the report with an
//! issue.
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fmt;
use std::sync::Arc;

use crate::diagnostics::{Cancellation, Cancelled, FormatError};
use crate::math::{Matrix3, Point3, Vector3};
use crate::molecule::{AtomId, BondOrder, Element, Molecule};
use crate::parallel::Run;
use crate::stereo::{StereoCarrier, StereoElement, StereoElementKind};
use crate::structure::{Structure, StructureEditError, StructureEditor, StructureError};
use crate::topology::edit::EditAtomId;
use crate::topology::{
    AtomSelection, ChainId, InstanceAtomId, InstanceBondId, SelectionError, Topology,
};
use crate::units::{Quantity, ANGSTROM, CANONICAL_LENGTH_UNIT};

use super::super::MmcifBlock;
use super::atom_site::optional;
use super::connectivity::{
    covalent_bonding_radius, within_bonding_distance, COVALENT_BOND_TOLERANCE_ANGSTROM,
};
use super::disulfides::CYSTEINES;
use super::types::{AssemblyOmissionReason, StructureInterpretIssue};
use crate::io::mmcif::category::MmcifCategory;
use crate::io::pdb::document::PdbBiomolecule;

/// The most operator uses the assemblies of one source hold together: an
/// operation that is the product of three operators uses three. Sizes are
/// counted before anything expands, so a damaged expression cannot ask for
/// an unbounded product. Real entries use a few thousand.
pub(super) const MAX_OPERATOR_USES: usize = 100_000;

/// How far a stated operator's rotation part may be from orthonormal, per
/// entry, and its determinant from 1. Sources state matrices to four to six
/// decimals; mirrors and scalings are far further off. An accepted rotation
/// is made exactly orthonormal, so products of many stay rigid.
const ROTATION_TOLERANCE: f64 = 1.0e-3;

/// Whether `rotation` is a proper rotation: RᵀR the identity and the
/// determinant 1, each entry within `tolerance`.
fn is_rotation(rotation: Matrix3, tolerance: f64) -> bool {
    let gram = rotation.transpose().multiply(rotation).columns();
    let identity = Matrix3::identity().columns();
    let orthonormal = gram.iter().zip(&identity).all(|(column, unit)| {
        (column.x - unit.x).abs() <= tolerance
            && (column.y - unit.y).abs() <= tolerance
            && (column.z - unit.z).abs() <= tolerance
    });
    orthonormal && (rotation.determinant() - 1.0).abs() <= tolerance
}

/// The rotation nearest `rotation`, one that [`is_rotation`] accepted: its
/// orthogonal polar factor, by Newton-Schulz steps X <- X(3I - XᵀX)/2. Four
/// steps take a matrix within the tolerance to machine precision, and the
/// identity stays exactly the identity.
fn orthonormalized(mut rotation: Matrix3) -> Matrix3 {
    let identity = Matrix3::identity().columns();
    for _ in 0..4 {
        let gram = rotation.transpose().multiply(rotation).columns();
        let [a, b, c] = [0, 1, 2].map(|i| (identity[i] * 3.0 - gram[i]) / 2.0);
        rotation = rotation.multiply(Matrix3::from_columns(a, b, c));
    }
    rotation
}

/// One biological assembly a source states: copies of its chains, each
/// placed by an operation, such as a virus capsid from 60 point-symmetry
/// operators or a dimer from a crystal's two-fold axis.
///
/// [`Self::to_structure`] builds it from the interpreted structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Assembly {
    pub(crate) id: String,
    pub(crate) details: Option<String>,
    pub(crate) oligomeric_details: Option<String>,
    pub(crate) oligomeric_count: Option<usize>,
    pub(crate) generators: Vec<AssemblyGenerator>,
}

impl Assembly {
    /// The assembly's identifier, such as `1`.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// How the source describes the assembly, such as `author_defined_assembly`
    /// or `complete icosahedral assembly`. A PDB biomolecule is
    /// `author_defined_assembly`, `software_defined_assembly` or
    /// `author_and_software_defined_assembly` by the units `REMARK 350`
    /// states for it.
    pub fn details(&self) -> Option<&str> {
        self.details.as_deref()
    }

    /// The oligomeric state as stated, such as `dimeric` or `60-meric`.
    pub fn oligomeric_details(&self) -> Option<&str> {
        self.oligomeric_details.as_deref()
    }

    /// The number of polymer chains the source says the assembly holds
    /// (mmCIF only).
    pub const fn oligomeric_count(&self) -> Option<usize> {
        self.oligomeric_count
    }

    /// The chains and the operations that place their copies, in source order.
    pub fn generators(&self) -> &[AssemblyGenerator] {
        &self.generators
    }

    /// The number of atoms [`Self::to_structure`] places from `topology`: the
    /// atoms of each generator's chains, once per operation.
    pub fn atom_count(&self, topology: &Topology) -> usize {
        let mut per_chain = BTreeMap::<&str, usize>::new();
        for chain in topology.chains() {
            *per_chain.entry(chain.label_id()).or_default() += chain
                .residues()
                .map(|residue| residue.atom_sites().len())
                .sum::<usize>();
        }
        self.generators
            .iter()
            .map(|generator| {
                let listed = generator
                    .chains
                    .iter()
                    .map(String::as_str)
                    .collect::<HashSet<_>>();
                let atoms: usize = per_chain
                    .iter()
                    .filter(|(chain, _)| listed.contains(*chain))
                    .map(|(_, count)| count)
                    .sum();
                atoms * generator.operations.len()
            })
            .sum()
    }

    /// Builds the assembly from `structure`, the one the report was
    /// interpreted with or one with the same chains and coordinates: for each
    /// generator, a copy of the structure's chains whose label identifiers it
    /// lists, moved by each operation. Listed chains the structure lacks are
    /// skipped.
    ///
    /// A copy keeps the chemistry, perception, hierarchy, occupancies, B
    /// factors and entity properties of its atoms; copies of one molecule share
    /// its definition. A copy under the identity keeps its chains'
    /// identifiers. Any other copy appends its operation's name to them, label
    /// and author alike, on its chains and atom sites, as wwPDB assembly files
    /// name copies: chain `A` under operator `2` becomes `A-2`, under the
    /// product of `12` and `60` `A-12-60`.
    /// Atoms come in generator, operation and source dense order.
    ///
    /// A bond between chains of two generators joins their copies under each
    /// operation both generators apply. A copy whose chain is bonded to a
    /// chain no copy under its operation holds loses the bond, and its atom
    /// with a fixed hydrogen count takes a hydrogen in its place, as a
    /// subset does; [`AssemblyReport::cut_bonds`] lists these. A copy of
    /// only some chains of a molecule loses the stereo elements that name
    /// atoms of its other chains, even when the copies are bonded again;
    /// [`AssemblyReport::lost_stereo`] lists these. Cysteine `SG` atoms of
    /// different copies within bonding distance form disulfides, as reading
    /// forms them within the asymmetric unit: when neither carries a hydrogen
    /// atom or a bond to another residue, nor lies within bonding distance of
    /// a metal it may coordinate, nearest pairs first;
    /// [`AssemblyReport::created_bonds`] lists them. Other links a source
    /// states to a symmetry mate stay import issues. The structure has no
    /// cell, since copies are not a lattice's content.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::NoAtoms`] when no listed chain has atoms in `structure`,
    /// [`AssemblyError::RepeatedCopy`] when two copies of one chain would take
    /// the same identifiers, and the structure errors of taking the chains
    /// and publishing the copies.
    pub fn to_structure(&self, structure: &Structure) -> Result<AssemblyStructure, AssemblyError> {
        self.to_structure_with(structure, Run::default())
    }

    /// [`Self::to_structure`] under a [`Run`], on the calling thread.
    ///
    /// # Errors
    ///
    /// As [`Self::to_structure`], and [`AssemblyError::Cancelled`] once
    /// `run`'s cancellation is requested.
    pub fn to_structure_with(
        &self,
        structure: &Structure,
        run: Run<'_>,
    ) -> Result<AssemblyStructure, AssemblyError> {
        let cancellation = run.cancellation();
        let cancelled = |source| AssemblyError::Cancelled {
            assembly: self.id.clone(),
            source,
        };
        let edit = |source| AssemblyError::Edit {
            assembly: self.id.clone(),
            source,
        };
        let topology = structure.shared_topology();
        let links = ChainLinks::of(&topology);
        let mut editor = StructureEditor::new();
        let mut copies = HashSet::<(ChainId, String)>::new();
        let mut placed = Vec::<PlacedCopy>::new();
        let mut lost_stereo = Vec::new();
        for (generator_index, generator) in self.generators.iter().enumerate() {
            Cancellation::poll(cancellation).map_err(cancelled)?;
            let listed = generator
                .chains
                .iter()
                .map(String::as_str)
                .collect::<HashSet<_>>();
            let chains = topology
                .chains()
                .filter(|chain| listed.contains(chain.label_id()))
                .map(|chain| (chain.id(), chain.label_id()))
                .collect::<Vec<_>>();
            if chains.is_empty() || generator.operations.is_empty() {
                continue;
            }
            let selection =
                AtomSelection::for_chains(Arc::clone(&topology), chains.iter().map(|(id, _)| *id))
                    .map_err(|source| AssemblyError::Selection {
                        assembly: self.id.clone(),
                        source,
                    })?;
            if selection.is_empty() {
                continue;
            }
            let whole = selection.len() == structure.atom_count();
            let part = if whole {
                structure.clone()
            } else {
                structure
                    .subset(&selection)
                    .map_err(|source| AssemblyError::Structure {
                        assembly: self.id.clone(),
                        source,
                    })?
            };
            let lost = if whole {
                Vec::new()
            } else {
                stereo_outside(&topology, &selection, part.topology())
            };
            // The atoms each copy is tracked by: the ends of bonds that leave
            // the generator's chains, and sulfurs that can form disulfides.
            let listed_ids = chains.iter().map(|(id, _)| *id).collect::<HashSet<_>>();
            let tracked = links.tracked(&listed_ids);
            let in_part = if whole {
                tracked.iter().map(|&atom| (atom, atom)).collect()
            } else {
                let tracked = tracked.iter().copied().collect::<HashSet<_>>();
                selection
                    .atom_ids()
                    .zip(part.topology().atom_ids().iter().copied())
                    .filter(|(atom, _)| tracked.contains(atom))
                    .collect::<HashMap<_, _>>()
            };
            for (operation_index, operation) in generator.operations.iter().enumerate() {
                Cancellation::poll(cancellation).map_err(cancelled)?;
                let suffix = (!operation.is_identity()).then(|| operation.name());
                for (chain, label) in &chains {
                    let named = copy_name(label, suffix.as_deref());
                    if !copies.insert((*chain, named.clone())) {
                        return Err(AssemblyError::RepeatedCopy {
                            assembly: self.id.clone(),
                            chain: named,
                        });
                    }
                }
                let mapping = editor
                    .append_moved_copy(
                        &part,
                        operation.rotation,
                        operation.translation,
                        &|label, author| {
                            (
                                copy_name(label, suffix.as_deref()),
                                author.map(|author| copy_name(author, suffix.as_deref())),
                            )
                        },
                    )
                    .map_err(edit)?;
                lost_stereo.extend(lost.iter().map(|atoms| AssemblyStereoLoss {
                    atoms: atoms.clone(),
                    operation: operation.name(),
                }));
                placed.push(PlacedCopy {
                    generator: generator_index,
                    operation: operation_index,
                    chains: listed_ids.clone(),
                    atoms: in_part
                        .iter()
                        .map(|(&atom, part_atom)| (atom, mapping.atoms[part_atom]))
                        .collect(),
                });
            }
        }
        if placed.is_empty() {
            return Err(AssemblyError::NoAtoms {
                assembly: self.id.clone(),
            });
        }
        Cancellation::poll(cancellation).map_err(cancelled)?;
        let cut_bonds = self
            .join_copies(&links, &placed, &mut editor)
            .map_err(edit)?;
        Cancellation::poll(cancellation).map_err(cancelled)?;
        let disulfides = bond_disulfides(&links, &placed, &mut editor).map_err(edit)?;
        let (structure, correspondence) = editor.finish_with_correspondence().map_err(edit)?;
        let created_bonds = disulfides
            .into_iter()
            .map(|atoms| AssemblyBond {
                atoms: atoms.map(|atom| {
                    correspondence
                        .atom(atom)
                        .expect("bonded atoms are published")
                        .0
                }),
            })
            .collect();
        Ok(AssemblyStructure {
            structure,
            report: AssemblyReport {
                cut_bonds,
                lost_stereo,
                created_bonds,
            },
        })
    }

    /// Bonds each copy to the copy under the same operation that holds the
    /// partner of a bond leaving its chains, and returns the bonds no copy
    /// holds a partner for.
    fn join_copies(
        &self,
        links: &ChainLinks,
        placed: &[PlacedCopy],
        editor: &mut StructureEditor,
    ) -> Result<Vec<AssemblyCut>, StructureEditError> {
        let mut cut = Vec::new();
        for copy in placed {
            let operation = &self.generators[copy.generator].operations[copy.operation];
            for link in &links.bonds {
                let (near, far) = match link.chains.map(|chain| copy.chains.contains(&chain)) {
                    [true, false] => (0, 1),
                    [false, true] => (1, 0),
                    _ => continue,
                };
                let partner = placed.iter().find(|other| {
                    other.chains.contains(&link.chains[far])
                        && self.generators[other.generator].operations[other.operation]
                            == *operation
                });
                match partner {
                    // Each bond is made once, from the copy of its first atom.
                    Some(other) if near == 0 => {
                        editor.add_bond(
                            copy.atoms[&link.atoms[0]],
                            other.atoms[&link.atoms[1]],
                            link.order,
                        )?;
                    }
                    Some(_) => {}
                    None => cut.push(AssemblyCut {
                        bond: link.bond,
                        operation: operation.name(),
                    }),
                }
            }
        }
        Ok(cut)
    }
}

/// The focus atoms of the stereo elements of `topology` whose first focus
/// atom `selection` holds and that name an atom outside it, so that `part`,
/// its subset, lost them.
fn stereo_outside(
    topology: &Topology,
    selection: &AtomSelection,
    part: &Topology,
) -> Vec<Vec<InstanceAtomId>> {
    let mut in_part = None::<HashMap<InstanceAtomId, InstanceAtomId>>;
    let mut lost = Vec::new();
    for instance in topology.instances() {
        let molecule = instance.molecule();
        let id = |atom| InstanceAtomId::new(instance.id(), atom);
        for (_, element) in molecule.stereo_elements() {
            let focus = stereo_focus(molecule, element);
            let carriers = stereo_carriers(element);
            let selected =
                |atoms: &[AtomId]| atoms.iter().all(|&atom| selection.contains(id(atom)));
            if !focus
                .first()
                .is_some_and(|&atom| selection.contains(id(atom)))
                || (selected(&focus) && selected(&carriers))
            {
                continue;
            }
            // A subset keeps an element whose cut carrier its capping hands
            // to an implicit hydrogen.
            if selected(&focus) {
                let map = in_part.get_or_insert_with(|| {
                    selection
                        .atom_ids()
                        .zip(part.atom_ids().iter().copied())
                        .collect()
                });
                let mut mapped = focus.iter().map(|&atom| map[&id(atom)]).collect::<Vec<_>>();
                mapped.sort_unstable();
                let kept = part.instance(mapped[0].instance()).is_some_and(|target| {
                    let molecule = target.molecule();
                    molecule.stereo_elements().any(|(_, kept)| {
                        let mut atoms = stereo_focus(molecule, kept)
                            .into_iter()
                            .map(|atom| InstanceAtomId::new(target.id(), atom))
                            .collect::<Vec<_>>();
                        atoms.sort_unstable();
                        atoms == mapped
                    })
                });
                if kept {
                    continue;
                }
            }
            lost.push(focus.into_iter().map(id).collect());
        }
    }
    lost
}

/// A tetrahedral center, or the two atoms of a double bond or axis.
fn stereo_focus(molecule: &Molecule, element: &StereoElement) -> Vec<AtomId> {
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => vec![stereo.center],
        StereoElementKind::DoubleBond(stereo) => vec![stereo.left, stereo.right],
        StereoElementKind::Axis(stereo) => molecule
            .bond(stereo.axis)
            .map(|bond| vec![bond.a(), bond.b()])
            .unwrap_or_default(),
    }
}

/// The atoms an element's carriers name.
fn stereo_carriers(element: &StereoElement) -> Vec<AtomId> {
    let carriers: &[StereoCarrier] = match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => &stereo.carriers,
        StereoElementKind::DoubleBond(stereo) => &[stereo.left_carrier, stereo.right_carrier],
        StereoElementKind::Axis(stereo) => &stereo.carriers,
    };
    carriers
        .iter()
        .filter_map(|carrier| match carrier {
            StereoCarrier::Atom(atom) => Some(*atom),
            _ => None,
        })
        .collect()
}

/// One copy an assembly placed, with the draft handles of its tracked atoms.
struct PlacedCopy {
    generator: usize,
    operation: usize,
    chains: HashSet<ChainId>,
    atoms: HashMap<InstanceAtomId, EditAtomId>,
}

/// A bond of the source between two chains.
struct ChainLink {
    bond: InstanceBondId,
    atoms: [InstanceAtomId; 2],
    chains: [ChainId; 2],
    order: BondOrder,
}

/// What building copies needs of the source's bonds: those between chains,
/// and the cysteine sulfurs free to form a disulfide.
struct ChainLinks {
    bonds: Vec<ChainLink>,
    sulfurs: Vec<(InstanceAtomId, ChainId)>,
}

impl ChainLinks {
    fn of(topology: &Topology) -> Self {
        let chain = |atom: crate::topology::AtomView<'_>| atom.chain().map(|chain| chain.id());
        let bonds = topology
            .bonds()
            .filter_map(|bond| {
                let atoms = bond.atoms();
                let [Some(left), Some(right)] = atoms.map(chain) else {
                    return None;
                };
                (left != right).then(|| ChainLink {
                    bond: bond.id(),
                    atoms: atoms.map(|atom| atom.id()),
                    chains: [left, right],
                    order: bond.bond().order(),
                })
            })
            .collect();
        let sulfur = Element::from_symbol("S").expect("sulfur");
        let sulfurs = topology
            .atoms()
            .filter(|atom| atom.atom().element == sulfur)
            .filter_map(|atom| {
                let site = atom.atom_site()?;
                let residue = site.residue();
                let cysteine = CYSTEINES
                    .iter()
                    .any(|name| residue.name().eq_ignore_ascii_case(name));
                let named = site
                    .metadata()
                    .label_atom_id()
                    .or(site.metadata().auth_atom_id());
                let free = atom.neighbors().all(|neighbor| {
                    neighbor.atom().element.atomic_number() != 1
                        && neighbor.residue().map(|other| other.id()) == Some(residue.id())
                });
                (cysteine && named == Some("SG") && free).then(|| (atom.id(), residue.chain().id()))
            })
            .collect();
        Self { bonds, sulfurs }
    }

    /// The atoms of `chains` a copy of them is tracked by.
    fn tracked(&self, chains: &HashSet<ChainId>) -> Vec<InstanceAtomId> {
        let mut atoms = self
            .bonds
            .iter()
            .flat_map(|link| link.atoms.into_iter().zip(link.chains))
            .chain(self.sulfurs.iter().copied())
            .filter(|(_, chain)| chains.contains(chain))
            .map(|(atom, _)| atom)
            .collect::<Vec<_>>();
        atoms.sort_unstable();
        atoms.dedup();
        atoms
    }
}

/// Bonds the free cysteine sulfurs of different copies that lie within
/// bonding distance, nearest pairs first, each sulfur once, and returns them.
fn bond_disulfides(
    links: &ChainLinks,
    placed: &[PlacedCopy],
    editor: &mut StructureEditor,
) -> Result<Vec<[EditAtomId; 2]>, StructureEditError> {
    let sulfur = Element::from_symbol("S").expect("sulfur");
    let mut sulfurs = Vec::<(usize, EditAtomId, Point3)>::new();
    for (index, copy) in placed.iter().enumerate() {
        for (atom, _) in &links.sulfurs {
            let Some(&handle) = copy.atoms.get(atom) else {
                continue;
            };
            let point = editor
                .position(handle)
                .expect("copied atoms are live")
                .value_in(ANGSTROM)
                .expect("positions are lengths");
            sulfurs.push((index, handle, point));
        }
    }
    // A sulfur within bonding distance of a metal may coordinate it, as a
    // coordination record states for the asymmetric unit.
    let metals = editor
        .topology()
        .atoms()
        .filter(|(_, atom)| atom.element.is_metal())
        .map(|(handle, atom)| {
            let point = editor
                .position(handle)
                .expect("live atoms have positions")
                .value_in(ANGSTROM)
                .expect("positions are lengths");
            (atom.element, point)
        })
        .collect::<Vec<_>>();
    sulfurs.retain(|&(_, _, point)| {
        !metals
            .iter()
            .any(|&metal| within_bonding_distance((sulfur, point), metal))
    });
    // Sorted along x, each sulfur meets only those within bonding reach on x.
    let reach = 2.0 * covalent_bonding_radius(sulfur).0 + COVALENT_BOND_TOLERANCE_ANGSTROM;
    sulfurs.sort_by(|a, b| a.2.x.total_cmp(&b.2.x).then(a.1.cmp(&b.1)));
    let mut pairs = Vec::new();
    for (index, &(left_copy, left, left_point)) in sulfurs.iter().enumerate() {
        for &(right_copy, right, right_point) in &sulfurs[index + 1..] {
            if right_point.x - left_point.x > reach {
                break;
            }
            if left_copy != right_copy
                && within_bonding_distance((sulfur, left_point), (sulfur, right_point))
            {
                pairs.push((
                    (left_point - right_point).norm(),
                    left.min(right),
                    left.max(right),
                ));
            }
        }
    }
    // Ties fall to copy order, so the result does not depend on the sweep.
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    let mut bonded = HashSet::new();
    let mut created = Vec::new();
    for (_, left, right) in pairs {
        if bonded.contains(&left) || bonded.contains(&right) {
            continue;
        }
        bonded.extend([left, right]);
        editor.add_bond(left, right, BondOrder::Single)?;
        created.push([left, right]);
    }
    Ok(created)
}

/// An assembly built from a structure, with the bonds the build changed.
#[derive(Debug, Clone)]
pub struct AssemblyStructure {
    structure: Structure,
    report: AssemblyReport,
}

impl AssemblyStructure {
    /// The assembly.
    pub fn structure(&self) -> &Structure {
        &self.structure
    }

    /// Takes the assembly.
    pub fn into_structure(self) -> Structure {
        self.structure
    }

    /// The bonds the build cut and created.
    pub fn report(&self) -> &AssemblyReport {
        &self.report
    }
}

/// The bonds building an assembly cut and created; see
/// [`Assembly::to_structure`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AssemblyReport {
    cut_bonds: Vec<AssemblyCut>,
    lost_stereo: Vec<AssemblyStereoLoss>,
    created_bonds: Vec<AssemblyBond>,
}

impl AssemblyReport {
    /// The source bonds a copy lost because no copy under its operation
    /// holds the partner, in copy order. Its atom with a fixed hydrogen count
    /// took a hydrogen in the bond's place.
    pub fn cut_bonds(&self) -> &[AssemblyCut] {
        &self.cut_bonds
    }

    /// The source stereo elements a copy lost because they name atoms of
    /// chains outside its generator, in copy order.
    pub fn lost_stereo(&self) -> &[AssemblyStereoLoss] {
        &self.lost_stereo
    }

    /// The disulfides formed between copies, nearest first.
    pub fn created_bonds(&self) -> &[AssemblyBond] {
        &self.created_bonds
    }
}

/// A source bond one copy of an assembly lost.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyCut {
    bond: InstanceBondId,
    operation: String,
}

impl AssemblyCut {
    /// The bond in the source structure.
    pub const fn bond(&self) -> InstanceBondId {
        self.bond
    }

    /// The name of the operation that placed the copy; see
    /// [`AssemblyOperation::name`].
    pub fn operation(&self) -> &str {
        &self.operation
    }
}

/// A source stereo element one copy of an assembly lost.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyStereoLoss {
    atoms: Vec<InstanceAtomId>,
    operation: String,
}

impl AssemblyStereoLoss {
    /// The element's focus in the source structure: a tetrahedral center,
    /// or the two atoms of a double bond or axis.
    pub fn atoms(&self) -> &[InstanceAtomId] {
        &self.atoms
    }

    /// The name of the operation that placed the copy; see
    /// [`AssemblyOperation::name`].
    pub fn operation(&self) -> &str {
        &self.operation
    }
}

/// A bond building an assembly created between copies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssemblyBond {
    atoms: [InstanceAtomId; 2],
}

impl AssemblyBond {
    /// The bonded atoms in the assembly.
    pub const fn atoms(&self) -> [InstanceAtomId; 2] {
        self.atoms
    }
}

/// A chain identifier with an operation's name appended, or as it is under
/// the identity.
fn copy_name(id: &str, suffix: Option<&str>) -> String {
    match suffix {
        Some(suffix) => format!("{id}-{suffix}"),
        None => id.to_owned(),
    }
}

/// The chains one part of an assembly copies and the operations that place
/// the copies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssemblyGenerator {
    pub(crate) chains: Vec<String>,
    pub(crate) operations: Vec<AssemblyOperation>,
}

impl AssemblyGenerator {
    /// The chains, by label identifier: mmCIF `label_asym_id`, or the PDB
    /// chain identifiers `REMARK 350` lists, which PDB chains take as labels.
    pub fn chains(&self) -> &[String] {
        &self.chains
    }

    /// The operations, in the order the source states them; each places one
    /// copy of the chains.
    pub fn operations(&self) -> &[AssemblyOperation] {
        &self.operations
    }
}

/// One operation of an assembly: a rotation, then a translation, from one
/// stated operator or the product of several.
#[derive(Debug, Clone)]
pub struct AssemblyOperation {
    pub(crate) operators: Vec<String>,
    pub(crate) rotation: Matrix3,
    /// In the canonical length unit.
    pub(crate) translation: Vector3,
}

impl AssemblyOperation {
    /// The identifiers of the operators it combines, in the order the source
    /// writes them; the rightmost applies first.
    pub fn operators(&self) -> &[String] {
        &self.operators
    }

    /// The operators' identifiers joined by dashes, such as `2` or `12-60`.
    pub fn name(&self) -> String {
        self.operators.join("-")
    }

    /// The rotation, applied before the translation. Each stated operator's
    /// matrix is made exactly orthonormal.
    pub const fn rotation(&self) -> Matrix3 {
        self.rotation
    }

    /// The translation.
    pub fn translation(&self) -> Quantity<Vector3> {
        Quantity::new(self.translation, CANONICAL_LENGTH_UNIT)
    }

    /// Whether it leaves every point where it is: an identity rotation and
    /// no translation, exactly as stated.
    pub fn is_identity(&self) -> bool {
        self.rotation == Matrix3::identity() && self.translation == Vector3::zero()
    }

    /// Whether it is finite and its rotation proper, within the tolerance
    /// reading allows a stated operator.
    pub(crate) fn is_rigid(&self) -> bool {
        let columns = self.rotation.columns();
        columns
            .iter()
            .chain([&self.translation])
            .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite())
            && is_rotation(self.rotation, ROTATION_TOLERANCE)
    }

    /// Applies `self` after the rotation and translation given, in place.
    fn apply_after(&self, rotation: &mut Matrix3, translation: &mut Vector3) {
        let moved = self.rotation.transform_vector(*translation);
        *rotation = self.rotation.multiply(*rotation);
        *translation = Vector3::new(
            moved.x + self.translation.x,
            moved.y + self.translation.y,
            moved.z + self.translation.z,
        );
    }
}

/// Operations compare entry by entry, bit for bit, so an operation equals
/// itself and its decoded copy.
impl PartialEq for AssemblyOperation {
    fn eq(&self, other: &Self) -> bool {
        let bits = |rotation: Matrix3, translation: Vector3| {
            let [a, b, c] = rotation.columns();
            [a, b, c, translation].map(|v| [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()])
        };
        self.operators == other.operators
            && bits(self.rotation, self.translation) == bits(other.rotation, other.translation)
    }
}

impl Eq for AssemblyOperation {}

/// Why an assembly could not be built from a structure.
#[non_exhaustive]
#[derive(Debug)]
pub enum AssemblyError {
    /// None of the chains the assembly lists has atoms in the structure, as
    /// when the structure is not the one the report was read with.
    NoAtoms {
        /// The assembly.
        assembly: String,
    },
    /// Two copies of one chain would take the same identifiers: the source
    /// lists the chain twice under one operation, or under two identity
    /// operations.
    RepeatedCopy {
        /// The assembly.
        assembly: String,
        /// The repeated chain identifier.
        chain: String,
    },
    /// The cancellation stopped the build.
    Cancelled {
        /// The assembly.
        assembly: String,
        /// The cancellation.
        source: Cancelled,
    },
    /// The structure's chains could not be selected.
    Selection {
        /// The assembly.
        assembly: String,
        /// The selection failure.
        source: SelectionError,
    },
    /// The structure's chains could not be taken apart from the rest.
    Structure {
        /// The assembly.
        assembly: String,
        /// The subset failure.
        source: StructureError,
    },
    /// The copies could not be appended or published.
    Edit {
        /// The assembly.
        assembly: String,
        /// The edit failure.
        source: StructureEditError,
    },
}

impl fmt::Display for AssemblyError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoAtoms { assembly } => write!(
                formatter,
                "Assembly {assembly:?} copies no atoms. The structure has none of its chains."
            ),
            Self::RepeatedCopy { assembly, chain } => write!(
                formatter,
                "Assembly {assembly:?} places two copies named chain {chain:?}."
            ),
            Self::Cancelled { assembly, .. }
            | Self::Selection { assembly, .. }
            | Self::Structure { assembly, .. }
            | Self::Edit { assembly, .. } => {
                write!(formatter, "Cannot build assembly {assembly:?}.")
            }
        }
    }
}

impl std::error::Error for AssemblyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::NoAtoms { .. } | Self::RepeatedCopy { .. } => None,
            Self::Cancelled { source, .. } => Some(source),
            Self::Selection { source, .. } => Some(source),
            Self::Structure { source, .. } => Some(source),
            Self::Edit { source, .. } => Some(source),
        }
    }
}

/// A stated operator's entries in ångströms: the matrix by rows and the
/// vector, or `None` when any is missing or not a number.
pub(super) type OperatorEntries = Option<([[f64; 3]; 3], [f64; 3])>;

/// The operation one stated operator gives, checked to be a rotation and
/// made exactly one.
pub(super) fn operator(
    id: &str,
    entries: OperatorEntries,
) -> Result<AssemblyOperation, AssemblyOmissionReason> {
    let Some((matrix, vector)) = entries.filter(|(matrix, vector)| {
        matrix
            .iter()
            .flatten()
            .chain(vector)
            .all(|value| value.is_finite())
    }) else {
        return Err(AssemblyOmissionReason::UnreadableOperator {
            operator: id.to_owned(),
        });
    };
    let rotation = Matrix3::from_columns(
        Vector3::new(matrix[0][0], matrix[1][0], matrix[2][0]),
        Vector3::new(matrix[0][1], matrix[1][1], matrix[2][1]),
        Vector3::new(matrix[0][2], matrix[1][2], matrix[2][2]),
    );
    if !is_rotation(rotation, ROTATION_TOLERANCE) {
        return Err(AssemblyOmissionReason::NotRotation {
            operator: id.to_owned(),
        });
    }
    let scale = ANGSTROM
        .conversion_factor_to(CANONICAL_LENGTH_UNIT)
        .expect("lengths convert");
    Ok(AssemblyOperation {
        operators: vec![id.to_owned()],
        rotation: orthonormalized(rotation),
        translation: Vector3::new(vector[0], vector[1], vector[2]) * scale,
    })
}

/// One group of an operator expression, its numeric ranges not yet
/// expanded.
#[derive(Debug)]
struct Group<'a> {
    /// Identifiers, and ranges such as `1-60` as their first and last.
    items: Vec<Result<&'a str, (usize, usize)>>,
    /// The identifiers it names.
    len: usize,
}

impl Group<'_> {
    /// The identifiers, ranges expanded.
    fn ids(&self) -> Vec<String> {
        let mut ids = Vec::with_capacity(self.len);
        for item in &self.items {
            match *item {
                Ok(id) => ids.push(id.to_owned()),
                Err((first, last)) => ids.extend((first..=last).map(|n| n.to_string())),
            }
        }
        ids
    }
}

/// The groups of operator identifiers an mmCIF operator expression names,
/// left to right: `1` or `1,2,5` is one group, `(1-60)` one group of 60, and
/// `(1-60)(61-88)` two groups whose product it names. `None` when the
/// expression cannot be read.
fn expression_groups(expression: &str) -> Option<Vec<Group<'_>>> {
    let text = expression.trim();
    if text.is_empty() {
        return None;
    }
    if !text.contains(['(', ')']) {
        return group_items(text).map(|group| vec![group]);
    }
    let mut groups = Vec::new();
    let mut rest = text;
    loop {
        rest = rest.trim_start();
        if rest.is_empty() {
            return Some(groups);
        }
        let (group, after) = rest.strip_prefix('(')?.split_once(')')?;
        groups.push(group_items(group)?);
        rest = after;
    }
}

/// One comma-separated group; numeric ranges such as `1-60` stay
/// unexpanded.
fn group_items(group: &str) -> Option<Group<'_>> {
    let mut items = Vec::new();
    let mut len = 0_usize;
    for item in group.split(',') {
        let item = item.trim();
        if item.is_empty() || item.contains(|c: char| c.is_whitespace() || "()".contains(c)) {
            return None;
        }
        let range = item.split_once('-').and_then(|(first, last)| {
            Some((first.parse::<usize>().ok()?, last.parse::<usize>().ok()?))
        });
        match range {
            Some((first, last)) => {
                if first > last {
                    return None;
                }
                len = len.checked_add(last - first)?.checked_add(1)?;
                items.push(Err((first, last)));
            }
            None => {
                len = len.checked_add(1)?;
                items.push(Ok(item));
            }
        }
    }
    Some(Group { items, len })
}

/// The operations an expression names, as products of the operators in
/// `operators`: the rightmost group's operator applies first. `budget` is
/// the operator uses left for the source, and they are taken from it.
fn expression_operations(
    expression: &str,
    operators: &BTreeMap<String, OperatorEntries>,
    budget: &mut usize,
) -> Result<Vec<AssemblyOperation>, AssemblyOmissionReason> {
    let groups =
        expression_groups(expression).ok_or_else(|| AssemblyOmissionReason::Expression {
            expression: expression.to_owned(),
        })?;
    // Every group names at least one identifier, so the expanded groups hold
    // no more identifiers than the uses counted here.
    let count = groups
        .iter()
        .try_fold(1_usize, |count, group| count.checked_mul(group.len));
    let uses = count.and_then(|count| count.checked_mul(groups.len()));
    let (Some(count), Some(uses)) = (count, uses) else {
        return Err(AssemblyOmissionReason::TooManyOperations);
    };
    if uses > *budget {
        return Err(AssemblyOmissionReason::TooManyOperations);
    }
    *budget -= uses;
    let groups = groups.iter().map(Group::ids).collect::<Vec<_>>();
    let mut resolved = BTreeMap::<&str, AssemblyOperation>::new();
    for id in groups.iter().flatten() {
        if resolved.contains_key(id.as_str()) {
            continue;
        }
        let entries = operators
            .get(id)
            .ok_or_else(|| AssemblyOmissionReason::UnknownOperator {
                operator: id.clone(),
            })?;
        resolved.insert(id, operator(id, *entries)?);
    }
    let mut operations = Vec::with_capacity(count);
    let mut picks = vec![0_usize; groups.len()];
    for _ in 0..count {
        let mut picked = groups
            .iter()
            .zip(&picks)
            .map(|(group, &pick)| &resolved[group[pick].as_str()])
            .rev();
        let first = picked.next().expect("an expression has a group");
        let (mut rotation, mut translation) = (first.rotation, first.translation);
        for outer in picked {
            outer.apply_after(&mut rotation, &mut translation);
        }
        operations.push(AssemblyOperation {
            operators: groups
                .iter()
                .zip(&picks)
                .map(|(group, &pick)| group[pick].clone())
                .collect(),
            rotation,
            translation,
        });
        // The next pick, the last group fastest.
        for (pick, group) in picks.iter_mut().zip(&groups).rev() {
            *pick += 1;
            if *pick < group.len() {
                break;
            }
            *pick = 0;
        }
    }
    Ok(operations)
}

/// The assemblies an mmCIF block states, in `_pdbx_struct_assembly` order,
/// then those only `_pdbx_struct_assembly_gen` names. Those that cannot be
/// read are left out with an issue.
pub(super) fn read_mmcif_assemblies(
    block: &MmcifBlock,
    issues: &mut Vec<StructureInterpretIssue>,
) -> Result<Vec<Assembly>, FormatError> {
    struct Stated {
        id: Option<String>,
        details: Option<String>,
        oligomeric_details: Option<String>,
        oligomeric_count: Option<usize>,
        line: Option<usize>,
        generators: Vec<(Vec<String>, String)>,
    }
    let text = |table: &MmcifCategory<'_>, row: usize, tag: &str| {
        optional(table, row, tag).map(str::to_owned)
    };
    let mut stated = Vec::<Stated>::new();
    if let Some(table) = block.category("_pdbx_struct_assembly")? {
        for row in 0..table.row_count() {
            stated.push(Stated {
                id: text(&table, row, "_pdbx_struct_assembly.id"),
                details: text(&table, row, "_pdbx_struct_assembly.details"),
                oligomeric_details: text(&table, row, "_pdbx_struct_assembly.oligomeric_details"),
                oligomeric_count: optional(&table, row, "_pdbx_struct_assembly.oligomeric_count")
                    .and_then(|count| count.parse().ok()),
                line: table.row_line(row),
                generators: Vec::new(),
            });
        }
    }
    if let Some(table) = block.category("_pdbx_struct_assembly_gen")? {
        for row in 0..table.row_count() {
            let id = text(&table, row, "_pdbx_struct_assembly_gen.assembly_id");
            let chains = optional(&table, row, "_pdbx_struct_assembly_gen.asym_id_list")
                .map(|list| {
                    list.split(',')
                        .map(str::trim)
                        .filter(|chain| !chain.is_empty())
                        .map(str::to_owned)
                        .collect()
                })
                .unwrap_or_default();
            let expression =
                text(&table, row, "_pdbx_struct_assembly_gen.oper_expression").unwrap_or_default();
            let index = match stated
                .iter()
                .position(|assembly| assembly.id.is_some() && assembly.id == id)
            {
                Some(index) => index,
                None => {
                    stated.push(Stated {
                        id,
                        details: None,
                        oligomeric_details: None,
                        oligomeric_count: None,
                        line: table.row_line(row),
                        generators: Vec::new(),
                    });
                    stated.len() - 1
                }
            };
            stated[index].generators.push((chains, expression));
        }
    }
    if stated.is_empty() {
        return Ok(Vec::new());
    }
    let mut operators = BTreeMap::<String, OperatorEntries>::new();
    if let Some(table) = block.category("_pdbx_struct_oper_list")? {
        for row in 0..table.row_count() {
            let Some(id) = optional(&table, row, "_pdbx_struct_oper_list.id") else {
                continue;
            };
            let number = |tag: String| optional(&table, row, &tag)?.parse::<f64>().ok();
            let mut matrix = [[0.0; 3]; 3];
            let mut vector = [0.0; 3];
            let mut complete = true;
            for (i, (row, shift)) in matrix.iter_mut().zip(&mut vector).enumerate() {
                for (j, entry) in row.iter_mut().enumerate() {
                    match number(format!(
                        "_pdbx_struct_oper_list.matrix[{}][{}]",
                        i + 1,
                        j + 1
                    )) {
                        Some(value) => *entry = value,
                        None => complete = false,
                    }
                }
                match number(format!("_pdbx_struct_oper_list.vector[{}]", i + 1)) {
                    Some(value) => *shift = value,
                    None => complete = false,
                }
            }
            operators.insert(id.to_owned(), complete.then_some((matrix, vector)));
        }
    }
    let mut assemblies = Vec::new();
    let mut budget = MAX_OPERATOR_USES;
    for assembly in stated {
        let omitted = |reason| StructureInterpretIssue::AssemblyOmitted {
            assembly: assembly.id.clone(),
            reason,
            line: assembly.line,
        };
        let Some(id) = assembly.id.clone() else {
            issues.push(omitted(AssemblyOmissionReason::Unidentified));
            continue;
        };
        if assembly.generators.is_empty()
            || assembly
                .generators
                .iter()
                .any(|(chains, _)| chains.is_empty())
        {
            issues.push(omitted(AssemblyOmissionReason::NoGenerators));
            continue;
        }
        let generators = assembly
            .generators
            .iter()
            .map(|(chains, expression)| {
                Ok(AssemblyGenerator {
                    chains: chains.clone(),
                    operations: expression_operations(expression, &operators, &mut budget)?,
                })
            })
            .collect::<Result<Vec<_>, AssemblyOmissionReason>>();
        match generators {
            Ok(generators) => assemblies.push(Assembly {
                id,
                details: assembly.details,
                oligomeric_details: assembly.oligomeric_details,
                oligomeric_count: assembly.oligomeric_count,
                generators,
            }),
            Err(reason) => issues.push(omitted(reason)),
        }
    }
    Ok(assemblies)
}

/// The assemblies of a PDB entry's `REMARK 350` biomolecules, in source
/// order. A biomolecule is `author_defined_assembly` when it states an
/// author-determined unit, `software_defined_assembly` when it states a
/// software-determined one, and `author_and_software_defined_assembly` when
/// both; its oligomeric details are the author's unit, else the software's.
/// Those that cannot be read are left out with an issue.
pub(super) fn read_pdb_assemblies(
    biomolecules: &[PdbBiomolecule],
    issues: &mut Vec<StructureInterpretIssue>,
) -> Vec<Assembly> {
    let mut assemblies = Vec::new();
    let mut budget = MAX_OPERATOR_USES;
    for biomolecule in biomolecules {
        let id = Some(biomolecule.id.clone()).filter(|id| !id.is_empty());
        let omitted = |reason| StructureInterpretIssue::AssemblyOmitted {
            assembly: id.clone(),
            reason,
            line: Some(biomolecule.location.line()),
        };
        let Some(id) = id.clone() else {
            issues.push(omitted(AssemblyOmissionReason::Unidentified));
            continue;
        };
        if biomolecule.generators.is_empty()
            || biomolecule
                .generators
                .iter()
                .any(|generator| generator.chains.is_empty() || generator.operators.is_empty())
        {
            issues.push(omitted(AssemblyOmissionReason::NoGenerators));
            continue;
        }
        let uses = biomolecule
            .generators
            .iter()
            .map(|generator| generator.operators.len())
            .sum::<usize>();
        if uses > budget {
            issues.push(omitted(AssemblyOmissionReason::TooManyOperations));
            continue;
        }
        budget -= uses;
        let generators = biomolecule
            .generators
            .iter()
            .map(|generator| {
                let operations = generator
                    .operators
                    .iter()
                    .map(|(serial, rows)| {
                        let entries = match rows {
                            [Some(x), Some(y), Some(z)] => Some((
                                [[x[0], x[1], x[2]], [y[0], y[1], y[2]], [z[0], z[1], z[2]]],
                                [x[3], y[3], z[3]],
                            )),
                            _ => None,
                        };
                        operator(serial, entries)
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(AssemblyGenerator {
                    chains: generator.chains.clone(),
                    operations,
                })
            })
            .collect::<Result<Vec<_>, AssemblyOmissionReason>>();
        let details = match (&biomolecule.author_unit, &biomolecule.software_unit) {
            (Some(_), Some(_)) => Some("author_and_software_defined_assembly"),
            (Some(_), None) => Some("author_defined_assembly"),
            (None, Some(_)) => Some("software_defined_assembly"),
            (None, None) => None,
        };
        match generators {
            Ok(generators) => assemblies.push(Assembly {
                id,
                details: details.map(str::to_owned),
                oligomeric_details: biomolecule
                    .author_unit
                    .clone()
                    .or_else(|| biomolecule.software_unit.clone()),
                oligomeric_count: None,
                generators,
            }),
            Err(reason) => issues.push(omitted(reason)),
        }
    }
    assemblies
}

#[cfg(test)]
mod tests;
