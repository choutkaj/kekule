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
use std::collections::{BTreeMap, HashSet};
use std::fmt;
use std::sync::Arc;

use crate::diagnostics::{Cancellation, Cancelled, FormatError};
use crate::math::{Matrix3, Vector3};
use crate::structure::{Structure, StructureEditError, StructureEditor, StructureError};
use crate::topology::{AtomSelection, ChainId, SelectionError, Topology};
use crate::units::{Quantity, ANGSTROM, CANONICAL_LENGTH_UNIT};

use super::super::MmcifBlock;
use super::atom_site::optional;
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
    /// Atoms come in generator, operation and source dense order. A generator's
    /// chains that are bonded to chains it does not list lose those bonds. The
    /// structure has no cell, since copies are not a lattice's content.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::NoAtoms`] when no listed chain has atoms in `structure`,
    /// [`AssemblyError::RepeatedCopy`] when two copies of one chain would take
    /// the same identifiers, [`AssemblyError::Cancelled`] when `cancellation`
    /// stops it, and the structure errors of taking the chains and publishing
    /// the copies.
    pub fn to_structure(
        &self,
        structure: &Structure,
        cancellation: Option<&Cancellation>,
    ) -> Result<Structure, AssemblyError> {
        let cancelled = |source| AssemblyError::Cancelled {
            assembly: self.id.clone(),
            source,
        };
        let topology = structure.shared_topology();
        let mut editor = StructureEditor::new();
        let mut copies = HashSet::<(ChainId, String)>::new();
        let mut copied = false;
        for generator in &self.generators {
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
                        source: Box::new(source),
                    })?;
            if selection.is_empty() {
                continue;
            }
            let part = if selection.len() == structure.atom_count() {
                structure.clone()
            } else {
                structure
                    .subset(&selection)
                    .map_err(|source| AssemblyError::Structure {
                        assembly: self.id.clone(),
                        source: Box::new(source),
                    })?
            };
            for operation in &generator.operations {
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
                editor
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
                    .map_err(|source| AssemblyError::Edit {
                        assembly: self.id.clone(),
                        source: Box::new(source),
                    })?;
                copied = true;
            }
        }
        if !copied {
            return Err(AssemblyError::NoAtoms {
                assembly: self.id.clone(),
            });
        }
        Cancellation::poll(cancellation).map_err(cancelled)?;
        editor.finish().map_err(|source| AssemblyError::Edit {
            assembly: self.id.clone(),
            source: Box::new(source),
        })
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
        source: Box<SelectionError>,
    },
    /// The structure's chains could not be taken apart from the rest.
    Structure {
        /// The assembly.
        assembly: String,
        /// The subset failure.
        source: Box<StructureError>,
    },
    /// The copies could not be appended or published.
    Edit {
        /// The assembly.
        assembly: String,
        /// The edit failure.
        source: Box<StructureEditError>,
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
            Self::Selection { source, .. } => Some(source.as_ref()),
            Self::Structure { source, .. } => Some(source.as_ref()),
            Self::Edit { source, .. } => Some(source.as_ref()),
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
mod tests {
    use super::*;

    fn entries(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> OperatorEntries {
        Some((matrix, vector))
    }

    const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    /// A quarter turn about z.
    const TURN: [[f64; 3]; 3] = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];

    #[test]
    fn oversized_expressions_are_refused_before_they_expand() {
        let operators = BTreeMap::from([("1".to_owned(), entries(IDENTITY, [0.0; 3]))]);
        let refused = |expression: &str, budget: usize| {
            let mut left = budget;
            let result = expression_operations(expression, &operators, &mut left);
            (result.err(), left)
        };
        // A product of ranges that would be 10^10 operations, and a list of a
        // hundred million identifiers, refused by their counts alone.
        let too_many = Some(AssemblyOmissionReason::TooManyOperations);
        assert_eq!(
            refused("(1-99999)(1-99999)", MAX_OPERATOR_USES),
            (too_many.clone(), MAX_OPERATOR_USES)
        );
        let list = "1-99999,".repeat(1_000);
        assert_eq!(
            refused(list.trim_end_matches(','), MAX_OPERATOR_USES).0,
            too_many
        );
        // One operation of 30,000 operators uses 30,000.
        let long = "(1)".repeat(30_000);
        assert_eq!(refused(&long, 29_999).0, too_many);
        assert_eq!(refused(&long, 30_000), (None, 0));
        // The budget is shared: what one expression takes, the next lacks.
        let mut budget = 5;
        assert!(expression_operations("(1,1,1)", &operators, &mut budget).is_ok());
        assert_eq!(budget, 2);
        assert_eq!(
            expression_operations("(1,1,1)", &operators, &mut budget).err(),
            too_many
        );
    }

    #[test]
    fn expressions_name_lists_ranges_and_products() {
        let names = |expression| {
            expression_groups(expression).map(|groups| {
                groups
                    .into_iter()
                    .map(|group| group.ids().join(" "))
                    .collect::<Vec<_>>()
            })
        };
        assert_eq!(names("1"), Some(vec!["1".to_owned()]));
        assert_eq!(names("1,2,5"), Some(vec!["1 2 5".to_owned()]));
        assert_eq!(names("(1-4)"), Some(vec!["1 2 3 4".to_owned()]));
        assert_eq!(
            names("(X0)(1-2)"),
            Some(vec!["X0".to_owned(), "1 2".to_owned()])
        );
        assert_eq!(
            names("(1,2)(3, 4)"),
            Some(vec!["1 2".to_owned(), "3 4".to_owned()])
        );
        for unreadable in ["", "(1", "1)", "(1)x", "1,,2", "(4-1)", "1 2"] {
            assert_eq!(names(unreadable), None, "{unreadable:?}");
        }
    }

    #[test]
    fn products_apply_the_rightmost_operator_first() {
        let operators = BTreeMap::from([
            ("1".to_owned(), entries(IDENTITY, [0.0; 3])),
            ("T".to_owned(), entries(IDENTITY, [10.0, 0.0, 0.0])),
            ("R".to_owned(), entries(TURN, [0.0; 3])),
        ]);
        let operations =
            expression_operations("(R)(T)", &operators, &mut MAX_OPERATOR_USES.clone())
                .expect("readable");
        assert_eq!(operations.len(), 1);
        let operation = &operations[0];
        assert_eq!(operation.operators(), ["R", "T"]);
        assert_eq!(operation.name(), "R-T");
        // Translated along x first, then turned onto y.
        let moved = operation
            .translation()
            .into_unit(ANGSTROM)
            .expect("a length")
            .into_value();
        assert!((moved.x - 0.0).abs() < 1e-12 && (moved.y - 10.0).abs() < 1e-12);
        let product =
            expression_operations("(1,R)(1,T)", &operators, &mut MAX_OPERATOR_USES.clone())
                .expect("readable");
        let names = product
            .iter()
            .map(AssemblyOperation::name)
            .collect::<Vec<_>>();
        assert_eq!(names, ["1-1", "1-T", "R-1", "R-T"]);
        assert!(product[0].is_identity() && !product[1].is_identity());
    }

    #[test]
    fn rounded_rotations_compose_without_drift() {
        // Within the tolerance, yet the product of 1,000 as stated would
        // stretch x by half.
        let stretched = [[1.0004, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let operators = BTreeMap::from([("S".to_owned(), entries(stretched, [0.0; 3]))]);
        let product = expression_operations(
            &"(S)".repeat(1000),
            &operators,
            &mut MAX_OPERATOR_USES.clone(),
        )
        .unwrap();
        assert_eq!(product.len(), 1);
        assert!(is_rotation(product[0].rotation, 1.0e-12));
        assert!(product[0].is_rigid());
        assert!(operator("I", entries(IDENTITY, [0.0; 3]))
            .unwrap()
            .is_identity());
    }

    #[test]
    fn operators_that_are_not_rotations_are_refused() {
        let mirror = [[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let scaled = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for matrix in [mirror, scaled] {
            assert_eq!(
                operator("2", entries(matrix, [0.0; 3])),
                Err(AssemblyOmissionReason::NotRotation {
                    operator: "2".to_owned()
                })
            );
        }
        assert_eq!(
            operator("3", None),
            Err(AssemblyOmissionReason::UnreadableOperator {
                operator: "3".to_owned()
            })
        );
        // Six decimals, as files state them.
        let rounded = [
            [0.309017, -0.809017, 0.5],
            [0.809017, 0.5, 0.309017],
            [-0.5, 0.309017, 0.809017],
        ];
        let operation = operator("4", entries(rounded, [1.0, 2.0, 3.0])).unwrap();
        assert!(is_rotation(operation.rotation, 1.0e-12));
        let operators = BTreeMap::from([("1".to_owned(), entries(IDENTITY, [0.0; 3]))]);
        assert_eq!(
            expression_operations("(1-2)", &operators, &mut MAX_OPERATOR_USES.clone()),
            Err(AssemblyOmissionReason::UnknownOperator {
                operator: "2".to_owned()
            })
        );
    }
}
