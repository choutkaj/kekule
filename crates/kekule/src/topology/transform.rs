//! Immutable topology transformations.
//!
//! Whole-instance filters preserve definitions and instances. [`Topology::subset`]
//! is the separate hierarchy-aware induced-graph operation and returns its own
//! narrow dense-state correspondence.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;
use std::sync::Arc;

use super::{
    AtomSelection, HierarchyError, InstanceAtomId, InstanceBondId, MoleculeInstanceId,
    SelectionError, Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuildError,
    TopologyBuilder,
};
use crate::core::{MoleculeError, MoleculePublicationError};

/// Retains complete molecule instances in source topology order.
///
/// Duplicate identifiers are treated as one membership request. A request
/// containing every source instance preserves the source `Arc<Topology>`. Empty
/// results and invalid identifiers are rejected before target construction.
///
/// # Examples
///
/// ```
/// use kekule::core::{Atom, Element, MoleculeEditor};
/// use kekule::topology::{transform, TopologyBuilder};
/// use std::sync::Arc;
///
/// let mut water_builder = MoleculeEditor::new();
/// water_builder.add_atom(Atom::new(Element::from_symbol("O").unwrap()))?;
/// let water = water_builder.finish()?;
///
/// let mut builder = TopologyBuilder::new();
/// let definition = builder.add_molecule_definition(&water)?;
/// let first = builder.add_instance(definition)?;
/// builder.add_instance(definition)?;
/// let source = Arc::new(builder.build()?);
///
/// let target = transform::retain_instances(&source, [first])?;
/// assert_eq!(target.instance_count(), 1);
/// assert_eq!(target.definition_count(), 1);
/// assert!(!Arc::ptr_eq(&source, &target));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn retain_instances(
    topology: &Arc<Topology>,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<Arc<Topology>, TopologyTransformError> {
    let retained = validate_instances(topology, instances)?;
    retain_normalized(topology, &retained)
}

/// Removes complete molecule instances while preserving filtered source order.
///
/// Duplicate identifiers are harmless. Removing no instances preserves the
/// source `Arc<Topology>`; removing every instance is rejected.
pub fn remove_instances(
    topology: &Arc<Topology>,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<Arc<Topology>, TopologyTransformError> {
    let removed = validate_instances(topology, instances)?;
    let retained = InstanceMembership {
        members: removed
            .members
            .into_iter()
            .map(|removed| !removed)
            .collect(),
        len: topology.instance_count() - removed.len,
    };
    retain_normalized(topology, &retained)
}

struct InstanceMembership {
    members: Vec<bool>,
    len: usize,
}

impl InstanceMembership {
    fn contains(&self, instance: MoleculeInstanceId) -> bool {
        self.members[instance.index()]
    }
}

fn validate_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<InstanceMembership, TopologyTransformError> {
    let mut normalized = InstanceMembership {
        members: vec![false; topology.instance_count()],
        len: 0,
    };
    for instance in instances {
        if topology.instance(instance).is_err() {
            return Err(TopologyTransformError::InvalidSourceInstance(instance));
        }
        if !normalized.members[instance.index()] {
            normalized.members[instance.index()] = true;
            normalized.len += 1;
        }
    }
    Ok(normalized)
}

fn retain_normalized(
    topology: &Arc<Topology>,
    retained: &InstanceMembership,
) -> Result<Arc<Topology>, TopologyTransformError> {
    if retained.len == 0 {
        return Err(TopologyTransformError::EmptyTargetTopology);
    }
    if retained.len == topology.instance_count() {
        return Ok(Arc::clone(topology));
    }

    let mut referenced_definitions = vec![false; topology.definition_count()];
    for (instance_id, instance) in topology.instances() {
        if retained.contains(instance_id) {
            referenced_definitions[instance.definition().index()] = true;
        }
    }
    let retained_definition_count = referenced_definitions
        .iter()
        .filter(|referenced| **referenced)
        .count();

    let mut builder = TopologyBuilder::new();
    builder.reserve_definitions(retained_definition_count)?;
    builder.reserve_instances(retained.len)?;

    let mut definition_targets = vec![None; topology.definition_count()];
    for (source_id, definition) in topology
        .definitions()
        .filter(|(id, _)| referenced_definitions[id.index()])
    {
        let target_id = builder.add_molecule_definition(definition.molecule())?;
        definition_targets[source_id.index()] = Some(target_id);
    }

    let mut atom_targets = BTreeMap::new();
    for (source_instance, instance) in topology
        .instances()
        .filter(|(id, _)| retained.contains(*id))
    {
        let target_definition = definition_targets[instance.definition().index()]
            .expect("retained instance has a retained definition");
        let target_instance = builder.add_instance(target_definition)?;
        for atom in topology
            .definition_for_instance(source_instance)
            .expect("retained instance references a live definition")
            .molecule()
            .atom_ids()
        {
            atom_targets.insert(
                InstanceAtomId::new(source_instance, atom),
                InstanceAtomId::new(target_instance, atom),
            );
        }
    }

    copy_filtered_hierarchy(topology, builder.hierarchy_mut(), &atom_targets)?;

    Ok(Arc::new(builder.build()?))
}

/// Failure to construct an immutable whole-instance topology subset.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyTransformError {
    /// A requested instance does not exist in the source topology.
    InvalidSourceInstance(MoleculeInstanceId),
    /// The requested membership would produce a topology with no instances.
    EmptyTargetTopology,
    /// The filtered target topology could not be constructed.
    TopologyBuild(TopologyBuildError),
    /// The retained hierarchy could not be reconstructed.
    Hierarchy(HierarchyError),
}

impl fmt::Display for TopologyTransformError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidSourceInstance(instance) => {
                write!(formatter, "invalid source molecule instance: {instance}")
            }
            Self::EmptyTargetTopology => {
                formatter.write_str("topology transformation would remove every instance")
            }
            Self::TopologyBuild(error) => {
                write!(formatter, "cannot build target topology: {error}")
            }
            Self::Hierarchy(error) => {
                write!(formatter, "cannot retain topology hierarchy: {error}")
            }
        }
    }
}

impl std::error::Error for TopologyTransformError {}

impl From<TopologyBuildError> for TopologyTransformError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}

impl From<HierarchyError> for TopologyTransformError {
    fn from(error: HierarchyError) -> Self {
        Self::Hierarchy(error)
    }
}

/// Result of hierarchy-aware structural subsetting.
#[derive(Debug, Clone)]
pub struct TopologySubset {
    topology: Arc<Topology>,
    correspondence: TopologySubsetCorrespondence,
}

impl TopologySubset {
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub const fn correspondence(&self) -> &TopologySubsetCorrespondence {
        &self.correspondence
    }

    pub fn into_parts(self) -> (Arc<Topology>, TopologySubsetCorrespondence) {
        (self.topology, self.correspondence)
    }
}

/// Source-to-target correspondence specific to [`Topology::subset`].
#[derive(Debug, Clone)]
pub struct TopologySubsetCorrespondence {
    source: Arc<Topology>,
    target: Arc<Topology>,
    atom_targets: BTreeMap<InstanceAtomId, InstanceAtomId>,
    bond_targets: BTreeMap<InstanceBondId, InstanceBondId>,
    source_atom_indices: Vec<TopologyAtomIndex>,
    source_bond_indices: Vec<TopologyBondIndex>,
}

impl TopologySubsetCorrespondence {
    pub fn source_topology(&self) -> &Topology {
        &self.source
    }

    pub fn target_topology(&self) -> &Topology {
        &self.target
    }

    pub fn target_atom(&self, source: InstanceAtomId) -> Option<InstanceAtomId> {
        self.atom_targets.get(&source).copied()
    }

    pub fn target_bond(&self, source: InstanceBondId) -> Option<InstanceBondId> {
        self.bond_targets.get(&source).copied()
    }

    /// Source dense atom index for every target atom in target dense order.
    pub fn source_atom_indices(&self) -> &[TopologyAtomIndex] {
        &self.source_atom_indices
    }

    /// Source dense bond index for every target bond in target dense order.
    pub fn source_bond_indices(&self) -> &[TopologyBondIndex] {
        &self.source_bond_indices
    }
}

impl Topology {
    /// Constructs the induced hierarchy-aware topology selected by `selection`.
    ///
    /// Every selected source molecule is partitioned into connected induced
    /// components. Target molecule and dense ordering follow source instance,
    /// source atom, and source bond order deterministically.
    pub fn subset(&self, selection: &AtomSelection) -> Result<TopologySubset, TopologySubsetError> {
        if !std::ptr::eq(self, selection.topology()) {
            return Err(SelectionError::TopologyMismatch.into());
        }
        let source = selection.shared_topology();
        if selection.indices().is_empty() {
            return Err(TopologySubsetError::EmptySelection);
        }
        let selected = selection
            .indices()
            .iter()
            .filter_map(|index| self.atom_id(*index))
            .collect::<BTreeSet<_>>();
        let mut builder = TopologyBuilder::new();
        let mut atom_targets = BTreeMap::new();
        let mut bond_targets = BTreeMap::new();

        for molecule_view in self.molecules() {
            let source_molecule = molecule_view.molecule();
            let selected_local = source_molecule
                .atom_ids()
                .filter(|atom| selected.contains(&InstanceAtomId::new(molecule_view.id(), *atom)))
                .collect::<BTreeSet<_>>();
            if selected_local.is_empty() {
                continue;
            }
            let mut visited = BTreeSet::new();
            for seed in source_molecule.atom_ids() {
                if !selected_local.contains(&seed) || !visited.insert(seed) {
                    continue;
                }
                let mut queue = VecDeque::from([seed]);
                let mut component = Vec::new();
                while let Some(atom) = queue.pop_front() {
                    component.push(atom);
                    for neighbor in source_molecule.neighbors(atom)? {
                        if selected_local.contains(&neighbor) && visited.insert(neighbor) {
                            queue.push_back(neighbor);
                        }
                    }
                }
                component.sort_unstable();
                let membership = component.iter().copied().collect::<BTreeSet<_>>();
                let mut editor = source_molecule.edit();
                let discarded = source_molecule
                    .atom_ids()
                    .filter(|atom| !membership.contains(atom))
                    .collect::<Vec<_>>();
                for atom in discarded {
                    editor.delete_atom(atom)?;
                }
                let target_molecule = editor.finish()?;
                let target_instance = builder.add_molecule(&target_molecule)?;
                for source_atom in membership.iter().copied() {
                    atom_targets.insert(
                        InstanceAtomId::new(molecule_view.id(), source_atom),
                        InstanceAtomId::new(target_instance, source_atom),
                    );
                }
                for (source_bond, bond) in source_molecule.bonds() {
                    if !membership.contains(&bond.a()) || !membership.contains(&bond.b()) {
                        continue;
                    }
                    bond_targets.insert(
                        InstanceBondId::new(molecule_view.id(), source_bond),
                        InstanceBondId::new(target_instance, source_bond),
                    );
                }
            }
        }

        copy_filtered_hierarchy(self, builder.hierarchy_mut(), &atom_targets)?;
        let target = Arc::new(builder.build()?);
        let atom_sources = atom_targets
            .iter()
            .map(|(source, target)| (*target, *source))
            .collect::<BTreeMap<_, _>>();
        let bond_sources = bond_targets
            .iter()
            .map(|(source, target)| (*target, *source))
            .collect::<BTreeMap<_, _>>();
        let source_atom_indices = target
            .atom_ids()
            .iter()
            .map(|target_atom| {
                let source_atom = atom_sources[target_atom];
                self.atom_index(source_atom)
                    .expect("source subset atom has dense index")
            })
            .collect();
        let source_bond_indices = target
            .bond_ids()
            .iter()
            .map(|target_bond| {
                let source_bond = bond_sources[target_bond];
                self.bond_index(source_bond)
                    .expect("source subset bond has dense index")
            })
            .collect();
        let correspondence = TopologySubsetCorrespondence {
            source,
            target: Arc::clone(&target),
            atom_targets,
            bond_targets,
            source_atom_indices,
            source_bond_indices,
        };
        Ok(TopologySubset {
            topology: target,
            correspondence,
        })
    }
}

fn copy_filtered_hierarchy(
    source: &Topology,
    target: &mut super::Hierarchy,
    atom_targets: &BTreeMap<InstanceAtomId, InstanceAtomId>,
) -> Result<(), HierarchyError> {
    target.props_mut().clone_from(source.hierarchy().props());
    for (_, source_chain) in source.hierarchy().chains() {
        let retained_residues = source_chain
            .residues()
            .iter()
            .filter_map(|id| source.hierarchy().residue(*id).ok())
            .filter(|residue| {
                residue.atom_sites().iter().any(|site| {
                    source
                        .hierarchy()
                        .atom_site(*site)
                        .is_ok_and(|site| atom_targets.contains_key(&site.atom()))
                })
            })
            .collect::<Vec<_>>();
        if retained_residues.is_empty() {
            continue;
        }
        let chain = target.add_chain(
            source_chain.label_id().to_owned(),
            source_chain.author_id().map(str::to_owned),
        )?;
        target
            .chain_props_mut(chain)?
            .clone_from(source_chain.props());
        for source_residue in retained_residues {
            let residue = target.add_residue(
                chain,
                source_residue.name().to_owned(),
                source_residue.label_seq_id(),
                source_residue.author_seq_id().map(str::to_owned),
                source_residue.insertion_code().map(str::to_owned),
            )?;
            target.set_residue_component_ids(
                residue,
                source_residue.label_comp_id().map(str::to_owned),
                source_residue.author_comp_id().map(str::to_owned),
            )?;
            target
                .residue_props_mut(residue)?
                .clone_from(source_residue.props());
            for source_site_id in source_residue.atom_sites() {
                let source_site = source.hierarchy().atom_site(*source_site_id)?;
                let Some(target_atom) = atom_targets.get(&source_site.atom()).copied() else {
                    continue;
                };
                let site =
                    target.add_atom_site(residue, target_atom, source_site.metadata().clone())?;
                target
                    .atom_site_props_mut(site)?
                    .clone_from(source_site.props());
            }
        }
    }
    Ok(())
}

/// Failure to construct a hierarchy-aware induced topology subset.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologySubsetError {
    EmptySelection,
    Selection(SelectionError),
    Molecule(MoleculeError),
    Publication(MoleculePublicationError),
    Hierarchy(HierarchyError),
    TopologyBuild(TopologyBuildError),
}

impl fmt::Display for TopologySubsetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySelection => formatter.write_str("topology subset selection is empty"),
            Self::Selection(error) => write!(formatter, "invalid subset selection: {error}"),
            Self::Molecule(error) => write!(formatter, "cannot construct subset molecule: {error}"),
            Self::Publication(error) => {
                write!(formatter, "cannot publish subset molecule: {error}")
            }
            Self::Hierarchy(error) => write!(formatter, "cannot filter subset hierarchy: {error}"),
            Self::TopologyBuild(error) => {
                write!(formatter, "cannot publish subset topology: {error}")
            }
        }
    }
}

impl std::error::Error for TopologySubsetError {}

impl From<SelectionError> for TopologySubsetError {
    fn from(error: SelectionError) -> Self {
        Self::Selection(error)
    }
}
impl From<MoleculeError> for TopologySubsetError {
    fn from(error: MoleculeError) -> Self {
        Self::Molecule(error)
    }
}
impl From<MoleculePublicationError> for TopologySubsetError {
    fn from(error: MoleculePublicationError) -> Self {
        Self::Publication(error)
    }
}
impl From<HierarchyError> for TopologySubsetError {
    fn from(error: HierarchyError) -> Self {
        Self::Hierarchy(error)
    }
}
impl From<TopologyBuildError> for TopologySubsetError {
    fn from(error: TopologyBuildError) -> Self {
        Self::TopologyBuild(error)
    }
}
