//! Immutable topology transformations.
//!
//! Whole-instance filters preserve definitions and instances. [`Topology::subset`]
//! is the separate hierarchy-aware induced-graph operation and returns its own
//! narrow dense-state correspondence.

use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;
use std::sync::Arc;

use super::{
    AtomSelection, HierarchyError, InstanceAtomId, InstanceBondId, MoleculeInstanceId, ResidueId,
    SelectionError, Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuildError,
    TopologyBuilder,
};
use crate::core::{
    AtomId, BondId, MoleculeEditor, MoleculeError, MoleculePublicationError, StereoCarrier,
    StereoElement, StereoElementKind, StereoGroup,
};
use crate::properties::PropertyError;

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
        builder.set_molecule_class(target_id, definition.class())?;
        definition_targets[source_id.index()] = Some(target_id);
    }

    let mut atom_targets = BTreeMap::new();
    let mut bond_targets = BTreeMap::new();
    let mut instance_sources = Vec::with_capacity(retained.len);
    for (source_instance, instance) in topology
        .instances()
        .filter(|(id, _)| retained.contains(*id))
    {
        let target_definition = definition_targets[instance.definition().index()]
            .expect("retained instance has a retained definition");
        let target_instance = builder.add_instance(target_definition)?;
        instance_sources.push(Some(source_instance.index()));
        let source_definition = topology
            .definition_for_instance(source_instance)
            .expect("retained instance references a live definition")
            .molecule();
        for atom in source_definition.atom_ids() {
            atom_targets.insert(
                InstanceAtomId::new(source_instance, atom),
                InstanceAtomId::new(target_instance, atom),
            );
        }
        for bond in source_definition.bond_ids() {
            bond_targets.insert(
                InstanceBondId::new(source_instance, bond),
                InstanceBondId::new(target_instance, bond),
            );
        }
    }

    let hierarchy_projection =
        copy_filtered_hierarchy(topology, builder.hierarchy_mut(), &atom_targets)?;
    for (target_index, source_index) in hierarchy_projection.residues.iter().copied().enumerate() {
        let class = topology
            .hierarchy()
            .residue(ResidueId::new(source_index as u32))
            .expect("retained residue projection references a live source residue")
            .class();
        builder.set_residue_class(ResidueId::new(target_index as u32), class)?;
    }

    let mut target = builder.build()?;
    target.properties = project_topology_properties(
        topology,
        &target,
        &instance_sources,
        &atom_targets,
        &bond_targets,
        &hierarchy_projection,
    )?;
    Ok(Arc::new(target))
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
    Property(PropertyError),
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
            Self::Property(error) => {
                write!(formatter, "cannot project topology properties: {error}")
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

impl From<PropertyError> for TopologyTransformError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
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
    /// Components use compact local atom and bond IDs; use the returned
    /// correspondence to translate source identities into target identities.
    /// Instances selecting the same atoms from a shared definition reuse the
    /// reconstructed definitions. Complete molecules and residues retain their
    /// classifications; changed entities are reclassified.
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
            .map(|index| {
                self.atom_id(*index)
                    .expect("validated atom selection references a live dense atom")
            })
            .collect::<BTreeSet<_>>();
        let mut builder = TopologyBuilder::new();
        let mut atom_targets = BTreeMap::new();
        let mut bond_targets = BTreeMap::new();
        let mut instance_sources = Vec::new();
        let mut definition_subsets = BTreeMap::new();

        for molecule_view in self.molecules() {
            let source_molecule = molecule_view.molecule();
            let selected_local = source_molecule
                .atom_ids()
                .filter(|atom| selected.contains(&InstanceAtomId::new(molecule_view.id(), *atom)))
                .collect::<BTreeSet<_>>();
            if selected_local.is_empty() {
                continue;
            }
            let whole_instance_selected = selected_local.len() == source_molecule.atom_count();
            let definitions =
                match definition_subsets.entry((molecule_view.definition_id(), selected_local)) {
                    std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        let definitions = build_subset_definitions(
                            source_molecule,
                            &entry.key().1,
                            molecule_view.class(),
                            &mut builder,
                        )?;
                        entry.insert(definitions)
                    }
                };
            for definition in definitions.iter() {
                let target_instance = builder.add_instance(definition.id)?;
                instance_sources
                    .push(whole_instance_selected.then_some(molecule_view.id().index()));
                for (index, &source_atom) in definition.source_atoms.iter().enumerate() {
                    atom_targets.insert(
                        InstanceAtomId::new(molecule_view.id(), source_atom),
                        InstanceAtomId::new(target_instance, AtomId::new(index as u32)),
                    );
                }
                for (index, &source_bond) in definition.source_bonds.iter().enumerate() {
                    bond_targets.insert(
                        InstanceBondId::new(molecule_view.id(), source_bond),
                        InstanceBondId::new(target_instance, BondId::new(index as u32)),
                    );
                }
            }
        }

        let hierarchy_projection =
            copy_filtered_hierarchy(self, builder.hierarchy_mut(), &atom_targets)?;
        for (target_index, &source_index) in hierarchy_projection.residues.iter().enumerate() {
            let residue = self
                .hierarchy()
                .residue(ResidueId::new(source_index as u32))?;
            let complete = residue.atom_sites().iter().all(|site| {
                let atom = self
                    .hierarchy()
                    .atom_site(*site)
                    .expect("published residue references a live atom site")
                    .atom();
                atom_targets.contains_key(&atom)
            });
            if complete {
                builder.set_residue_class(ResidueId::new(target_index as u32), residue.class())?;
            }
        }
        let mut target = builder.build()?;
        target.properties = project_topology_properties(
            self,
            &target,
            &instance_sources,
            &atom_targets,
            &bond_targets,
            &hierarchy_projection,
        )?;
        let target = Arc::new(target);
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

struct SubsetDefinition {
    id: super::MoleculeDefinitionId,
    source_atoms: Vec<AtomId>,
    source_bonds: Vec<BondId>,
}

fn build_subset_definitions(
    source_molecule: &crate::core::Molecule,
    selected_local: &BTreeSet<AtomId>,
    class: super::MoleculeClass,
    builder: &mut TopologyBuilder,
) -> Result<Vec<SubsetDefinition>, TopologySubsetError> {
    let whole_instance_selected = selected_local.len() == source_molecule.atom_count();
    let mut definitions = Vec::new();
    let mut visited = BTreeSet::new();
    let mut components = Vec::new();
    let mut local_atoms = BTreeMap::new();
    let mut local_bonds = BTreeMap::new();
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
        let mut editor = MoleculeEditor::new();
        for &atom in &component {
            let target_atom = editor.add_atom(source_molecule.atom(atom)?.clone())?;
            local_atoms.insert(atom, (components.len(), target_atom));
        }
        components.push(SubsetComponent {
            editor,
            source_atoms: component,
            source_bonds: Vec::new(),
        });
    }

    // Scan source bonds and stereo once for all components, rather than
    // scanning or cloning the whole source molecule for every fragment.
    for (source_bond, bond) in source_molecule.bonds() {
        let (Some(&(component, a)), Some(&(other_component, b))) =
            (local_atoms.get(&bond.a()), local_atoms.get(&bond.b()))
        else {
            continue;
        };
        debug_assert_eq!(component, other_component);
        let target_bond = components[component].editor.add_bond(a, b, bond.order)?;
        components[component].source_bonds.push(source_bond);
        local_bonds.insert(source_bond, (component, target_bond));
    }
    let mut local_stereo = BTreeMap::new();
    for (source_element, element) in source_molecule.stereo_elements() {
        if let Some((component, element)) = remap_subset_stereo(element, &local_atoms, &local_bonds)
        {
            let target_element = components[component].editor.add_stereo_element(element)?;
            local_stereo.insert(source_element, (component, target_element));
        }
    }
    for (_, group) in source_molecule.stereo_groups() {
        let mut members_by_component = BTreeMap::<_, Vec<_>>::new();
        for member in &group.members {
            if let Some(&(component, target_element)) = local_stereo.get(member) {
                members_by_component
                    .entry(component)
                    .or_default()
                    .push(target_element);
            }
        }
        for (component, members) in members_by_component {
            components[component].editor.add_stereo_group(StereoGroup {
                kind: group.kind,
                members,
            })?;
        }
    }

    for mut component in components {
        let properties = &mut component.editor.working_mut().properties;
        *properties.atoms_mut() = source_molecule.atom_properties().select_indices(
            &component
                .source_atoms
                .iter()
                .map(|atom| atom.index())
                .collect::<Vec<_>>(),
        )?;
        *properties.bonds_mut() = source_molecule.bond_properties().select_indices(
            &component
                .source_bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<Vec<_>>(),
        )?;
        if whole_instance_selected {
            for (key, value) in source_molecule.properties().iter() {
                properties.insert(key.clone(), value.clone())?;
            }
        }
        let target_molecule = component.editor.finish()?;
        let definition = builder.add_molecule_definition_owned(target_molecule)?;
        if whole_instance_selected {
            builder.set_molecule_class(definition, class)?;
        }
        definitions.push(SubsetDefinition {
            id: definition,
            source_atoms: component.source_atoms,
            source_bonds: component.source_bonds,
        });
    }
    Ok(definitions)
}

#[derive(Debug, Default)]
struct SubsetComponent {
    editor: MoleculeEditor,
    source_atoms: Vec<AtomId>,
    source_bonds: Vec<BondId>,
}

// Retain an assertion only when all its represented references survive in the
// same component, matching atom/bond deletion's stereo-pruning semantics.
fn remap_subset_stereo(
    element: &StereoElement,
    atoms: &BTreeMap<AtomId, (usize, AtomId)>,
    bonds: &BTreeMap<BondId, (usize, BondId)>,
) -> Option<(usize, StereoElement)> {
    let component = match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => atoms.get(&stereo.center)?.0,
        StereoElementKind::DoubleBond(stereo) => bonds.get(&stereo.bond)?.0,
        StereoElementKind::Axis(stereo) => bonds.get(&stereo.axis)?.0,
    };
    let atom = |source| {
        atoms
            .get(&source)
            .and_then(|&(owner, target)| (owner == component).then_some(target))
    };
    let bond = |source| {
        bonds
            .get(&source)
            .and_then(|&(owner, target)| (owner == component).then_some(target))
    };
    let carrier = |source| match source {
        StereoCarrier::Atom(source) => atom(source).map(StereoCarrier::Atom),
        other => Some(other),
    };
    let mut kind = element.kind.clone();
    match &mut kind {
        StereoElementKind::Tetrahedral(stereo) => {
            stereo.center = atom(stereo.center)?;
            for value in &mut stereo.carriers {
                *value = carrier(*value)?;
            }
        }
        StereoElementKind::DoubleBond(stereo) => {
            stereo.bond = bond(stereo.bond)?;
            stereo.left = atom(stereo.left)?;
            stereo.right = atom(stereo.right)?;
            stereo.left_carrier = carrier(stereo.left_carrier)?;
            stereo.right_carrier = carrier(stereo.right_carrier)?;
        }
        StereoElementKind::Axis(stereo) => {
            stereo.axis = bond(stereo.axis)?;
            for value in &mut stereo.carriers {
                *value = carrier(*value)?;
            }
        }
    }
    Some((component, StereoElement::new(kind)))
}

#[derive(Debug, Default)]
struct HierarchyProjection {
    chains: Vec<usize>,
    residues: Vec<usize>,
    atom_sites: Vec<usize>,
}

fn copy_filtered_hierarchy(
    source: &Topology,
    target: &mut super::Hierarchy,
    atom_targets: &BTreeMap<InstanceAtomId, InstanceAtomId>,
) -> Result<HierarchyProjection, HierarchyError> {
    let mut projection = HierarchyProjection::default();
    for (source_chain_id, source_chain) in source.hierarchy().chains() {
        let retained_residues = source_chain
            .residues()
            .iter()
            .map(|id| {
                source
                    .hierarchy()
                    .residue(*id)
                    .expect("published hierarchy chain references a live residue")
            })
            .filter(|residue| {
                residue.atom_sites().iter().any(|site| {
                    source
                        .hierarchy()
                        .atom_site(*site)
                        .map(|site| atom_targets.contains_key(&site.atom()))
                        .expect("published hierarchy residue references a live atom site")
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
        projection.chains.push(source_chain_id.index());
        for source_residue in retained_residues {
            let residue = target.add_residue(
                chain,
                source_residue.name().to_owned(),
                source_residue.label_seq_id(),
                source_residue.author_seq_id().map(str::to_owned),
                source_residue.insertion_code().map(str::to_owned),
            )?;
            projection.residues.push(source_residue.id().index());
            target.set_residue_component_ids(
                residue,
                source_residue.label_comp_id().map(str::to_owned),
                source_residue.author_comp_id().map(str::to_owned),
            )?;
            for source_site_id in source_residue.atom_sites() {
                let source_site = source.hierarchy().atom_site(*source_site_id)?;
                let Some(target_atom) = atom_targets.get(&source_site.atom()).copied() else {
                    continue;
                };
                target.add_atom_site(residue, target_atom, source_site.metadata().clone())?;
                projection.atom_sites.push(source_site.id().index());
            }
        }
    }
    Ok(projection)
}

fn project_topology_properties(
    source: &Topology,
    target: &Topology,
    instance_sources: &[Option<usize>],
    atom_targets: &BTreeMap<InstanceAtomId, InstanceAtomId>,
    bond_targets: &BTreeMap<InstanceBondId, InstanceBondId>,
    hierarchy: &HierarchyProjection,
) -> Result<crate::properties::Properties, PropertyError> {
    let atom_sources = atom_targets
        .iter()
        .map(|(source, target)| (*target, *source))
        .collect::<BTreeMap<_, _>>();
    let bond_sources = bond_targets
        .iter()
        .map(|(source, target)| (*target, *source))
        .collect::<BTreeMap<_, _>>();
    let atom_indices = target
        .atom_ids()
        .iter()
        .map(|target_atom| {
            source
                .atom_index(atom_sources[target_atom])
                .expect("projected source atom has a dense index")
                .index()
        })
        .collect::<Vec<_>>();
    let bond_indices = target
        .bond_ids()
        .iter()
        .map(|target_bond| {
            source
                .bond_index(bond_sources[target_bond])
                .expect("projected source bond has a dense index")
                .index()
        })
        .collect::<Vec<_>>();
    source.properties().project_topology(
        instance_sources,
        &atom_indices,
        &bond_indices,
        &hierarchy.chains,
        &hierarchy.residues,
        &hierarchy.atom_sites,
    )
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
    Property(PropertyError),
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
            Self::Property(error) => write!(formatter, "cannot project subset properties: {error}"),
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

impl From<PropertyError> for TopologySubsetError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}
