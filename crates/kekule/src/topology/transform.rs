//! Hierarchy-aware induced topology subsets.
//!
//! [`Topology::subset`] returns the subset together with its narrow
//! dense-state correspondence.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Arc;

use super::{
    AtomSelection, HierarchyError, InstanceAtomId, InstanceBondId, ResidueId, SelectionError,
    Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuildError, TopologyBuilder,
};
use crate::molecule::{AtomId, BondId, MoleculeError, MoleculePublicationError};
use crate::properties::PropertyError;

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
    /// components. Target instances follow source instance order, then source
    /// atom order within each instance. Target dense atom order keeps the
    /// selected atoms' source dense order, so
    /// [`TopologySubsetCorrespondence::source_atom_indices`] is ascending.
    /// Components use compact local atom and bond IDs; use the returned
    /// correspondence to translate source identities into target identities.
    /// Instances selecting the same atoms from a shared definition reuse the
    /// reconstructed definitions. Complete molecules and residues retain their
    /// classifications; changed entities are reclassified. The selection may
    /// come from any snapshot sharing this layout; the subset and its
    /// correspondence always use the `source` snapshot, which they retain.
    pub fn subset(
        source: impl Into<Arc<Self>>,
        selection: &AtomSelection,
    ) -> Result<TopologySubset, TopologySubsetError> {
        let source: Arc<Self> = source.into();
        let this = &source;
        if !this.shares_layout(selection.topology()) {
            return Err(SelectionError::TopologyMismatch.into());
        }
        if selection.indices().is_empty() {
            return Err(TopologySubsetError::EmptySelection);
        }
        let selected = selection
            .indices()
            .iter()
            .map(|index| {
                this.atom_id(*index)
                    .expect("validated atom selection references a live dense atom")
            })
            .collect::<BTreeSet<_>>();
        let mut builder = TopologyBuilder::new();
        let mut atom_targets = BTreeMap::new();
        let mut bond_targets = BTreeMap::new();
        let mut instance_sources = Vec::new();
        let mut definition_subsets = BTreeMap::new();

        for molecule_view in this.instances() {
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
                        let definitions = super::components::build_component_definitions(
                            source_molecule,
                            &entry.key().1,
                            Some((
                                molecule_view.class(),
                                this.molecule_class_overrides
                                    .contains_key(&molecule_view.definition_id()),
                            )),
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

        // Selected atoms keep their relative source dense order.
        builder.set_atom_order(selection.indices().iter().map(|&index| {
            atom_targets[&this
                .atom_id(index)
                .expect("validated atom selection references a live dense atom")]
        }))?;
        let hierarchy_projection =
            copy_filtered_hierarchy(this, builder.hierarchy_mut().storage_mut(), &atom_targets)?;
        for (target_index, &source_index) in hierarchy_projection.residues.iter().enumerate() {
            let residue = this
                .hierarchy()
                .require_residue(ResidueId::new(source_index as u32))?;
            let complete = residue.atom_sites().iter().all(|site| {
                let atom = this
                    .hierarchy()
                    .atom_site(*site)
                    .expect("published residue references a live atom site")
                    .atom();
                atom_targets.contains_key(&atom)
            });
            if complete {
                builder.preserve_residue_class(
                    ResidueId::new(target_index as u32),
                    residue.class(),
                    residue.explicit_class.is_some(),
                )?;
            }
        }
        let mut target = builder.build()?;
        target.properties = Arc::new(project_topology_properties(
            this,
            &target,
            &instance_sources,
            &atom_targets,
            &bond_targets,
            &hierarchy_projection,
        )?);
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
                this.atom_index(source_atom)
                    .expect("source subset atom has dense index")
            })
            .collect();
        let source_bond_indices = target
            .bond_ids()
            .iter()
            .map(|target_bond| {
                let source_bond = bond_sources[target_bond];
                this.bond_index(source_bond)
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
                let source_site = source.hierarchy().require_atom_site(*source_site_id)?;
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
) -> Result<crate::properties::TopologyProperties, PropertyError> {
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
    source.properties().project(
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
    /// The selection contains no atoms.
    EmptySelection,
    /// The selection does not match the topology; see the source.
    Selection(SelectionError),
    /// Extracting a selected molecule fragment failed; see the source.
    Molecule(MoleculeError),
    /// A selected molecule fragment could not be published; see the source.
    Publication(MoleculePublicationError),
    /// Building the subset hierarchy failed; see the source.
    Hierarchy(HierarchyError),
    /// Publishing the subset topology failed; see the source.
    TopologyBuild(TopologyBuildError),
    /// Projecting property rows onto the subset failed; see the source.
    Property(PropertyError),
}

impl fmt::Display for TopologySubsetError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySelection => formatter.write_str("topology subset selection is empty"),
            Self::Selection(_) => formatter.write_str("invalid subset selection"),
            Self::Molecule(_) => formatter.write_str("cannot construct subset molecule"),
            Self::Publication(_) => formatter.write_str("cannot publish subset molecule"),
            Self::Hierarchy(_) => formatter.write_str("cannot filter subset hierarchy"),
            Self::TopologyBuild(_) => formatter.write_str("cannot publish subset topology"),
            Self::Property(_) => formatter.write_str("cannot project subset properties"),
        }
    }
}

impl std::error::Error for TopologySubsetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Selection(error) => Some(error),
            Self::Molecule(error) => Some(error),
            Self::Publication(error) => Some(error),
            Self::Hierarchy(error) => Some(error),
            Self::TopologyBuild(error) => Some(error),
            Self::Property(error) => Some(error),
            Self::EmptySelection => None,
        }
    }
}

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

impl From<super::components::ComponentBuildError> for TopologySubsetError {
    fn from(error: super::components::ComponentBuildError) -> Self {
        use super::components::ComponentBuildError;
        match error {
            ComponentBuildError::Molecule(e) => Self::Molecule(e),
            ComponentBuildError::Publication(e) => Self::Publication(e),
            ComponentBuildError::Property(e) => Self::Property(e),
            ComponentBuildError::Topology(e) => Self::TopologyBuild(e),
        }
    }
}
