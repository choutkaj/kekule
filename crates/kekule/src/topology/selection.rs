use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::core::Element;
use crate::substructure::QueryMatch;

use super::{
    AtomSiteId, AtomSiteView, ChainId, ChainView, InstanceAtomId, MoleculeClass,
    MoleculeDefinitionId, MoleculeInstanceId, ResidueClass, ResidueId, ResidueView, Topology,
    TopologyAtomIndex,
};

/// A topology-bound, sorted, unique dense atom selection.
#[derive(Debug, Clone)]
pub struct AtomSelection {
    topology: Arc<Topology>,
    indices: Vec<TopologyAtomIndex>,
}

impl PartialEq for AtomSelection {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.indices == other.indices
    }
}

impl Eq for AtomSelection {}

impl AtomSelection {
    /// Selects every atom in authoritative dense order, sharing this exact topology.
    ///
    /// This is infallible because the topology has already validated its layout.
    pub fn all(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            indices: (0..topology.atom_count())
                .map(|index| TopologyAtomIndex::new(index as u32))
                .collect(),
        }
    }

    pub fn from_atoms(
        topology: &Arc<Topology>,
        atoms: impl IntoIterator<Item = InstanceAtomId>,
    ) -> Result<Self, SelectionError> {
        let mut indices = atoms
            .into_iter()
            .map(|atom| {
                topology
                    .atom_index(atom)
                    .ok_or(SelectionError::InvalidAtomId(atom))
            })
            .collect::<Result<Vec<_>, _>>()?;
        indices.sort_unstable();
        indices.dedup();
        Ok(Self {
            topology: Arc::clone(topology),
            indices,
        })
    }

    pub fn from_indices(
        topology: &Arc<Topology>,
        indices: impl IntoIterator<Item = TopologyAtomIndex>,
    ) -> Result<Self, SelectionError> {
        let atoms = indices
            .into_iter()
            .map(|index| {
                topology
                    .atom_id(index)
                    .ok_or(SelectionError::InvalidAtomIndex(index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(topology, atoms)
    }

    pub fn for_instances(
        topology: &Arc<Topology>,
        instances: impl IntoIterator<Item = MoleculeInstanceId>,
    ) -> Result<Self, SelectionError> {
        let instances = instances.into_iter().collect::<BTreeSet<_>>();
        for instance in &instances {
            topology
                .instance(*instance)
                .map_err(|_| SelectionError::InvalidMoleculeInstanceId(*instance))?;
        }
        Self::from_atoms(
            topology,
            topology
                .atom_ids()
                .iter()
                .copied()
                .filter(|atom| instances.contains(&atom.molecule())),
        )
    }

    pub fn for_definitions(
        topology: &Arc<Topology>,
        definitions: impl IntoIterator<Item = MoleculeDefinitionId>,
    ) -> Result<Self, SelectionError> {
        let definitions = definitions.into_iter().collect::<BTreeSet<_>>();
        for definition in &definitions {
            topology
                .definition(*definition)
                .map_err(|_| SelectionError::InvalidMoleculeDefinitionId(*definition))?;
        }
        let instances = topology
            .instances()
            .filter(|(_, instance)| definitions.contains(&instance.definition()))
            .map(|(id, _)| id);
        Self::for_instances(topology, instances)
    }

    /// Selects all atoms in molecule instances having one of `classes`.
    pub fn for_molecule_classes(
        topology: &Arc<Topology>,
        classes: impl IntoIterator<Item = MoleculeClass>,
    ) -> Result<Self, SelectionError> {
        let classes = classes.into_iter().collect::<BTreeSet<_>>();
        Self::for_instances(
            topology,
            topology
                .molecules()
                .filter(|molecule| classes.contains(&molecule.class()))
                .map(|molecule| molecule.id()),
        )
    }

    pub fn for_elements(
        topology: &Arc<Topology>,
        elements: impl IntoIterator<Item = Element>,
    ) -> Result<Self, SelectionError> {
        let elements = elements.into_iter().collect::<BTreeSet<_>>();
        Self::from_atoms(
            topology,
            topology
                .atoms()
                .filter(|(_, atom)| elements.contains(&atom.element))
                .map(|(id, _)| id),
        )
    }

    pub fn for_atom_sites(
        topology: &Arc<Topology>,
        atom_sites: impl IntoIterator<Item = AtomSiteId>,
    ) -> Result<Self, SelectionError> {
        let atom_sites = atom_sites.into_iter().collect::<BTreeSet<_>>();
        let atoms = atom_sites
            .into_iter()
            .map(|site| {
                topology
                    .atom_for_site(site)
                    .map_err(|_| SelectionError::InvalidAtomSiteId(site))
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_atoms(topology, atoms)
    }

    pub fn for_residues(
        topology: &Arc<Topology>,
        residues: impl IntoIterator<Item = ResidueId>,
    ) -> Result<Self, SelectionError> {
        let residues = residues.into_iter().collect::<BTreeSet<_>>();
        for residue in &residues {
            topology
                .residue(*residue)
                .map_err(|_| SelectionError::InvalidResidueId(*residue))?;
        }
        Self::for_atom_sites(
            topology,
            topology
                .atom_sites()
                .filter(|site| residues.contains(&site.residue().id()))
                .map(AtomSiteView::id),
        )
    }

    /// Selects all atoms in residues having one of `classes`.
    pub fn for_residue_classes(
        topology: &Arc<Topology>,
        classes: impl IntoIterator<Item = ResidueClass>,
    ) -> Result<Self, SelectionError> {
        let classes = classes.into_iter().collect::<BTreeSet<_>>();
        Self::for_residues(
            topology,
            topology
                .residues()
                .filter(|residue| classes.contains(&residue.class()))
                .map(ResidueView::id),
        )
    }

    pub fn for_chains(
        topology: &Arc<Topology>,
        chains: impl IntoIterator<Item = ChainId>,
    ) -> Result<Self, SelectionError> {
        let chains = chains.into_iter().collect::<BTreeSet<_>>();
        for chain in &chains {
            topology
                .chain(*chain)
                .map_err(|_| SelectionError::InvalidChainId(*chain))?;
        }
        Self::for_residues(
            topology,
            topology
                .residues()
                .filter(|residue| chains.contains(&residue.chain().id()))
                .map(ResidueView::id),
        )
    }

    pub fn for_chain_label(topology: &Arc<Topology>, label: &str) -> Result<Self, SelectionError> {
        Self::for_chains(
            topology,
            topology
                .chains()
                .filter(|chain| chain.label_id() == label)
                .map(ChainView::id),
        )
    }

    pub fn from_query_matches(
        topology: &Arc<Topology>,
        instance: MoleculeInstanceId,
        matches: &[QueryMatch],
    ) -> Result<Self, SelectionError> {
        topology
            .instance(instance)
            .map_err(|_| SelectionError::InvalidMoleculeInstanceId(instance))?;
        Self::from_atoms(
            topology,
            matches.iter().flat_map(|query_match| {
                query_match
                    .atoms()
                    .iter()
                    .copied()
                    .map(move |atom| InstanceAtomId::new(instance, atom))
            }),
        )
    }

    pub fn ensure_compatible(&self, topology: &Arc<Topology>) -> Result<(), SelectionError> {
        if !Arc::ptr_eq(&self.topology, topology) {
            return Err(SelectionError::TopologyMismatch);
        }
        Ok(())
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub(crate) fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn indices(&self) -> &[TopologyAtomIndex] {
        &self.indices
    }

    pub fn semantic_ids(
        &self,
        topology: &Arc<Topology>,
    ) -> Result<Vec<InstanceAtomId>, SelectionError> {
        self.ensure_compatible(topology)?;
        Ok(self
            .indices
            .iter()
            .map(|index| {
                topology
                    .atom_id(*index)
                    .expect("selection contains validated dense indices")
            })
            .collect())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum SelectionError {
    TopologyMismatch,
    InvalidMoleculeDefinitionId(MoleculeDefinitionId),
    InvalidMoleculeInstanceId(MoleculeInstanceId),
    InvalidAtomId(InstanceAtomId),
    InvalidChainId(ChainId),
    InvalidResidueId(ResidueId),
    InvalidAtomSiteId(AtomSiteId),
    InvalidAtomIndex(TopologyAtomIndex),
}

impl fmt::Display for SelectionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("atom selection belongs to a different topology")
            }
            Self::InvalidMoleculeDefinitionId(id) => {
                write!(formatter, "invalid selected molecule definition: {id}")
            }
            Self::InvalidMoleculeInstanceId(id) => {
                write!(formatter, "invalid selected molecule instance: {id}")
            }
            Self::InvalidAtomId(id) => write!(formatter, "invalid selected atom: {id}"),
            Self::InvalidChainId(id) => write!(formatter, "invalid selected chain: {id}"),
            Self::InvalidResidueId(id) => write!(formatter, "invalid selected residue: {id}"),
            Self::InvalidAtomSiteId(id) => write!(formatter, "invalid selected atom site: {id}"),
            Self::InvalidAtomIndex(index) => write!(formatter, "invalid selected {index}"),
        }
    }
}

impl std::error::Error for SelectionError {}
