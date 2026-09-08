use super::*;
use crate::topology::{TopologyAtomIndex, TopologyBondIndex};

/// One published edit and its transaction-specific identity correspondence.
#[derive(Debug, Clone)]
pub struct TopologyEdit {
    topology: Arc<Topology>,
    correspondence: TopologyEditCorrespondence,
}
impl TopologyEdit {
    pub fn topology(&self) -> &Topology {
        &self.topology
    }
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }
    pub fn correspondence(&self) -> &TopologyEditCorrespondence {
        &self.correspondence
    }
    pub fn into_parts(self) -> (Arc<Topology>, TopologyEditCorrespondence) {
        (self.topology, self.correspondence)
    }
}

/// Actual source/draft-to-output correspondence for one structural transaction.
///
/// Deleted entities return `None`; new entities have no source dense index.
/// Instance correspondence is one-to-many for splits and may share a target after
/// merges. This does not transfer bindings of selections or prepared calculations.
#[derive(Debug, Clone)]
pub struct TopologyEditCorrespondence {
    source: Option<Arc<Topology>>,
    target: Arc<Topology>,
    atoms: BTreeMap<EditAtomId, InstanceAtomId>,
    bonds: BTreeMap<EditBondId, InstanceBondId>,
    source_atoms: BTreeMap<InstanceAtomId, InstanceAtomId>,
    source_bonds: BTreeMap<InstanceBondId, InstanceBondId>,
    instances: BTreeMap<MoleculeInstanceId, Vec<MoleculeInstanceId>>,
    chains: BTreeMap<EditChainId, ChainId>,
    residues: BTreeMap<EditResidueId, ResidueId>,
    sites: BTreeMap<EditAtomSiteId, AtomSiteId>,
    source_chains: BTreeMap<ChainId, ChainId>,
    source_residues: BTreeMap<ResidueId, ResidueId>,
    source_sites: BTreeMap<AtomSiteId, AtomSiteId>,
    source_atom_indices: Vec<Option<TopologyAtomIndex>>,
    source_bond_indices: Vec<Option<TopologyBondIndex>>,
    pub(crate) atom_slots: Vec<usize>,
    pub(crate) bond_slots: Vec<usize>,
    append_tokens: Vec<Arc<()>>,
}
impl TopologyEditCorrespondence {
    pub(crate) fn contains_append(&self, token: &Arc<()>) -> bool {
        self.append_tokens
            .iter()
            .any(|known| Arc::ptr_eq(known, token))
    }

    pub fn source_topology(&self) -> Option<&Topology> {
        self.source.as_deref()
    }
    pub fn target_topology(&self) -> &Topology {
        &self.target
    }
    pub fn atom(&self, handle: EditAtomId) -> Option<InstanceAtomId> {
        self.atoms.get(&handle).copied()
    }
    pub fn bond(&self, handle: EditBondId) -> Option<InstanceBondId> {
        self.bonds.get(&handle).copied()
    }
    pub fn chain(&self, handle: EditChainId) -> Option<ChainId> {
        self.chains.get(&handle).copied()
    }
    pub fn residue(&self, handle: EditResidueId) -> Option<ResidueId> {
        self.residues.get(&handle).copied()
    }
    pub fn atom_site(&self, handle: EditAtomSiteId) -> Option<AtomSiteId> {
        self.sites.get(&handle).copied()
    }
    pub fn target_atom(&self, source: InstanceAtomId) -> Option<InstanceAtomId> {
        self.source_atoms.get(&source).copied()
    }
    pub fn target_bond(&self, source: InstanceBondId) -> Option<InstanceBondId> {
        self.source_bonds.get(&source).copied()
    }
    pub fn target_chain(&self, source: ChainId) -> Option<ChainId> {
        self.source_chains.get(&source).copied()
    }
    pub fn target_residue(&self, source: ResidueId) -> Option<ResidueId> {
        self.source_residues.get(&source).copied()
    }
    pub fn target_atom_site(&self, source: AtomSiteId) -> Option<AtomSiteId> {
        self.source_sites.get(&source).copied()
    }
    pub fn target_instances(&self, source: MoleculeInstanceId) -> Option<&[MoleculeInstanceId]> {
        self.instances.get(&source).map(Vec::as_slice)
    }
    /// Source dense rows in output atom order; added atoms have `None`.
    pub fn source_atom_indices(&self) -> &[Option<TopologyAtomIndex>] {
        &self.source_atom_indices
    }
    pub fn source_bond_indices(&self) -> &[Option<TopologyBondIndex>] {
        &self.source_bond_indices
    }
}

#[derive(Default)]
struct Targets {
    atoms: BTreeMap<EditAtomId, InstanceAtomId>,
    bonds: BTreeMap<EditBondId, InstanceBondId>,
    chains: BTreeMap<EditChainId, ChainId>,
    residues: BTreeMap<EditResidueId, ResidueId>,
    sites: BTreeMap<EditAtomSiteId, AtomSiteId>,
    atom_slots: Vec<usize>,
    bond_slots: Vec<usize>,
}

impl TopologyEditor {
    /// Checks full component publication, classification, hierarchy and properties.
    /// This evaluates a snapshot and leaves all draft state and handles unchanged.
    pub fn validate(&self) -> Result<(), TopologyEditError> {
        self.publish().map(|_| ())
    }
    pub fn finish(self) -> Result<Arc<Topology>, TopologyEditError> {
        self.publish().map(|r| r.topology)
    }
    pub fn finish_with_correspondence(self) -> Result<TopologyEdit, TopologyEditError> {
        self.publish()
    }
    pub fn try_finish(self) -> Result<Arc<Topology>, TopologyFinishError> {
        match self.publish() {
            Ok(result) => Ok(result.topology),
            Err(error) => Err(TopologyFinishError {
                error: Box::new(error),
                editor: Box::new(self),
            }),
        }
    }
    pub fn try_finish_with_correspondence(self) -> Result<TopologyEdit, TopologyFinishError> {
        match self.publish() {
            Ok(result) => Ok(result),
            Err(error) => Err(TopologyFinishError {
                error: Box::new(error),
                editor: Box::new(self),
            }),
        }
    }

    fn publish(&self) -> Result<TopologyEdit, TopologyEditError> {
        if self.is_empty() {
            return Err(TopologyEditError::EmptyTopology);
        }
        if self.revision == 0 {
            if let Some(source) = &self.source {
                let targets = Targets {
                    atoms: self.source_atoms.iter().map(|(&s, &h)| (h, s)).collect(),
                    bonds: self.source_bonds.iter().map(|(&s, &h)| (h, s)).collect(),
                    chains: self
                        .hierarchy
                        .source_chains
                        .iter()
                        .map(|(&s, &h)| (h, s))
                        .collect(),
                    residues: self
                        .hierarchy
                        .source_residues
                        .iter()
                        .map(|(&s, &h)| (h, s))
                        .collect(),
                    sites: self
                        .hierarchy
                        .source_sites
                        .iter()
                        .map(|(&s, &h)| (h, s))
                        .collect(),
                    atom_slots: (0..source.atom_count()).collect(),
                    bond_slots: (0..source.bond_count()).collect(),
                };
                return Ok(self.result(Arc::clone(source), targets));
            }
        }
        let mut builder = TopologyBuilder::new();
        let mut targets = Targets::default();
        let mut source_definitions = BTreeMap::new();
        let mut added_definitions = BTreeMap::new();
        // Preserve definition order/reuse for untouched source occurrences.
        let retained = self
            .groups
            .iter()
            .flatten()
            .filter_map(|g| match g.chemistry {
                GroupChemistry::Source(id) => Some(id),
                _ => None,
            })
            .collect::<BTreeSet<_>>();
        if let Some(source) = &self.source {
            for id in retained {
                let definition = source.definition(id).unwrap();
                let target = builder.add_molecule_definition(definition.molecule())?;
                builder.set_molecule_class(target, definition.class())?;
                source_definitions.insert(id, target);
            }
        }
        let mut instance_sources = Vec::new();
        for (group_index, group) in self
            .groups
            .iter()
            .enumerate()
            .filter_map(|(i, g)| g.as_ref().map(|g| (i, g)))
        {
            if group.atoms.is_empty() {
                continue;
            }
            let molecule = self.molecule(group_index);
            match &group.chemistry {
                GroupChemistry::Source(id) => {
                    let instance = builder.add_instance(source_definitions[id])?;
                    instance_sources.push(group.instance_slot);
                    self.record_component(
                        group,
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids().map(|b| (b, b)),
                        instance,
                        &mut targets,
                    );
                }
                GroupChemistry::Added(owned) => {
                    // Only explicit shared definitions are reused, never chemical guesses.
                    let definition = if let Some(&id) = added_definitions.get(&Arc::as_ptr(owned)) {
                        id
                    } else {
                        let id = builder.add_molecule_definition(molecule)?;
                        if let Some(class) = group.class {
                            builder.set_molecule_class(id, class)?;
                        }
                        added_definitions.insert(Arc::as_ptr(owned), id);
                        id
                    };
                    let instance = builder.add_instance(definition)?;
                    instance_sources.push(group.instance_slot);
                    self.record_component(
                        group,
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids().map(|b| (b, b)),
                        instance,
                        &mut targets,
                    );
                }
                GroupChemistry::Draft(draft) if draft.is_connected() => {
                    let published = draft.as_ref().clone().finish()?;
                    let definition = builder.add_molecule_definition_owned(published)?;
                    if let Some(class) =
                        self.component_class(&group.atoms.values().copied().collect(), group.class)
                    {
                        builder.set_molecule_class(definition, class)?;
                    }
                    let instance = builder.add_instance(definition)?;
                    instance_sources
                        .push((!group.changed).then_some(group.instance_slot).flatten());
                    self.record_component(
                        group,
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids().map(|b| (b, b)),
                        instance,
                        &mut targets,
                    );
                }
                GroupChemistry::Draft(_) => {
                    let components = super::super::components::build_component_definitions(
                        molecule,
                        &molecule.atom_ids().collect(),
                        None,
                        &mut builder,
                    )?;
                    for component in components {
                        let handles = component
                            .source_atoms
                            .iter()
                            .map(|id| group.atoms[id])
                            .collect();
                        if let Some(class) = self.component_class(&handles, None) {
                            builder.set_molecule_class(component.id, class)?;
                        }
                        let instance = builder.add_instance(component.id)?;
                        // No old occurrence annotation is inherited across a split.
                        instance_sources.push(None);
                        self.record_component(
                            group,
                            component
                                .source_atoms
                                .into_iter()
                                .enumerate()
                                .map(|(i, s)| (s, AtomId::new(i as u32))),
                            component
                                .source_bonds
                                .into_iter()
                                .enumerate()
                                .map(|(i, s)| (s, BondId::new(i as u32))),
                            instance,
                            &mut targets,
                        );
                    }
                }
            }
        }
        let mut chain_slots = Vec::new();
        let mut residue_slots = Vec::new();
        let mut site_slots = Vec::new();
        for (&handle, chain) in &self.hierarchy.chains {
            let id = builder
                .hierarchy_mut()
                .add_chain(chain.label.clone(), chain.author.clone())?;
            targets.chains.insert(handle, id);
            chain_slots.push(chain.slot);
        }
        for (&handle, residue) in &self.hierarchy.residues {
            let chain = targets.chains[&residue.chain];
            let id = builder.hierarchy_mut().add_residue(
                chain,
                residue.name.clone(),
                residue.label_seq,
                residue.author_seq.clone(),
                residue.insertion.clone(),
            )?;
            builder.hierarchy_mut().set_residue_component_ids(
                id,
                residue.label_comp.clone(),
                residue.author_comp.clone(),
            )?;
            if let Some(class) = residue.class {
                builder.set_residue_class(id, class)?;
            }
            targets.residues.insert(handle, id);
            residue_slots.push(residue.slot);
        }
        for (&handle, site) in &self.hierarchy.sites {
            let atom = targets.atoms[&site.atom];
            let id = builder.hierarchy_mut().add_atom_site(
                targets.residues[&site.residue],
                atom,
                site.metadata.clone(),
            )?;
            targets.sites.insert(handle, id);
            site_slots.push(site.slot);
        }
        let mut properties = self.properties.project_topology(
            &instance_sources,
            &targets.atom_slots,
            &targets.bond_slots,
            &chain_slots,
            &residue_slots,
            &site_slots,
        )?;
        for (key, value) in self.properties.iter() {
            properties.insert(key.clone(), value.clone())?;
        }
        builder.install_properties(properties);
        let target = Arc::new(builder.build()?);
        Ok(self.result(target, targets))
    }

    fn record_component(
        &self,
        group: &Group,
        atoms: impl Iterator<Item = (AtomId, AtomId)>,
        bonds: impl Iterator<Item = (BondId, BondId)>,
        instance: MoleculeInstanceId,
        targets: &mut Targets,
    ) {
        for (source, local) in atoms {
            let handle = group.atoms[&source];
            targets
                .atoms
                .insert(handle, InstanceAtomId::new(instance, local));
            targets.atom_slots.push(self.atoms[&handle].slot);
        }
        for (source, local) in bonds {
            let handle = group.bonds[&source];
            targets
                .bonds
                .insert(handle, InstanceBondId::new(instance, local));
            targets.bond_slots.push(self.bonds[&handle].slot);
        }
    }
    fn result(&self, target: Arc<Topology>, targets: Targets) -> TopologyEdit {
        let source_atoms = self
            .source_atoms
            .iter()
            .filter_map(|(&s, h)| targets.atoms.get(h).map(|&t| (s, t)))
            .collect::<BTreeMap<_, _>>();
        let source_bonds = self
            .source_bonds
            .iter()
            .filter_map(|(&s, h)| targets.bonds.get(h).map(|&t| (s, t)))
            .collect::<BTreeMap<_, _>>();
        let atom_sources = source_atoms
            .iter()
            .map(|(&s, &t)| (t, s))
            .collect::<BTreeMap<_, _>>();
        let bond_sources = source_bonds
            .iter()
            .map(|(&s, &t)| (t, s))
            .collect::<BTreeMap<_, _>>();
        let mut instances = BTreeMap::<_, BTreeSet<_>>::new();
        if let Some(source) = &self.source {
            for (id, _) in source.instances() {
                instances.insert(id, BTreeSet::new());
            }
        }
        for (&source, &target) in &source_atoms {
            instances
                .get_mut(&source.molecule())
                .unwrap()
                .insert(target.molecule());
        }
        let source_atom_indices = target
            .atom_ids()
            .iter()
            .map(|id| {
                atom_sources
                    .get(id)
                    .and_then(|id| self.source.as_ref().unwrap().atom_index(*id))
            })
            .collect();
        let source_bond_indices = target
            .bond_ids()
            .iter()
            .map(|id| {
                bond_sources
                    .get(id)
                    .and_then(|id| self.source.as_ref().unwrap().bond_index(*id))
            })
            .collect();
        let source_chains = self
            .hierarchy
            .source_chains
            .iter()
            .filter_map(|(&s, h)| targets.chains.get(h).map(|&t| (s, t)))
            .collect();
        let source_residues = self
            .hierarchy
            .source_residues
            .iter()
            .filter_map(|(&s, h)| targets.residues.get(h).map(|&t| (s, t)))
            .collect();
        let source_sites = self
            .hierarchy
            .source_sites
            .iter()
            .filter_map(|(&s, h)| targets.sites.get(h).map(|&t| (s, t)))
            .collect();
        let correspondence = TopologyEditCorrespondence {
            source: self.source.clone(),
            target: Arc::clone(&target),
            atoms: targets.atoms,
            bonds: targets.bonds,
            source_atoms,
            source_bonds,
            instances: instances
                .into_iter()
                .map(|(s, t)| (s, t.into_iter().collect()))
                .collect(),
            chains: targets.chains,
            residues: targets.residues,
            sites: targets.sites,
            source_chains,
            source_residues,
            source_sites,
            source_atom_indices,
            source_bond_indices,
            atom_slots: targets.atom_slots,
            bond_slots: targets.bond_slots,
            append_tokens: self.append_tokens.clone(),
        };
        TopologyEdit {
            topology: target,
            correspondence,
        }
    }
}
