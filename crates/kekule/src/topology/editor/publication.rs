use super::*;

/// Internal row projection needed to publish a model's realization arrays.
pub(crate) struct TopologyPublication {
    pub topology: Arc<Topology>,
    pub atom_slots: Vec<usize>,
    pub bond_slots: Vec<usize>,
}

#[derive(Default)]
struct Targets {
    atoms: BTreeMap<EditAtomId, InstanceAtomId>,
    chains: BTreeMap<EditChainId, ChainId>,
    residues: BTreeMap<EditResidueId, ResidueId>,
    atom_slots: Vec<usize>,
    bond_slots: Vec<usize>,
}

impl TopologyEditor {
    /// Checks full component publication, classification, hierarchy and properties.
    /// This evaluates a snapshot and leaves all draft state and handles unchanged.
    pub fn validate(&self) -> Result<(), TopologyEditError> {
        self.publish().map(|_| ())
    }
    /// Publishes an immutable topology. Editing handles are draft-only.
    pub fn finish(self) -> Result<Arc<Topology>, TopologyEditError> {
        self.publish().map(|r| r.topology)
    }
    /// Publishes a topology, returning the draft with any publication error.
    pub fn try_finish(self) -> Result<Arc<Topology>, TopologyFinishError> {
        match self.publish() {
            Ok(result) => Ok(result.topology),
            Err(error) => Err(TopologyFinishError {
                error: Box::new(error),
                editor: Box::new(self),
            }),
        }
    }
    pub(crate) fn publish(&self) -> Result<TopologyPublication, TopologyEditError> {
        if self.is_empty() {
            return Err(TopologyEditError::EmptyTopology);
        }
        if self.revision == 0 {
            if let Some(source) = &self.source {
                return Ok(TopologyPublication {
                    topology: Arc::clone(source),
                    atom_slots: (0..source.atom_count()).collect(),
                    bond_slots: (0..source.bond_count()).collect(),
                });
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
                        molecule.bond_ids(),
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
                        molecule.bond_ids(),
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
                        molecule.bond_ids(),
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
                            component.source_bonds.into_iter(),
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
        for site in self.hierarchy.sites.values() {
            let atom = targets.atoms[&site.atom];
            builder.hierarchy_mut().add_atom_site(
                targets.residues[&site.residue],
                atom,
                site.metadata.clone(),
            )?;
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
        Ok(TopologyPublication {
            topology: target,
            atom_slots: targets.atom_slots,
            bond_slots: targets.bond_slots,
        })
    }

    fn record_component(
        &self,
        group: &Group,
        atoms: impl Iterator<Item = (AtomId, AtomId)>,
        bonds: impl Iterator<Item = BondId>,
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
        for source in bonds {
            let handle = group.bonds[&source];
            targets.bond_slots.push(self.bonds[&handle].slot);
        }
    }
}
