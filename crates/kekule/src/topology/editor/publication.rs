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
    atom_order: Vec<(usize, InstanceAtomId)>,
    bond_slots: Vec<usize>,
}

/// Correspondence for one successful editor publication. Deleted or foreign
/// handles have no target. Dense indices refer to the returned topology only.
#[derive(Debug, Clone, Default)]
pub struct EditCorrespondence {
    atoms: BTreeMap<EditAtomId, (InstanceAtomId, super::super::TopologyAtomIndex)>,
    bonds: BTreeMap<EditBondId, (InstanceBondId, super::super::TopologyBondIndex)>,
}
impl EditCorrespondence {
    pub fn atom(
        &self,
        id: EditAtomId,
    ) -> Option<(InstanceAtomId, super::super::TopologyAtomIndex)> {
        self.atoms.get(&id).copied()
    }
    pub fn bond(
        &self,
        id: EditBondId,
    ) -> Option<(InstanceBondId, super::super::TopologyBondIndex)> {
        self.bonds.get(&id).copied()
    }
}

impl TopologyEditor {
    /// Publishes a topology with mappings for every surviving editing handle.
    pub fn finish_with_correspondence(
        self,
    ) -> Result<(Arc<Topology>, EditCorrespondence), TopologyEditError> {
        let (published, correspondence) = self.into_mapped_publication()?;
        Ok((published.topology, correspondence))
    }

    pub(crate) fn into_mapped_publication(
        self,
    ) -> Result<(TopologyPublication, EditCorrespondence), TopologyEditError> {
        let atoms: BTreeMap<_, _> = self.atoms.iter().map(|(&id, loc)| (loc.slot, id)).collect();
        let bonds: BTreeMap<_, _> = self.bonds.iter().map(|(&id, loc)| (loc.slot, id)).collect();
        let published = self.into_publication()?;
        let mut correspondence = EditCorrespondence::default();
        for (slot, &target) in published
            .atom_slots
            .iter()
            .zip(published.topology.atom_ids())
        {
            correspondence.atoms.insert(
                atoms[slot],
                (target, published.topology.atom_index(target).unwrap()),
            );
        }
        for (slot, &target) in published
            .bond_slots
            .iter()
            .zip(published.topology.bond_ids())
        {
            correspondence.bonds.insert(
                bonds[slot],
                (target, published.topology.bond_index(target).unwrap()),
            );
        }
        Ok((published, correspondence))
    }

    /// Checks full component publication, classification, hierarchy and properties.
    /// This evaluates a snapshot and leaves all draft state and handles unchanged.
    pub fn validate(&self) -> Result<(), TopologyEditError> {
        self.publish().map(|_| ())
    }
    /// Publishes an immutable topology, moving owned molecular drafts and any
    /// uniquely owned source definitions. Shared source definitions are cloned.
    /// Inferred classes are refreshed only when their hierarchy or chemical
    /// evidence changes. Editing handles are draft-only.
    pub fn finish(self) -> Result<Arc<Topology>, TopologyEditError> {
        self.into_publication().map(|r| r.topology)
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
        self.clone().into_publication()
    }

    pub(crate) fn into_publication(mut self) -> Result<TopologyPublication, TopologyEditError> {
        if self.is_empty() {
            return Err(TopologyEditError::EmptyTopology);
        }
        if self.revision == 0 {
            if let Some(source) = self.source.take() {
                return Ok(TopologyPublication {
                    atom_slots: (0..source.atom_count()).collect(),
                    bond_slots: (0..source.bond_count()).collect(),
                    topology: source,
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
        let invalidated = self
            .groups
            .iter()
            .flatten()
            .filter_map(|group| match group.chemistry {
                GroupChemistry::Source(id) if group.class.is_none() => Some(id),
                _ => None,
            })
            .collect::<BTreeSet<_>>();
        let invalidated_added = self
            .groups
            .iter()
            .flatten()
            .filter_map(|group| match &group.chemistry {
                GroupChemistry::Added(molecule) if group.class.is_none() => {
                    Some(Arc::as_ptr(molecule))
                }
                _ => None,
            })
            .collect::<BTreeSet<_>>();
        if let Some(source) = self.source.take() {
            let (definitions, overrides) = match Arc::try_unwrap(source) {
                Ok(source) => {
                    let overrides = source.layout.molecule_class_overrides.clone();
                    (source.definitions, overrides)
                }
                Err(source) => (
                    source
                        .definitions
                        .iter()
                        .filter(|d| retained.contains(&d.id()))
                        .cloned()
                        .collect(),
                    source.layout.molecule_class_overrides.clone(),
                ),
            };
            for definition in definitions {
                let id = definition.id();
                if !retained.contains(&id) {
                    continue;
                }
                let class = definition.class();
                let target = builder.add_molecule_definition(definition.molecule)?;
                if !invalidated.contains(&id) || overrides.contains_key(&id) {
                    builder.preserve_molecule_class(target, class, overrides.contains_key(&id))?;
                }
                source_definitions.insert(id, target);
            }
        }
        let mut instance_sources = Vec::new();
        for group in std::mem::take(&mut self.groups).into_iter().flatten() {
            if group.atoms.is_empty() {
                continue;
            }
            match group.chemistry {
                GroupChemistry::Source(id) => {
                    let definition = source_definitions[&id];
                    let instance = builder.add_instance(definition)?;
                    let molecule = builder.definition(definition)?.molecule();
                    instance_sources.push(group.instance_slot);
                    self.record_component(
                        (&group.atoms, &group.bonds),
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids(),
                        instance,
                        &mut targets,
                    );
                }
                GroupChemistry::Added(owned) => {
                    // Only explicit shared definitions are reused, never chemical guesses.
                    let identity = Arc::as_ptr(&owned);
                    let definition = if let Some(&id) = added_definitions.get(&identity) {
                        id
                    } else {
                        let molecule =
                            Arc::try_unwrap(owned).unwrap_or_else(|shared| (*shared).clone());
                        let id = builder.add_molecule_definition(molecule)?;
                        if let Some(class) = group.class {
                            if !invalidated_added.contains(&identity) || group.class_explicit {
                                builder.preserve_molecule_class(id, class, group.class_explicit)?;
                            }
                        }
                        added_definitions.insert(identity, id);
                        id
                    };
                    let instance = builder.add_instance(definition)?;
                    let molecule = builder.definition(definition)?.molecule();
                    instance_sources.push(group.instance_slot);
                    self.record_component(
                        (&group.atoms, &group.bonds),
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids(),
                        instance,
                        &mut targets,
                    );
                }
                GroupChemistry::Draft(draft) if draft.is_connected() => {
                    // Publication renumbers deleted draft slots densely.
                    let (published, ids) = (*draft).finish_with_correspondence()?;
                    let definition = builder.add_molecule_definition(published)?;
                    let handles = group.atoms.values().copied().collect();
                    if let Some(class) = self.component_class(&handles, group.class) {
                        let explicit = group.class_explicit
                            || self.molecule_classes.keys().any(|id| handles.contains(id));
                        builder.preserve_molecule_class(definition, class, explicit)?;
                    }
                    let instance = builder.add_instance(definition)?;
                    instance_sources
                        .push((!group.changed).then_some(group.instance_slot).flatten());
                    self.record_component(
                        (&group.atoms, &group.bonds),
                        ids.atoms().iter().map(|(&draft, &local)| (draft, local)),
                        ids.bonds().keys().copied(),
                        instance,
                        &mut targets,
                    );
                }
                GroupChemistry::Draft(draft) => {
                    let molecule = draft.working();
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
                            (&group.atoms, &group.bonds),
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
        // Publish atoms in slot order: survivors keep their source dense
        // order, and added atoms follow in creation order.
        targets.atom_order.sort_unstable_by_key(|&(slot, _)| slot);
        builder.set_atom_order(targets.atom_order.iter().map(|&(_, atom)| atom))?;
        let atom_slots = targets
            .atom_order
            .iter()
            .map(|&(slot, _)| slot)
            .collect::<Vec<_>>();
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
                builder.preserve_residue_class(id, class, residue.class_explicit)?;
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
        let mut properties = self.properties.project(
            &instance_sources,
            &atom_slots,
            &targets.bond_slots,
            &chain_slots,
            &residue_slots,
            &site_slots,
        )?;
        for (key, value) in self.properties.owner().iter() {
            properties.owner_mut().insert(key.clone(), value.clone())?;
        }
        builder.install_properties(properties);
        let target = Arc::new(builder.build()?);
        Ok(TopologyPublication {
            topology: target,
            atom_slots,
            bond_slots: targets.bond_slots,
        })
    }

    fn record_component(
        &self,
        handles: (&BTreeMap<AtomId, EditAtomId>, &BTreeMap<BondId, EditBondId>),
        atoms: impl Iterator<Item = (AtomId, AtomId)>,
        bonds: impl Iterator<Item = BondId>,
        instance: MoleculeInstanceId,
        targets: &mut Targets,
    ) {
        for (source, local) in atoms {
            let handle = handles.0[&source];
            let target = InstanceAtomId::new(instance, local);
            targets.atoms.insert(handle, target);
            targets.atom_order.push((self.atoms[&handle].slot, target));
        }
        for source in bonds {
            let handle = handles.1[&source];
            targets.bond_slots.push(self.bonds[&handle].slot);
        }
    }
}
