use super::*;
use crate::topology::{ResidueLabels, TopologyAtomIndex, TopologyBondIndex};
use crate::DraftError;

/// Internal row projection needed to publish a model's structure arrays.
pub(crate) struct TopologyPublication {
    pub topology: Arc<Topology>,
    pub atom_slots: Vec<usize>,
    pub bond_slots: Vec<usize>,
    /// Chain and residue handles in published ID order.
    chains: Vec<EditChainId>,
    residues: Vec<EditResidueId>,
}

struct Targets {
    atoms: DenseMap<EditAtomId, InstanceAtomId>,
    chains: DenseMap<EditChainId, ChainId>,
    residues: DenseMap<EditResidueId, ResidueId>,
    // Published atoms by draft slot.
    atom_order: Vec<Option<InstanceAtomId>>,
    bond_slots: Vec<usize>,
}

/// One dense index space's correspondence between an editor's source
/// topology and the topology it published.
///
/// [`Self::to_published`] has one entry per source index, in source order:
/// the published index, or `None` for a deleted entity. [`Self::to_source`]
/// has one entry per published index: the source index, or `None` for an
/// added entity. A host remaps a column of per-entity values in one pass:
///
/// ```
/// use kekule::{smiles, topology::{Topology, edit::TopologyEditor}};
///
/// let source = Topology::from_molecules(smiles::to_molecules("CCO.O")?)?;
/// let mut editor = TopologyEditor::from_topology(source.clone());
/// let carbon = editor.atom_handle(source.atom_ids()[0]).expect("source atom");
/// editor.delete_atom(carbon)?;
/// let (published, correspondence) = editor.finish_with_correspondence()?;
///
/// let colors = ["red", "green", "blue", "cyan"];
/// let atoms = correspondence.atoms();
/// let remapped = atoms
///     .to_source()
///     .iter()
///     .map(|source| source.map(|index| colors[index.index()]))
///     .collect::<Vec<_>>();
/// assert_eq!(remapped, [Some("green"), Some("blue"), Some("cyan")]);
/// assert_eq!(atoms.to_published()[0], None);
/// assert_eq!(atoms.to_published().len(), source.atom_count());
/// assert_eq!(atoms.to_source().len(), published.atom_count());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DenseCorrespondence<I> {
    published: Vec<Option<I>>,
    source: Vec<Option<I>>,
}

impl<I> Default for DenseCorrespondence<I> {
    fn default() -> Self {
        Self {
            published: Vec::new(),
            source: Vec::new(),
        }
    }
}

impl<I: Copy> DenseCorrespondence<I> {
    /// The published index of every source index, in source order; `None`
    /// for a deleted entity. Empty for an editor without a source topology.
    pub fn to_published(&self) -> &[Option<I>] {
        &self.published
    }

    /// The source index of every published index, in published order;
    /// `None` for an added entity.
    pub fn to_source(&self) -> &[Option<I>] {
        &self.source
    }
}

macro_rules! dense_lookups {
    ($($index:ty),*) => {$(
        impl DenseCorrespondence<$index> {
            /// The published index of a source index, or `None` when the
            /// entity was deleted or `source` is out of range.
            pub fn published(&self, source: $index) -> Option<$index> {
                self.published.get(source.index()).copied().flatten()
            }

            /// The source index of a published index, or `None` when the
            /// entity was added or `published` is out of range.
            pub fn source(&self, published: $index) -> Option<$index> {
                self.source.get(published.index()).copied().flatten()
            }

            /// Builds both directions from the source index of each
            /// published index.
            fn from_sources(source_count: usize, source: Vec<Option<$index>>) -> Self {
                let mut published = vec![None; source_count];
                for (index, value) in source.iter().enumerate() {
                    if let Some(value) = value {
                        published[value.index()] = Some(<$index>::new(index as u32));
                    }
                }
                Self { published, source }
            }
        }
    )*};
}
dense_lookups!(TopologyAtomIndex, TopologyBondIndex, ChainId, ResidueId);

/// Correspondence for one successful editor publication. Deleted or foreign
/// handles have no target. Dense indices refer to the returned topology only.
///
/// Besides the handle lookups, it maps the source topology's dense atom and
/// bond indices, residues, and chains to the published ones and back
/// ([`DenseCorrespondence`]), so hosts remap per-entity columns in one pass.
/// The source is the topology the editor was started from; entities appended
/// or added since have no source, and a draft started empty has an empty
/// source.
#[derive(Debug, Clone, Default)]
pub struct EditCorrespondence {
    atoms: DenseMap<EditAtomId, (InstanceAtomId, TopologyAtomIndex)>,
    bonds: DenseMap<EditBondId, (InstanceBondId, TopologyBondIndex)>,
    dense_atoms: DenseCorrespondence<TopologyAtomIndex>,
    dense_bonds: DenseCorrespondence<TopologyBondIndex>,
    residues: DenseCorrespondence<ResidueId>,
    chains: DenseCorrespondence<ChainId>,
}
impl EditCorrespondence {
    /// The published qualified ID and dense index of a draft atom, or `None` when
    /// the handle was deleted or is foreign.
    pub fn atom(&self, id: EditAtomId) -> Option<(InstanceAtomId, TopologyAtomIndex)> {
        self.atoms.get(&id).copied()
    }
    /// The published qualified ID and dense index of a draft bond, or `None` when
    /// the handle was deleted or is foreign.
    pub fn bond(&self, id: EditBondId) -> Option<(InstanceBondId, TopologyBondIndex)> {
        self.bonds.get(&id).copied()
    }
    /// Source to published dense atom indices and back.
    pub fn atoms(&self) -> &DenseCorrespondence<TopologyAtomIndex> {
        &self.dense_atoms
    }
    /// Source to published dense bond indices and back.
    pub fn bonds(&self) -> &DenseCorrespondence<TopologyBondIndex> {
        &self.dense_bonds
    }
    /// Source to published hierarchy residues and back. A residue the edit
    /// emptied is deleted.
    pub fn residues(&self) -> &DenseCorrespondence<ResidueId> {
        &self.residues
    }
    /// Source to published hierarchy chains and back. A chain the edit
    /// emptied is deleted.
    pub fn chains(&self) -> &DenseCorrespondence<ChainId> {
        &self.chains
    }
}

impl TopologyEditor {
    /// Publishes a topology with mappings for every surviving editing handle.
    ///
    /// # Errors
    ///
    /// As for [`Self::finish`].
    pub fn finish_with_correspondence(
        self,
    ) -> Result<(Arc<Topology>, EditCorrespondence), TopologyEditError> {
        let (published, correspondence) = self.into_mapped_publication()?;
        Ok((published.topology, correspondence))
    }

    pub(crate) fn into_mapped_publication(
        self,
    ) -> Result<(TopologyPublication, EditCorrespondence), TopologyEditError> {
        let mut atom_handles = vec![None; self.atom_slot_count()];
        for (&id, loc) in &self.atoms {
            atom_handles[loc.slot] = Some(id);
        }
        let mut bond_handles = vec![None; self.bond_slot_count()];
        for (&id, loc) in &self.bonds {
            bond_handles[loc.slot] = Some(id);
        }
        let blocks = self.source_blocks;
        let mut correspondence = EditCorrespondence {
            atoms: DenseMap::block_of(&self.atoms),
            bonds: DenseMap::block_of(&self.bonds),
            ..EditCorrespondence::default()
        };
        let published = self.into_publication()?;
        let topology = &published.topology;
        let mut atom_sources = Vec::with_capacity(topology.atom_count());
        for (index, (&slot, &target)) in published
            .atom_slots
            .iter()
            .zip(topology.atom_ids())
            .enumerate()
        {
            let handle = atom_handles[slot].expect("published slots are live");
            let index = TopologyAtomIndex::new(index as u32);
            correspondence.atoms.insert(handle, (target, index));
            atom_sources.push(blocks.atoms.source(handle).map(TopologyAtomIndex::new));
        }
        let mut bond_sources = Vec::with_capacity(topology.bond_count());
        for (index, (&slot, &target)) in published
            .bond_slots
            .iter()
            .zip(topology.bond_ids())
            .enumerate()
        {
            let handle = bond_handles[slot].expect("published slots are live");
            let index = TopologyBondIndex::new(index as u32);
            correspondence.bonds.insert(handle, (target, index));
            bond_sources.push(blocks.bonds.source(handle).map(TopologyBondIndex::new));
        }
        correspondence.dense_atoms =
            DenseCorrespondence::<TopologyAtomIndex>::from_sources(blocks.atoms.len, atom_sources);
        correspondence.dense_bonds =
            DenseCorrespondence::<TopologyBondIndex>::from_sources(blocks.bonds.len, bond_sources);
        correspondence.residues = DenseCorrespondence::<ResidueId>::from_sources(
            blocks.residues.len,
            published
                .residues
                .iter()
                .map(|&handle| blocks.residues.source(handle).map(ResidueId::new))
                .collect(),
        );
        correspondence.chains = DenseCorrespondence::<ChainId>::from_sources(
            blocks.chains.len,
            published
                .chains
                .iter()
                .map(|&handle| blocks.chains.source(handle).map(ChainId::new))
                .collect(),
        );
        Ok((published, correspondence))
    }

    /// Checks full component publication, classification, hierarchy and properties.
    /// This evaluates a snapshot and leaves all draft state and handles unchanged.
    pub fn validate(&self) -> Result<(), TopologyEditError> {
        self.publish().map(|_| ())
    }
    /// Publishes an immutable topology, moving owned molecular drafts. Source
    /// definitions no edit touched are shared with the source snapshot, with
    /// their perception, rather than copied.
    /// Inferred classes are refreshed only when their hierarchy or chemical
    /// evidence changes. Editing handles are draft-only.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::EmptyTopology`] when no atoms remain, or a
    /// [`TopologyEditError::Molecule`], [`TopologyEditError::Publication`],
    /// [`TopologyEditError::Topology`], or [`TopologyEditError::Property`]
    /// error when a component, the topology, or its properties cannot be
    /// published.
    pub fn finish(self) -> Result<Arc<Topology>, TopologyEditError> {
        self.into_publication().map(|r| r.topology)
    }
    /// Publishes a topology, returning the draft with any publication error.
    pub fn try_finish(self) -> Result<Arc<Topology>, DraftError<Self, TopologyEditError>> {
        match self.publish() {
            Ok(result) => Ok(result.topology),
            Err(error) => Err(DraftError::new(error, self)),
        }
    }
    pub(crate) fn publish(&self) -> Result<TopologyPublication, TopologyEditError> {
        self.clone().into_publication()
    }

    pub(crate) fn into_publication(mut self) -> Result<TopologyPublication, TopologyEditError> {
        if self.is_empty() {
            return Err(TopologyEditError::EmptyTopology);
        }
        let owner_unchanged = self
            .source
            .as_ref()
            .is_some_and(|source| source.properties().owner() == self.properties.owner());
        if self.revision == 0 && owner_unchanged {
            if let Some(source) = self.source.take() {
                return Ok(TopologyPublication {
                    atom_slots: (0..source.atom_count()).collect(),
                    bond_slots: (0..source.bond_count()).collect(),
                    chains: self.hierarchy.chains.keys().copied().collect(),
                    residues: self.hierarchy.residues.keys().copied().collect(),
                    topology: source,
                });
            }
        }
        let mut builder = TopologyBuilder::new();
        let mut targets = Targets {
            atoms: DenseMap::block_of(&self.atoms),
            chains: DenseMap::block_of(&self.hierarchy.chains),
            residues: DenseMap::block_of(&self.hierarchy.residues),
            atom_order: vec![None; self.atom_slot_count()],
            bond_slots: Vec::with_capacity(self.bonds.len()),
        };
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
            // Untouched definitions share the source molecule, whether or not
            // the source snapshot is still shared.
            let overrides = &source.molecule_class_overrides;
            for definition in &source.definitions {
                let id = definition.id();
                if !retained.contains(&id) {
                    continue;
                }
                let class = definition.class();
                let target = builder.add_shared_definition(
                    Arc::clone(&definition.molecule),
                    definition.default_perception,
                )?;
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
                    let molecule = builder.require_definition(definition)?.molecule();
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
                        let id = builder.add_shared_definition(owned, false)?;
                        if let Some(class) = group.class {
                            if !invalidated_added.contains(&identity) || group.class_explicit {
                                builder.preserve_molecule_class(id, class, group.class_explicit)?;
                            }
                        }
                        added_definitions.insert(identity, id);
                        id
                    };
                    let instance = builder.add_instance(definition)?;
                    let molecule = builder.require_definition(definition)?.molecule();
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
        let atom_slots = (0..targets.atom_order.len())
            .filter(|&slot| targets.atom_order[slot].is_some())
            .collect::<Vec<_>>();
        builder.set_atom_order(targets.atom_order.iter().flatten().copied())?;
        let mut chain_slots = Vec::with_capacity(self.hierarchy.chains.len());
        let mut residue_slots = Vec::with_capacity(self.hierarchy.residues.len());
        let mut site_slots = Vec::with_capacity(self.hierarchy.sites.len());
        for (&handle, chain) in &self.hierarchy.chains {
            let id = builder
                .hierarchy_mut()
                .add_chain(chain.label.clone(), chain.author.clone())?;
            targets.chains.insert(handle, id);
            chain_slots.push(chain.slot);
        }
        let mut residue_sites = DenseMap::block_of(&self.hierarchy.residues);
        for site in self.hierarchy.sites.values() {
            let count = residue_sites.get(&site.residue).copied().unwrap_or(0);
            residue_sites.insert(site.residue, count + 1);
        }
        for (&handle, residue) in &self.hierarchy.residues {
            let chain = targets.chains[&residue.chain];
            let id = builder.hierarchy_mut().add_shared_residue(
                chain,
                ResidueLabels {
                    name: Arc::clone(&residue.name),
                    label_comp_id: residue.label_comp.clone(),
                    author_comp_id: residue.author_comp.clone(),
                    label_seq_id: residue.label_seq,
                    author_seq_id: residue.author_seq.clone(),
                    insertion_code: residue.insertion.clone(),
                },
                residue_sites.get(&handle).copied().unwrap_or(0),
            )?;
            if let Some(class) = residue.class {
                builder.preserve_residue_class(id, class, residue.class_explicit)?;
            }
            targets.residues.insert(handle, id);
            residue_slots.push(residue.slot);
        }
        for site in self.hierarchy.sites.values() {
            let atom = targets.atoms[&site.atom];
            builder.hierarchy_mut().add_shared_atom_site(
                targets.residues[&site.residue],
                atom,
                Arc::clone(&site.metadata),
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
            chains: self.hierarchy.chains.keys().copied().collect(),
            residues: self.hierarchy.residues.keys().copied().collect(),
        })
    }

    fn record_component(
        &self,
        handles: (&DenseMap<AtomId, EditAtomId>, &DenseMap<BondId, EditBondId>),
        atoms: impl Iterator<Item = (AtomId, AtomId)>,
        bonds: impl Iterator<Item = BondId>,
        instance: MoleculeInstanceId,
        targets: &mut Targets,
    ) {
        for (source, local) in atoms {
            let handle = handles.0[&source];
            let target = InstanceAtomId::new(instance, local);
            targets.atoms.insert(handle, target);
            targets.atom_order[self.atoms[&handle].slot] = Some(target);
        }
        for source in bonds {
            let handle = handles.1[&source];
            targets.bond_slots.push(self.bonds[&handle].slot);
        }
    }
}
