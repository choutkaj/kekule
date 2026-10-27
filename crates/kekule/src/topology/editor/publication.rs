use super::overlay::{Entry, Run};
use super::*;
use crate::topology::dense_order::{DenseOrder, ABSENT};
use crate::topology::hierarchy::{Chain, HierarchyColumns};
use crate::topology::hierarchy_index::HierarchyIndex;
use crate::topology::shared_rows::SharedRows;
use crate::topology::{Hierarchy, MoleculeInstance, TopologyAtomIndex, TopologyBondIndex};
use crate::DraftError;

/// Internal row projection needed to publish a model's structure arrays.
pub(crate) struct TopologyPublication {
    pub topology: Arc<Topology>,
    pub atom_slots: Vec<usize>,
    pub bond_slots: Vec<usize>,
    /// Chain and residue handles in published ID order.
    pub(super) chains: Vec<EditChainId>,
    pub(super) residues: Vec<EditResidueId>,
    /// Source to published dense atoms and bonds and back.
    pub(super) atoms: DenseCorrespondence<TopologyAtomIndex>,
    pub(super) bonds: DenseCorrespondence<TopologyBondIndex>,
}

/// The instances publication stages, with the atoms of the occurrences an
/// edit touched and the bonds of every instance.
#[derive(Default)]
struct StagedInstances {
    instances: Vec<MoleculeInstance>,
    /// The source property row of each instance.
    sources: Vec<Option<usize>>,
    /// Each touched atom's slot, handle and published ID.
    touched: Vec<(usize, EditAtomId, InstanceAtomId)>,
    /// Bonds in published dense order, instance by instance.
    bonds: Vec<InstanceBondId>,
    bond_slots: Vec<usize>,
}

impl StagedInstances {
    fn add(
        &mut self,
        definition: MoleculeDefinitionId,
        source: Option<usize>,
    ) -> MoleculeInstanceId {
        let id = MoleculeInstanceId::new(
            u32::try_from(self.instances.len()).expect("published instances fit their IDs"),
        );
        self.instances.push(MoleculeInstance { id, definition });
        self.sources.push(source);
        id
    }
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
        }
    )*};
}
dense_lookups!(TopologyAtomIndex, TopologyBondIndex, ChainId, ResidueId);

macro_rules! dense_from_sources {
    ($($index:ty),*) => {$(
        impl DenseCorrespondence<$index> {
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
dense_from_sources!(ChainId, ResidueId);

impl<I> DenseCorrespondence<I> {
    fn from_maps(published: Vec<Option<I>>, source: Vec<Option<I>>) -> Self {
        Self { published, source }
    }
}

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
    // Source handles are one block in source dense order, so the dense maps
    // resolve them; added handles are stored.
    atom_block: HandleBlock,
    bond_block: HandleBlock,
    atom_order: Arc<DenseOrder<InstanceAtomId>>,
    bond_order: Arc<DenseOrder<InstanceBondId>>,
    added_atoms: DenseMap<EditAtomId, TopologyAtomIndex>,
    added_bonds: DenseMap<EditBondId, TopologyBondIndex>,
    dense_atoms: DenseCorrespondence<TopologyAtomIndex>,
    dense_bonds: DenseCorrespondence<TopologyBondIndex>,
    residues: DenseCorrespondence<ResidueId>,
    chains: DenseCorrespondence<ChainId>,
}
impl EditCorrespondence {
    /// The published qualified ID and dense index of a draft atom, or `None` when
    /// the handle was deleted or is foreign.
    pub fn atom(&self, id: EditAtomId) -> Option<(InstanceAtomId, TopologyAtomIndex)> {
        let index = match self.atom_block.source(id) {
            Some(source) => self.dense_atoms.published(TopologyAtomIndex::new(source))?,
            None => *self.added_atoms.get(&id)?,
        };
        Some((self.atom_order.order()[index.index()], index))
    }
    /// The published qualified ID and dense index of a draft bond, or `None` when
    /// the handle was deleted or is foreign.
    pub fn bond(&self, id: EditBondId) -> Option<(InstanceBondId, TopologyBondIndex)> {
        let index = match self.bond_block.source(id) {
            Some(source) => self.dense_bonds.published(TopologyBondIndex::new(source))?,
            None => *self.added_bonds.get(&id)?,
        };
        Some((self.bond_order.order()[index.index()], index))
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

/// Stores the handle of every published entity without a source: an added
/// one, whose handle `added` gives by slot in slot order.
fn record_added<K: DenseKey + Copy, I: Copy>(
    sources: &[Option<I>],
    slots: &[usize],
    added: &[(usize, K)],
    index: impl Fn(u32) -> I,
    stored: &mut DenseMap<K, I>,
) {
    for (published, _) in sources.iter().enumerate().filter(|(_, s)| s.is_none()) {
        let slot = slots[published];
        let position = added
            .binary_search_by_key(&slot, |&(slot, _)| slot)
            .expect("published slots are live");
        stored.insert(added[position].1, index(published as u32));
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
        let (atom_block, bond_block) = (self.atoms.block(), self.bonds.block());
        let chain_block = self.hierarchy.chains.block();
        let residue_block = self.hierarchy.residues.block();
        // Added atoms and bonds in handle order, which is also slot order.
        let added_atoms = self
            .atoms
            .added()
            .iter()
            .map(|(&handle, loc)| (loc.slot, handle))
            .collect::<Vec<_>>();
        let added_bonds = self
            .bonds
            .added()
            .iter()
            .map(|(&handle, loc)| (loc.slot, handle))
            .collect::<Vec<_>>();
        let mut published = self.into_publication()?;
        let topology = &published.topology;
        let mut correspondence = EditCorrespondence {
            atom_block,
            bond_block,
            atom_order: Arc::clone(&topology.atom_order),
            bond_order: Arc::clone(&topology.bond_order),
            dense_atoms: std::mem::take(&mut published.atoms),
            dense_bonds: std::mem::take(&mut published.bonds),
            ..EditCorrespondence::default()
        };
        record_added(
            correspondence.dense_atoms.to_source(),
            &published.atom_slots,
            &added_atoms,
            TopologyAtomIndex::new,
            &mut correspondence.added_atoms,
        );
        record_added(
            correspondence.dense_bonds.to_source(),
            &published.bond_slots,
            &added_bonds,
            TopologyBondIndex::new,
            &mut correspondence.added_bonds,
        );
        correspondence.residues = DenseCorrespondence::<ResidueId>::from_sources(
            residue_block.len,
            published
                .residues
                .iter()
                .map(|&handle| residue_block.source(handle).map(ResidueId::new))
                .collect(),
        );
        correspondence.chains = DenseCorrespondence::<ChainId>::from_sources(
            chain_block.len,
            published
                .chains
                .iter()
                .map(|&handle| chain_block.source(handle).map(ChainId::new))
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
    /// Publication costs what the edit touched plus a few passes over plain
    /// per-atom arrays: occurrences, residues and atom sites no edit touched
    /// are carried over from the source by index, and only the touched ones
    /// are classified.
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
                let atom_map = (0..source.atom_count() as u32)
                    .map(|index| Some(TopologyAtomIndex::new(index)))
                    .collect::<Vec<_>>();
                let bond_map = (0..source.bond_count() as u32)
                    .map(|index| Some(TopologyBondIndex::new(index)))
                    .collect::<Vec<_>>();
                return Ok(TopologyPublication {
                    atom_slots: (0..source.atom_count()).collect(),
                    bond_slots: (0..source.bond_count()).collect(),
                    chains: self.hierarchy.chains.keys().collect(),
                    residues: self.hierarchy.residues.keys().collect(),
                    atoms: DenseCorrespondence::from_maps(atom_map.clone(), atom_map),
                    bonds: DenseCorrespondence::from_maps(bond_map.clone(), bond_map),
                    topology: source,
                });
            }
        }
        let source = self.source.clone();
        let source = source.as_deref();
        let sourced = source.is_some() && !self.cleared;
        let mut builder = TopologyBuilder::new();
        let groups = std::mem::take(&mut self.groups);
        let definition_count = source.map_or(0, Topology::definition_count);
        let source_definition = |group: usize| {
            source
                .and_then(|source| source.instance_record(MoleculeInstanceId::new(group as u32)))
                .expect("source groups name source instances")
                .definition()
        };
        // Preserve definition order/reuse for untouched source occurrences.
        let mut retained = vec![false; definition_count];
        let mut invalidated = vec![false; definition_count];
        let mut invalidated_added = BTreeSet::new();
        for (index, group) in groups.iter().enumerate() {
            match group {
                GroupSlot::Source => retained[source_definition(index).index()] = true,
                GroupSlot::Live(group) => match &group.chemistry {
                    GroupChemistry::Source(id) => {
                        retained[id.index()] = true;
                        invalidated[id.index()] |= group.class.is_none();
                    }
                    GroupChemistry::Added(molecule) if group.class.is_none() => {
                        invalidated_added.insert(Arc::as_ptr(molecule));
                    }
                    _ => {}
                },
                GroupSlot::Merged => {}
            }
        }
        let mut source_definitions = vec![None; definition_count];
        if let Some(source) = source {
            // Untouched definitions share the source molecule, whether or not
            // the source snapshot is still shared.
            let overrides = &source.molecule_class_overrides;
            let kept = source
                .definitions
                .iter()
                .filter(|definition| retained[definition.id().index()])
                .collect::<Vec<_>>();
            for (target, definition) in kept.iter().enumerate() {
                source_definitions[definition.id().index()] =
                    Some(MoleculeDefinitionId::new(target as u32));
            }
            builder.extend_shared_definitions(kept.into_iter().map(|definition| {
                let id = definition.id();
                let explicit = overrides.contains_key(&id);
                let class = (!invalidated[id.index()] || explicit)
                    .then_some((definition.class(), explicit));
                (
                    Arc::clone(&definition.molecule),
                    definition.default_perception,
                    class,
                )
            }))?;
        }
        let mut added_definitions = BTreeMap::new();
        let mut staged = StagedInstances::default();
        let mut plain = vec![None; source.map_or(0, Topology::instance_count)];
        for (index, group) in groups.into_iter().enumerate() {
            let group = match group {
                GroupSlot::Merged => continue,
                GroupSlot::Source => {
                    // An untouched occurrence keeps its local IDs, and its bonds
                    // their source slots.
                    let source = source.expect("source groups have a source");
                    let definition = source_definitions[source_definition(index).index()]
                        .expect("retained source definition");
                    let instance = staged.add(definition, Some(index));
                    plain[index] = Some(instance);
                    for (local, &dense) in
                        source.bond_order.instance_lookup(index).iter().enumerate()
                    {
                        if dense != ABSENT {
                            staged
                                .bonds
                                .push(InstanceBondId::new(instance, BondId::new(local as u32)));
                            staged.bond_slots.push(dense as usize);
                        }
                    }
                    continue;
                }
                GroupSlot::Live(group) => *group,
            };
            if group.atoms.is_empty() {
                continue;
            }
            let handles = (&group.atoms, &group.bonds);
            match group.chemistry {
                GroupChemistry::Source(id) => {
                    let definition =
                        source_definitions[id.index()].expect("retained source definition");
                    let instance = staged.add(definition, group.instance_slot);
                    let molecule = builder.require_definition(definition)?.molecule();
                    self.record_component(
                        handles,
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids().map(|b| (b, b)),
                        instance,
                        &mut staged,
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
                    let instance = staged.add(definition, group.instance_slot);
                    let molecule = builder.require_definition(definition)?.molecule();
                    self.record_component(
                        handles,
                        molecule.atom_ids().map(|a| (a, a)),
                        molecule.bond_ids().map(|b| (b, b)),
                        instance,
                        &mut staged,
                    );
                }
                GroupChemistry::Draft(draft) if draft.is_connected() => {
                    // Publication renumbers deleted draft slots densely.
                    let (published, ids) = (*draft).finish_with_correspondence()?;
                    let definition = builder.add_molecule_definition(published)?;
                    let atoms = group.atoms.values().copied().collect();
                    if let Some(class) = self.component_class(&atoms, group.class) {
                        let explicit = group.class_explicit
                            || self.molecule_classes.keys().any(|id| atoms.contains(id));
                        builder.preserve_molecule_class(definition, class, explicit)?;
                    }
                    let instance = staged.add(
                        definition,
                        (!group.changed).then_some(group.instance_slot).flatten(),
                    );
                    self.record_component(
                        handles,
                        ids.atoms().iter().map(|(&draft, &local)| (draft, local)),
                        ids.bonds().iter().map(|(&draft, &local)| (draft, local)),
                        instance,
                        &mut staged,
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
                        let atoms = component
                            .source_atoms
                            .iter()
                            .map(|id| group.atoms[id])
                            .collect();
                        if let Some(class) = self.component_class(&atoms, None) {
                            builder.set_molecule_class(component.id, class)?;
                        }
                        // No old occurrence annotation is inherited across a split.
                        let instance = staged.add(component.id, None);
                        self.record_component(
                            handles,
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
                            &mut staged,
                        );
                    }
                }
            }
        }
        // Publish atoms in slot order: survivors keep their source dense
        // order, and added atoms follow in creation order. Touched atoms are
        // in slot order too.
        staged.touched.sort_unstable_by_key(|&(slot, _, _)| slot);
        let atom_count = self.atoms.len();
        let mut order = Vec::with_capacity(atom_count);
        let mut atom_slots = Vec::with_capacity(atom_count);
        let mut atom_published = vec![None; self.atoms.block().len];
        let mut atom_sources = Vec::with_capacity(atom_count);
        let mut next = staged.touched.iter().peekable();
        if let (Some(source), true) = (source, sourced) {
            for (slot, id) in source.atom_ids().iter().enumerate() {
                let target = if let Some(instance) = plain[id.instance().index()] {
                    InstanceAtomId::new(instance, id.atom())
                } else if let Some(&(_, _, target)) = next.next_if(|touched| touched.0 == slot) {
                    target
                } else {
                    continue;
                };
                atom_published[slot] = Some(TopologyAtomIndex::new(order.len() as u32));
                atom_sources.push(Some(TopologyAtomIndex::new(slot as u32)));
                order.push(target);
                atom_slots.push(slot);
            }
        }
        for &(slot, _, target) in next {
            atom_sources.push(None);
            order.push(target);
            atom_slots.push(slot);
        }
        let slot_counts = |slots: fn(&Molecule) -> usize| {
            staged
                .instances
                .iter()
                .map(|instance| {
                    slots(
                        builder
                            .definition(instance.definition())
                            .expect("staged definition")
                            .molecule(),
                    )
                })
                .collect::<Vec<_>>()
        };
        let atom_order =
            DenseOrder::from_order(order, slot_counts(|m| m.graph().atom_slot_count()));
        let bond_order = DenseOrder::from_order(
            std::mem::take(&mut staged.bonds),
            slot_counts(|m| m.graph().bond_slot_count()),
        );
        let bond_block = self.bonds.block();
        let mut bond_published = vec![None; bond_block.len];
        let bond_sources = staged
            .bond_slots
            .iter()
            .enumerate()
            .map(|(published, &slot)| {
                (sourced && slot < bond_block.len).then(|| {
                    bond_published[slot] = Some(TopologyBondIndex::new(published as u32));
                    TopologyBondIndex::new(slot as u32)
                })
            })
            .collect();
        let published =
            self.publish_hierarchy(source, (&atom_published, &atom_slots, atom_order.order()));
        let index = HierarchyIndex::from_site_atoms(
            published.site_dense,
            atom_order.len(),
            staged.instances.len(),
        );
        let mut properties = self.properties.project(
            &staged.sources,
            &atom_slots,
            &staged.bond_slots,
            &published.chain_slots,
            &published.residue_slots,
            &published.site_slots,
        )?;
        for (key, value) in self.properties.owner().iter() {
            properties.owner_mut().insert(key.clone(), value.clone())?;
        }
        builder.install_properties(properties);
        builder.install_instances(staged.instances, atom_order, bond_order);
        let target =
            Arc::new(builder.publish_prepared(published.hierarchy, index, &published.infer)?);
        Ok(TopologyPublication {
            topology: target,
            atom_slots,
            bond_slots: staged.bond_slots,
            chains: published.chains,
            residues: published.residues,
            atoms: DenseCorrespondence::from_maps(atom_published, atom_sources),
            bonds: DenseCorrespondence::from_maps(bond_published, bond_sources),
        })
    }

    /// The published hierarchy: chains, residues and atom sites in handle
    /// order, each residue's sites in site handle order and each chain's
    /// residues in residue handle order. Unchanged source nodes carry their
    /// labels and metadata over by key. `atoms` gives each live source atom's
    /// published index, each published atom's slot, and the published atoms.
    fn publish_hierarchy(
        &self,
        source: Option<&Topology>,
        atoms: (&[Option<TopologyAtomIndex>], &[usize], &[InstanceAtomId]),
    ) -> PublishedHierarchy {
        let (atom_published, atom_slots, order) = atoms;
        let hierarchy = &self.hierarchy;
        let source_hierarchy = source.map(Topology::hierarchy);
        let mut labels =
            source_hierarchy.map_or_else(SharedRows::default, |source| source.label_rows().share());
        let mut metadata = source_hierarchy
            .map_or_else(SharedRows::default, |source| source.metadata_rows().share());
        let source_of = || source_hierarchy.expect("source nodes have a source");

        let chain_block = hierarchy.chains.block();
        let mut chain_ids = vec![u32::MAX; chain_block.len];
        let mut added_chains = BTreeMap::new();
        let mut chains = Vec::with_capacity(hierarchy.chains.len());
        let mut chain_handles = Vec::with_capacity(hierarchy.chains.len());
        let mut chain_slots = Vec::with_capacity(hierarchy.chains.len());
        for (handle, entry) in hierarchy.chains.entries() {
            let id = chains.len() as u32;
            let (label_id, author_id, slot) = match entry {
                Entry::Source(index) => {
                    let chain = source_of()
                        .chain(ChainId::new(index as u32))
                        .expect("source chain");
                    chain_ids[index] = id;
                    (chain.label_id.clone(), chain.author_id.clone(), index)
                }
                Entry::Stored(chain) => {
                    match chain_block.source(handle) {
                        Some(index) => chain_ids[index as usize] = id,
                        None => {
                            added_chains.insert(handle, id);
                        }
                    }
                    (chain.label.clone(), chain.author.clone(), chain.slot)
                }
            };
            chains.push(Chain {
                id: ChainId::new(id),
                label_id,
                author_id,
                residues: Vec::new(),
            });
            chain_handles.push(handle);
            chain_slots.push(slot);
        }
        let chain_of = |handle: EditChainId| match chain_block.source(handle) {
            Some(index) => chain_ids[index as usize],
            None => added_chains[&handle],
        };

        let residue_block = hierarchy.residues.block();
        let residue_count = hierarchy.residues.len();
        let mut residue_ids = vec![u32::MAX; residue_block.len];
        let mut added_residues = BTreeMap::new();
        let mut residue_chains = Vec::with_capacity(residue_count);
        let mut residue_classes = Vec::with_capacity(residue_count);
        let mut residue_explicit = Vec::with_capacity(residue_count);
        let mut residue_labels = Vec::with_capacity(residue_count);
        let mut residue_handles = Vec::with_capacity(residue_count);
        let mut residue_slots = Vec::with_capacity(residue_count);
        let mut infer = Vec::new();
        for (handle, entry) in hierarchy.residues.entries() {
            let id = residue_chains.len() as u32;
            let (chain, class, explicit, key, slot) = match entry {
                Entry::Source(index) => {
                    let residue = source_of()
                        .residue(ResidueId::new(index as u32))
                        .expect("source residue");
                    residue_ids[index] = id;
                    (
                        chain_ids[residue.chain().index()],
                        Some(residue.class()),
                        residue.explicit_class(),
                        residue.label_key(),
                        index,
                    )
                }
                Entry::Stored(residue) => {
                    match residue_block.source(handle) {
                        Some(index) => residue_ids[index as usize] = id,
                        None => {
                            added_residues.insert(handle, id);
                        }
                    }
                    (
                        chain_of(residue.chain),
                        residue.class,
                        residue.class.filter(|_| residue.class_explicit),
                        labels.key(&residue.labels),
                        residue.slot,
                    )
                }
            };
            chains[chain as usize].residues.push(ResidueId::new(id));
            residue_chains.push(ChainId::new(chain));
            residue_classes.push(class.unwrap_or(ResidueClass::Other));
            if class.is_none() {
                infer.push(ResidueId::new(id));
            }
            residue_explicit.push(explicit);
            residue_labels.push(key);
            residue_handles.push(handle);
            residue_slots.push(slot);
        }
        let residue_of = |handle: EditResidueId| match residue_block.source(handle) {
            Some(index) => residue_ids[index as usize],
            None => added_residues[&handle],
        };

        let site_count = hierarchy.sites.len();
        let mut site_residues = Vec::with_capacity(site_count);
        let mut site_atoms = Vec::with_capacity(site_count);
        let mut site_dense = Vec::with_capacity(site_count);
        let mut site_metadata = Vec::with_capacity(site_count);
        let mut site_slots = Vec::with_capacity(site_count);
        for run in hierarchy.sites.runs() {
            match run {
                Run::Source(range) => {
                    // Unchanged source sites are copied a run at a time.
                    let source = source.expect("source sites have a source");
                    let site = source_of();
                    site_residues.extend(
                        site.site_residue_ids()[range.clone()]
                            .iter()
                            .map(|residue| ResidueId::new(residue_ids[residue.index()])),
                    );
                    for atom in &source.site_atoms()[range.clone()] {
                        let dense = atom_published[atom.index()].expect("site atoms are live");
                        site_dense.push(dense);
                        site_atoms.push(order[dense.index()]);
                    }
                    site_metadata.extend_from_slice(&site.site_metadata_keys()[range.clone()]);
                    site_slots.extend(range);
                }
                Run::Stored(_, site) => {
                    // Published atoms are in slot order.
                    let atom_slot = self.atom_location(site.atom).expect("live site atom").slot;
                    let dense = atom_slots
                        .binary_search(&atom_slot)
                        .expect("site atoms are published");
                    site_residues.push(ResidueId::new(residue_of(site.residue)));
                    site_atoms.push(order[dense]);
                    site_dense.push(TopologyAtomIndex::new(dense as u32));
                    site_metadata.push(metadata.key(&site.metadata));
                    site_slots.push(site.slot);
                }
            }
        }
        // Each residue lists its sites in site order.
        let mut residue_site_starts = vec![0u32; residue_chains.len() + 1];
        for residue in &site_residues {
            residue_site_starts[residue.index() + 1] += 1;
        }
        for residue in 0..residue_chains.len() {
            residue_site_starts[residue + 1] += residue_site_starts[residue];
        }
        let mut next = residue_site_starts.clone();
        let mut residue_site_items = vec![AtomSiteId::new(0); site_residues.len()];
        for (site, residue) in site_residues.iter().enumerate() {
            residue_site_items[next[residue.index()] as usize] = AtomSiteId::new(site as u32);
            next[residue.index()] += 1;
        }
        PublishedHierarchy {
            hierarchy: Hierarchy::from_columns(HierarchyColumns {
                chains,
                residue_chains,
                residue_classes,
                residue_explicit,
                residue_labels,
                labels,
                residue_site_starts,
                residue_site_items,
                site_residues,
                site_atoms,
                site_metadata,
                metadata,
            }),
            site_dense,
            infer,
            chains: chain_handles,
            residues: residue_handles,
            chain_slots,
            residue_slots,
            site_slots,
        }
    }

    fn record_component(
        &self,
        handles: (&DenseMap<AtomId, EditAtomId>, &DenseMap<BondId, EditBondId>),
        atoms: impl Iterator<Item = (AtomId, AtomId)>,
        bonds: impl Iterator<Item = (BondId, BondId)>,
        instance: MoleculeInstanceId,
        staged: &mut StagedInstances,
    ) {
        for (source, local) in atoms {
            let handle = handles.0[&source];
            let slot = self.atom_location(handle).expect("live atom").slot;
            staged
                .touched
                .push((slot, handle, InstanceAtomId::new(instance, local)));
        }
        for (source, local) in bonds {
            let handle = handles.1[&source];
            staged.bonds.push(InstanceBondId::new(instance, local));
            staged
                .bond_slots
                .push(self.bond_location(handle).expect("live bond").slot);
        }
    }
}

/// A published hierarchy with each site's dense atom, the residues to
/// classify, and the draft handles and slots of its nodes in published order.
struct PublishedHierarchy {
    hierarchy: Hierarchy,
    site_dense: Vec<TopologyAtomIndex>,
    infer: Vec<ResidueId>,
    chains: Vec<EditChainId>,
    residues: Vec<EditResidueId>,
    chain_slots: Vec<usize>,
    residue_slots: Vec<usize>,
    site_slots: Vec<usize>,
}
