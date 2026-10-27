use super::*;

/// Coordinate-free chain metadata in a draft. Use editor methods to change it.
#[derive(Debug, Clone)]
pub struct EditChain {
    pub(super) label: String,
    pub(super) author: Option<String>,
    pub(super) slot: usize,
}
impl EditChain {
    /// The label chain identifier.
    pub fn label_id(&self) -> &str {
        &self.label
    }
    /// The author chain identifier, or `None`.
    pub fn author_id(&self) -> Option<&str> {
        self.author.as_deref()
    }
}
/// Residue metadata whose chain and atom-site identities are editing handles.
#[derive(Debug, Clone)]
pub struct EditResidue {
    pub(super) chain: EditChainId,
    // Shared with the source and published snapshots that keep them.
    pub(super) labels: SharedRow<ResidueLabels>,
    pub(super) class: Option<ResidueClass>,
    pub(super) class_explicit: bool,
    pub(super) slot: usize,
}
impl EditResidue {
    /// The handle of the chain containing this residue.
    pub fn chain(&self) -> EditChainId {
        self.chain
    }
    /// The residue (component) name.
    pub fn name(&self) -> &str {
        &self.labels.name
    }
    /// The label component identifier, or `None`.
    pub fn label_comp_id(&self) -> Option<&str> {
        self.labels.label_comp_id.as_deref()
    }
    /// The author component identifier, or `None`.
    pub fn author_comp_id(&self) -> Option<&str> {
        self.labels.author_comp_id.as_deref()
    }
    /// The label sequence number, or `None`.
    pub fn label_seq_id(&self) -> Option<i32> {
        self.labels.label_seq_id
    }
    /// The author sequence identifier, or `None`.
    pub fn author_seq_id(&self) -> Option<&str> {
        self.labels.author_seq_id.as_deref()
    }
    /// The insertion code, or `None`.
    pub fn insertion_code(&self) -> Option<&str> {
        self.labels.insertion_code.as_deref()
    }
    /// Returns an explicit assignment, if present. `None` retains an unchanged
    /// inferred class or requests inference when its evidence has changed.
    pub fn class_override(&self) -> Option<ResidueClass> {
        self.class.filter(|_| self.class_explicit)
    }
}
/// Atom-site organization referring to a stable editing atom handle.
#[derive(Debug, Clone)]
pub struct EditAtomSite {
    pub(super) residue: EditResidueId,
    pub(super) atom: EditAtomId,
    // Shared with the source and published snapshots that keep it.
    pub(super) metadata: SharedRow<AtomSiteMetadata>,
    pub(super) slot: usize,
}
impl EditAtomSite {
    /// The handle of the residue containing this site.
    pub fn residue(&self) -> EditResidueId {
        self.residue
    }
    /// The handle of the atom placed by this site.
    pub fn atom(&self) -> EditAtomId {
        self.atom
    }
    /// Source atom-site labels.
    pub fn metadata(&self) -> &AtomSiteMetadata {
        &self.metadata
    }
}

/// The draft hierarchy: source chains, residues and atom sites are read from
/// the source until an edit changes them.
#[derive(Debug, Clone, Default)]
pub(super) struct EditHierarchy {
    pub(super) chains: Overlay<EditChainId, EditChain>,
    pub(super) residues: Overlay<EditResidueId, EditResidue>,
    pub(super) sites: Overlay<EditAtomSiteId, EditAtomSite>,
    // The site of each atom whose site changed, `None` for a removed one;
    // other source atoms keep their source site.
    changed_atom_sites: BTreeMap<EditAtomId, Option<EditAtomSiteId>>,
    added_atom_sites: DenseMap<EditAtomId, EditAtomSiteId>,
    // Whether the draft was cleared, so no source atom keeps a source site.
    cleared: bool,
    // The atoms of residues whose membership changed; other source residues
    // keep their source atoms, read through `residue_atom_list`.
    pub(super) residue_atoms: DenseMap<EditResidueId, BTreeSet<EditAtomId>>,
}

impl EditHierarchy {
    pub(super) fn clear(&mut self) {
        self.chains.clear();
        self.residues.clear();
        self.sites.clear();
        self.changed_atom_sites.clear();
        self.added_atom_sites.clear();
        self.cleared = true;
        self.residue_atoms.clear();
    }
}

/// Reads unchanged source hierarchy nodes as draft nodes.
#[derive(Clone, Copy)]
pub(super) struct SourceHierarchy<'a> {
    source: Option<&'a Topology>,
    atoms: HandleBlock,
    chains: HandleBlock,
    residues: HandleBlock,
}

impl<'a> SourceHierarchy<'a> {
    pub(super) fn new(
        source: Option<&'a Topology>,
        atoms: HandleBlock,
        hierarchy: &EditHierarchy,
    ) -> Self {
        Self {
            source,
            atoms,
            chains: hierarchy.chains.block(),
            residues: hierarchy.residues.block(),
        }
    }

    fn source(self) -> &'a Topology {
        self.source.expect("source nodes have a source")
    }

    pub(super) fn chain(self, index: usize) -> EditChain {
        let chain = self
            .source()
            .hierarchy()
            .chain(ChainId::new(index as u32))
            .expect("source chain");
        EditChain {
            label: chain.label_id().to_owned(),
            author: chain.author_id().map(str::to_owned),
            slot: index,
        }
    }

    pub(super) fn residue(self, index: usize) -> EditResidue {
        let residue = self
            .source()
            .hierarchy()
            .residue(ResidueId::new(index as u32))
            .expect("source residue");
        EditResidue {
            chain: self.chains.handle(residue.chain().index()),
            labels: residue.shared_labels(),
            class: Some(residue.class()),
            class_explicit: residue.explicit_class().is_some(),
            slot: index,
        }
    }

    pub(super) fn site(self, index: usize) -> EditAtomSite {
        let source = self.source();
        let site = source
            .hierarchy()
            .atom_site(AtomSiteId::new(index as u32))
            .expect("source site");
        EditAtomSite {
            residue: self.residues.handle(site.residue().index()),
            atom: self.atoms.handle(source.site_atoms()[index].index()),
            metadata: site.shared_metadata(),
            slot: index,
        }
    }
}

impl TopologyEditor {
    /// The reader of unchanged source hierarchy nodes.
    pub(super) fn source_hierarchy(&self) -> SourceHierarchy<'_> {
        SourceHierarchy {
            source: self.source.as_deref(),
            atoms: self.atoms.block(),
            chains: self.hierarchy.chains.block(),
            residues: self.hierarchy.residues.block(),
        }
    }
    /// Live draft chains in handle order.
    pub fn chains(&self) -> impl ExactSizeIterator<Item = (EditChainId, EditChain)> + '_ {
        let source = self.source_hierarchy();
        self.hierarchy.chains.iter(move |index| source.chain(index))
    }
    /// Live draft residues in handle order.
    pub fn residues(&self) -> impl ExactSizeIterator<Item = (EditResidueId, EditResidue)> + '_ {
        let source = self.source_hierarchy();
        self.hierarchy
            .residues
            .iter(move |index| source.residue(index))
    }
    /// Live draft atom sites in handle order.
    pub fn atom_sites(&self) -> impl ExactSizeIterator<Item = (EditAtomSiteId, EditAtomSite)> + '_ {
        let source = self.source_hierarchy();
        self.hierarchy.sites.iter(move |index| source.site(index))
    }
    /// One draft chain, or `None` when `id` names no live chain of this draft.
    pub fn chain(&self, id: EditChainId) -> Option<EditChain> {
        let source = self.source_hierarchy();
        self.hierarchy.chains.get(&id, |index| source.chain(index))
    }
    /// One draft residue, or `None` when `id` names no live residue of this draft.
    pub fn residue(&self, id: EditResidueId) -> Option<EditResidue> {
        let source = self.source_hierarchy();
        self.hierarchy
            .residues
            .get(&id, |index| source.residue(index))
    }
    /// One draft atom site, or `None` when `id` names no live site of this draft.
    pub fn atom_site(&self, id: EditAtomSiteId) -> Option<EditAtomSite> {
        let source = self.source_hierarchy();
        self.hierarchy.sites.get(&id, |index| source.site(index))
    }
    /// The draft handle of a source chain, or `None` when the chain is not in
    /// the source or was deleted.
    pub fn chain_handle(&self, source: ChainId) -> Option<EditChainId> {
        let block = self.hierarchy.chains.block();
        (source.index() < block.len)
            .then(|| block.handle(source.index()))
            .filter(|id| self.hierarchy.chains.contains_key(id))
    }
    /// The draft handle of a source residue, or `None` when the residue is not
    /// in the source or was deleted.
    pub fn residue_handle(&self, source: ResidueId) -> Option<EditResidueId> {
        let block = self.hierarchy.residues.block();
        (source.index() < block.len)
            .then(|| block.handle(source.index()))
            .filter(|id| self.hierarchy.residues.contains_key(id))
    }
    pub(crate) fn require_chain(&self, id: EditChainId) -> Result<EditChain, TopologyEditError> {
        self.chain(id).ok_or(TopologyEditError::InvalidChain(id))
    }
    pub(crate) fn require_residue(
        &self,
        id: EditResidueId,
    ) -> Result<EditResidue, TopologyEditError> {
        self.residue(id)
            .ok_or(TopologyEditError::InvalidResidue(id))
    }
    pub(crate) fn require_atom_site(
        &self,
        id: EditAtomSiteId,
    ) -> Result<EditAtomSite, TopologyEditError> {
        self.atom_site(id)
            .ok_or(TopologyEditError::InvalidAtomSite(id))
    }

    /// The draft handle of a source atom site.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidSourceAtomSite`] when the site is not in the
    /// source or was deleted.
    pub fn atom_site_handle(
        &self,
        source: AtomSiteId,
    ) -> Result<EditAtomSiteId, TopologyEditError> {
        let block = self.hierarchy.sites.block();
        (source.index() < block.len)
            .then(|| block.handle(source.index()))
            .filter(|id| self.hierarchy.sites.contains_key(id))
            .ok_or(TopologyEditError::InvalidSourceAtomSite(source))
    }
    /// The site placing `atom`, or `None` when the atom has no site.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtom`] when `atom` is deleted or foreign.
    pub fn atom_site_for_atom(
        &self,
        atom: EditAtomId,
    ) -> Result<Option<EditAtomSiteId>, TopologyEditError> {
        self.atom_location(atom)?;
        Ok(self.site_of_atom(atom))
    }
    /// The site placing a live atom, or `None`.
    pub(super) fn site_of_atom(&self, atom: EditAtomId) -> Option<EditAtomSiteId> {
        let hierarchy = &self.hierarchy;
        let Some(index) = self.atoms.block().source(atom) else {
            return hierarchy.added_atom_sites.get(&atom).copied();
        };
        if let Some(&site) = hierarchy.changed_atom_sites.get(&atom) {
            return site;
        }
        if hierarchy.cleared {
            return None;
        }
        self.source_ref()
            .hierarchy_index
            .atom_site(crate::topology::TopologyAtomIndex::new(index))
            .map(|site| hierarchy.sites.block().handle(site.index()))
    }
    /// Records the site placing a live atom.
    fn set_site_of_atom(&mut self, atom: EditAtomId, site: Option<EditAtomSiteId>) {
        let hierarchy = &mut self.hierarchy;
        if self.atoms.block().source(atom).is_some() {
            hierarchy.changed_atom_sites.insert(atom, site);
        } else if let Some(site) = site {
            hierarchy.added_atom_sites.insert(atom, site);
        } else {
            hierarchy.added_atom_sites.remove(&atom);
        }
    }
    /// Adds an empty chain with a label and optional author identifier.
    ///
    /// # Errors
    ///
    /// Currently never fails.
    pub fn add_chain(
        &mut self,
        label: impl Into<String>,
        author: Option<String>,
    ) -> Result<EditChainId, TopologyEditError> {
        let id = EditChainId::new();
        let slot = self.properties.chains().len();
        self.hierarchy.chains.insert(
            id,
            EditChain {
                label: label.into(),
                author,
                slot,
            },
        );
        self.properties.chains_mut().resize_missing(slot + 1);
        self.changed();
        Ok(id)
    }
    /// Appends a residue to `chain`. `name` also becomes the label component ID;
    /// the class is inferred at publication unless assigned.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidChain`] when `chain` is deleted or foreign.
    pub fn add_residue(
        &mut self,
        chain: EditChainId,
        name: impl Into<String>,
        label_seq: Option<i32>,
        author_seq: Option<String>,
        insertion: Option<String>,
    ) -> Result<EditResidueId, TopologyEditError> {
        self.require_chain(chain)?;
        let id = EditResidueId::new();
        let name = Arc::<str>::from(name.into());
        let slot = self.properties.residues().len();
        self.hierarchy.residues.insert(
            id,
            EditResidue {
                chain,
                labels: SharedRow::new(ResidueLabels {
                    label_comp_id: Some(Arc::clone(&name)),
                    name,
                    author_comp_id: None,
                    label_seq_id: label_seq,
                    author_seq_id: author_seq.map(Arc::from),
                    insertion_code: insertion.map(Arc::from),
                }),
                class: None,
                class_explicit: false,
                slot,
            },
        );
        self.hierarchy.residue_atoms.insert(id, BTreeSet::new());
        self.properties.residues_mut().resize_missing(slot + 1);
        self.changed();
        Ok(id)
    }
    /// Places `atom` in `residue`. Invalidates inferred classes the placement
    /// affects.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidResidue`] or [`TopologyEditError::InvalidAtom`]
    /// for deleted or foreign handles, or
    /// [`TopologyEditError::DuplicateAtomPlacement`] when `atom` already has a site.
    pub fn add_atom_site(
        &mut self,
        residue: EditResidueId,
        atom: EditAtomId,
        metadata: AtomSiteMetadata,
    ) -> Result<EditAtomSiteId, TopologyEditError> {
        self.require_residue(residue)?;
        self.atom_location(atom)?;
        if self.site_of_atom(atom).is_some() {
            return Err(TopologyEditError::DuplicateAtomPlacement(atom));
        }
        let id = EditAtomSiteId::new();
        let slot = self.properties.atom_sites().len();
        self.hierarchy.sites.insert(
            id,
            EditAtomSite {
                residue,
                atom,
                metadata: SharedRow::new(metadata),
                slot,
            },
        );
        self.set_site_of_atom(atom, Some(id));
        self.residue_atoms_mut(residue).insert(atom);
        self.invalidate_molecule_class_for_atom(atom);
        self.invalidate_residue_class(residue, true);
        self.properties.atom_sites_mut().resize_missing(slot + 1);
        self.changed();
        Ok(id)
    }
    /// The stored chain `id`, which the caller checked is live.
    fn chain_mut(&mut self, id: EditChainId) -> &mut EditChain {
        let source = SourceHierarchy {
            source: self.source.as_deref(),
            ..self.source_hierarchy_blocks()
        };
        self.hierarchy
            .chains
            .get_mut(&id, |index| source.chain(index))
            .expect("live chain")
    }
    /// The stored residue `id`, which the caller checked is live.
    pub(super) fn residue_mut(&mut self, id: EditResidueId) -> &mut EditResidue {
        let source = SourceHierarchy {
            source: self.source.as_deref(),
            ..self.source_hierarchy_blocks()
        };
        self.hierarchy
            .residues
            .get_mut(&id, |index| source.residue(index))
            .expect("live residue")
    }
    /// The stored site `id`, which the caller checked is live.
    fn site_mut(&mut self, id: EditAtomSiteId) -> &mut EditAtomSite {
        let source = SourceHierarchy {
            source: self.source.as_deref(),
            ..self.source_hierarchy_blocks()
        };
        self.hierarchy
            .sites
            .get_mut(&id, |index| source.site(index))
            .expect("live atom site")
    }
    /// The handle blocks of [`Self::source_hierarchy`], without the source.
    fn source_hierarchy_blocks(&self) -> SourceHierarchy<'static> {
        SourceHierarchy {
            source: None,
            atoms: self.atoms.block(),
            chains: self.hierarchy.chains.block(),
            residues: self.hierarchy.residues.block(),
        }
    }
    /// Replaces a chain's label and author identifiers.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidChain`] when `id` is deleted or foreign.
    pub fn set_chain_identifiers(
        &mut self,
        id: EditChainId,
        label: impl Into<String>,
        author: Option<String>,
    ) -> Result<(), TopologyEditError> {
        let chain = self.require_chain(id)?;
        let label = label.into();
        if chain.label == label && chain.author == author {
            return Ok(());
        }
        let chain = self.chain_mut(id);
        chain.label = label;
        chain.author = author;
        self.changed();
        Ok(())
    }
    /// Replaces a residue's label and author component identifiers, invalidating
    /// its inferred class.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidResidue`] when `id` is deleted or foreign.
    pub fn set_residue_component_ids(
        &mut self,
        id: EditResidueId,
        label: Option<String>,
        author: Option<String>,
    ) -> Result<(), TopologyEditError> {
        let residue = self.require_residue(id)?;
        if residue.label_comp_id() == label.as_deref()
            && residue.author_comp_id() == author.as_deref()
        {
            return Ok(());
        }
        self.residue_mut(id).labels = SharedRow::new(ResidueLabels {
            label_comp_id: label.map(Arc::from),
            author_comp_id: author.map(Arc::from),
            ..(*residue.labels).clone()
        });
        self.invalidate_residue_class(id, false);
        self.changed();
        Ok(())
    }
    /// Explicitly assigns a residue's class, overriding inference at publication.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidResidue`] when `id` is deleted or foreign.
    pub fn set_residue_class(
        &mut self,
        id: EditResidueId,
        class: ResidueClass,
    ) -> Result<(), TopologyEditError> {
        let residue = self.require_residue(id)?;
        let previous = residue.class;
        if residue.class == Some(class) && residue.class_explicit {
            return Ok(());
        }
        let residue = self.residue_mut(id);
        residue.class = Some(class);
        residue.class_explicit = true;
        if previous != Some(class) {
            self.invalidate_molecule_classes_for_residue(id);
        }
        self.revision += 1;
        Ok(())
    }
    /// Clears one residue's class, explicit or retained, so publication
    /// infers it afresh.
    pub fn clear_residue_class(&mut self, id: EditResidueId) -> Result<(), TopologyEditError> {
        if self.require_residue(id)?.class.is_none() {
            return Ok(());
        }
        let residue = self.residue_mut(id);
        residue.class = None;
        residue.class_explicit = false;
        self.invalidate_molecule_classes_for_residue(id);
        self.revision += 1;
        Ok(())
    }
    /// Replaces one atom site's source metadata. Changed atom names invalidate
    /// the inferred class of the atom's molecule.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtomSite`] when `id` is deleted or foreign.
    pub fn set_atom_site_metadata(
        &mut self,
        id: EditAtomSiteId,
        metadata: AtomSiteMetadata,
    ) -> Result<(), TopologyEditError> {
        let site = self.require_atom_site(id)?;
        if *site.metadata == metadata {
            return Ok(());
        }
        let names_changed = site.metadata.label_atom_id != metadata.label_atom_id
            || site.metadata.auth_atom_id != metadata.auth_atom_id;
        self.site_mut(id).metadata = SharedRow::new(metadata);
        if names_changed {
            self.invalidate_molecule_class_for_atom(site.atom);
        }
        self.changed();
        Ok(())
    }
    /// Moves an atom site to another residue.
    ///
    /// # Errors
    ///
    /// [`TopologyEditError::InvalidAtomSite`] or [`TopologyEditError::InvalidResidue`]
    /// for deleted or foreign handles.
    pub fn set_atom_site_residue(
        &mut self,
        id: EditAtomSiteId,
        residue: EditResidueId,
    ) -> Result<(), TopologyEditError> {
        self.require_residue(residue)?;
        let site = self.require_atom_site(id)?;
        let previous = site.residue;
        if previous == residue {
            return Ok(());
        }
        self.site_mut(id).residue = residue;
        let atom = site.atom;
        self.residue_atoms_mut(previous).remove(&atom);
        self.residue_atoms_mut(residue).insert(atom);
        self.invalidate_molecule_class_for_atom(atom);
        self.invalidate_residue_class(residue, true);
        self.invalidate_residue_class(previous, true);
        self.prune_empty_residue(previous);
        self.changed();
        Ok(())
    }
    /// Removes organization only; the chemical atom remains in the system.
    pub fn delete_atom_site(
        &mut self,
        id: EditAtomSiteId,
    ) -> Result<EditAtomSite, TopologyEditError> {
        self.require_atom_site(id)?;
        let source = self.source_hierarchy_blocks();
        let source = SourceHierarchy {
            source: self.source.as_deref(),
            ..source
        };
        let site = self
            .hierarchy
            .sites
            .remove(&id, |index| source.site(index))
            .unwrap();
        self.set_site_of_atom(site.atom, None);
        self.residue_atoms_mut(site.residue).remove(&site.atom);
        self.properties.atom_sites_mut().clear_index(site.slot);
        self.invalidate_residue_class(site.residue, true);
        self.invalidate_molecule_class_for_atom(site.atom);
        self.prune_empty_residue(site.residue);
        self.changed();
        Ok(site)
    }
    /// Removes a residue and its sites, retaining all chemical atoms.
    pub fn delete_residue(&mut self, id: EditResidueId) -> Result<EditResidue, TopologyEditError> {
        let previous = self.require_residue(id)?;
        let mut sites = self
            .residue_atom_list(id)
            .into_iter()
            .filter_map(|atom| self.site_of_atom(atom))
            .collect::<Vec<_>>();
        sites.sort_unstable();
        for site in sites {
            self.delete_atom_site(site)?;
        }
        // Empty residues can exist when deliberately added during assembly.
        if self.hierarchy.residues.contains_key(&id) {
            self.prune_empty_residue(id);
            self.changed();
        }
        Ok(previous)
    }
    /// Removes a chain and its organization, retaining all chemical atoms.
    pub fn delete_chain(&mut self, id: EditChainId) -> Result<EditChain, TopologyEditError> {
        let previous = self.require_chain(id)?;
        for residue in self.chain_residues(id) {
            self.delete_residue(residue)?;
        }
        let source = SourceHierarchy {
            source: self.source.as_deref(),
            ..self.source_hierarchy_blocks()
        };
        if self
            .hierarchy
            .chains
            .remove(&id, |index| source.chain(index))
            .is_some()
        {
            self.properties.chains_mut().clear_index(previous.slot);
            self.changed();
        }
        Ok(previous)
    }
    /// The live residues of a live chain, in handle order.
    fn chain_residues(&self, chain: EditChainId) -> Vec<EditResidueId> {
        // A residue keeps the chain it was added to.
        let residues = &self.hierarchy.residues;
        let mut listed = match (
            self.source.as_deref(),
            self.hierarchy.chains.block().source(chain),
        ) {
            (Some(source), Some(index)) => source
                .hierarchy()
                .chain(ChainId::new(index))
                .expect("source chain")
                .residues()
                .iter()
                .map(|residue| residues.block().handle::<EditResidueId>(residue.index()))
                .filter(|residue| residues.contains_key(residue))
                .collect(),
            _ => Vec::new(),
        };
        listed.sort_unstable();
        listed.extend(
            residues
                .added()
                .iter()
                .filter(|(_, residue)| residue.chain == chain)
                .map(|(&id, _)| id),
        );
        listed
    }

    pub(super) fn invalidate_residue_for_atom(&mut self, atom: EditAtomId) {
        if let Some(residue) = self
            .site_of_atom(atom)
            .and_then(|id| self.atom_site(id))
            .map(|site| site.residue)
        {
            self.invalidate_residue_class(residue, true);
        }
    }

    fn invalidate_residue_class(&mut self, id: EditResidueId, composition_changed: bool) {
        let current = self.require_residue(id).expect("live residue");
        let had_class = current.class.is_some();
        if (composition_changed || !current.class_explicit)
            && (current.class.is_some() || current.class_explicit)
        {
            let residue = self.residue_mut(id);
            residue.class = None;
            residue.class_explicit = false;
        }
        // Previously attached atoms were invalidated when this residue first
        // lost its cached class. Newly attached/detached atoms are handled by
        // the site operation, so assembly does not repeatedly scan a growing residue.
        if had_class {
            self.invalidate_molecule_classes_for_residue(id);
        }
    }

    fn invalidate_molecule_classes_for_residue(&mut self, residue: EditResidueId) {
        let atoms = self.residue_atom_list(residue);
        for atom in atoms {
            self.invalidate_molecule_class_for_atom(atom);
        }
    }

    fn invalidate_molecule_class_for_atom(&mut self, atom: EditAtomId) {
        if let Ok(loc) = self.atom_location(atom) {
            let (class, explicit) = self.group_class(loc.group);
            if !explicit && class.is_some() {
                self.group_mut(loc.group).class = None;
            }
        }
    }
    /// The atoms of one live residue, copied from the source on first use so
    /// that a membership change starts from the current atoms.
    fn residue_atoms_mut(&mut self, residue: EditResidueId) -> &mut BTreeSet<EditAtomId> {
        if !self.hierarchy.residue_atoms.contains_key(&residue) {
            let atoms = self.source_residue_atoms(residue).collect();
            self.hierarchy.residue_atoms.insert(residue, atoms);
        }
        self.hierarchy
            .residue_atoms
            .get_mut(&residue)
            .expect("materialized residue atoms")
    }
    /// The atoms of one live residue in handle order.
    pub(super) fn residue_atom_list(&self, residue: EditResidueId) -> Vec<EditAtomId> {
        match self.hierarchy.residue_atoms.get(&residue) {
            Some(atoms) => atoms.iter().copied().collect(),
            None => {
                let mut atoms = self.source_residue_atoms(residue).collect::<Vec<_>>();
                atoms.sort_unstable();
                atoms
            }
        }
    }
    /// The source atoms of a source residue whose membership is unchanged.
    fn source_residue_atoms(
        &self,
        residue: EditResidueId,
    ) -> impl Iterator<Item = EditAtomId> + '_ {
        let source = self
            .source
            .as_deref()
            .zip(self.hierarchy.residues.block().source(residue));
        let atoms = self.atoms.block();
        source.into_iter().flat_map(move |(source, id)| {
            source
                .hierarchy()
                .residue(ResidueId::new(id))
                .expect("source residue")
                .atom_sites()
                .iter()
                .map(move |&site| atoms.handle(source.site_atoms()[site.index()].index()))
        })
    }
    pub(super) fn remove_site_for_atom(&mut self, atom: EditAtomId) {
        if let Some(site) = self.site_of_atom(atom) {
            self.delete_atom_site(site).expect("live hierarchy site");
        }
    }
    fn prune_empty_residue(&mut self, residue: EditResidueId) {
        if !self.residue_atom_list(residue).is_empty() {
            return;
        }
        self.hierarchy.residue_atoms.remove(&residue);
        let source = SourceHierarchy {
            source: self.source.as_deref(),
            ..self.source_hierarchy_blocks()
        };
        let removed = self
            .hierarchy
            .residues
            .remove(&residue, |index| source.residue(index))
            .unwrap();
        self.properties.residues_mut().clear_index(removed.slot);
        if self.chain_residues(removed.chain).is_empty() {
            let chain = self
                .hierarchy
                .chains
                .remove(&removed.chain, |index| source.chain(index))
                .unwrap();
            self.properties.chains_mut().clear_index(chain.slot);
        }
    }
    pub(super) fn import_hierarchy(&mut self, source: &Topology) {
        // A fresh draft takes one handle block per node kind, in source ID
        // order, and reads the source nodes until an edit changes them.
        let hierarchy = source.hierarchy();
        self.hierarchy = EditHierarchy {
            chains: Overlay::new(HandleBlock::reserve(hierarchy.chains().len())),
            residues: Overlay::new(HandleBlock::reserve(hierarchy.residues().len())),
            sites: Overlay::new(HandleBlock::reserve(hierarchy.atom_sites().len())),
            ..EditHierarchy::default()
        };
    }
}
