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
    pub(super) name: String,
    pub(super) label_comp: Option<String>,
    pub(super) author_comp: Option<String>,
    pub(super) label_seq: Option<i32>,
    pub(super) author_seq: Option<String>,
    pub(super) insertion: Option<String>,
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
        &self.name
    }
    /// The label component identifier, or `None`.
    pub fn label_comp_id(&self) -> Option<&str> {
        self.label_comp.as_deref()
    }
    /// The author component identifier, or `None`.
    pub fn author_comp_id(&self) -> Option<&str> {
        self.author_comp.as_deref()
    }
    /// The label sequence number, or `None`.
    pub fn label_seq_id(&self) -> Option<i32> {
        self.label_seq
    }
    /// The author sequence identifier, or `None`.
    pub fn author_seq_id(&self) -> Option<&str> {
        self.author_seq.as_deref()
    }
    /// The insertion code, or `None`.
    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion.as_deref()
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
    pub(super) metadata: AtomSiteMetadata,
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

#[derive(Debug, Clone, Default)]
pub(super) struct EditHierarchy {
    pub(super) chains: BTreeMap<EditChainId, EditChain>,
    pub(super) residues: BTreeMap<EditResidueId, EditResidue>,
    pub(super) sites: BTreeMap<EditAtomSiteId, EditAtomSite>,
    pub(super) atom_sites: BTreeMap<EditAtomId, EditAtomSiteId>,
    pub(super) residue_atoms: BTreeMap<EditResidueId, BTreeSet<EditAtomId>>,
    pub(super) source_chains: BTreeMap<ChainId, EditChainId>,
    pub(super) source_residues: BTreeMap<ResidueId, EditResidueId>,
    pub(super) source_sites: BTreeMap<AtomSiteId, EditAtomSiteId>,
}

impl TopologyEditor {
    /// Live draft chains in handle order.
    pub fn chains(&self) -> impl ExactSizeIterator<Item = (EditChainId, &EditChain)> {
        self.hierarchy.chains.iter().map(|(&id, c)| (id, c))
    }
    /// Live draft residues in handle order.
    pub fn residues(&self) -> impl ExactSizeIterator<Item = (EditResidueId, &EditResidue)> {
        self.hierarchy.residues.iter().map(|(&id, r)| (id, r))
    }
    /// Live draft atom sites in handle order.
    pub fn atom_sites(&self) -> impl ExactSizeIterator<Item = (EditAtomSiteId, &EditAtomSite)> {
        self.hierarchy.sites.iter().map(|(&id, s)| (id, s))
    }
    /// One draft chain, or `None` when `id` names no live chain of this draft.
    pub fn chain(&self, id: EditChainId) -> Option<&EditChain> {
        self.hierarchy.chains.get(&id)
    }
    /// One draft residue, or `None` when `id` names no live residue of this draft.
    pub fn residue(&self, id: EditResidueId) -> Option<&EditResidue> {
        self.hierarchy.residues.get(&id)
    }
    /// One draft atom site, or `None` when `id` names no live site of this draft.
    pub fn atom_site(&self, id: EditAtomSiteId) -> Option<&EditAtomSite> {
        self.hierarchy.sites.get(&id)
    }
    /// The draft handle of a source chain, or `None` when the chain is not in
    /// the source or was deleted.
    pub fn chain_handle(&self, source: ChainId) -> Option<EditChainId> {
        self.hierarchy
            .source_chains
            .get(&source)
            .copied()
            .filter(|id| self.hierarchy.chains.contains_key(id))
    }
    /// The draft handle of a source residue, or `None` when the residue is not
    /// in the source or was deleted.
    pub fn residue_handle(&self, source: ResidueId) -> Option<EditResidueId> {
        self.hierarchy
            .source_residues
            .get(&source)
            .copied()
            .filter(|id| self.hierarchy.residues.contains_key(id))
    }
    pub(crate) fn require_chain(&self, id: EditChainId) -> Result<&EditChain, TopologyEditError> {
        self.chain(id).ok_or(TopologyEditError::InvalidChain(id))
    }
    pub(crate) fn require_residue(
        &self,
        id: EditResidueId,
    ) -> Result<&EditResidue, TopologyEditError> {
        self.residue(id)
            .ok_or(TopologyEditError::InvalidResidue(id))
    }
    pub(crate) fn require_atom_site(
        &self,
        id: EditAtomSiteId,
    ) -> Result<&EditAtomSite, TopologyEditError> {
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
        self.hierarchy
            .source_sites
            .get(&source)
            .copied()
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
        Ok(self.hierarchy.atom_sites.get(&atom).copied())
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
        let name = name.into();
        let slot = self.properties.residues().len();
        self.hierarchy.residues.insert(
            id,
            EditResidue {
                chain,
                label_comp: Some(name.clone()),
                name,
                author_comp: None,
                label_seq,
                author_seq,
                insertion,
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
        if self.hierarchy.atom_sites.contains_key(&atom) {
            return Err(TopologyEditError::DuplicateAtomPlacement(atom));
        }
        let id = EditAtomSiteId::new();
        let slot = self.properties.atom_sites().len();
        self.hierarchy.sites.insert(
            id,
            EditAtomSite {
                residue,
                atom,
                metadata,
                slot,
            },
        );
        self.hierarchy.atom_sites.insert(atom, id);
        self.hierarchy
            .residue_atoms
            .get_mut(&residue)
            .unwrap()
            .insert(atom);
        self.invalidate_molecule_class_for_atom(atom);
        self.invalidate_residue_class(residue, true);
        self.properties.atom_sites_mut().resize_missing(slot + 1);
        self.changed();
        Ok(id)
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
        self.require_chain(id)?;
        let label = label.into();
        let chain = self.hierarchy.chains.get_mut(&id).unwrap();
        if chain.label == label && chain.author == author {
            return Ok(());
        }
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
        self.require_residue(id)?;
        let residue = self.hierarchy.residues.get_mut(&id).unwrap();
        if residue.label_comp == label && residue.author_comp == author {
            return Ok(());
        }
        residue.label_comp = label;
        residue.author_comp = author;
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
        let previous = self.require_residue(id)?.class;
        if self.require_residue(id)?.class == Some(class)
            && self.require_residue(id)?.class_explicit
        {
            return Ok(());
        }
        self.hierarchy.residues.get_mut(&id).unwrap().class = Some(class);
        self.hierarchy.residues.get_mut(&id).unwrap().class_explicit = true;
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
        let residue = self.hierarchy.residues.get_mut(&id).unwrap();
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
        if self.require_atom_site(id)?.metadata == metadata {
            return Ok(());
        }
        let previous = &self.require_atom_site(id)?.metadata;
        let names_changed = previous.label_atom_id != metadata.label_atom_id
            || previous.auth_atom_id != metadata.auth_atom_id;
        self.hierarchy.sites.get_mut(&id).unwrap().metadata = metadata;
        if names_changed {
            self.invalidate_molecule_class_for_atom(self.hierarchy.sites[&id].atom);
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
        let previous = self.require_atom_site(id)?.residue;
        if previous == residue {
            return Ok(());
        }
        self.hierarchy.sites.get_mut(&id).unwrap().residue = residue;
        let atom = self.hierarchy.sites[&id].atom;
        self.hierarchy
            .residue_atoms
            .get_mut(&previous)
            .unwrap()
            .remove(&atom);
        self.hierarchy
            .residue_atoms
            .get_mut(&residue)
            .unwrap()
            .insert(atom);
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
        let site = self.hierarchy.sites.remove(&id).unwrap();
        self.hierarchy.atom_sites.remove(&site.atom);
        self.hierarchy
            .residue_atoms
            .get_mut(&site.residue)
            .unwrap()
            .remove(&site.atom);
        self.properties.atom_sites_mut().clear_index(site.slot);
        self.invalidate_residue_class(site.residue, true);
        self.invalidate_molecule_class_for_atom(site.atom);
        self.prune_empty_residue(site.residue);
        self.changed();
        Ok(site)
    }
    /// Removes a residue and its sites, retaining all chemical atoms.
    pub fn delete_residue(&mut self, id: EditResidueId) -> Result<EditResidue, TopologyEditError> {
        let previous = self.require_residue(id)?.clone();
        let sites = self
            .atom_sites()
            .filter(|(_, s)| s.residue == id)
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
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
        let previous = self.require_chain(id)?.clone();
        let residues = self
            .residues()
            .filter(|(_, r)| r.chain == id)
            .map(|(id, _)| id)
            .collect::<Vec<_>>();
        for residue in residues {
            self.delete_residue(residue)?;
        }
        if self.hierarchy.chains.remove(&id).is_some() {
            self.properties.chains_mut().clear_index(previous.slot);
            self.changed();
        }
        Ok(previous)
    }

    pub(super) fn invalidate_residue_for_atom(&mut self, atom: EditAtomId) {
        if let Some(residue) = self
            .hierarchy
            .atom_sites
            .get(&atom)
            .and_then(|id| self.hierarchy.sites.get(id))
            .map(|site| site.residue)
        {
            self.invalidate_residue_class(residue, true);
        }
    }

    fn invalidate_residue_class(&mut self, id: EditResidueId, composition_changed: bool) {
        let residue = self.hierarchy.residues.get_mut(&id).unwrap();
        let had_class = residue.class.is_some();
        if composition_changed || !residue.class_explicit {
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
        let atoms = self.hierarchy.residue_atoms[&residue]
            .iter()
            .copied()
            .collect::<Vec<_>>();
        for atom in atoms {
            self.invalidate_molecule_class_for_atom(atom);
        }
    }

    fn invalidate_molecule_class_for_atom(&mut self, atom: EditAtomId) {
        if let Some(loc) = self.atoms.get(&atom) {
            let group = self.groups[loc.group].as_mut().unwrap();
            if !group.class_explicit {
                group.class = None;
            }
        }
    }
    pub(super) fn remove_site_for_atom(&mut self, atom: EditAtomId) {
        if let Some(site) = self.hierarchy.atom_sites.get(&atom).copied() {
            self.delete_atom_site(site).expect("live hierarchy site");
        }
    }
    fn prune_empty_residue(&mut self, residue: EditResidueId) {
        if !self.hierarchy.residue_atoms[&residue].is_empty() {
            return;
        }
        self.hierarchy.residue_atoms.remove(&residue);
        let removed = self.hierarchy.residues.remove(&residue).unwrap();
        self.properties.residues_mut().clear_index(removed.slot);
        if !self
            .hierarchy
            .residues
            .values()
            .any(|r| r.chain == removed.chain)
        {
            let chain = self.hierarchy.chains.remove(&removed.chain).unwrap();
            self.properties.chains_mut().clear_index(chain.slot);
        }
    }
    pub(super) fn import_hierarchy(&mut self, source: &Topology) {
        for (source_id, chain) in source.hierarchy().chains() {
            let id = EditChainId::new();
            self.hierarchy.chains.insert(
                id,
                EditChain {
                    label: chain.label_id().to_owned(),
                    author: chain.author_id().map(str::to_owned),
                    slot: source_id.index(),
                },
            );
            self.hierarchy.source_chains.insert(source_id, id);
        }
        for (source_id, residue) in source.hierarchy().residues() {
            let id = EditResidueId::new();
            self.hierarchy.residues.insert(
                id,
                EditResidue {
                    chain: self.hierarchy.source_chains[&residue.chain()],
                    name: residue.name().to_owned(),
                    label_comp: residue.label_comp_id().map(str::to_owned),
                    author_comp: residue.author_comp_id().map(str::to_owned),
                    label_seq: residue.label_seq_id(),
                    author_seq: residue.author_seq_id().map(str::to_owned),
                    insertion: residue.insertion_code().map(str::to_owned),
                    class: Some(residue.class()),
                    class_explicit: source
                        .hierarchy()
                        .residue(source_id)
                        .is_some_and(|residue| residue.explicit_class.is_some()),
                    slot: source_id.index(),
                },
            );
            self.hierarchy.source_residues.insert(source_id, id);
            self.hierarchy.residue_atoms.insert(id, BTreeSet::new());
        }
        for (source_id, site) in source.hierarchy().atom_sites() {
            let id = EditAtomSiteId::new();
            let atom = self.source_atoms[&site.atom()];
            self.hierarchy.sites.insert(
                id,
                EditAtomSite {
                    residue: self.hierarchy.source_residues[&site.residue()],
                    atom,
                    metadata: site.metadata().clone(),
                    slot: source_id.index(),
                },
            );
            self.hierarchy.source_sites.insert(source_id, id);
            self.hierarchy.atom_sites.insert(atom, id);
            self.hierarchy
                .residue_atoms
                .get_mut(&self.hierarchy.source_residues[&site.residue()])
                .unwrap()
                .insert(atom);
        }
    }
}
