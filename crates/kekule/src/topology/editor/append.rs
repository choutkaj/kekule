//! Import one complete topology into an already staged model edit.
use super::*;

/// The label and author identifiers an appended chain takes from its own.
pub(crate) type ChainIds<'a> = dyn Fn(&str, Option<&str>) -> (String, Option<String>) + 'a;

#[derive(Debug, Clone)]
pub(crate) struct AppendMapping {
    pub source: Arc<Topology>,
    pub atoms: BTreeMap<InstanceAtomId, EditAtomId>,
    pub bonds: BTreeMap<InstanceBondId, EditBondId>,
    pub chains: BTreeMap<ChainId, EditChainId>,
    pub residues: BTreeMap<ResidueId, EditResidueId>,
    pub sites: BTreeMap<AtomSiteId, EditAtomSiteId>,
}

impl TopologyEditor {
    /// The coordinating model operation owns rollback. Source identity is local
    /// to this import, separate from the editor's original source correspondence.
    /// `chain_ids` gives each source chain's label and author identifiers in
    /// the draft from its own.
    pub(crate) fn append_topology(
        &mut self,
        source: Arc<Topology>,
        chain_ids: &ChainIds<'_>,
    ) -> Result<AppendMapping, TopologyEditError> {
        let mut mapping = AppendMapping {
            source: Arc::clone(&source),
            atoms: BTreeMap::new(),
            bonds: BTreeMap::new(),
            chains: BTreeMap::new(),
            residues: BTreeMap::new(),
            sites: BTreeMap::new(),
        };
        // An occurrence shares an identical definition the draft holds, from
        // its source or an earlier append. Otherwise its definition is copied
        // once per import, so occurrences of one source definition share the
        // copy. The copy keeps the definition's perception, and so its profile.
        let mut copies = BTreeMap::<MoleculeDefinitionId, Arc<Molecule>>::new();
        let instance_start = self.properties.molecule_instances().len();
        self.properties
            .molecule_instances_mut()
            .resize_missing(instance_start + source.instance_count());
        // Imported atoms keep the source dense order after existing rows.
        let bases = (self.properties.atoms().len(), self.properties.bonds().len());
        let handles = GroupIdentity::source_handles(&source);
        for (instance, value) in source.instance_records() {
            let definition = source
                .definition(value.definition())
                .expect("instances name source definitions");
            let molecule = definition.molecule();
            let labels = sharing::instance_labels(&source, instance);
            let appended = sharing::Appended {
                molecule,
                default_perception: definition.default_perception,
                labels: &labels,
                class: definition.class(),
                class_explicit: source
                    .molecule_class_overrides
                    .contains_key(&definition.id()),
            };
            let chemistry = match self.shared.find(self.source.as_deref(), &appended) {
                Some(shared) => shared,
                None => {
                    let copy = Arc::clone(
                        copies
                            .entry(definition.id())
                            .or_insert_with(|| Arc::new(molecule.clone())),
                    );
                    self.shared.insert(Arc::clone(&copy), &appended);
                    GroupChemistry::Added(copy, definition.default_perception)
                }
            };
            let identity = GroupIdentity::dense(&source, instance, handles, bases);
            let group = self.register_group(
                chemistry,
                molecule,
                Some(appended.class),
                appended.class_explicit,
                identity,
            );
            let added = self.group_handles(group);
            self.group_mut(group).instance_slot = Some(instance_start + instance.index());
            self.properties
                .resize_atoms(self.properties.atoms().len() + molecule.atom_count());
            self.properties
                .resize_bonds(self.properties.bonds().len() + molecule.bond_count());
            mapping.atoms.extend(
                added
                    .atoms
                    .into_iter()
                    .map(|(atom, handle)| (InstanceAtomId::new(instance, atom), handle)),
            );
            mapping.bonds.extend(
                added
                    .bonds
                    .into_iter()
                    .map(|(bond, handle)| (InstanceBondId::new(instance, bond), handle)),
            );
        }
        // Hierarchy is imported independently of molecular connectedness. Even
        // equal source labels remain distinct chains, residues and atom sites.
        // A renamed chain's sites take its new identifiers where they named
        // the old ones, so writers name the chain as it is now.
        let mut renamed = std::collections::HashMap::new();
        for (id, chain) in source.hierarchy().chains() {
            let (label, author) = chain_ids(chain.label_id(), chain.author_id());
            if label != chain.label_id() || author.as_deref() != chain.author_id() {
                renamed.insert(
                    id,
                    (
                        (chain.label_id().to_owned(), label.clone()),
                        (chain.author_id().map(str::to_owned), author.clone()),
                    ),
                );
            }
            let handle = self.add_chain(label, author)?;
            mapping.chains.insert(id, handle);
        }
        for (id, residue) in source.hierarchy().residues() {
            let handle = self.add_residue(
                mapping.chains[&residue.chain()],
                residue.name(),
                residue.label_seq_id(),
                residue.author_seq_id().map(str::to_owned),
                residue.insertion_code().map(str::to_owned),
            )?;
            self.set_residue_component_ids(
                handle,
                residue.label_comp_id().map(str::to_owned),
                residue.author_comp_id().map(str::to_owned),
            )?;
            mapping.residues.insert(id, handle);
        }
        for (id, site) in source.hierarchy().atom_sites() {
            let mut metadata = site.metadata().clone();
            let chain = source
                .hierarchy()
                .residue(site.residue())
                .map(|r| r.chain());
            if let Some(((old_label, new_label), (old_author, new_author))) =
                chain.and_then(|chain| renamed.get(&chain))
            {
                if metadata.label_asym_id.as_ref() == Some(old_label) {
                    metadata.label_asym_id = Some(new_label.clone());
                }
                if metadata.auth_asym_id.is_some() && metadata.auth_asym_id == *old_author {
                    metadata.auth_asym_id = new_author.clone();
                }
            }
            let handle = self.add_atom_site(
                mapping.residues[&site.residue()],
                mapping.atoms[&site.atom()],
                metadata,
            )?;
            mapping.sites.insert(id, handle);
        }
        // Insertion invalidates classification while the residue is incomplete.
        for (id, residue) in source.hierarchy().residues() {
            let target = self.residue_mut(mapping.residues[&id]);
            target.class = Some(residue.class());
            target.class_explicit = source
                .hierarchy()
                .residue(id)
                .is_some_and(|residue| residue.explicit_class().is_some());
        }
        for (instance, value) in source.instance_records() {
            let definition = source.definition(value.definition()).unwrap();
            let local = definition.molecule().atom_ids().next().unwrap();
            let handle = mapping.atoms[&InstanceAtomId::new(instance, local)];
            let group = self.atom_location(handle)?.group;
            let group = self.group_mut(group);
            group.class = Some(definition.class());
            group.class_explicit = source
                .molecule_class_overrides
                .contains_key(&value.definition());
        }

        let atom_rows = source
            .atom_ids()
            .iter()
            .map(|id| self.atom_slot(mapping.atoms[id]))
            .collect::<Result<Vec<_>, _>>()?;
        let bond_rows = source
            .bond_ids()
            .iter()
            .map(|id| self.bond_slot(mapping.bonds[id]))
            .collect::<Result<Vec<_>, _>>()?;
        let chain_rows = source
            .hierarchy()
            .chains()
            .map(|(id, _)| {
                self.require_chain(mapping.chains[&id])
                    .map(|chain| chain.slot)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let residue_rows = source
            .hierarchy()
            .residues()
            .map(|(id, _)| {
                self.require_residue(mapping.residues[&id])
                    .map(|residue| residue.slot)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let site_rows = source
            .hierarchy()
            .atom_sites()
            .map(|(id, _)| {
                self.require_atom_site(mapping.sites[&id])
                    .map(|site| site.slot)
            })
            .collect::<Result<Vec<_>, _>>()?;
        self.properties
            .molecule_instances_mut()
            .copy_rows_from(
                source.properties().molecule_instances().raw(),
                &(instance_start..instance_start + source.instance_count()).collect::<Vec<_>>(),
            )
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::MoleculeInstances,
                error: Box::new(error),
            })?;
        self.properties
            .atoms_mut()
            .copy_rows_from(source.properties().atoms().raw(), &atom_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::TopologyAtoms,
                error: Box::new(error),
            })?;
        self.properties
            .bonds_mut()
            .copy_rows_from(source.properties().bonds().raw(), &bond_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::TopologyBonds,
                error: Box::new(error),
            })?;
        self.properties
            .chains_mut()
            .copy_rows_from(source.properties().chains().raw(), &chain_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::Chains,
                error: Box::new(error),
            })?;
        self.properties
            .residues_mut()
            .copy_rows_from(source.properties().residues().raw(), &residue_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::Residues,
                error: Box::new(error),
            })?;
        self.properties
            .atom_sites_mut()
            .copy_rows_from(source.properties().atom_sites().raw(), &site_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: crate::properties::PropertyDomain::AtomSites,
                error: Box::new(error),
            })?;
        self.changed();
        Ok(mapping)
    }
}
