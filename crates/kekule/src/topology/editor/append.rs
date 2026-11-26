//! Import one complete topology into an already staged structure edit.
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
    /// Whether [`Self::append_topology`] takes every annotation of `source`:
    /// each table of `source` has a row per entity it copies, and no column of
    /// the draft clashes with one of the same key in `source`.
    pub(crate) fn accepts_properties_of(&self, source: &Topology) -> bool {
        let rows = [
            source.instance_count(),
            source.atom_ids().len(),
            source.bond_ids().len(),
            source.chains().len(),
            source.residues().len(),
            source.atom_sites().len(),
        ];
        source.properties().validate_dimensions(rows).is_ok()
            && self.properties.accepts_rows_of(source.properties())
    }

    /// The coordinating structure operation owns rollback. Source identity is
    /// local to this import, separate from the editor's original source
    /// correspondence. `chain_ids` gives each source chain's label and author
    /// identifiers in the draft from its own.
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
        // An instance shares an identical definition the draft holds, from
        // its source or an earlier append. Otherwise its definition is copied
        // once per import, so instances of one source definition share the
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
        for chain in source.chains() {
            let id = chain.id();
            let (label, author) = chain_ids(chain.label_id(), chain.author_id());
            if label != chain.label_id() || author.as_deref() != chain.author_id() {
                renamed.insert(
                    id,
                    (
                        (chain.label_id(), chain.author_id()),
                        (label.clone(), author.clone()),
                    ),
                );
            }
            let handle = self.add_chain(label, author)?;
            mapping.chains.insert(id, handle);
        }
        // Residues and sites share the source's label and metadata rows.
        let hierarchy = source.hierarchy();
        for residue in source.residues() {
            let handle = self.insert_residue(
                mapping.chains[&residue.chain().id()],
                hierarchy.shared_residue_labels(residue.id()),
            )?;
            mapping.residues.insert(residue.id(), handle);
        }
        let mut retargeted = Vec::new();
        for site in source.atom_sites() {
            let handle = self.insert_atom_site(
                mapping.residues[&site.residue().id()],
                mapping.atoms[&site.atom().id()],
                hierarchy.shared_site_metadata(site.id()),
            )?;
            if let Some((old, (label, author))) = renamed.get(&site.residue().chain().id()) {
                if let Some(metadata) =
                    retargeted_site_ids(site.metadata(), *old, (label, author.as_deref()))
                {
                    retargeted.push((handle, metadata));
                }
            }
            mapping.sites.insert(site.id(), handle);
        }
        self.set_site_rows(retargeted);
        // Insertion invalidates classification while the residue is incomplete.
        for residue in source.residues() {
            let assignment = residue.class_assignment();
            let target = self.residue_mut(mapping.residues[&residue.id()]);
            target.class = Some(assignment.class());
            target.class_explicit = assignment.is_explicit();
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
            .chains()
            .map(|chain| {
                self.require_chain(mapping.chains[&chain.id()])
                    .map(|chain| chain.slot)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let residue_rows = source
            .residues()
            .map(|residue| {
                self.require_residue(mapping.residues[&residue.id()])
                    .map(|residue| residue.slot)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let site_rows = source
            .atom_sites()
            .map(|site| {
                self.require_atom_site(mapping.sites[&site.id()])
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
