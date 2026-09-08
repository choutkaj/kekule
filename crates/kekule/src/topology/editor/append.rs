//! Import one complete topology into an already staged model edit.
use super::*;

#[derive(Debug, Clone)]
pub(crate) struct AppendMapping {
    pub source: Arc<Topology>,
    pub token: Arc<()>,
    pub atoms: BTreeMap<InstanceAtomId, EditAtomId>,
    pub bonds: BTreeMap<InstanceBondId, EditBondId>,
    pub chains: BTreeMap<ChainId, EditChainId>,
    pub residues: BTreeMap<ResidueId, EditResidueId>,
    pub sites: BTreeMap<AtomSiteId, EditAtomSiteId>,
}

impl TopologyEditor {
    /// The coordinating model operation owns rollback. Source identity is local
    /// to this import, separate from the editor's original source correspondence.
    pub(crate) fn append_topology(
        &mut self,
        source: Arc<Topology>,
    ) -> Result<AppendMapping, TopologyEditError> {
        let mut mapping = AppendMapping {
            source: Arc::clone(&source),
            token: Arc::new(()),
            atoms: BTreeMap::new(),
            bonds: BTreeMap::new(),
            chains: BTreeMap::new(),
            residues: BTreeMap::new(),
            sites: BTreeMap::new(),
        };
        let definitions = source
            .definitions()
            .map(|(id, definition)| (id, Arc::new(definition.molecule().clone())))
            .collect::<BTreeMap<_, _>>();
        let instance_start = self.properties.molecule_instances().len();
        self.properties
            .molecule_instances_mut()
            .resize_missing(instance_start + source.instance_count());
        for (instance, value) in source.instances() {
            let molecule = &definitions[&value.definition()];
            let group_index = self.groups.len();
            let added = self.register_group(
                GroupChemistry::Added(Arc::clone(molecule)),
                molecule,
                None,
                Some(source.definition(value.definition()).unwrap().class()),
            );
            self.groups[group_index].as_mut().unwrap().instance_slot =
                Some(instance_start + instance.index());
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
        for (id, chain) in source.hierarchy().chains() {
            let handle = self.add_chain(chain.label_id(), chain.author_id().map(str::to_owned))?;
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
            let handle = self.add_atom_site(
                mapping.residues[&site.residue()],
                mapping.atoms[&site.atom()],
                site.metadata().clone(),
            )?;
            mapping.sites.insert(id, handle);
        }
        // Insertion invalidates classification while the residue is incomplete.
        for (id, residue) in source.hierarchy().residues() {
            self.set_residue_class(mapping.residues[&id], residue.class())?;
        }

        let atom_rows = source
            .atom_ids()
            .iter()
            .map(|id| self.atoms[&mapping.atoms[id]].slot)
            .collect::<Vec<_>>();
        let bond_rows = source
            .bond_ids()
            .iter()
            .map(|id| self.bonds[&mapping.bonds[id]].slot)
            .collect::<Vec<_>>();
        let chain_rows = source
            .hierarchy()
            .chains()
            .map(|(id, _)| self.hierarchy.chains[&mapping.chains[&id]].slot)
            .collect::<Vec<_>>();
        let residue_rows = source
            .hierarchy()
            .residues()
            .map(|(id, _)| self.hierarchy.residues[&mapping.residues[&id]].slot)
            .collect::<Vec<_>>();
        let site_rows = source
            .hierarchy()
            .atom_sites()
            .map(|(id, _)| self.hierarchy.sites[&mapping.sites[&id]].slot)
            .collect::<Vec<_>>();
        self.properties
            .molecule_instances_mut()
            .copy_rows_from(
                source.molecule_instance_properties(),
                &(instance_start..instance_start + source.instance_count()).collect::<Vec<_>>(),
            )
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: "molecule instance",
                error: Box::new(error),
            })?;
        self.properties
            .atoms_mut()
            .copy_rows_from(source.atom_properties(), &atom_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: "topology atom",
                error: Box::new(error),
            })?;
        self.properties
            .bonds_mut()
            .copy_rows_from(source.bond_properties(), &bond_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: "topology bond",
                error: Box::new(error),
            })?;
        self.properties
            .chains_mut()
            .copy_rows_from(source.chain_properties(), &chain_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: "chain",
                error: Box::new(error),
            })?;
        self.properties
            .residues_mut()
            .copy_rows_from(source.residue_properties(), &residue_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: "residue",
                error: Box::new(error),
            })?;
        self.properties
            .atom_sites_mut()
            .copy_rows_from(source.atom_site_properties(), &site_rows)
            .map_err(|error| TopologyEditError::AppendProperty {
                domain: "atom site",
                error: Box::new(error),
            })?;
        self.append_tokens.push(Arc::clone(&mapping.token));
        self.changed();
        Ok(mapping)
    }
}
