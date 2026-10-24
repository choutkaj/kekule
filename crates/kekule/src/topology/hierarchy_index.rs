//! Hierarchy lookups derived once per published topology.

use super::dense_order::DenseOrder;
use super::{AtomSiteId, ChainId, Hierarchy, InstanceAtomId, ResidueId, TopologyAtomIndex};

const NO_SITE: u32 = u32::MAX;

/// Lookups between a published hierarchy, dense atoms, and instances, so
/// walks cost what they return.
///
/// Derived at publication from the hierarchy and dense atom order, and shared
/// by every snapshot of that layout. Residue sites are in residue order;
/// instance lists are in ascending ID order. Compact copies of the residue
/// lists keep walks off the larger hierarchy nodes.
#[derive(Debug, Clone)]
pub(super) struct HierarchyIndex {
    /// The dense atom of each site, by site ID.
    site_atoms: Vec<TopologyAtomIndex>,
    /// The residue of each site, by site ID.
    site_residues: Vec<ResidueId>,
    /// The site of each dense atom, or [`NO_SITE`].
    atom_sites: Vec<u32>,
    residue_sites: Groups<AtomSiteId>,
    instance_sites: Groups<AtomSiteId>,
    instance_residues: Groups<ResidueId>,
    instance_chains: Groups<ChainId>,
}

impl HierarchyIndex {
    /// Indexes `hierarchy`, whose sites all name atoms of `atoms`.
    pub(super) fn new(
        hierarchy: &Hierarchy,
        atoms: &DenseOrder<InstanceAtomId>,
        instance_count: usize,
    ) -> Self {
        let mut site_atoms = Vec::with_capacity(hierarchy.atom_sites().len());
        let mut site_residues = Vec::with_capacity(hierarchy.atom_sites().len());
        let mut site_instances = Vec::with_capacity(hierarchy.atom_sites().len());
        let mut atom_sites = vec![NO_SITE; atoms.len()];
        for (id, site) in hierarchy.atom_sites() {
            let index = atoms
                .index(site.atom())
                .expect("published atom sites name topology atoms");
            site_atoms.push(TopologyAtomIndex::new(index as u32));
            site_residues.push(site.residue());
            site_instances.push(site.atom().instance().index());
            atom_sites[index] = id.raw();
        }
        let residue_sites = Groups::from_lists(
            hierarchy
                .residues()
                .map(|(_, residue)| residue.atom_sites()),
        );
        // A residue or chain reaches the instances of its sites.
        let site_instances = &site_instances;
        let residue_instances = |residue: usize| {
            residue_sites
                .get(residue)
                .iter()
                .map(move |site| site_instances[site.index()])
        };
        Self {
            instance_sites: Groups::new(
                instance_count,
                site_instances
                    .iter()
                    .enumerate()
                    .map(|(site, &instance)| (instance, AtomSiteId::new(site as u32))),
            ),
            instance_residues: Groups::new(
                instance_count,
                (0..hierarchy.residues().len()).flat_map(|residue| {
                    residue_instances(residue)
                        .map(move |instance| (instance, ResidueId::new(residue as u32)))
                }),
            ),
            instance_chains: Groups::new(
                instance_count,
                hierarchy.chains().flat_map(|(id, chain)| {
                    chain
                        .residues()
                        .iter()
                        .flat_map(|residue| residue_instances(residue.index()))
                        .map(move |instance| (instance, id))
                }),
            ),
            site_atoms,
            site_residues,
            atom_sites,
            residue_sites,
        }
    }

    /// The dense atom of a valid site.
    pub(super) fn site_atom(&self, site: AtomSiteId) -> TopologyAtomIndex {
        self.site_atoms[site.index()]
    }

    /// The residue of a valid site.
    pub(super) fn site_residue(&self, site: AtomSiteId) -> ResidueId {
        self.site_residues[site.index()]
    }

    /// The sites of a valid residue, in residue order.
    pub(super) fn residue_sites(&self, residue: ResidueId) -> &[AtomSiteId] {
        self.residue_sites.get(residue.index())
    }

    /// The site of a valid dense atom, or `None` when it has none.
    pub(super) fn atom_site(&self, atom: TopologyAtomIndex) -> Option<AtomSiteId> {
        let site = self.atom_sites[atom.index()];
        (site != NO_SITE).then(|| AtomSiteId::new(site))
    }

    /// The sites of a valid instance's atoms.
    pub(super) fn instance_sites(&self, instance: usize) -> &[AtomSiteId] {
        self.instance_sites.get(instance)
    }

    /// The residues with a site of a valid instance's atoms.
    pub(super) fn instance_residues(&self, instance: usize) -> &[ResidueId] {
        self.instance_residues.get(instance)
    }

    /// The chains with a site of a valid instance's atoms.
    pub(super) fn instance_chains(&self, instance: usize) -> &[ChainId] {
        self.instance_chains.get(instance)
    }
}

/// Items of numbered groups, stored contiguously: group `g` owns
/// `items[starts[g]..starts[g + 1]]`.
#[derive(Debug, Clone)]
struct Groups<T> {
    starts: Vec<usize>,
    items: Vec<T>,
}

impl<T: Copy + PartialEq> Groups<T> {
    /// Groups `(group, item)` pairs given in ascending item order, keeping
    /// each item once per group, so every group stays ascending.
    fn new(group_count: usize, members: impl Iterator<Item = (usize, T)>) -> Self {
        let mut last = vec![None; group_count];
        let members = members
            .filter(|&(group, item)| last[group].replace(item) != Some(item))
            .collect::<Vec<_>>();
        let mut starts = vec![0; group_count + 1];
        for &(group, _) in &members {
            starts[group + 1] += 1;
        }
        for group in 0..group_count {
            starts[group + 1] += starts[group];
        }
        let Some(&(_, first)) = members.first() else {
            return Self {
                starts,
                items: Vec::new(),
            };
        };
        let mut next = starts.clone();
        let mut items = vec![first; members.len()];
        for (group, item) in members {
            items[next[group]] = item;
            next[group] += 1;
        }
        Self { starts, items }
    }

    /// Groups each list in turn, keeping its order.
    fn from_lists<'a>(lists: impl ExactSizeIterator<Item = &'a [T]>) -> Self
    where
        T: 'a,
    {
        let mut starts = Vec::with_capacity(lists.len() + 1);
        starts.push(0);
        let mut items = Vec::new();
        for list in lists {
            items.extend_from_slice(list);
            starts.push(items.len());
        }
        Self { starts, items }
    }

    fn get(&self, group: usize) -> &[T] {
        &self.items[self.starts[group]..self.starts[group + 1]]
    }
}
