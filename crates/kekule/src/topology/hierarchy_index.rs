//! Hierarchy lookups derived once per published topology.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::dense_order::DenseOrder;
use super::{AtomSiteId, ChainId, Hierarchy, InstanceAtomId, ResidueId, TopologyAtomIndex};

const NO_SITE: u32 = u32::MAX;

/// Lookups between a published hierarchy, dense atoms, and instances, so
/// walks cost what they return.
///
/// Derived at publication from the hierarchy and dense atom order, and shared
/// by every snapshot of that layout. The per-instance lists are built on first
/// use, so publishing an edit costs two passes over plain per-site and
/// per-atom arrays. Instance lists are in ascending ID order.
#[derive(Debug, Clone)]
pub(super) struct HierarchyIndex {
    /// The dense atom of each site, by site ID.
    site_atoms: Vec<TopologyAtomIndex>,
    /// The site of each dense atom, or [`NO_SITE`].
    atom_sites: Vec<u32>,
    instance_count: usize,
    /// The sites, residues and chains of each instance, built on first use.
    instances: OnceLock<InstanceGroups>,
    /// The sites' atom names, numbered on first use.
    names: OnceLock<SiteNames>,
}

/// The sites, residues and chains with an atom of each instance.
#[derive(Debug, Clone)]
struct InstanceGroups {
    sites: Groups<AtomSiteId>,
    residues: Groups<ResidueId>,
    chains: Groups<ChainId>,
}

impl HierarchyIndex {
    /// Indexes `hierarchy`, whose sites all name atoms of `atoms`.
    pub(super) fn new(
        hierarchy: &Hierarchy,
        atoms: &DenseOrder<InstanceAtomId>,
        instance_count: usize,
    ) -> Self {
        let site_atoms = hierarchy
            .site_atom_ids()
            .iter()
            .map(|&atom| {
                let index = atoms
                    .index(atom)
                    .expect("published atom sites name topology atoms");
                TopologyAtomIndex::new(index as u32)
            })
            .collect();
        Self::from_site_atoms(site_atoms, atoms.len(), instance_count)
    }

    /// Indexes a hierarchy from the dense atom of each of its sites.
    pub(super) fn from_site_atoms(
        site_atoms: Vec<TopologyAtomIndex>,
        atom_count: usize,
        instance_count: usize,
    ) -> Self {
        let mut atom_sites = vec![NO_SITE; atom_count];
        for (site, atom) in site_atoms.iter().enumerate() {
            atom_sites[atom.index()] = site as u32;
        }
        Self {
            site_atoms,
            atom_sites,
            instance_count,
            instances: OnceLock::new(),
            names: OnceLock::new(),
        }
    }

    /// The per-instance lists of `hierarchy`, the hierarchy this indexes.
    fn instances(&self, hierarchy: &Hierarchy) -> &InstanceGroups {
        self.instances.get_or_init(|| {
            let instance_count = self.instance_count;
            let site_instances = hierarchy
                .site_atom_ids()
                .iter()
                .map(|atom| atom.instance().index())
                .collect::<Vec<_>>();
            // A residue or chain reaches the instances of its sites.
            let site_instances = &site_instances;
            let residue_instances = |residue: ResidueId| {
                hierarchy
                    .residue(residue)
                    .expect("chains list hierarchy residues")
                    .atom_sites()
                    .iter()
                    .map(move |site| site_instances[site.index()])
            };
            InstanceGroups {
                sites: Groups::new(
                    instance_count,
                    site_instances
                        .iter()
                        .enumerate()
                        .map(|(site, &instance)| (instance, AtomSiteId::new(site as u32))),
                ),
                residues: Groups::new(
                    instance_count,
                    hierarchy.residues().flat_map(|(residue, _)| {
                        residue_instances(residue).map(move |instance| (instance, residue))
                    }),
                ),
                chains: Groups::new(
                    instance_count,
                    hierarchy.chains().flat_map(|(id, chain)| {
                        chain
                            .residues()
                            .iter()
                            .flat_map(|&residue| residue_instances(residue))
                            .map(move |instance| (instance, id))
                    }),
                ),
            }
        })
    }

    /// The dense atom of each site, by site ID.
    pub(super) fn site_atoms(&self) -> &[TopologyAtomIndex] {
        &self.site_atoms
    }

    /// The atom names of `hierarchy`, the hierarchy this indexes, numbered
    /// on first use.
    pub(super) fn names(&self, hierarchy: &Hierarchy) -> &SiteNames {
        self.names.get_or_init(|| SiteNames::new(hierarchy))
    }

    /// The dense atom of a valid site.
    pub(super) fn site_atom(&self, site: AtomSiteId) -> TopologyAtomIndex {
        self.site_atoms[site.index()]
    }

    /// The site of a valid dense atom, or `None` when it has none.
    pub(super) fn atom_site(&self, atom: TopologyAtomIndex) -> Option<AtomSiteId> {
        let site = self.atom_sites[atom.index()];
        (site != NO_SITE).then(|| AtomSiteId::new(site))
    }

    /// The sites of a valid instance's atoms in `hierarchy`, the hierarchy
    /// this indexes.
    pub(super) fn instance_sites(&self, hierarchy: &Hierarchy, instance: usize) -> &[AtomSiteId] {
        self.instances(hierarchy).sites.get(instance)
    }

    /// The residues with a site of a valid instance's atoms.
    pub(super) fn instance_residues(&self, hierarchy: &Hierarchy, instance: usize) -> &[ResidueId] {
        self.instances(hierarchy).residues.get(instance)
    }

    /// The chains with a site of a valid instance's atoms.
    pub(super) fn instance_chains(&self, hierarchy: &Hierarchy, instance: usize) -> &[ChainId] {
        self.instances(hierarchy).chains.get(instance)
    }
}

/// The atom names of a hierarchy's sites as numbers into a table of the
/// distinct names, so a test of names decides each distinct name once.
#[derive(Debug, Clone)]
pub(crate) struct SiteNames {
    /// The distinct names; `None` stands for a site without one.
    table: Vec<Option<Box<str>>>,
    /// Each site's author atom name, else its label atom name, by site ID.
    names: Vec<u32>,
    /// Each site's label atom name, by site ID.
    labels: Vec<u32>,
}

impl SiteNames {
    fn new(hierarchy: &Hierarchy) -> Self {
        fn number<'h>(
            numbers: &mut HashMap<Option<&'h str>, u32>,
            table: &mut Vec<Option<Box<str>>>,
            name: Option<&'h str>,
        ) -> u32 {
            *numbers.entry(name).or_insert_with(|| {
                table.push(name.map(Box::from));
                (table.len() - 1) as u32
            })
        }
        let (mut table, mut numbers) = (Vec::new(), HashMap::new());
        let mut names = Vec::with_capacity(hierarchy.atom_sites().len());
        let mut labels = Vec::with_capacity(hierarchy.atom_sites().len());
        for (_, site) in hierarchy.atom_sites() {
            let label = site.metadata().label_atom_id();
            let name = site.metadata().auth_atom_id().or(label);
            names.push(number(&mut numbers, &mut table, name));
            labels.push(number(&mut numbers, &mut table, label));
        }
        Self {
            table,
            names,
            labels,
        }
    }

    /// The distinct names.
    pub(crate) fn table(&self) -> &[Option<Box<str>>] {
        &self.table
    }

    /// Each site's author atom name, else its label atom name, as a position
    /// in [`Self::table`], by site ID.
    pub(crate) fn names(&self) -> &[u32] {
        &self.names
    }

    /// Each site's label atom name as a position in [`Self::table`], by site
    /// ID.
    pub(crate) fn labels(&self) -> &[u32] {
        &self.labels
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

    fn get(&self, group: usize) -> &[T] {
        &self.items[self.starts[group]..self.starts[group + 1]]
    }
}
