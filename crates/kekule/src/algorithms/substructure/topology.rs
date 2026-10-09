use super::*;
use crate::topology::{InstanceAtomId, MoleculeInstanceId, Topology};
use std::ops::ControlFlow;
use std::sync::Arc;

/// Full query mapping bound to an immutable topology snapshot.
#[derive(Debug, Clone)]
pub struct TopologyQueryMatch {
    topology: Arc<Topology>,
    atoms: Vec<InstanceAtomId>,
}
impl TopologyQueryMatch {
    pub fn topology(&self) -> &Arc<Topology> {
        &self.topology
    }
    pub fn atoms(&self) -> &[InstanceAtomId] {
        &self.atoms
    }
    pub fn atom(&self, query_atom: QueryAtomId) -> Option<InstanceAtomId> {
        self.atoms.get(query_atom.index()).copied()
    }
    pub fn tagged_atoms(
        &self,
        query: &TaggedQuery<'_>,
    ) -> Result<Vec<InstanceAtomId>, TaggedMatchError> {
        if self.atoms.len() != query.query().atom_count() {
            return Err(TaggedMatchError::MappingSize);
        }
        Ok(query
            .tags()
            .iter()
            .map(|(_, a)| self.atoms[a.index()])
            .collect())
    }
}

/// Prepared system target. Per-atom chemistry facts are shared by reusable
/// molecule definition; occurrence identities and connectivity remain distinct.
#[derive(Debug)]
pub struct PreparedTopologyTarget<'a> {
    topology: &'a Arc<Topology>,
    data: engine::TargetData<'a>,
    ids: Vec<InstanceAtomId>,
    definitions: Vec<(engine::TargetData<'a>, Vec<MoleculeInstanceId>)>,
}
impl<'a> PreparedTopologyTarget<'a> {
    pub fn new(topology: &'a Arc<Topology>) -> Self {
        let data = engine::TargetData::new(topology.molecules().map(|m| m.molecule()));
        let ids = topology
            .molecules()
            .flat_map(|m| m.atoms().map(|atom| atom.id()))
            .collect();
        let definitions = topology
            .definitions()
            .map(|d| {
                let instances = topology
                    .molecules()
                    .filter(|m| m.definition_id() == d.id())
                    .map(|m| m.id())
                    .collect();
                (engine::TargetData::new([d.molecule()]), instances)
            })
            .collect();
        Self {
            definitions,
            topology,
            data,
            ids,
        }
    }
    pub fn visit_matches(
        &self,
        query: &QueryGraph,
        visitor: impl FnMut(&TopologyQueryMatch) -> ControlFlow<()>,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        self.visit_matches_with_options(query, SubstructureMatchOptions::default(), visitor)
    }

    /// Dot-separated fragments can map to the same or different molecule instances.
    pub fn visit_matches_with_options(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
        mut visitor: impl FnMut(&TopologyQueryMatch) -> ControlFlow<()>,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        if query.is_component_local() {
            validate_options(options)?;
            let mut count = 0usize;
            let mut work = SubstructureMatchWork {
                query_atoms: query.atom_count(),
                target_atoms: self.ids.len(),
                ..Default::default()
            };
            for (data, instances) in &self.definitions {
                let completion = engine::visit_with_work(
                    data,
                    query,
                    SubstructureMatchOptions {
                        max_matches: usize::MAX,
                        ..options
                    },
                    true,
                    &mut |indices| {
                        for &instance in instances {
                            count = count.saturating_add(1);
                            if count > options.max_matches {
                                return false;
                            }
                            let matched = TopologyQueryMatch {
                                topology: self.topology.clone(),
                                atoms: indices
                                    .iter()
                                    .map(|&i| InstanceAtomId::new(instance, data.atoms[i].local))
                                    .collect(),
                            };
                            if visitor(&matched).is_break() {
                                return false;
                            }
                        }
                        true
                    },
                    &mut work,
                )?;
                if count > options.max_matches {
                    work.matches = count;
                    return Err(SubstructureMatchError::ResourceLimit {
                        resource: "matches",
                        observed: count,
                        limit: options.max_matches,
                        work,
                    });
                }
                if completion == MatchCompletion::Stopped {
                    return Ok(completion);
                }
            }
            return Ok(MatchCompletion::Complete);
        }
        engine::visit(&self.data, query, options, true, &mut |indices| {
            visitor(&TopologyQueryMatch {
                topology: self.topology.clone(),
                atoms: indices.iter().map(|&i| self.ids[i]).collect(),
            })
            .is_continue()
        })
    }
    pub fn find_matches(
        &self,
        query: &QueryGraph,
    ) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
        self.find_matches_with_options(query, SubstructureMatchOptions::default())
    }

    /// Every match; exceeding the match cap is an error, never a partial list.
    pub fn find_matches_with_options(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
    ) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
        let mut matches = Vec::new();
        self.visit_matches_with_options(query, options, |m| {
            matches.push(m.clone());
            ControlFlow::Continue(())
        })?;
        Ok(matches)
    }
}

/// Every topology match with default options.
pub fn find_topology_matches(
    topology: &Arc<Topology>,
    query: &QueryGraph,
) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
    PreparedTopologyTarget::new(topology).find_matches(query)
}

/// Every topology match; exceeding the match cap is an error.
pub fn find_topology_matches_with_options(
    topology: &Arc<Topology>,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
    PreparedTopologyTarget::new(topology).find_matches_with_options(query, options)
}
