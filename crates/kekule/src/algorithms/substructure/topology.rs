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
            .flat_map(|m| m.atoms().map(|(id, _)| id))
            .collect();
        let definitions = topology
            .definitions()
            .map(|(_, d)| {
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
    /// Dot-separated fragments can map to the same or different molecule instances.
    pub fn visit_matches(
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
    pub fn find_matches_complete(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
    ) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
        let mut matches = Vec::new();
        self.visit_matches(query, options, |m| {
            matches.push(m.clone());
            ControlFlow::Continue(())
        })?;
        Ok(matches)
    }
}

pub fn find_topology_substructure_matches_complete(
    topology: &Arc<Topology>,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
    PreparedTopologyTarget::new(topology).find_matches_complete(query, options)
}
