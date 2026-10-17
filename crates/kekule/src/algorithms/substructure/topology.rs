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
    /// The topology snapshot the match was found in.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// Clones the handle to the matched snapshot, without copying topology data.
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }
    /// Matched atoms in query-atom order.
    pub fn atoms(&self) -> &[InstanceAtomId] {
        &self.atoms
    }
    /// The atom matched by `query_atom`, or `None` when the query has no such
    /// atom.
    pub fn atom(&self, query_atom: QueryAtomId) -> Option<InstanceAtomId> {
        self.atoms.get(query_atom.index()).copied()
    }
    /// The matched atoms of `query`'s tags, in tag order; `query` must be the
    /// matched query.
    ///
    /// # Errors
    ///
    /// [`TaggedMatchError::MappingSize`] when the match does not have the
    /// query's atom count.
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
    /// Prepares every molecule instance of `topology` for repeated queries.
    /// No perception is run.
    pub fn new(topology: &'a Arc<Topology>) -> Self {
        let data = engine::TargetData::new(topology.instances().map(|m| m.molecule()));
        let ids = topology
            .instances()
            .flat_map(|m| m.atoms().map(|atom| atom.id()))
            .collect();
        let definitions = topology
            .definitions()
            .map(|d| {
                let instances = topology
                    .instances()
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
    /// The first match in search order, if any.
    ///
    /// Stops at the first embedding, so the match cap never applies.
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::MissingPerception`] or
    /// [`SubstructureMatchError::IncompatiblePerception`] when the query needs
    /// target perception that is not installed or uses another model, and
    /// [`SubstructureMatchError::LimitExceeded`] when the query exceeds
    /// [`SubstructureMatchOptions::max_query_atoms`] or the search exceeds a
    /// bound of `options`.
    pub fn first(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
    ) -> Result<Option<TopologyQueryMatch>, SubstructureMatchError> {
        let mut first = None;
        self.visit(query, options, |matched| {
            first = Some(matched.clone());
            ControlFlow::Break(())
        })?;
        Ok(first)
    }

    /// Streams matches until the visitor breaks. Dot-separated fragments can
    /// map to the same or different molecule instances. Delivered matches
    /// remain provisional until `Complete` is returned; a resource error may
    /// follow earlier callbacks.
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::MissingPerception`] or
    /// [`SubstructureMatchError::IncompatiblePerception`] when the query needs
    /// target perception that is not installed or uses another model, and
    /// [`SubstructureMatchError::LimitExceeded`] when the query exceeds
    /// [`SubstructureMatchOptions::max_query_atoms`] or the search exceeds a
    /// bound of `options`, including
    /// [`SubstructureMatchOptions::max_matches`].
    pub fn visit(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
        mut visitor: impl FnMut(&TopologyQueryMatch) -> ControlFlow<()>,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        if query.is_component_local() {
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
                    options
                        .with_max_matches(usize::MAX)
                        .expect("usize::MAX is a positive match bound"),
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
                    return Err(SubstructureMatchError::LimitExceeded {
                        limit: crate::diagnostics::LimitExceeded::new(
                            crate::diagnostics::Resource::Matches,
                            options.max_matches,
                            count,
                        ),
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
    /// Every match; exceeding the match cap is an error, never a partial list.
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::MissingPerception`] or
    /// [`SubstructureMatchError::IncompatiblePerception`] when the query needs
    /// target perception that is not installed or uses another model, and
    /// [`SubstructureMatchError::LimitExceeded`] when the query exceeds
    /// [`SubstructureMatchOptions::max_query_atoms`] or the search exceeds a
    /// bound of `options`, including
    /// [`SubstructureMatchOptions::max_matches`].
    pub fn all(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
    ) -> Result<Vec<TopologyQueryMatch>, SubstructureMatchError> {
        let mut matches = Vec::new();
        self.visit(query, options, |m| {
            matches.push(m.clone());
            ControlFlow::Continue(())
        })?;
        Ok(matches)
    }
}
