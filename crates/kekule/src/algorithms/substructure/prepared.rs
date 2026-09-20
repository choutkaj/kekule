use super::*;
use std::collections::BTreeSet;
use std::ops::ControlFlow;

/// Reusable facts borrowed from one exact molecule and its installed perception.
/// The borrow prevents chemical edits while this prepared target is live.
/// A parsed or built `QueryGraph` is already a reusable compiled query.
#[derive(Debug)]
pub struct PreparedTarget<'a> {
    data: engine::TargetData<'a>,
}
impl<'a> PreparedTarget<'a> {
    pub fn new(target: &'a Molecule) -> Self {
        Self {
            data: engine::TargetData::new([target]),
        }
    }

    pub(super) fn visit(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
        complete: bool,
        visitor: &mut dyn FnMut(&[AtomId]) -> bool,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        engine::visit(&self.data, query, options, complete, &mut |indices| {
            let atoms = indices
                .iter()
                .map(|&i| self.data.atoms[i].local)
                .collect::<Vec<_>>();
            visitor(&atoms)
        })
    }
    pub fn find_matches_complete(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
    ) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
        let mut matches = Vec::new();
        self.visit(query, options, true, &mut |atoms| {
            matches.push(QueryMatch {
                atoms: atoms.to_vec(),
            });
            true
        })?;
        Ok(matches)
    }
    pub fn visit_matches(
        &self,
        query: &QueryGraph,
        options: SubstructureMatchOptions,
        mut visitor: impl FnMut(&QueryMatch) -> ControlFlow<()>,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        self.visit(query, options, true, &mut |atoms| {
            visitor(&QueryMatch {
                atoms: atoms.to_vec(),
            })
            .is_continue()
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TaggedMatchError {
    MissingTags,
    ZeroTag,
    DuplicateTag(u32),
    MappingSize,
    Match(SubstructureMatchError),
}
impl fmt::Display for TaggedMatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTags => f.write_str("tagged projection requires at least one atom tag"),
            Self::ZeroTag => f.write_str("tagged projection requires positive tags"),
            Self::DuplicateTag(t) => write!(f, "duplicate atom tag {t}"),
            Self::MappingSize => f.write_str("mapping does not have this query's atom count"),
            Self::Match(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for TaggedMatchError {}

/// Validated positive, unique top-level tags in ascending numeric order.
/// Tags inside recursive atom predicates describe local subqueries and are not outputs.
#[derive(Debug)]
pub struct TaggedQuery<'a> {
    query: &'a QueryGraph,
    tags: Vec<(u32, QueryAtomId)>,
}
impl<'a> TaggedQuery<'a> {
    pub fn new(query: &'a QueryGraph) -> Result<Self, TaggedMatchError> {
        let mut tags = query.tagged_atoms().collect::<Vec<_>>();
        tags.sort_unstable();
        if tags.is_empty() {
            return Err(TaggedMatchError::MissingTags);
        }
        if tags[0].0 == 0 {
            return Err(TaggedMatchError::ZeroTag);
        }
        for pair in tags.windows(2) {
            if pair[0].0 == pair[1].0 {
                return Err(TaggedMatchError::DuplicateTag(pair[0].0));
            }
        }
        Ok(Self { query, tags })
    }
    pub fn query(&self) -> &QueryGraph {
        self.query
    }
    pub fn tags(&self) -> &[(u32, QueryAtomId)] {
        &self.tags
    }
    /// Projects a full mapping of this query; callers must supply its matching query.
    pub fn project(&self, matched: &QueryMatch) -> Result<Vec<AtomId>, TaggedMatchError> {
        if matched.atoms.len() != self.query.atom_count() {
            return Err(TaggedMatchError::MappingSize);
        }
        Ok(self
            .tags
            .iter()
            .map(|(_, a)| matched.atoms[a.index()])
            .collect())
    }
    /// Complete enumeration without atom-set deduplication. Optional deduplication
    /// collapses only identical ordered tagged tuples, never tag permutations.
    pub fn find_matches(
        &self,
        target: &PreparedTarget<'_>,
        mut options: SubstructureMatchOptions,
        deduplicate_tuples: bool,
    ) -> Result<Vec<Vec<AtomId>>, TaggedMatchError> {
        options.uniquify = false;
        let mut result = Vec::new();
        let mut seen = BTreeSet::new();
        target
            .visit(self.query, options, true, &mut |atoms| {
                let tuple = self
                    .tags
                    .iter()
                    .map(|(_, a)| atoms[a.index()])
                    .collect::<Vec<_>>();
                if !deduplicate_tuples || seen.insert(tuple.clone()) {
                    result.push(tuple);
                }
                true
            })
            .map_err(TaggedMatchError::Match)?;
        Ok(result)
    }
}
