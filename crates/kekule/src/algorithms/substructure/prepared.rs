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
    /// Prepares `target` for repeated queries. No perception is run.
    pub fn new(target: &'a Molecule) -> Self {
        Self {
            data: engine::TargetData::new([target]),
        }
    }

    pub(super) fn visit_atoms(
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
    ) -> Result<Option<QueryMatch>, SubstructureMatchError> {
        let mut first = None;
        self.visit(query, options, |matched| {
            first = Some(matched.clone());
            ControlFlow::Break(())
        })?;
        Ok(first)
    }

    /// Every match, retaining every embedding unless `uniquify` is requested.
    ///
    /// The match cap is a resource bound: exceeding it returns an error, never
    /// a partial list. To stop early, use [`Self::visit`].
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
    ) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
        let mut matches = Vec::new();
        self.visit_atoms(query, options, true, &mut |atoms| {
            matches.push(QueryMatch {
                atoms: atoms.to_vec(),
            });
            true
        })?;
        Ok(matches)
    }

    /// Streams matches until the visitor breaks. Delivered matches remain
    /// provisional until `Complete` is returned; a resource error may follow
    /// earlier callbacks.
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
        mut visitor: impl FnMut(&QueryMatch) -> ControlFlow<()>,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        self.visit_atoms(query, options, true, &mut |atoms| {
            visitor(&QueryMatch {
                atoms: atoms.to_vec(),
            })
            .is_continue()
        })
    }
}

/// Why a tagged query cannot be prepared or matched.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum TaggedMatchError {
    /// The query has no tagged atom.
    MissingTags,
    /// A tag is zero; tags are positive.
    ZeroTag,
    /// Two atoms carry this tag.
    DuplicateTag(u32),
    /// A mapping does not have the query's atom count.
    MappingSize,
    /// Matching failed; see the source.
    Match(SubstructureMatchError),
}
impl fmt::Display for TaggedMatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingTags => f.write_str("tagged projection requires at least one atom tag"),
            Self::ZeroTag => f.write_str("tagged projection requires positive tags"),
            Self::DuplicateTag(t) => write!(f, "duplicate atom tag {t}"),
            Self::MappingSize => f.write_str("mapping does not have this query's atom count"),
            Self::Match(_) => f.write_str("tagged substructure search failed"),
        }
    }
}
impl std::error::Error for TaggedMatchError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Match(error) => Some(error),
            _ => None,
        }
    }
}

/// Which tagged tuples [`TaggedQuery::find_matches`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TupleMode {
    /// Every tuple of every embedding, including repeated tuples.
    All,
    /// Identical ordered tuples collapse to one; tag permutations remain
    /// distinct.
    Distinct,
}

/// Validated positive, unique top-level tags in ascending numeric order.
/// Tags inside recursive atom predicates describe local subqueries and are not outputs.
#[derive(Debug)]
pub struct TaggedQuery<'a> {
    query: &'a QueryGraph,
    tags: Vec<(u32, QueryAtomId)>,
}
impl<'a> TaggedQuery<'a> {
    /// Collects and validates the top-level tags of `query`.
    ///
    /// # Errors
    ///
    /// [`TaggedMatchError::MissingTags`] when no top-level atom is tagged,
    /// [`TaggedMatchError::ZeroTag`] when a tag is zero, and
    /// [`TaggedMatchError::DuplicateTag`] when two atoms share a tag.
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
    /// The tagged query.
    pub fn query(&self) -> &QueryGraph {
        self.query
    }
    /// Every tag with its query atom, in ascending tag order.
    pub fn tags(&self) -> &[(u32, QueryAtomId)] {
        &self.tags
    }
    /// Projects a full mapping of this query onto its tagged atoms, in tag
    /// order; callers must supply a match of this query.
    ///
    /// # Errors
    ///
    /// [`TaggedMatchError::MappingSize`] when the match does not have this
    /// query's atom count.
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
    /// Complete enumeration without atom-set deduplication.
    /// [`TupleMode::Distinct`] collapses only identical ordered tagged tuples,
    /// never tag permutations. Each tuple lists target atoms in tag order.
    ///
    /// # Errors
    ///
    /// [`TaggedMatchError::Match`] wrapping the errors of
    /// [`PreparedTarget::all`].
    pub fn find_matches(
        &self,
        target: &PreparedTarget<'_>,
        mut options: SubstructureMatchOptions,
        tuples: TupleMode,
    ) -> Result<Vec<Vec<AtomId>>, TaggedMatchError> {
        options.uniquify = false;
        let mut result = Vec::new();
        let mut seen = BTreeSet::new();
        target
            .visit_atoms(self.query, options, true, &mut |atoms| {
                let tuple = self
                    .tags
                    .iter()
                    .map(|(_, a)| atoms[a.index()])
                    .collect::<Vec<_>>();
                if tuples == TupleMode::All || seen.insert(tuple.clone()) {
                    result.push(tuple);
                }
                true
            })
            .map_err(TaggedMatchError::Match)?;
        Ok(result)
    }
}
