use std::fmt;

use crate::diagnostics::{LimitExceeded, Resource};

/// Why a query stereo constraint or group is invalid.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum QueryStereoIssue {
    /// An atom's carrier frame names a different tetrahedral focus.
    FrameFocusMismatch,
    /// Two constraints share one focus.
    DuplicateFocus,
    /// A tetrahedral predicate has no carrier frame.
    MissingCarrierFrame,
    /// A stereo group has no members.
    EmptyGroup,
    /// A stereo group lists a member twice.
    DuplicateGroupMember,
    /// A member already belongs to another group.
    MemberAlreadyGrouped,
    /// A group member has no tetrahedral carrier frame.
    GroupMemberWithoutFrame,
    /// Tetrahedral carriers do not cover every query neighbor, or exceed
    /// four.
    IncompleteTetrahedralCarriers,
    /// A tetrahedral carrier is repeated.
    DuplicateTetrahedralCarrier,
    /// Double-bond carriers are not distinct substituents at their
    /// endpoints.
    InvalidDoubleBondCarriers,
}

impl fmt::Display for QueryStereoIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::FrameFocusMismatch => "atom frame must have the same tetrahedral focus",
            Self::DuplicateFocus => "duplicate query stereo focus",
            Self::MissingCarrierFrame => "tetrahedral predicate requires a carrier frame",
            Self::EmptyGroup => "empty query stereo group",
            Self::DuplicateGroupMember => "duplicate query stereo group member",
            Self::MemberAlreadyGrouped => "query stereo member already grouped",
            Self::GroupMemberWithoutFrame => {
                "query stereo group member requires a tetrahedral carrier frame"
            }
            Self::IncompleteTetrahedralCarriers => {
                "tetrahedral carriers must cover all query neighbors, at most four"
            }
            Self::DuplicateTetrahedralCarrier => "tetrahedral carriers must be distinct neighbors",
            Self::InvalidDoubleBondCarriers => {
                "double-bond carriers must be distinct substituents at the corresponding endpoints"
            }
        })
    }
}

/// Maximum nesting depth of recursive SMARTS queries in one query graph.
pub const MAX_QUERY_RECURSION_DEPTH: usize = 32;
/// Maximum atoms of one query graph, counting every recursive query.
pub const MAX_QUERY_TOTAL_ATOMS: usize = 4_096;
/// Maximum expression nodes of one query graph, counting every recursive
/// query.
pub const MAX_QUERY_TOTAL_EXPRESSION_NODES: usize = 65_536;

use super::{
    AtomExpression, AtomPredicate, BondExpression, QueryStereoConstraint, QueryStereoGroup,
};

fixed_u32_id!(QueryAtomId, "qa");
fixed_u32_id!(QueryBondId, "qb");

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryAtom {
    expression: AtomExpression,
    tag: Option<u32>,
    stereo_frame: Option<QueryStereoConstraint>,
}

impl QueryAtom {
    /// SMARTS atom-map label. This is output metadata, not a target predicate.
    pub const fn tag(&self) -> Option<u32> {
        self.tag
    }
    pub(crate) fn stereo_frame(&self) -> Option<&QueryStereoConstraint> {
        self.stereo_frame.as_ref()
    }

    pub const fn expression(&self) -> &AtomExpression {
        &self.expression
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryBond {
    a: QueryAtomId,
    b: QueryAtomId,
    expression: BondExpression,
}

impl QueryBond {
    pub const fn a(&self) -> QueryAtomId {
        self.a
    }

    pub const fn b(&self) -> QueryAtomId {
        self.b
    }

    pub const fn endpoints(&self) -> (QueryAtomId, QueryAtomId) {
        (self.a, self.b)
    }

    pub const fn expression(&self) -> &BondExpression {
        &self.expression
    }

    pub(crate) const fn other_atom(&self, atom: QueryAtomId) -> QueryAtomId {
        if self.a.0 == atom.0 {
            self.b
        } else {
            self.a
        }
    }

    fn connects(&self, a: QueryAtomId, b: QueryAtomId) -> bool {
        (self.a == a && self.b == b) || (self.a == b && self.b == a)
    }
}

/// An immutable graph whose vertices and edges carry boolean query expressions.
///
/// Query graphs describe matching predicates rather than represented molecular
/// chemistry. Construct one programmatically with [`QueryGraphBuilder`] or parse
/// SMARTS through [`crate::smarts::parse_str`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryGraph {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    adjacency: Vec<Vec<QueryBondId>>,
    stereo: Vec<QueryStereoConstraint>,
    stereo_groups: Vec<QueryStereoGroup>,
}

impl QueryGraph {
    pub fn stereo_groups(&self) -> &[QueryStereoGroup] {
        &self.stereo_groups
    }

    /// A detached builder preserving all query predicates and stereo relationships.
    pub fn to_builder(&self) -> QueryGraphBuilder {
        QueryGraphBuilder {
            atoms: self.atoms.clone(),
            bonds: self.bonds.clone(),
            adjacency: self.adjacency.clone(),
            stereo: self.stereo.clone(),
            stereo_groups: self.stereo_groups.clone(),
        }
    }
    /// Query atoms carrying a label, in query order (duplicates are retained).
    pub fn tagged_atoms(&self) -> impl Iterator<Item = (u32, QueryAtomId)> + '_ {
        self.atom_ids()
            .filter_map(|id| self.atoms[id.index()].tag.map(|tag| (tag, id)))
    }

    pub(crate) fn validate_complexity(&self) -> Result<(), QueryGraphError> {
        let mut pending = vec![(self, 0usize)];
        let mut atoms = 0usize;
        let mut nodes = 0usize;
        while let Some((query, depth)) = pending.pop() {
            LimitExceeded::check(Resource::RecursionDepth, MAX_QUERY_RECURSION_DEPTH, depth)?;
            atoms = atoms.saturating_add(query.atom_count());
            for atom in &query.atoms {
                nodes = nodes.saturating_add(atom.expression.node_count());
                for p in atom.expression.predicates() {
                    if let AtomPredicate::Recursive(q) = p {
                        pending.push((q, depth + 1));
                    }
                }
            }
            for bond in &query.bonds {
                nodes = nodes.saturating_add(bond.expression.node_count());
            }
            LimitExceeded::check(Resource::Atoms, MAX_QUERY_TOTAL_ATOMS, atoms)?;
            LimitExceeded::check(
                Resource::ExpressionNodes,
                MAX_QUERY_TOTAL_EXPRESSION_NODES,
                nodes,
            )?;
        }
        Ok(())
    }

    /// Whether all outer query atoms belong to one connected component.
    pub fn is_connected(&self) -> bool {
        let mut seen = vec![false; self.atom_count()];
        let mut stack = vec![QueryAtomId::new(0)];
        seen[0] = true;
        while let Some(a) = stack.pop() {
            for b in self.neighbors(a).expect("valid adjacency") {
                if !seen[b.index()] {
                    seen[b.index()] = true;
                    stack.push(b);
                }
            }
        }
        seen.into_iter().all(|s| s)
    }
    pub(crate) fn is_component_local(&self) -> bool {
        self.is_connected()
            && self.atoms.iter().all(|a| {
                a.expression.predicates().iter().all(|p| match p {
                    AtomPredicate::Recursive(q) => q.is_component_local(),
                    _ => true,
                })
            })
    }

    pub fn stereo_constraints(&self) -> &[QueryStereoConstraint] {
        &self.stereo
    }

    pub fn builder() -> QueryGraphBuilder {
        QueryGraphBuilder::new()
    }

    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn atom(&self, id: QueryAtomId) -> Result<&QueryAtom, QueryGraphError> {
        self.atoms
            .get(id.index())
            .ok_or(QueryGraphError::InvalidAtomId(id))
    }

    pub fn bond(&self, id: QueryBondId) -> Result<&QueryBond, QueryGraphError> {
        self.bonds
            .get(id.index())
            .ok_or(QueryGraphError::InvalidBondId(id))
    }

    pub fn atom_ids(&self) -> impl Iterator<Item = QueryAtomId> + '_ {
        (0..=u32::MAX).take(self.atoms.len()).map(QueryAtomId::new)
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = QueryBondId> + '_ {
        (0..=u32::MAX).take(self.bonds.len()).map(QueryBondId::new)
    }

    pub fn incident_bonds(
        &self,
        atom: QueryAtomId,
    ) -> Result<impl Iterator<Item = (QueryBondId, &QueryBond)> + '_, QueryGraphError> {
        self.atom(atom)?;
        Ok(self.adjacency[atom.index()]
            .iter()
            .map(|bond_id| (*bond_id, &self.bonds[bond_id.index()])))
    }

    pub fn neighbors(
        &self,
        atom: QueryAtomId,
    ) -> Result<impl Iterator<Item = QueryAtomId> + '_, QueryGraphError> {
        Ok(self
            .incident_bonds(atom)?
            .map(move |(_, bond)| bond.other_atom(atom)))
    }

    pub fn bond_between(
        &self,
        a: QueryAtomId,
        b: QueryAtomId,
    ) -> Result<Option<QueryBondId>, QueryGraphError> {
        self.atom(a)?;
        self.atom(b)?;
        Ok(self.adjacency[a.index()]
            .iter()
            .copied()
            .find(|bond_id| self.bonds[bond_id.index()].connects(a, b)))
    }
}

/// Transactional builder for an expression-bearing [`QueryGraph`].
///
/// Atoms must be added before bonds that reference them. [`Self::build`]
/// validates the final adjacency and rejects an empty query.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct QueryGraphBuilder {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    adjacency: Vec<Vec<QueryBondId>>,
    stereo: Vec<QueryStereoConstraint>,
    stereo_groups: Vec<QueryStereoGroup>,
}

impl QueryGraphBuilder {
    pub const fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: Vec::new(),
            stereo: Vec::new(),
            stereo_groups: Vec::new(),
        }
    }

    pub fn with_capacity(atoms: usize, bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(atoms),
            bonds: Vec::with_capacity(bonds),
            adjacency: Vec::with_capacity(atoms),
            stereo: Vec::new(),
            stereo_groups: Vec::new(),
        }
    }

    pub fn add_atom(&mut self, expression: AtomExpression) -> Result<QueryAtomId, QueryGraphError> {
        let raw = crate::molecule::checked_raw_id(self.atoms.len()).map_err(|_| {
            QueryGraphError::LimitExceeded(LimitExceeded::first_beyond(
                Resource::Atoms,
                u32::MAX as usize,
            ))
        })?;
        let id = QueryAtomId::new(raw);
        self.atoms.push(QueryAtom {
            expression,
            tag: None,
            stereo_frame: None,
        });
        self.adjacency.push(Vec::new());
        Ok(id)
    }

    /// Sets an optional atom-map label without imposing tagged-projection rules.
    pub fn set_atom_tag(
        &mut self,
        atom: QueryAtomId,
        tag: Option<u32>,
    ) -> Result<(), QueryGraphError> {
        self.validate_atom(atom)?;
        self.atoms[atom.index()].tag = tag;
        Ok(())
    }

    /// Sets the local carrier frame used by Boolean tetrahedral predicates.
    pub fn set_atom_stereo_frame(
        &mut self,
        atom: QueryAtomId,
        mut frame: QueryStereoConstraint,
    ) -> Result<(), QueryGraphError> {
        self.validate_atom(atom)?;
        if !matches!(&frame, QueryStereoConstraint::Tetrahedral { center, .. } if *center == atom) {
            return Err(QueryGraphError::InvalidStereo(
                QueryStereoIssue::FrameFocusMismatch,
            ));
        }
        frame.validate(self.atoms.len(), &self.bonds, &self.adjacency)?;
        frame.canonicalize();
        self.atoms[atom.index()].stereo_frame = Some(frame);
        Ok(())
    }

    pub fn add_bond(
        &mut self,
        a: QueryAtomId,
        b: QueryAtomId,
        expression: BondExpression,
    ) -> Result<QueryBondId, QueryGraphError> {
        self.validate_atom(a)?;
        self.validate_atom(b)?;
        if a == b {
            return Err(QueryGraphError::SelfBond(a));
        }
        if self.adjacency[a.index()]
            .iter()
            .any(|bond_id| self.bonds[bond_id.index()].connects(a, b))
        {
            return Err(QueryGraphError::DuplicateBond { a, b });
        }
        let raw = crate::molecule::checked_raw_id(self.bonds.len()).map_err(|_| {
            QueryGraphError::LimitExceeded(LimitExceeded::first_beyond(
                Resource::Bonds,
                u32::MAX as usize,
            ))
        })?;
        let id = QueryBondId::new(raw);
        self.bonds.push(QueryBond { a, b, expression });
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        Ok(id)
    }

    /// Adds a checked local stereo constraint after its query bonds exist.
    /// Later graph additions are revalidated at publication.
    pub fn add_stereo_constraint(
        &mut self,
        mut constraint: QueryStereoConstraint,
    ) -> Result<(), QueryGraphError> {
        constraint.validate(self.atoms.len(), &self.bonds, &self.adjacency)?;
        if self
            .stereo
            .iter()
            .any(|existing| existing.focus() == constraint.focus())
        {
            return Err(QueryGraphError::InvalidStereo(
                QueryStereoIssue::DuplicateFocus,
            ));
        }
        constraint.canonicalize();
        self.stereo.push(constraint);
        self.stereo.sort_by_key(QueryStereoConstraint::focus);
        Ok(())
    }

    pub fn build(self) -> Result<QueryGraph, QueryGraphError> {
        if self.atoms.is_empty() {
            return Err(QueryGraphError::EmptyGraph);
        }
        for constraint in &self.stereo {
            constraint.validate(self.atoms.len(), &self.bonds, &self.adjacency)?;
        }
        for atom in &self.atoms {
            if let Some(frame) = &atom.stereo_frame {
                frame.validate(self.atoms.len(), &self.bonds, &self.adjacency)?;
            }
            if atom
                .expression
                .contains_predicate(|p| matches!(p, AtomPredicate::Tetrahedral(_)))
                && atom.stereo_frame.is_none()
            {
                return Err(QueryGraphError::InvalidStereo(
                    QueryStereoIssue::MissingCarrierFrame,
                ));
            }
        }
        let graph = QueryGraph {
            atoms: self.atoms,
            bonds: self.bonds,
            adjacency: self.adjacency,
            stereo: self.stereo,
            stereo_groups: self.stereo_groups,
        };
        graph.validate_complexity()?;
        Ok(graph)
    }

    /// Adds one enhanced stereo relationship. No member may occur in two groups.
    pub fn add_stereo_group(&mut self, mut group: QueryStereoGroup) -> Result<(), QueryGraphError> {
        if group.members.is_empty() {
            return Err(QueryGraphError::InvalidStereo(QueryStereoIssue::EmptyGroup));
        }
        group.members.sort_unstable();
        if group.members.windows(2).any(|p| p[0] == p[1]) {
            return Err(QueryGraphError::InvalidStereo(
                QueryStereoIssue::DuplicateGroupMember,
            ));
        }
        for &member in &group.members {
            self.validate_atom(member)?;
            if self
                .stereo_groups
                .iter()
                .any(|g| g.members.contains(&member))
            {
                return Err(QueryGraphError::InvalidStereo(
                    QueryStereoIssue::MemberAlreadyGrouped,
                ));
            }
            let atom = &self.atoms[member.index()];
            let has_frame = atom.stereo_frame.is_some()
                && atom
                    .expression
                    .contains_predicate(|p| matches!(p, AtomPredicate::Tetrahedral(_)));
            let has_constraint = self.stereo.iter().any(|s| matches!(s, QueryStereoConstraint::Tetrahedral { center, .. } if *center == member));
            if !has_frame && !has_constraint {
                return Err(QueryGraphError::InvalidStereo(
                    QueryStereoIssue::GroupMemberWithoutFrame,
                ));
            }
        }
        self.stereo_groups.push(group);
        Ok(())
    }

    fn validate_atom(&self, id: QueryAtomId) -> Result<(), QueryGraphError> {
        if id.index() < self.atoms.len() {
            Ok(())
        } else {
            Err(QueryGraphError::InvalidAtomId(id))
        }
    }
}

/// A query graph edit or build was rejected.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryGraphError {
    /// The graph has no atoms.
    EmptyGraph,
    /// An atom reference is not in the graph.
    InvalidAtomId(QueryAtomId),
    /// A bond reference is not in the graph.
    InvalidBondId(QueryBondId),
    /// A stereo constraint or group is invalid.
    InvalidStereo(QueryStereoIssue),
    /// A bond would connect an atom to itself.
    SelfBond(QueryAtomId),
    /// The two atoms are already bonded.
    DuplicateBond {
        /// The first atom.
        a: QueryAtomId,
        /// The second atom.
        b: QueryAtomId,
    },
    /// The query exceeds a built-in bound: [`MAX_QUERY_RECURSION_DEPTH`],
    /// [`MAX_QUERY_TOTAL_ATOMS`], [`MAX_QUERY_TOTAL_EXPRESSION_NODES`], or the
    /// identifier range.
    LimitExceeded(LimitExceeded),
}

impl From<LimitExceeded> for QueryGraphError {
    fn from(limit: LimitExceeded) -> Self {
        Self::LimitExceeded(limit)
    }
}

impl fmt::Display for QueryGraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGraph => f.write_str("query graph must contain at least one atom"),
            Self::InvalidAtomId(id) => write!(f, "invalid query atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid query bond id: {id}"),
            Self::InvalidStereo(issue) => write!(f, "invalid query stereo: {issue}"),
            Self::SelfBond(id) => write!(f, "cannot create a query bond from {id} to itself"),
            Self::DuplicateBond { a, b } => write!(f, "duplicate query bond between {a} and {b}"),
            Self::LimitExceeded(limit) => write!(f, "query graph {limit}"),
        }
    }
}

impl std::error::Error for QueryGraphError {}
