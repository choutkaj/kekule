use std::fmt;

use super::{AtomExpression, BondExpression};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QueryAtomId(u32);

impl QueryAtomId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for QueryAtomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "qa{}", self.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct QueryBondId(u32);

impl QueryBondId {
    pub const fn new(raw: u32) -> Self {
        Self(raw)
    }

    pub const fn raw(self) -> u32 {
        self.0
    }

    pub const fn index(self) -> usize {
        self.0 as usize
    }
}

impl fmt::Display for QueryBondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "qb{}", self.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryAtom {
    expression: AtomExpression,
}

impl QueryAtom {
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryGraph {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    adjacency: Vec<Vec<QueryBondId>>,
}

impl QueryGraph {
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
        (0..self.atoms.len()).map(|index| QueryAtomId::new(index as u32))
    }

    pub fn bond_ids(&self) -> impl Iterator<Item = QueryBondId> + '_ {
        (0..self.bonds.len()).map(|index| QueryBondId::new(index as u32))
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

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct QueryGraphBuilder {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    adjacency: Vec<Vec<QueryBondId>>,
}

impl QueryGraphBuilder {
    pub const fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: Vec::new(),
        }
    }

    pub fn with_capacity(atoms: usize, bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(atoms),
            bonds: Vec::with_capacity(bonds),
            adjacency: Vec::with_capacity(atoms),
        }
    }

    pub fn add_atom(&mut self, expression: AtomExpression) -> Result<QueryAtomId, QueryGraphError> {
        let raw = u32::try_from(self.atoms.len()).map_err(|_| QueryGraphError::ResourceLimit {
            resource: "atoms",
            limit: u32::MAX as usize,
        })?;
        let id = QueryAtomId::new(raw);
        self.atoms.push(QueryAtom { expression });
        self.adjacency.push(Vec::new());
        Ok(id)
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
        let raw = u32::try_from(self.bonds.len()).map_err(|_| QueryGraphError::ResourceLimit {
            resource: "bonds",
            limit: u32::MAX as usize,
        })?;
        let id = QueryBondId::new(raw);
        self.bonds.push(QueryBond { a, b, expression });
        self.adjacency[a.index()].push(id);
        self.adjacency[b.index()].push(id);
        Ok(id)
    }

    pub fn build(self) -> Result<QueryGraph, QueryGraphError> {
        if self.atoms.is_empty() {
            return Err(QueryGraphError::EmptyGraph);
        }
        Ok(QueryGraph {
            atoms: self.atoms,
            bonds: self.bonds,
            adjacency: self.adjacency,
        })
    }

    fn validate_atom(&self, id: QueryAtomId) -> Result<(), QueryGraphError> {
        if id.index() < self.atoms.len() {
            Ok(())
        } else {
            Err(QueryGraphError::InvalidAtomId(id))
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryGraphError {
    EmptyGraph,
    InvalidAtomId(QueryAtomId),
    InvalidBondId(QueryBondId),
    SelfBond(QueryAtomId),
    DuplicateBond {
        a: QueryAtomId,
        b: QueryAtomId,
    },
    ResourceLimit {
        resource: &'static str,
        limit: usize,
    },
}

impl fmt::Display for QueryGraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyGraph => f.write_str("query graph must contain at least one atom"),
            Self::InvalidAtomId(id) => write!(f, "invalid query atom id: {id}"),
            Self::InvalidBondId(id) => write!(f, "invalid query bond id: {id}"),
            Self::SelfBond(id) => write!(f, "cannot create a query bond from {id} to itself"),
            Self::DuplicateBond { a, b } => write!(f, "duplicate query bond between {a} and {b}"),
            Self::ResourceLimit { resource, limit } => {
                write!(f, "query graph {resource} limit exceeded: limit {limit}")
            }
        }
    }
}

impl std::error::Error for QueryGraphError {}
