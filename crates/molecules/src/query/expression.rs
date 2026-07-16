use std::fmt;

use crate::core::{BondOrder, Element};

/// Maximum size of one programmatically constructed query expression.
pub const MAX_QUERY_EXPRESSION_NODES: usize = 4_096;
/// Maximum nesting depth of one programmatically constructed query expression.
pub const MAX_QUERY_EXPRESSION_DEPTH: usize = 64;

/// A syntax-independent boolean expression over query primitives.
///
/// Fields are intentionally private so every value respects the documented
/// node and depth bounds. Parsers are frontends for this representation; they
/// do not add syntax-specific state to concrete atoms or bonds.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryExpression<P> {
    root: ExpressionNode<P>,
    node_count: usize,
    depth: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ExpressionNode<P> {
    Constant(bool),
    Predicate(P),
    Not(Box<ExpressionNode<P>>),
    And(Vec<ExpressionNode<P>>),
    Or(Vec<ExpressionNode<P>>),
}

impl<P> QueryExpression<P> {
    pub const fn always() -> Self {
        Self {
            root: ExpressionNode::Constant(true),
            node_count: 1,
            depth: 1,
        }
    }

    pub const fn never() -> Self {
        Self {
            root: ExpressionNode::Constant(false),
            node_count: 1,
            depth: 1,
        }
    }

    pub const fn predicate(predicate: P) -> Self {
        Self {
            root: ExpressionNode::Predicate(predicate),
            node_count: 1,
            depth: 1,
        }
    }

    pub fn negate(self) -> Result<Self, QueryExpressionError> {
        let root = match self.root {
            ExpressionNode::Constant(value) => ExpressionNode::Constant(!value),
            ExpressionNode::Not(inner) => *inner,
            other => ExpressionNode::Not(Box::new(other)),
        };
        Self::from_root(root)
    }

    pub fn all(expressions: impl IntoIterator<Item = Self>) -> Result<Self, QueryExpressionError> {
        let mut terms = Vec::new();
        for expression in expressions {
            match expression.root {
                ExpressionNode::Constant(true) => {}
                ExpressionNode::Constant(false) => return Ok(Self::never()),
                ExpressionNode::And(children) => terms.extend(children),
                other => terms.push(other),
            }
        }
        match terms.len() {
            0 => Ok(Self::always()),
            1 => Self::from_root(terms.pop().expect("one expression term")),
            _ => Self::from_root(ExpressionNode::And(terms)),
        }
    }

    pub fn any(expressions: impl IntoIterator<Item = Self>) -> Result<Self, QueryExpressionError> {
        let mut terms = Vec::new();
        for expression in expressions {
            match expression.root {
                ExpressionNode::Constant(false) => {}
                ExpressionNode::Constant(true) => return Ok(Self::always()),
                ExpressionNode::Or(children) => terms.extend(children),
                other => terms.push(other),
            }
        }
        match terms.len() {
            0 => Ok(Self::never()),
            1 => Self::from_root(terms.pop().expect("one expression term")),
            _ => Self::from_root(ExpressionNode::Or(terms)),
        }
    }

    pub const fn node_count(&self) -> usize {
        self.node_count
    }

    pub const fn depth(&self) -> usize {
        self.depth
    }

    /// Evaluate the boolean expression using caller-provided primitive logic.
    pub fn evaluate_with(&self, mut predicate: impl FnMut(&P) -> bool) -> bool {
        evaluate_node(&self.root, &mut predicate)
    }

    /// Return whether any primitive in the expression satisfies a predicate.
    ///
    /// Negation does not suppress traversal; this is intended for capability
    /// and prerequisite discovery, not logical evaluation.
    pub fn contains_predicate(&self, mut predicate: impl FnMut(&P) -> bool) -> bool {
        any_predicate(&self.root, &mut predicate)
    }

    fn from_root(root: ExpressionNode<P>) -> Result<Self, QueryExpressionError> {
        let (node_count, depth) = measure(&root);
        if node_count > MAX_QUERY_EXPRESSION_NODES {
            return Err(QueryExpressionError::ResourceLimit {
                resource: "expression nodes",
                observed: node_count,
                limit: MAX_QUERY_EXPRESSION_NODES,
            });
        }
        if depth > MAX_QUERY_EXPRESSION_DEPTH {
            return Err(QueryExpressionError::ResourceLimit {
                resource: "expression depth",
                observed: depth,
                limit: MAX_QUERY_EXPRESSION_DEPTH,
            });
        }
        Ok(Self {
            root,
            node_count,
            depth,
        })
    }
}

fn measure<P>(node: &ExpressionNode<P>) -> (usize, usize) {
    match node {
        ExpressionNode::Constant(_) | ExpressionNode::Predicate(_) => (1, 1),
        ExpressionNode::Not(child) => {
            let (nodes, depth) = measure(child);
            (nodes.saturating_add(1), depth.saturating_add(1))
        }
        ExpressionNode::And(children) | ExpressionNode::Or(children) => {
            let mut nodes = 1usize;
            let mut depth = 0usize;
            for child in children {
                let (child_nodes, child_depth) = measure(child);
                nodes = nodes.saturating_add(child_nodes);
                depth = depth.max(child_depth);
            }
            (nodes, depth.saturating_add(1))
        }
    }
}

fn evaluate_node<P>(node: &ExpressionNode<P>, predicate: &mut impl FnMut(&P) -> bool) -> bool {
    match node {
        ExpressionNode::Constant(value) => *value,
        ExpressionNode::Predicate(value) => predicate(value),
        ExpressionNode::Not(child) => !evaluate_node(child, predicate),
        ExpressionNode::And(children) => {
            children.iter().all(|child| evaluate_node(child, predicate))
        }
        ExpressionNode::Or(children) => {
            children.iter().any(|child| evaluate_node(child, predicate))
        }
    }
}

fn any_predicate<P>(node: &ExpressionNode<P>, predicate: &mut impl FnMut(&P) -> bool) -> bool {
    match node {
        ExpressionNode::Constant(_) => false,
        ExpressionNode::Predicate(value) => predicate(value),
        ExpressionNode::Not(child) => any_predicate(child, predicate),
        ExpressionNode::And(children) | ExpressionNode::Or(children) => {
            children.iter().any(|child| any_predicate(child, predicate))
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum AtomPredicate {
    Element(Element),
    Isotope(u16),
    FormalCharge(i8),
    Aromatic(bool),
    Degree(u8),
    TotalHydrogens(u8),
    RingMembership(bool),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum BondPredicate {
    Order(BondOrder),
    Aromatic(bool),
    RingMembership(bool),
}

pub type AtomExpression = QueryExpression<AtomPredicate>;
pub type BondExpression = QueryExpression<BondPredicate>;

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryExpressionError {
    ResourceLimit {
        resource: &'static str,
        observed: usize,
        limit: usize,
    },
}

impl fmt::Display for QueryExpressionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                observed,
                limit,
            } => write!(
                f,
                "query {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
        }
    }
}

impl std::error::Error for QueryExpressionError {}
