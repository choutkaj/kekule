use std::cmp::Reverse;
use std::collections::BTreeSet;
use std::fmt;

use crate::core::{AtomId, BondId, Molecule};
use crate::query::{
    AtomPredicate, BondPredicate, QueryAtomId, QueryBond, QueryExpression, QueryGraph,
};

/// Absolute query-size ceiling for the recursive bounded matcher.
pub const MAX_SUBSTRUCTURE_QUERY_ATOMS: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubstructureMatchOptions {
    /// Maximum returned matches; reaching this cap stops successfully.
    pub max_matches: usize,
    /// Maximum candidate assignments visited by backtracking.
    pub max_search_states: usize,
    /// Maximum query size accepted by the recursive search.
    pub max_query_atoms: usize,
    /// Maximum query-atom by target-atom compatibility matrix size.
    pub max_candidate_pairs: usize,
    /// Collapse query-automorphism duplicates by target atom set.
    pub uniquify: bool,
}

impl Default for SubstructureMatchOptions {
    fn default() -> Self {
        Self {
            max_matches: 1_000,
            max_search_states: 1_000_000,
            max_query_atoms: MAX_SUBSTRUCTURE_QUERY_ATOMS,
            max_candidate_pairs: 2_000_000,
            uniquify: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SubstructureMatchWork {
    pub query_atoms: usize,
    pub target_atoms: usize,
    pub candidate_pairs: usize,
    pub search_states: usize,
    pub matches: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryMatch {
    atoms: Vec<AtomId>,
}

impl QueryMatch {
    pub fn atom(&self, query_atom: QueryAtomId) -> Option<AtomId> {
        self.atoms.get(query_atom.index()).copied()
    }

    /// Target atoms in query-atom order.
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryPerception {
    Valence,
    RingMembership,
    Aromaticity,
}

impl fmt::Display for QueryPerception {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Valence => f.write_str("valence"),
            Self::RingMembership => f.write_str("ring membership"),
            Self::Aromaticity => f.write_str("aromaticity"),
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstructureMatchError {
    InvalidOptions(&'static str),
    MissingPerception(QueryPerception),
    ResourceLimit {
        resource: &'static str,
        observed: usize,
        limit: usize,
        work: SubstructureMatchWork,
    },
}

impl SubstructureMatchError {
    pub fn work(&self) -> Option<SubstructureMatchWork> {
        match self {
            Self::ResourceLimit { work, .. } => Some(*work),
            Self::InvalidOptions(_) | Self::MissingPerception(_) => None,
        }
    }
}

impl fmt::Display for SubstructureMatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => write!(f, "invalid substructure options: {message}"),
            Self::MissingPerception(perception) => write!(
                f,
                "substructure query requires current {perception} perception on the target"
            ),
            Self::ResourceLimit {
                resource,
                observed,
                limit,
                ..
            } => write!(
                f,
                "substructure {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
        }
    }
}

impl std::error::Error for SubstructureMatchError {}

pub fn find_substructure_match(
    target: &Molecule,
    query: &QueryGraph,
) -> Result<Option<QueryMatch>, SubstructureMatchError> {
    let options = SubstructureMatchOptions {
        max_matches: 1,
        uniquify: false,
        ..SubstructureMatchOptions::default()
    };
    Ok(
        find_substructure_matches_with_options(target, query, options)?
            .into_iter()
            .next(),
    )
}

pub fn find_substructure_matches(
    target: &Molecule,
    query: &QueryGraph,
) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
    find_substructure_matches_with_options(target, query, SubstructureMatchOptions::default())
}

pub fn find_substructure_matches_with_options(
    target: &Molecule,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
    validate_options(options)?;
    let target_atoms = target.atom_ids().collect::<Vec<_>>();
    let mut work = SubstructureMatchWork {
        query_atoms: query.atom_count(),
        target_atoms: target_atoms.len(),
        ..SubstructureMatchWork::default()
    };
    if query.atom_count() > options.max_query_atoms {
        return Err(SubstructureMatchError::ResourceLimit {
            resource: "query atoms",
            observed: query.atom_count(),
            limit: options.max_query_atoms,
            work,
        });
    }
    if query.atom_count() > target_atoms.len() {
        return Ok(Vec::new());
    }
    require_perception(target, query)?;

    let candidate_pairs = query.atom_count().saturating_mul(target_atoms.len());
    work.candidate_pairs = candidate_pairs;
    if candidate_pairs > options.max_candidate_pairs {
        return Err(SubstructureMatchError::ResourceLimit {
            resource: "candidate pairs",
            observed: candidate_pairs,
            limit: options.max_candidate_pairs,
            work,
        });
    }

    let mut candidates = Vec::with_capacity(query.atom_count());
    for query_atom in query.atom_ids() {
        let expression = query
            .atom(query_atom)
            .expect("query atom ids are internally valid")
            .expression();
        let mut compatible = Vec::new();
        for target_atom in target_atoms.iter().copied() {
            if atom_matches(target, target_atom, expression) {
                compatible.push(target_atom);
            }
        }
        if compatible.is_empty() {
            return Ok(Vec::new());
        }
        candidates.push(compatible);
    }

    let mut search = Search {
        target,
        query,
        options,
        candidates,
        mapping: vec![None; query.atom_count()],
        used: BTreeSet::new(),
        unique_atom_sets: BTreeSet::new(),
        matches: Vec::new(),
        work,
    };
    search.visit()?;
    Ok(search.matches)
}

fn validate_options(options: SubstructureMatchOptions) -> Result<(), SubstructureMatchError> {
    for (name, value) in [
        ("max_matches must be greater than zero", options.max_matches),
        (
            "max_search_states must be greater than zero",
            options.max_search_states,
        ),
        (
            "max_query_atoms must be greater than zero",
            options.max_query_atoms,
        ),
        (
            "max_candidate_pairs must be greater than zero",
            options.max_candidate_pairs,
        ),
    ] {
        if value == 0 {
            return Err(SubstructureMatchError::InvalidOptions(name));
        }
    }
    if options.max_query_atoms > MAX_SUBSTRUCTURE_QUERY_ATOMS {
        return Err(SubstructureMatchError::InvalidOptions(
            "max_query_atoms exceeds the stack-safe matcher ceiling",
        ));
    }
    Ok(())
}

fn require_perception(target: &Molecule, query: &QueryGraph) -> Result<(), SubstructureMatchError> {
    let needs_valence = query.atom_ids().any(|id| {
        query
            .atom(id)
            .expect("query atom ids are internally valid")
            .expression()
            .contains_predicate(|predicate| matches!(predicate, AtomPredicate::TotalHydrogens(_)))
    });
    if needs_valence && !target.perception().has_valence() {
        return Err(SubstructureMatchError::MissingPerception(
            QueryPerception::Valence,
        ));
    }

    let needs_rings = query.atom_ids().any(|id| {
        query
            .atom(id)
            .expect("query atom ids are internally valid")
            .expression()
            .contains_predicate(|predicate| matches!(predicate, AtomPredicate::RingMembership(_)))
    }) || query.bond_ids().any(|id| {
        query
            .bond(id)
            .expect("query bond ids are internally valid")
            .expression()
            .contains_predicate(|predicate| matches!(predicate, BondPredicate::RingMembership(_)))
    });
    if needs_rings && !target.perception().has_rings() {
        return Err(SubstructureMatchError::MissingPerception(
            QueryPerception::RingMembership,
        ));
    }

    let needs_aromaticity = query.atom_ids().any(|id| {
        query
            .atom(id)
            .expect("query atom ids are internally valid")
            .expression()
            .contains_predicate(|predicate| matches!(predicate, AtomPredicate::Aromatic(_)))
    }) || query.bond_ids().any(|id| {
        query
            .bond(id)
            .expect("query bond ids are internally valid")
            .expression()
            .contains_predicate(|predicate| matches!(predicate, BondPredicate::Aromatic(_)))
    });
    if needs_aromaticity && !target.perception().has_aromaticity() {
        return Err(SubstructureMatchError::MissingPerception(
            QueryPerception::Aromaticity,
        ));
    }
    Ok(())
}

fn atom_matches(
    target: &Molecule,
    target_atom: AtomId,
    expression: &QueryExpression<AtomPredicate>,
) -> bool {
    let atom = target
        .atom(target_atom)
        .expect("target atom ids are internally valid");
    expression.evaluate_with(|predicate| match predicate {
        AtomPredicate::Element(element) => atom.element == *element,
        AtomPredicate::Isotope(isotope) => atom.isotope == Some(*isotope),
        AtomPredicate::FormalCharge(charge) => atom.formal_charge == *charge,
        AtomPredicate::Aromatic(aromatic) => {
            target.atom_is_aromatic(target_atom).ok().flatten() == Some(*aromatic)
        }
        AtomPredicate::Degree(degree) => {
            target
                .neighbors(target_atom)
                .expect("target atom adjacency is internally valid")
                .count()
                == usize::from(*degree)
        }
        AtomPredicate::TotalHydrogens(hydrogens) => {
            total_hydrogens(target, target_atom) == usize::from(*hydrogens)
        }
        AtomPredicate::RingMembership(in_ring) => target
            .ring_membership()
            .is_some_and(|membership| membership.atom_in_ring(target_atom) == *in_ring),
    })
}

fn total_hydrogens(target: &Molecule, target_atom: AtomId) -> usize {
    let atom = target
        .atom(target_atom)
        .expect("target atom ids are internally valid");
    let graph_hydrogens = target
        .neighbors(target_atom)
        .expect("target atom adjacency is internally valid")
        .filter(|neighbor| {
            target
                .atom(*neighbor)
                .is_ok_and(|atom| atom.element.atomic_number() == 1)
        })
        .count();
    usize::from(atom.explicit_hydrogens)
        + usize::from(
            target
                .implicit_hydrogens(target_atom)
                .ok()
                .flatten()
                .unwrap_or(0),
        )
        + graph_hydrogens
}

fn bond_matches(target: &Molecule, target_bond: BondId, query_bond: &QueryBond) -> bool {
    let bond = target
        .bond(target_bond)
        .expect("target bond ids are internally valid");
    query_bond
        .expression()
        .evaluate_with(|predicate| match predicate {
            BondPredicate::Order(order) => bond.order == *order,
            BondPredicate::Aromatic(aromatic) => {
                target.bond_is_aromatic(target_bond).ok().flatten() == Some(*aromatic)
            }
            BondPredicate::RingMembership(in_ring) => target
                .ring_membership()
                .is_some_and(|membership| membership.bond_in_ring(target_bond) == *in_ring),
        })
}

struct Search<'a> {
    target: &'a Molecule,
    query: &'a QueryGraph,
    options: SubstructureMatchOptions,
    candidates: Vec<Vec<AtomId>>,
    mapping: Vec<Option<AtomId>>,
    used: BTreeSet<AtomId>,
    unique_atom_sets: BTreeSet<Vec<AtomId>>,
    matches: Vec<QueryMatch>,
    work: SubstructureMatchWork,
}

impl Search<'_> {
    fn visit(&mut self) -> Result<bool, SubstructureMatchError> {
        if self.mapping.iter().all(Option::is_some) {
            self.record_match();
            return Ok(self.matches.len() >= self.options.max_matches);
        }

        let query_atom = self.select_next_atom();
        let candidates = self.candidates[query_atom.index()].clone();
        for target_atom in candidates {
            self.work.search_states = self.work.search_states.saturating_add(1);
            if self.work.search_states > self.options.max_search_states {
                return Err(SubstructureMatchError::ResourceLimit {
                    resource: "search states",
                    observed: self.work.search_states,
                    limit: self.options.max_search_states,
                    work: self.work,
                });
            }
            if self.used.contains(&target_atom) || !self.feasible(query_atom, target_atom) {
                continue;
            }
            self.mapping[query_atom.index()] = Some(target_atom);
            self.used.insert(target_atom);
            if self.visit()? {
                return Ok(true);
            }
            self.used.remove(&target_atom);
            self.mapping[query_atom.index()] = None;
        }
        Ok(false)
    }

    fn select_next_atom(&self) -> QueryAtomId {
        self.query
            .atom_ids()
            .filter(|id| self.mapping[id.index()].is_none())
            .max_by_key(|id| {
                let mapped_neighbors = self
                    .query
                    .neighbors(*id)
                    .expect("query adjacency is internally valid")
                    .filter(|neighbor| self.mapping[neighbor.index()].is_some())
                    .count();
                let available_candidates = self.candidates[id.index()]
                    .iter()
                    .filter(|candidate| !self.used.contains(candidate))
                    .count();
                let degree = self
                    .query
                    .neighbors(*id)
                    .expect("query adjacency is internally valid")
                    .count();
                (
                    mapped_neighbors,
                    Reverse(available_candidates),
                    degree,
                    Reverse(id.index()),
                )
            })
            .expect("an incomplete mapping has an unmapped query atom")
    }

    fn feasible(&self, query_atom: QueryAtomId, target_atom: AtomId) -> bool {
        let query_degree = self
            .query
            .neighbors(query_atom)
            .expect("query adjacency is internally valid")
            .count();
        let target_degree = self
            .target
            .neighbors(target_atom)
            .expect("target adjacency is internally valid")
            .count();
        if query_degree > target_degree {
            return false;
        }

        for (_, query_bond) in self
            .query
            .incident_bonds(query_atom)
            .expect("query adjacency is internally valid")
        {
            let query_neighbor = query_bond.other_atom(query_atom);
            if let Some(target_neighbor) = self.mapping[query_neighbor.index()] {
                let Ok(Some(target_bond)) = self.target.bond_between(target_atom, target_neighbor)
                else {
                    return false;
                };
                if !bond_matches(self.target, target_bond, query_bond) {
                    return false;
                }
            } else if !self.has_forward_candidate(target_atom, query_neighbor, query_bond) {
                return false;
            }
        }
        true
    }

    fn has_forward_candidate(
        &self,
        target_atom: AtomId,
        query_neighbor: QueryAtomId,
        query_bond: &QueryBond,
    ) -> bool {
        self.target
            .incident_bonds(target_atom)
            .expect("target adjacency is internally valid")
            .any(|(target_bond, bond)| {
                let target_neighbor = if bond.a() == target_atom {
                    bond.b()
                } else {
                    bond.a()
                };
                !self.used.contains(&target_neighbor)
                    && self.candidates[query_neighbor.index()].contains(&target_neighbor)
                    && bond_matches(self.target, target_bond, query_bond)
            })
    }

    fn record_match(&mut self) {
        let atoms = self
            .mapping
            .iter()
            .map(|atom| atom.expect("complete mapping contains every query atom"))
            .collect::<Vec<_>>();
        if self.options.uniquify {
            let mut atom_set = atoms.clone();
            atom_set.sort_unstable();
            if !self.unique_atom_sets.insert(atom_set) {
                return;
            }
        }
        self.matches.push(QueryMatch { atoms });
        self.work.matches = self.matches.len();
    }
}
