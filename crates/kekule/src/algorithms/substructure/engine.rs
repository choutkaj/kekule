use super::*;
use crate::core::{BondId, DoubleBondOrientation, RingBasisModel, StereoElementKind};
use crate::query::{AtomPredicate, BondPredicate, QueryBond, QueryStereoConstraint};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

#[derive(Debug)]
pub(super) struct AtomFacts {
    degree: usize,
    hydrogens: usize,
    nongraph_hydrogens: usize,
    valence: usize,
    ring_count: usize,
    ring_size: usize,
    ring_bonds: usize,
}

#[derive(Debug)]
pub(super) struct TargetAtom<'a> {
    pub molecule: &'a Molecule,
    pub local: AtomId,
    pub occurrence: usize,
    pub neighbors: Vec<(usize, BondId)>,
    facts: Arc<AtomFacts>,
}

#[derive(Debug)]
pub(super) struct TargetData<'a> {
    pub atoms: Vec<TargetAtom<'a>>,
}

impl<'a> TargetData<'a> {
    pub fn new(molecules: impl IntoIterator<Item = &'a Molecule>) -> Self {
        let mut atoms = Vec::new();
        let mut facts_cache = BTreeMap::new();
        for (occurrence, molecule) in molecules.into_iter().enumerate() {
            let offset = atoms.len();
            let indices: BTreeMap<_, _> = molecule
                .atom_ids()
                .enumerate()
                .map(|(i, a)| (a, offset + i))
                .collect();
            for (local, atom) in molecule.atoms() {
                let facts = facts_cache
                    .entry((molecule as *const Molecule as usize, local))
                    .or_insert_with(|| {
                        // Unperceived targets retain the existing partial-facts
                        // policy: only the specified contribution is known.
                        let nongraph_hydrogens = molecule
                            .implicit_hydrogens(local)
                            .expect("valid atom")
                            .unwrap_or_else(|| usize::from(atom.hydrogens.specified_count()));
                        let neighbors = molecule
                            .neighbors(local)
                            .expect("valid atom")
                            .collect::<Vec<_>>();
                        let rings = molecule
                            .ring_set()
                            .map(|s| {
                                s.rings()
                                    .iter()
                                    .filter(|r| r.atoms.contains(&local))
                                    .collect::<Vec<_>>()
                            })
                            .unwrap_or_default();
                        Arc::new(AtomFacts {
                            degree: neighbors.len(),
                            nongraph_hydrogens,
                            hydrogens: nongraph_hydrogens
                                + molecule.explicit_hydrogens(local).expect("valid atom"),
                            valence: crate::algorithms::valence::explicit_valence(molecule, local)
                                + nongraph_hydrogens,
                            ring_count: rings.len(),
                            ring_size: rings.iter().map(|r| r.atoms.len()).min().unwrap_or(0),
                            ring_bonds: molecule
                                .incident_bonds(local)
                                .expect("valid atom")
                                .filter(|(b, _)| {
                                    molecule
                                        .ring_membership()
                                        .is_some_and(|r| r.bond_in_ring(*b))
                                })
                                .count(),
                        })
                    })
                    .clone();
                let neighbors = molecule
                    .incident_bonds(local)
                    .expect("valid atom")
                    .map(|(id, b)| (indices[&if b.a() == local { b.b() } else { b.a() }], id))
                    .collect();
                atoms.push(TargetAtom {
                    occurrence,
                    molecule,
                    local,
                    neighbors,
                    facts,
                });
            }
        }
        Self { atoms }
    }
}

pub(super) fn visit(
    data: &TargetData<'_>,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
    complete: bool,
    visitor: &mut dyn FnMut(&[usize]) -> bool,
) -> Result<MatchCompletion, SubstructureMatchError> {
    let mut work = SubstructureMatchWork {
        query_atoms: query.atom_count(),
        target_atoms: data.atoms.len(),
        ..Default::default()
    };
    visit_with_work(data, query, options, complete, visitor, &mut work)
}

pub(super) fn visit_with_work(
    data: &TargetData<'_>,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
    complete: bool,
    visitor: &mut dyn FnMut(&[usize]) -> bool,
    work: &mut SubstructureMatchWork,
) -> Result<MatchCompletion, SubstructureMatchError> {
    validate_options(options)?;
    let mut context = Context {
        data,
        options,
        work: *work,
        candidates: BTreeMap::new(),
        recursive: BTreeMap::new(),
    };
    context.validate(query)?;
    let result = context.run(query, None, complete, visitor);
    *work = context.work;
    result
}

struct Context<'a, 't> {
    data: &'a TargetData<'t>,
    options: SubstructureMatchOptions,
    work: SubstructureMatchWork,
    candidates: BTreeMap<usize, Arc<Vec<Vec<usize>>>>,
    recursive: BTreeMap<(usize, usize), bool>,
}

impl Context<'_, '_> {
    fn limit(
        &self,
        resource: &'static str,
        observed: usize,
        limit: usize,
    ) -> SubstructureMatchError {
        SubstructureMatchError::ResourceLimit {
            resource,
            observed,
            limit,
            work: self.work,
        }
    }
    fn charge(&mut self) -> Result<(), SubstructureMatchError> {
        self.work.search_states = self.work.search_states.saturating_add(1);
        if self.work.search_states > self.options.max_search_states {
            return Err(self.limit(
                "search states",
                self.work.search_states,
                self.options.max_search_states,
            ));
        }
        Ok(())
    }
    fn validate(&self, query: &QueryGraph) -> Result<(), SubstructureMatchError> {
        if query.atom_count() > self.options.max_query_atoms {
            return Err(self.limit(
                "query atoms",
                query.atom_count(),
                self.options.max_query_atoms,
            ));
        }
        let mut molecules = BTreeMap::new();
        for a in &self.data.atoms {
            molecules.insert(a.molecule as *const Molecule as usize, a.molecule);
        }
        for atom in query.atom_ids() {
            for p in query
                .atom(atom)
                .expect("query atom")
                .expression()
                .predicates()
            {
                if let AtomPredicate::Recursive(q) = p {
                    self.validate(q)?;
                }
                for target in molecules.values() {
                    require_atom_perception(target, p)?;
                }
            }
        }
        for bond in query.bond_ids() {
            for p in query
                .bond(bond)
                .expect("query bond")
                .expression()
                .predicates()
            {
                for target in molecules.values() {
                    let required = match p {
                        BondPredicate::Aromatic(_) => Some((
                            target.perception().has_aromaticity(),
                            QueryPerception::Aromaticity,
                        )),
                        BondPredicate::RingMembership(_) => Some((
                            target.perception().has_rings(),
                            QueryPerception::RingMembership,
                        )),
                        _ => None,
                    };
                    if let Some((false, kind)) = required {
                        return Err(SubstructureMatchError::MissingPerception(kind));
                    }
                }
            }
        }
        Ok(())
    }
    fn candidates(
        &mut self,
        query: &QueryGraph,
    ) -> Result<Arc<Vec<Vec<usize>>>, SubstructureMatchError> {
        let key = query as *const QueryGraph as usize;
        if let Some(c) = self.candidates.get(&key) {
            return Ok(c.clone());
        }
        self.work.candidate_pairs = self
            .work
            .candidate_pairs
            .saturating_add(query.atom_count().saturating_mul(self.data.atoms.len()));
        if self.work.candidate_pairs > self.options.max_candidate_pairs {
            return Err(self.limit(
                "candidate pairs",
                self.work.candidate_pairs,
                self.options.max_candidate_pairs,
            ));
        }
        let mut candidates = Vec::new();
        for qa in query.atom_ids() {
            let expression = query.atom(qa).expect("query atom").expression();
            let degree = query.neighbors(qa).expect("query adjacency").count();
            let mut atoms = Vec::new();
            for a in 0..self.data.atoms.len() {
                self.charge()?;
                if self.data.atoms[a].neighbors.len() < degree {
                    continue;
                }
                if expression.possible_with(|p| {
                    if matches!(p, AtomPredicate::Tetrahedral(_)) {
                        Ok(None)
                    } else {
                        self.predicate(p, a).map(Some)
                    }
                })? {
                    atoms.push(a);
                }
            }
            candidates.push(atoms);
        }
        let candidates = Arc::new(candidates);
        self.candidates.insert(key, candidates.clone());
        Ok(candidates)
    }
    fn predicate(
        &mut self,
        p: &AtomPredicate,
        index: usize,
    ) -> Result<bool, SubstructureMatchError> {
        if let AtomPredicate::Recursive(q) = p {
            let key = (&**q as *const QueryGraph as usize, index);
            if let Some(found) = self.recursive.get(&key) {
                return Ok(*found);
            }
            let mut found = false;
            self.run(q, Some(index), false, &mut |_| {
                found = true;
                false
            })?;
            self.recursive.insert(key, found);
            return Ok(found);
        }
        let target = &self.data.atoms[index];
        let m = target.molecule;
        let id = target.local;
        let atom = m.atom(id).expect("target atom");
        let f = &target.facts;
        Ok(match p {
            AtomPredicate::Element(e) => atom.element == *e,
            AtomPredicate::Isotope(i) => atom.isotope.unwrap_or(0) == *i,
            AtomPredicate::FormalCharge(c) => atom.formal_charge == *c,
            AtomPredicate::Aromatic(a) => m.perception().atom_is_aromatic(id) == Some(*a),
            AtomPredicate::Degree(n) => f.degree == usize::from(*n),
            AtomPredicate::TotalConnectivity(n) => {
                f.degree + f.nongraph_hydrogens == usize::from(*n)
            }
            AtomPredicate::TotalHydrogens(n) => f.hydrogens == usize::from(*n),
            AtomPredicate::ImplicitHydrogens(n) => n.map_or(f.nongraph_hydrogens > 0, |n| {
                f.nongraph_hydrogens == usize::from(n)
            }),
            AtomPredicate::TotalValence(n) => f.valence == usize::from(*n),
            AtomPredicate::RingMembership(r) => m
                .ring_membership()
                .is_some_and(|s| s.atom_in_ring(id) == *r),
            AtomPredicate::RingBondCount(n) => f.ring_bonds == usize::from(*n),
            AtomPredicate::RingCount(n) => f.ring_count == usize::from(*n),
            AtomPredicate::SmallestRingSize(n) => f.ring_size == usize::from(*n),
            AtomPredicate::Recursive(_) | AtomPredicate::Tetrahedral(_) => {
                unreachable!("context dependent predicate")
            }
        })
    }
    fn run(
        &mut self,
        query: &QueryGraph,
        anchor: Option<usize>,
        complete: bool,
        visitor: &mut dyn FnMut(&[usize]) -> bool,
    ) -> Result<MatchCompletion, SubstructureMatchError> {
        if query.atom_count() > self.data.atoms.len() {
            return Ok(MatchCompletion::Complete);
        }
        let candidates = self.candidates(query)?;
        if candidates.iter().any(Vec::is_empty) {
            return Ok(MatchCompletion::Complete);
        }
        let mut state = Search {
            query,
            candidates,
            mapping: vec![None; query.atom_count()],
            used: vec![false; self.data.atoms.len()],
            unique: BTreeSet::new(),
            count: 0,
            complete,
            recursive: anchor.is_some(),
            visitor,
        };
        if let Some(a) = anchor {
            if !state.candidates[0].contains(&a) {
                return Ok(MatchCompletion::Complete);
            }
            state.mapping[0] = Some(a);
            state.used[a] = true;
        }
        if self.search(&mut state)? {
            Ok(MatchCompletion::Stopped)
        } else {
            Ok(MatchCompletion::Complete)
        }
    }
    fn search(&mut self, state: &mut Search<'_>) -> Result<bool, SubstructureMatchError> {
        if state.mapping.iter().all(Option::is_some) {
            return self.record(state);
        }
        let qa = state
            .query
            .atom_ids()
            .filter(|a| state.mapping[a.index()].is_none())
            .max_by_key(|a| {
                let neighbors = state
                    .query
                    .neighbors(*a)
                    .expect("query adjacency")
                    .collect::<Vec<_>>();
                (
                    neighbors
                        .iter()
                        .filter(|n| state.mapping[n.index()].is_some())
                        .count(),
                    Reverse(state.candidates[a.index()].len()),
                    neighbors.len(),
                    Reverse(a.index()),
                )
            })
            .expect("unmapped atom");
        let candidates = state.candidates.clone();
        for &target in &candidates[qa.index()] {
            self.charge()?;
            if state.used[target] || !self.feasible(state, qa, target) {
                continue;
            }
            state.mapping[qa.index()] = Some(target);
            state.used[target] = true;
            if self.search(state)? {
                return Ok(true);
            }
            state.mapping[qa.index()] = None;
            state.used[target] = false;
        }
        Ok(false)
    }
    fn feasible(&self, state: &Search<'_>, qa: QueryAtomId, target: usize) -> bool {
        for (_, bond) in state.query.incident_bonds(qa).expect("query adjacency") {
            let other = bond.other_atom(qa);
            if !self.data.atoms[target]
                .neighbors
                .iter()
                .any(|&(neighbor, id)| {
                    if let Some(mapped) = state.mapping[other.index()] {
                        neighbor == mapped
                            && bond_matches(self.data.atoms[target].molecule, id, bond)
                    } else {
                        !state.used[neighbor]
                            && state.candidates[other.index()]
                                .binary_search(&neighbor)
                                .is_ok()
                            && bond_matches(self.data.atoms[target].molecule, id, bond)
                    }
                })
            {
                return false;
            }
        }
        true
    }
    fn directional_matches(
        &mut self,
        query: &QueryGraph,
        mapping: &[usize],
        locals: &[AtomId],
    ) -> Result<bool, SubstructureMatchError> {
        let mut variables = Vec::new();
        for id in query.bond_ids() {
            let bond = query.bond(id).expect("query bond");
            if !bond
                .expression()
                .contains_predicate(|p| matches!(p, BondPredicate::Direction(_)))
            {
                continue;
            }
            let target = &self.data.atoms[mapping[bond.a().index()]];
            let target_bond = target
                .molecule
                .bond_between(target.local, locals[bond.b().index()])
                .expect("valid atoms")
                .expect("mapped bond");
            let mut domain = Vec::new();
            for direction in [None, Some(false), Some(true)] {
                if bond.expression().evaluate_with(|p| match p {
                    BondPredicate::Direction(up) => direction == Some(*up),
                    _ => bond_predicate(target.molecule, target_bond, p),
                }) {
                    domain.push(direction);
                }
            }
            if domain.is_empty() {
                return Ok(false);
            }
            variables.push((id, domain));
        }
        if variables.is_empty() {
            return Ok(true);
        }
        let mut relations = Vec::new();
        for double in query.bond_ids() {
            let bond = query.bond(double).expect("query bond");
            if !bond.expression().contains_predicate(|p| {
                matches!(p, BondPredicate::Order(crate::core::BondOrder::Double))
            }) {
                continue;
            }
            for (i, (left, _)) in variables.iter().enumerate() {
                let l = query.bond(*left).expect("query bond");
                if *left == double || (l.a() != bond.a() && l.b() != bond.a()) {
                    continue;
                }
                for (j, (right, _)) in variables.iter().enumerate() {
                    let r = query.bond(*right).expect("query bond");
                    if *right == double || (r.a() != bond.b() && r.b() != bond.b()) || i == j {
                        continue;
                    }
                    let target = self.data.atoms[mapping[bond.a().index()]].molecule;
                    let together = stereo::matches_constraint(
                        target,
                        query,
                        locals,
                        &QueryStereoConstraint::DoubleBond {
                            bond: double,
                            left_carrier: l.other_atom(bond.a()),
                            right_carrier: r.other_atom(bond.b()),
                            orientation: DoubleBondOrientation::Together,
                        },
                    );
                    let opposite = stereo::matches_constraint(
                        target,
                        query,
                        locals,
                        &QueryStereoConstraint::DoubleBond {
                            bond: double,
                            left_carrier: l.other_atom(bond.a()),
                            right_carrier: r.other_atom(bond.b()),
                            orientation: DoubleBondOrientation::Opposite,
                        },
                    );
                    relations.push((
                        i,
                        j,
                        l.a() != bond.a(),
                        r.a() != bond.b(),
                        together,
                        opposite,
                    ));
                }
            }
        }
        // Iterative finite-domain search keeps stack use independent of bond count.
        // Each attempted assignment consumes the same shared search budget.
        let mut assigned = Vec::new();
        let mut next = vec![0usize; variables.len()];
        loop {
            let depth = assigned.len();
            if depth == variables.len() {
                return Ok(true);
            }
            if next[depth] == variables[depth].1.len() {
                next[depth] = 0;
                if assigned.pop().is_none() {
                    return Ok(false);
                }
                continue;
            }
            self.charge()?;
            let value = variables[depth].1[next[depth]];
            next[depth] += 1;
            assigned.push(value);
            let valid = relations
                .iter()
                .all(|&(i, j, flip_i, flip_j, together, opposite)| {
                    if i >= assigned.len() || j >= assigned.len() {
                        return true;
                    }
                    match (assigned[i], assigned[j]) {
                        (Some(a), Some(b)) => {
                            if (a ^ flip_i) == (b ^ flip_j) {
                                together
                            } else {
                                opposite
                            }
                        }
                        _ => true,
                    }
                });
            if !valid {
                assigned.pop();
            }
        }
    }

    fn record(&mut self, state: &mut Search<'_>) -> Result<bool, SubstructureMatchError> {
        let mapping = state
            .mapping
            .iter()
            .map(|a| a.expect("complete mapping"))
            .collect::<Vec<_>>();
        let locals = mapping
            .iter()
            .map(|&a| self.data.atoms[a].local)
            .collect::<Vec<_>>();
        let enhanced = self.options.use_enhanced_stereo || !state.query.stereo_groups().is_empty();
        let mut relationships = stereo::GroupMatches::default();
        for constraint in state.query.stereo_constraints() {
            let focus = match constraint {
                QueryStereoConstraint::Tetrahedral { center, .. } => center.index(),
                QueryStereoConstraint::DoubleBond { bond, .. } => {
                    state.query.bond(*bond).expect("query bond").a().index()
                }
            };
            let matches = stereo::matches_constraint(
                self.data.atoms[mapping[focus]].molecule,
                state.query,
                &locals,
                constraint,
            );
            if enhanced && matches!(constraint, QueryStereoConstraint::Tetrahedral { .. }) {
                let target = &self.data.atoms[mapping[focus]];
                let mut opposite = constraint.clone();
                if let QueryStereoConstraint::Tetrahedral { orientation, .. } = &mut opposite {
                    *orientation = orientation.inverted();
                }
                let inverse_matches =
                    stereo::matches_constraint(target.molecule, state.query, &locals, &opposite);
                if (!matches && !inverse_matches)
                    || !relationships.add(
                        target.molecule,
                        target.local,
                        target.occurrence,
                        state.query,
                        crate::query::QueryAtomId::new(focus as u32),
                        (matches != inverse_matches).then_some(!matches),
                    )
                {
                    return Ok(false);
                }
            } else if !matches {
                return Ok(false);
            }
        }
        if !self.directional_matches(state.query, &mapping, &locals)? {
            return Ok(false);
        }
        for qa in state.query.atom_ids() {
            let atom = state.query.atom(qa).expect("query atom");
            if !atom
                .expression()
                .contains_predicate(|p| matches!(p, AtomPredicate::Tetrahedral(_)))
            {
                continue;
            }
            let target = &self.data.atoms[mapping[qa.index()]];
            let grouped = state
                .query
                .stereo_groups()
                .iter()
                .any(|g| g.members.contains(&qa));
            let specified=target.molecule.stereo_elements().any(|(_,s)|matches!(&s.kind,StereoElementKind::Tetrahedral(t) if t.center==target.local && t.orientation.is_some()));
            if grouped && !specified {
                return Ok(false);
            }
            let parity = stereo::matches_constraint(
                target.molecule,
                state.query,
                &locals,
                atom.stereo_frame().expect("validated frame"),
            );
            let matched = atom.expression().try_evaluate_with(|p| match p {
                AtomPredicate::Tetrahedral(same) => Ok(specified && parity == *same),
                _ => self.predicate(p, mapping[qa.index()]),
            })?;
            if enhanced && specified {
                let flexible = matches!(atom.stereo_frame(), Some(QueryStereoConstraint::Tetrahedral { carriers, .. }) if carriers.len() < 3);
                let inverted = atom.expression().try_evaluate_with(|p| match p {
                    AtomPredicate::Tetrahedral(same) => {
                        Ok((if flexible { parity } else { !parity }) == *same)
                    }
                    _ => self.predicate(p, mapping[qa.index()]),
                })?;
                if (!matched && !inverted)
                    || ((grouped || matched != inverted)
                        && !relationships.add(
                            target.molecule,
                            target.local,
                            target.occurrence,
                            state.query,
                            qa,
                            (matched != inverted).then_some(!matched),
                        ))
                {
                    return Ok(false);
                }
            } else if !matched {
                return Ok(false);
            }
        }
        if self.options.uniquify && !state.recursive {
            let mut key = mapping.clone();
            key.sort_unstable();
            if !state.unique.insert(key) {
                return Ok(false);
            }
        }
        state.count += 1;
        if !state.recursive {
            self.work.matches = state.count;
            if state.complete && state.count > self.options.max_matches {
                return Err(self.limit("matches", state.count, self.options.max_matches));
            }
        }
        if !(state.visitor)(&mapping) {
            return Ok(true);
        }
        Ok(!state.complete && state.count >= self.options.max_matches)
    }
}

struct Search<'a> {
    query: &'a QueryGraph,
    candidates: Arc<Vec<Vec<usize>>>,
    mapping: Vec<Option<usize>>,
    used: Vec<bool>,
    unique: BTreeSet<Vec<usize>>,
    count: usize,
    complete: bool,
    recursive: bool,
    visitor: &'a mut dyn FnMut(&[usize]) -> bool,
}

fn require_atom_perception(m: &Molecule, p: &AtomPredicate) -> Result<(), SubstructureMatchError> {
    let requirement = match p {
        AtomPredicate::TotalHydrogens(_)
        | AtomPredicate::ImplicitHydrogens(_)
        | AtomPredicate::TotalConnectivity(_)
        | AtomPredicate::TotalValence(_) => {
            Some((m.perception().has_valence(), QueryPerception::Valence))
        }
        AtomPredicate::Aromatic(_) => Some((
            m.perception().has_aromaticity(),
            QueryPerception::Aromaticity,
        )),
        AtomPredicate::RingMembership(_) | AtomPredicate::RingBondCount(_) => {
            Some((m.perception().has_rings(), QueryPerception::RingMembership))
        }
        AtomPredicate::RingCount(_) | AtomPredicate::SmallestRingSize(_) => {
            if m.ring_set().is_some()
                && m.perception().ring_basis_model() != Some(RingBasisModel::FiguerasSssrLike)
            {
                return Err(SubstructureMatchError::IncompatiblePerception(
                    QueryPerception::RingBasis,
                ));
            }
            Some((m.ring_set().is_some(), QueryPerception::RingBasis))
        }
        _ => None,
    };
    if let Some((false, kind)) = requirement {
        return Err(SubstructureMatchError::MissingPerception(kind));
    }
    Ok(())
}

fn bond_matches(target: &Molecule, id: BondId, query: &QueryBond) -> bool {
    query
        .expression()
        .possible_with::<()>(|p| match p {
            BondPredicate::Direction(_) => Ok(None),
            _ => Ok(Some(bond_predicate(target, id, p))),
        })
        .expect("infallible bond predicate")
}
fn bond_predicate(target: &Molecule, id: BondId, p: &BondPredicate) -> bool {
    let bond = target.bond(id).expect("target bond");
    match p {
        BondPredicate::Direction(_) => unreachable!("direction requires a complete mapping"),
        BondPredicate::Order(o) => bond.order == *o,
        BondPredicate::Aromatic(a) => target.perception().bond_is_aromatic(id) == Some(*a),
        BondPredicate::RingMembership(r) => target
            .ring_membership()
            .is_some_and(|s| s.bond_in_ring(id) == *r),
    }
}
